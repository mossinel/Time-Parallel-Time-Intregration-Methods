function [U,u,err_inf,err_2,Y]=MGRIT_Lorentz_F_Delta_theta2(sigma,r,b,T,MF,N,u0,K,algo,graph,check)
%%%% Lorentz MGRIT F-relaxation, \Delta, \theta
%%%% F relaxation
%%%% sigma,r,b: parameters
%%%% T: final time 
%%%% MF: F and G time steps for the forward euler 
%%%% N: divisions
%%%% u0:initial guess
%%%% K: parareal iterations  
%%%% algo: string which denote "FE" or "BE" method of integration
%%%% graph: string to produce plots
%%%% check: string to check if the jacobian is correct

% equation

ff=@(t,x) [sigma*(x(2)-x(1)) r*x(1)-x(2)-x(1)*x(3) x(1)*x(2)-b*x(3)];
df=@(t,x)[-sigma   sigma  0; r-x(3) -1  -x(1);  x(2)  x(1) -b];
if strcmp(algo,"FE")
    [~,u]=ForwardEuler(ff,[0 T],u0,MF*N); % fine solution
    theta=(MF+1)/(2*MF); %FE
else
    [~,u]=Backward_Euler_nonlinear(df,ff,[0 T],u0,MF*N); 
    theta=(MF-1)/(2*MF); %BE
end

%% definition of the variables

% time division
dT=T/N; TT=0:dT:T;      % coarse 
dt=dT/MF;               % fine

% coarse grid
U{1}=zeros(N+1,length(u0)); 
U{1}(1,:)=u0; % cell 

% fine grid
Y{1}=zeros(MF*N+1,length(u0)); 
Y{1}(1,:)=u0; % cell 

tau=zeros(N+1, length(u0));
Delta{1}=zeros(length(u0) ,length(u0));
%% method
% initial guess 
for n=1:N 
    if strcmp(algo,"FE")
        Go(n+1,:)=SForwardEuler(ff,TT(n),TT(n+1),U{1}(n,:),1);
    else
        Go(n+1,:)=SBackwardEuler(df,ff,TT(n),TT(n+1),U{1}(n,:),1);
    end
    U{1}(n+1,:)=Go(n+1,:);  
end
% iteration
if strcmp(algo,"FE")
% FORWARD EULER
    for k=1:K 
        for n=1:N
            [~,Y{k}((n-1)*MF+1:n*MF,:)]=ForwardEuler(ff,[TT(n),TT(n+1)-dt],U{k}(n,:),MF-1);
            Y{k}(n*MF+1,:)=U{k}(n+1,:);
        end

        U{k+1}(1,:)=u0;
    
        for n=1:N
            % phi^m, last step
            Fn{n+1,:}=SForwardEuler(ff,TT(n+1)-dt, TT(n+1),Y{k}(n*MF,:),1);
            
            % Dphi^m
            DFn{n+1}=SForwardEuler_Delta(df,ff,TT(n),TT(n+1),U{k}(n,:),MF,Y{k}((n-1)*MF+1:n*MF,:));%
            % Delta
            Delta{n+1}=DFn{n+1}-SEuler_Delta_theta(df,ff,TT(n), TT(n+1),U{k}(n,:),theta,check);
          % check
            if strcmp(check,"Y")         
                func=@(tt,uu) uu+dt*theta*ff(tt,uu);
                [flag(n), difference(n)]=check_jacobian(DFn{n+1},func,MF,U{k}(n,:),TT(n),TT(n+1));
                if max(difference)>=1e-5
                    error("!!wrong jacobian!!")
                end
            end
            % tau
            tau(n+1,:)=Fn{n+1,:}-SThetaEuler(df,ff,TT(n),TT(n+1),U{k}(n,:),1,theta)-(Delta{n+1}*U{k}(n,:)')';        
        
         end
    
        for n=1:N
            U{k+1}(n+1,:)=SThetaEuler(df,ff,TT(n),TT(n+1),U{k+1}(n,:),1,theta)+(Delta{n+1}*U{k+1}(n,:)')'+tau(n+1,:);
        end
    end
    %  F relaxation
    for n=1:N
        [~,Y{K+1}((n-1)*MF+1:n*MF,:)]=ForwardEuler(ff,[TT(n),TT(n+1)-dt],U{end}(n,:),MF-1);
        Y{K+1}(n*MF+1,:)=U{end}(n+1,:);
    end
else
% BACKWARD EULER
    for k=1:K 
    
        for n=1:N
            [~,Y{k}((n-1)*MF+1:n*MF,:)]=Backward_Euler_nonlinear(df,ff,[TT(n),TT(n+1)-dt],U{k}(n,:),MF-1);
            Y{k}(n*MF+1,:)=U{k}(n+1,:);

        end
    
        U{k+1}(1,:)=u0;
    
        for n=1:N
            % phi^m, last step
            Fn{n+1,:}=SBackwardEuler(df,ff,TT(n+1)-dt, TT(n+1),Y{k}(n*MF,:),1);
            % Dphi^m
            DFn{n+1}=SBackwardEuler_Delta(df,ff,TT(n),TT(n+1),U{k}(n,:),MF,[Y{k}((n-1)*MF+1:n*MF,:);Fn{n+1,:}]);
            
            % Delta
            Delta{n+1}=DFn{n+1}-SEuler_Delta_theta(df,ff,TT(n), TT(n+1),U{k}(n,:),theta,check);
            if strcmp(check,"Y")         
                func=@(tt,uu) uu+dt*theta*ff(tt,uu);
                [flag(n), difference(n)]=check_jacobian(DFn{n+1},func,MF,U{k}(n,:),TT(n),TT(n+1));
                if max(difference)>=1e-5
                    error("!!wrong jacobian!!")
                end
            end
            % tau
            tau(n+1,:)=Fn{n+1,:}-SThetaEuler(df,ff,TT(n),TT(n+1),U{k}(n,:),1,theta)-(Delta{n+1}*U{k}(n,:)')';        
        end
    
        for n=1:N
            U{k+1}(n+1,:)=SThetaEuler(df,ff,TT(n),TT(n+1),U{k+1}(n,:),1,theta)+(Delta{n+1}*U{k+1}(n,:)')'+tau(n+1,:);
        end
    end
    %  F relaxation
    for n=1:N
        [~,Y{K+1}((n-1)*MF+1:n*MF,:)]=Backward_Euler_nonlinear(df,ff,[TT(n),TT(n+1)-dt],U{end}(n,:),MF-1);
        Y{K+1}(n*MF+1,:)=U{end}(n+1,:);
    end
end
%% error

for k=1:K%+1
    err_inf(k)=max(max(abs(u(1:MF:end,:)-U{k}))); % compute error
    err_2(k)=sqrt(sum(sum(dT*(u(1:MF:end,:)-U{k}).^2)));
    err_fine_inf(k)=max(max(abs(u-Y{k}))); % compute error
    err_fine_2(k)=sqrt(sum(sum(dt*(u-Y{k}).^2)));
end

%% plot
if strcmp(graph,"Y")
    figure
    plot3(u(:,1),u(:,2),u(:,3),'-b'... % solution and
    ,U{K}(:,1),U{K}(:,2),U{K}(:,3),'or'); % parareal iterate
    axis([-20 30 -30 40 -10 60]); view([-13,8]);
    xlabel('x'); ylabel('y'); zlabel('z');
    grid on
    t_str=["Lorenz T = "+T+", k ="+K];
    title( t_str)
    %
    figure
    plot3(u(:,1),u(:,2),u(:,3),'-b'... % solution and
    ,Y{K}(:,1),Y{K}(:,2),Y{K}(:,3),'or'); % parareal iterate
    axis([-20 30 -30 40 -10 60]); view([-13,8]);
    xlabel('x'); ylabel('y'); zlabel('z');
    grid on
    t_str=["Lorenz fine T = "+T+", k ="+K];
    title( t_str)
    %
    figure
    semilogy([1:K],err_inf)
    hold on
    semilogy([1:K],err_2)
    legend("norm L^{\infty}","L^2")
    title(["error coarse, T="+T])
    %
    figure
    semilogy([1:K],err_fine_inf)
    hold on
    semilogy([1:K],err_fine_2)
    legend("norm L^{\infty}","L^2")
    title(["error fine, T="+T])
    %    
    figure
    subplot(3,1,1)
    plot([0:T/(length(u)-1):T],u(:,1),'r',[0:T/(length(u)-1):T],u(:,2),'g',[0:T/(length(u)-1):T],u(:,3), ...
    'b',[0:T/(length(U{1})-1):T],U{1}(:,1),'r-*',[0:T/(length(U{2})-1):T],U{2}(:,2),'g-*',[0:T/(length(U{1})-1):T],U{2}(:,3),'b-*')

    subplot(3,1,2)
    plot([0:T/(length(u)-1):T],u(:,1),'r',[0:T/(length(u)-1):T],u(:,2),'g',[0:T/(length(u)-1):T],u(:,3), ...
    'b',[0:T/(length(U{K/2})-1):T],U{K/2}(:,1),'r-*',[0:T/(length(U{K/2})-1):T],U{K/2}(:,2),'g-*',[0:T/(length(U{K/2})-1):T],U{K/2}(:,3),'b-*')

    subplot(3,1,3)
    plot([0:T/(length(u)-1):T],u(:,1),'r',[0:T/(length(u)-1):T],u(:,2),'g',[0:T/(length(u)-1):T],u(:,3), ...
    'b',[0:T/(length(U{K})-1):T],U{K}(:,1),'r-*',[0:T/(length(U{K})-1):T],U{K}(:,2),'g-*',[0:T/(length(U{K})-1):T],U{K}(:,3),'b-*')

end
end
