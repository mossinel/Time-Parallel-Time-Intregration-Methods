function [U,u,err_inf,err_2,Y]=MGRIT_Lorentz_F(sigma,r,b,T,MF,N,u0,K,algo,graph)
%%%% Lorentz MGRIT F-relaxation 
%%%% F relaxation
%%%% sigma,r,b: parameters
%%%% T: final time 
%%%% MF: F and G time steps for the forward euler 
%%%% N: divisions
%%%% u0:initial guess
%%%% K: parareal iterations  
%%%% algo: string which denote "FE" or "BE" method of integration
%%%% graph: string to produce plots

% equation

ff=@(t,x) [sigma*(x(2)-x(1)) r*x(1)-x(2)-x(1)*x(3) x(1)*x(2)-b*x(3)];
df=@(t,x)[-sigma   sigma  0; r-x(3) -1  -x(1);  x(2)  x(1) -b];
if strcmp(algo,"FE")
    [~,u]=ForwardEuler(ff,[0 T],u0,MF*N); % fine solution
else
    [~,u]=Backward_Euler_nonlinear(df,ff,[0 T],u0,MF*N); 
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

%F{1,:}=zeros(MF-1,1);
tau=zeros(N+1, length(u0));

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
            % tau
            tau(n+1,:)=Fn{n+1,:}-SForwardEuler(ff,TT(n),TT(n+1),U{k}(n,:),1);        
      end
    
        for n=1:N
            U{k+1}(n+1,:)=SForwardEuler(ff,TT(n),TT(n+1),U{k+1}(n,:),1)+tau(n+1,:);
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
            % tau
            tau(n+1,:)=Fn{n+1,:}-SBackwardEuler(df,ff,TT(n),TT(n+1),U{k}(n,:),1);        
        end
    
        for n=1:N
            U{k+1}(n+1,:)=SBackwardEuler(df,ff,TT(n),TT(n+1),U{k+1}(n,:),1)+tau(n+1,:);
        end
    end
    %  F relaxation
    for n=1:N
        [~,Y{K+1}((n-1)*MF+1:n*MF,:)]=Backward_Euler_nonlinear(df,ff,[TT(n),TT(n+1)-dt],U{end}(n,:),MF-1);
        Y{K+1}(n*MF+1,:)=U{end}(n+1,:);
    end
end
%% error
fun = @(t,x,y,z) ff(t, [x,y,z]); %

for k=1:K%+1
    err_inf(k)=max(max(abs(u(1:MF:end,:)-U{k}))); % compute error
%    err2(k)=norm(u(1:MF:end,:)-U{k},2);%to change
    err_2(k)=sqrt(sum(sum(dT*(u(1:MF:end,:)-U{k}).^2)));
    err_fine_inf(k)=max(max(abs(u-Y{k}))); % compute error
%    err2(k)=norm(u(1:MF:end,:)-U{k},2);%to change
    err_fine_2(k)=sqrt(sum(sum(dt*(u-Y{k}).^2)));
%     end
end
%% plot

if strcmp(graph,"Y")
    figure
    plot3(u(:,1),u(:,2),u(:,3),'-b'... % solution and
    ,U{K}(:,1),U{K}(:,2),U{K}(:,3),'or','Linewidth',2); % parareal iterate
    axis([-20 30 -30 40 -10 60]); view([-13,8]);
    xlabel('x'); ylabel('y'); zlabel('z');
    grid on
    %t_str=["Lorenz T = "+T+", k ="+K];
    %title( t_str)
    legend("Lorenz","Parareal")
    set(gca, 'FontSize', 15);
    %
    figure
    plot3(u(:,1),u(:,2),u(:,3),'-b'... % solution and
    ,Y{K}(:,1),Y{K}(:,2),Y{K}(:,3),'or','Linewidth',2); % parareal iterate
    axis([-20 30 -30 40 -10 60]); view([-13,8]);
    xlabel('x'); ylabel('y'); zlabel('z');
    grid on
    %t_str=["Lorenz fine T = "+T+", k ="+K];
    %title( t_str)
    %
    legend("Lorenz","Parareal")
    set(gca, 'FontSize', 10);
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
    'b',[0:T/(length(U{1})-1):T],U{1}(:,1),'r-*',[0:T/(length(U{2})-1):T],U{2}(:,2),'g-*',[0:T/(length(U{1})-1):T],U{2}(:,3),'b-*','Linewidth',2)

    subplot(3,1,2)
    plot([0:T/(length(u)-1):T],u(:,1),'r',[0:T/(length(u)-1):T],u(:,2),'g',[0:T/(length(u)-1):T],u(:,3), ...
    'b',[0:T/(length(U{K/2})-1):T],U{K/2}(:,1),'r-*',[0:T/(length(U{K/2})-1):T],U{K/2}(:,2),'g-*',[0:T/(length(U{K/2})-1):T],U{K/2}(:,3),'b-*','Linewidth',2)

    subplot(3,1,3)
    plot([0:T/(length(u)-1):T],u(:,1),'r',[0:T/(length(u)-1):T],u(:,2),'g',[0:T/(length(u)-1):T],u(:,3), ...
    'b',[0:T/(length(U{K})-1):T],U{K}(:,1),'r-*',[0:T/(length(U{K})-1):T],U{K}(:,2),'g-*',[0:T/(length(U{K})-1):T],U{K}(:,3),'b-*','Linewidth',2)

end
end




