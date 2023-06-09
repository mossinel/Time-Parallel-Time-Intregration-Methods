function [U,u,err_inf,err_2]=MGRIT_Lorentz_FCF_version(sigma,r ,b, T,MF,MG,N,u0,K)
%%%% Lorentz MGRIT F-relaxation 
%%%% la: parameter
%%%% T: final time 
%%%% MF,MG: F and G time steps for the forward euler 
%%%% N: divisions
%%%% u0:initial guess
%%%% K: parareal iterations  

% equation

ff=@(t,x) [sigma*(x(2)-x(1)) r*x(1)-x(2)-x(1)*x(3) x(1)*x(2)-b*x(3)];
[t,u]=ForwardEuler(ff,[0 T],u0,MF*N); % fine solution

%% method
dT=T/N; TT=0:dT:T; % coarse time mesh
dt=dT/MF; 


% M_=speye(N+1)-diag(ones(N,1),-1);
% coarse grid
U{1}=zeros(N+1,length(u0)); 
U{1}(1,:)=u0; % cell 

V{1}=zeros(N*MF+1,length(u0));
V{1}(1,:)=u0;
H{1}=zeros(N*MF+1,length(u0));
H{1}(1,:)=u0;

f=zeros(N+1,length(u0));
f(1,:)=ff(0,u0);

g=zeros(N+1,length(u0)); % given tearm

tau=zeros(N+1, length(u0));

for n=1:N % initial guess 
    Go(n+1,:)=SForwardEuler(ff,TT(n),TT(n+1),U{1}(n,:),1);
    U{1}(n+1,:)=Go(n+1,:);  
end

for k=1:K 
    
    V{k+1}(1,:)=u0;
    U{k+1}(1,:)=u0;
    
    for n=1:N
        
        % phi^m
        Fn{n+1,:}=SForwardEuler(ff,TT(n), TT(n+1),U{k}(n,:),MF);
        % tau
        tau(n+1,:)=Fn{n+1,:}-SForwardEuler(ff,TT(n),TT(n+1),U{k}(n,:),1);
    end
    
    
    for n=1:N
        U{k+1}(n+1,:)=SForwardEuler(ff,TT(n),TT(n+1),U{k+1}(n,:),1)+tau(n+1,:);
    end

 % FCF relaxation
 temp=zeros(N-1,3);
 for n=1:N-1
        
        [r,H{k+1}((n-1)*MF+1:(n+1)*MF,:)]=ForwardEuler(ff,[TT(n), TT(n+2)-dt],U{k+1}(n,:),2*MF-1);
         
        temp(n,:)=H{k+1}(MF*n+1,:);
        
 end
 for n=1:N-1
        
        U{k+1}(n+1,:)=temp(n,:);       
        
 end
end
%% error computation
DT=T/N;
for k=1:K%+1
    err_inf(k)=max(max(abs(u(1:MF:end,:)-U{k}))); % compute error
    err_2(k)=sqrt(sum(sum(DT*(u(1:MF:end,:)-U{k}).^2)));
end

%% plot
figure
plot3(u(:,1),u(:,2),u(:,3),'-b'... % solution and
,U{K}(:,1),U{K}(:,2),U{K}(:,3),'or'); % parareal iterate
axis([-20 30 -30 40 -10 60]); view([-13,8]);
xlabel('x'); ylabel('y'); zlabel('z');
grid on
%t_str=["Lorenz T = "+T+", k ="+K];
%title( t_str)
legend("u","MGRIT FCF-relaxation")
set(gca, 'FontSize', 10);
%%

figure
    semilogy([1:K],err_inf)
   hold on
    
   semilogy([1:K],err_2)
   legend("norm L^{\infty}","L^2")
    title(["error, T="+T])
    
    figure
    subplot(3,1,1)
    plot([0:T/(length(u)-1):T],u(:,1),'r',[0:T/(length(u)-1):T],u(:,2),'g',[0:T/(length(u)-1):T],u(:,3), ...
    'b',[0:T/(length(U{1})-1):T],U{1}(:,1),'r--',[0:T/(length(U{2})-1):T],U{2}(:,2),'g--',[0:T/(length(U{1})-1):T],U{2}(:,3),'b--')
    
    subplot(3,1,2)
    plot([0:T/(length(u)-1):T],u(:,1),'r',[0:T/(length(u)-1):T],u(:,2),'g',[0:T/(length(u)-1):T],u(:,3), ...
    'b',[0:T/(length(U{K/2})-1):T],U{K/2}(:,1),'r-*',[0:T/(length(U{K/2})-1):T],U{K/2}(:,2),'g-*',[0:T/(length(U{K/2})-1):T],U{K/2}(:,3),'b-*')

    subplot(3,1,3)
    plot([0:T/(length(u)-1):T],u(:,1),'r',[0:T/(length(u)-1):T],u(:,2),'g',[0:T/(length(u)-1):T],u(:,3), ...
    'b',[0:T/(length(U{K})-1):T],U{K}(:,1),'r-*',[0:T/(length(U{K})-1):T],U{K}(:,2),'g-*',[0:T/(length(U{K})-1):T],U{K}(:,3),'b-*')


end
