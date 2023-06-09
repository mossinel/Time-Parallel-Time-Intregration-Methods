function [U,u,err,err2]=Lorentz_Parareal_bound(sigma, r , b, T,MF,MG,N,u0,K,B,graph,fine_bdd)
%%%% Lorentz system, truncated
%%%% sigma, r, b: parameters
%%%% T: final time 
%%%% MF,MG: F and G time steps
%%%% N: divisions
%%%% u0:initial guess
%%%% K: parareal iterations 
%%%% B: bound to avoid Lorentz explosion
%%%% graph : flag to see trajectories of each component
%%%% fine_bdd : flag to have also the fine integrator with bounded input
%%%% FORWARD EULER

% equation
f=@(t,x) [sigma*(x(2)-x(1)) r*x(1)-x(2)-x(1)*x(3) x(1)*x(2)-b*x(3)];
fc=@(t,x)f(t,min(B,norm(x,inf))*x/norm(x,inf));
% solver (forward euler)
if strcmp(fine_bdd,'Y')
    F=@(t0,t1,u0) SForwardEuler(fc,t0,t1,u0,MF); % fine solver F
else
    F=@(t0,t1,u0) SForwardEuler(f,t0,t1,u0,MF);
end
G=@(t0,t1,u0) SForwardEuler(fc,t0,t1,u0,MG); % coarse solver G

% call to parareal and "exact" solution
U=Parareal(F,G,T,u0,N,K); % solve with parareal
[~,u]=ForwardEuler(f,[0 T],u0,MF*N); % fine solution

%% error computation
DT=T/N;
for k=1:K+1
    err(k)=max(max(abs(u(1:MF:end,:)-U{k}))); % compute error
    err2(k)=sqrt(sum(sum(DT*(u(1:MF:end,:)-U{k}).^2)));
end
 

%% plot
figure
plot3(u(:,1),u(:,2),u(:,3),'-b'... % solution and
,U{K}(:,1),U{K}(:,2),U{K}(:,3),'or','Linewidth',2); % parareal iterate
axis([-20 30 -30 40 -10 60]); view([-13,8]);
xlabel('x'); ylabel('y'); zlabel('z');
grid on
legend(["solution T=4 T_{\lambda}"],["Parareal K="+K],'Fontsize',15)
set(gca, 'FontSize', 15);
%t_str=["Lorenz T = "+T+", k ="+K];
%title( t_str)
if strcmp(graph,'Y')
    figure
    semilogy([0:K],err,'Linewidth',2)
    hold on
    
    semilogy([0:K],err2,'Linewidth',2)
    legend("norm L^{\infty}","L^2",'Fontsize',15)
    title(["error, T="+T])
    xlabel('iteration')
    set(gca, 'FontSize', 15);
    
    figure
    subplot(3,1,1)
    plot([0:T/(length(u)-1):T],u(:,1),'r',[0:T/(length(u)-1):T],u(:,2),'g',[0:T/(length(u)-1):T],u(:,3), ...
    'b',[0:T/(length(U{1})-1):T],U{1}(:,1),'ro-',[0:T/(length(U{1})-1):T],U{1}(:,2),'go-',[0:T/(length(U{1})-1):T],U{1}(:,3),'bo-','Linewidth',2)
    xlabel('time') 
    ylabel('x')
    set(gca, 'FontSize', 10);

    subplot(3,1,2)
    plot([0:T/(length(u)-1):T],u(:,1),'r',[0:T/(length(u)-1):T],u(:,2),'g',[0:T/(length(u)-1):T],u(:,3), ...
        'b',[0:T/(length(U{K/2})-1):T],U{K/2}(:,1),'ro-',[0:T/(length(U{K/2})-1):T],U{K/2}(:,2),'go-',[0:T/(length(U{K/2})-1):T],U{K/2}(:,3),'bo-','Linewidth',2)
     xlabel('time') 
    ylabel('x')
    set(gca, 'FontSize', 10);

    subplot(3,1,3)
    plot([0:T/(length(u)-1):T],u(:,1),'r',[0:T/(length(u)-1):T],u(:,2),'g',[0:T/(length(u)-1):T],u(:,3), ...
    'b',[0:T/(length(U{K})-1):T],U{K}(:,1),'ro-',[0:T/(length(U{K})-1):T],U{K}(:,2),'go-',[0:T/(length(U{K})-1):T],U{K}(:,3),'bo-','Linewidth',2)
     xlabel('time') 
    ylabel('x')
    set(gca, 'FontSize', 10);
end

end