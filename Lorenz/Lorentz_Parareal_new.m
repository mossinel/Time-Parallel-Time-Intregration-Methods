function [U,u,err,err2]=Lorentz_Parareal_new(sigma, r , b, T,MF,MG,N,u0,K)
%%%% Lorentz system 
%%%% sigma, r, b: parameters
%%%% T: final time 
%%%% MF,MG: F and G time steps
%%%% N: divisions
%%%% u0:initial guess
%%%% K: parareal iterations  
%%%% FORWARD EULER

% equation
f=@(t,x) [sigma*(x(2)-x(1)) r*x(1)-x(2)-x(1)*x(3) x(1)*x(2)-b*x(3)];

% solver (forward euler)
F=@(t0,t1,u0) SForwardEuler(f,t0,t1,u0,MF); % fine solver F
G=@(t0,t1,u0) SForwardEuler(f,t0,t1,u0,MG); % coarse solver G

% call to parareal and "exact" solution
U=Parareal(F,G,T,u0,N,K); % solve with parareal
[t,u]=ForwardEuler(f,[0 T],u0,MF*N); % fine solution

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
legend(["solution T="+T ],["Parareal K="+K],'Fontsize',15)
%t_str=["Lorenz T = "+T+", k ="+K];
%title( t_str)
set(gca, 'FontSize', 15);

end