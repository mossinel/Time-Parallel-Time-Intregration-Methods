function [t,u]=Dahlquist_theta(la,tspan,u0,N,theta,m)
% solves Dahlquist’s test equation using theta method ("Multigrid 
% reduction in time for chaotic dynamical systems")

dt=(tspan(2)-tspan(1))/N;
t=(tspan(1):dt:tspan(2))';
u(1,1)=u0;
for n=1:N
u(n+1,:)=u(n,:)*(1+theta*dt*la)/(1-(1-theta)*dt*la);
end