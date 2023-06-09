function [t,u]=DahlquistFE(la,tspan,u0,N)
% DAHLQUISTFE solves Dahlquist’s test equation using Forward Euler
% u=DahlquistFE(la,tspan,u0,N); solves Dahlquist’s test equation on
% the time interval tspan using N steps of Forward Euler starting
% with u0 and giving the result in u
dt=(tspan(2)-tspan(1))/N;
t=(tspan(1):dt:tspan(2))';
u(1,1)=u0;
for n=1:N
u(n+1,1)=u(n,1)*(1+dt*la);
end