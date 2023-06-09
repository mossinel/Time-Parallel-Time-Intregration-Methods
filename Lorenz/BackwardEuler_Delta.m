function [t,u_delta]=BackwardEuler_Delta(df,f,tspan,u0,N,u)
dt=(tspan(2)-tspan(1))/N;
t=[tspan(1):dt:tspan(2)];

u_delta=eye(length(u0),length(u0));
Dphi=@(tt,u) eye(length(u0))-dt*df(tt,u);
for n=1:N
    u_delta=Dphi(t(n+1),u(n+1,:))\u_delta;    
end


