function [t,u_delta]=ForwardEuler_Delta2(df,f,tspan,u0,N)

dt=(tspan(2)-tspan(1))/N;
[t,u]=ForwardEuler(f,tspan,u0,N);
u_delta=eye(length(u0),length(u0));
Dphi=@(tt,u) eye(length(u0))+dt*df(tt,u);
for n=1:N
    u_delta=Dphi(t(n),u(n,:))*u_delta;
end