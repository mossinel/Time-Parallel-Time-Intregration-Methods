function u=SBackwardEuler_Delta(df,f,t0,t1,u0,n,Y)
% Function wrapper
[~,u]=BackwardEuler_Delta(df,f,[t0 t1],u0,n,Y);
%u=u(end,:);