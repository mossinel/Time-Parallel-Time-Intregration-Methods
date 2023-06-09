function u=SBackwardEuler(df,f,t0,t1,u0,n)
% Function wrapper
[t,u]=Backward_Euler_nonlinear(df,f,[t0 t1],u0,n);
u=u(end,:);