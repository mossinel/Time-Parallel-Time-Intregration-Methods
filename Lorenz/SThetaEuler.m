function u=SThetaEuler(df,f,t0,t1,u0,n,theta)
% Function wrapper
[~,u]=theta_nonlinear(df,f,[t0 t1],u0,n,theta);
u=u(end,:);