function u=SEuler_Delta_theta(df,f,t0,t1,u0,theta,check)
% Function wrapper
[t,u]=Euler_Delta_theta(df,f,[t0 t1],u0,theta,check);
%u=u{end};

