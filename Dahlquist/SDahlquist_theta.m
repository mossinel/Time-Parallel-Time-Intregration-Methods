function u=SDahlquist_theta(la,t0,t1,u0,n,theta,m)
[t,u]=Dahlquist_theta(la,[t0 t1],u0,n,theta,m);
 u=u(end);
