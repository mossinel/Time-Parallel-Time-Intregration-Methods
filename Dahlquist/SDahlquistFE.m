function u=SDahlquistFE(la,t0,t1,u0,n)
[t,u]=DahlquistFE(la,[t0 t1],u0,n);
u=u(end);
