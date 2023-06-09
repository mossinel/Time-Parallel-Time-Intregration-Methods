function [flag, difference]=check_jacobian(J,func,m,u,t1,t2)
%%%% J: jacobian to evaluate
%%%% func: function to derivate
%%%% m: number of compositions
%%%% u: where to evaluate
%%%% t: time vector
dt=(t2-t1)/m;
dim=size(J); % Lorentz ->3x3
D=zeros(dim(1),dim(2));
for i=1:dim(2)
    e=zeros(1,length(u));
    e(1,i)=sqrt(eps);
    for j=1:dim(1)
        der1=func(t1,u+e);
        der2=func(t1,u-e);
        for k=2:m
            der1=func(t1+dt,der1);
            der2=func(t1+dt,der2);
        end
        D(j,i)=((der1(j)-der2(j))/(2*sqrt(eps)));
    end
end
difference=norm(D-J,inf);
if  difference<1e-6
    flag=1;
else
    flag=0;
end
end