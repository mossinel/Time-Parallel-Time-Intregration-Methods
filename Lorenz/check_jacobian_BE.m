function [flag, difference]=check_jacobian_BE(J,f,df,m,u,t1,t2,theta)
%%%% J: jacobian to evaluate
%%%% f: Lorentz
%%%% df: JAcobian of lorentz
%%%% m: number of compositions
%%%% u: where to evaluate
%%%% t: time vector

dim=size(J); % Lorentz ->3x3
D=zeros(dim(1),dim(2));
for i=1:dim(2)
    e=zeros(1,length(u));
    e(1,i)=sqrt(eps);
    for j=1:dim(1)
        if theta==0
            [~,der1]=Backward_Euler_nonlinear(df,f,[t1,t2], u+e, m);
            der1=der1(end,:);
            [~,der2]=Backward_Euler_nonlinear(df,f,[t1,t2], u-e, m);
            der2=der2(end,:);
        else
            [~, der1] = theta_nonlinear(df,f,[t1,t2], u+e, m,theta);
            der1=der1(end,:);
            [~, der2] = theta_nonlinear(df,f,[t1,t2], u-e, m,theta);
            der2=der2(end,:);
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