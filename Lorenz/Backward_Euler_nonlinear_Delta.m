function [t, U] = Backward_Euler_nonlinear_Delta(df,f,tspan, u0, N)

% Initialization
dt = (tspan(2) - tspan(1)) / N;
t = linspace(tspan(1), tspan(2), N+1);
U = zeros(N+1, 3);
U(1, :) = u0;

%  Backward Euler
for n = 1:N
    tn = t(n);
    Xn = U(n, :)';
    Xn1 = Xn;
    fn1 = f(tn+1, Xn1)';     
    J = eye(3) - dt * df(tn+1, Xn1);
    for k = 1:20
        fn = f(tn+1, Xn1)';
        Xn1 = Xn + J \ (dt * fn);
        if norm(fn - fn1) < 1e-10
            break;
        end
        fn1 = fn;
    end
    U(n+1, :) = Xn1';
end
u_delta{1}=zeros(length(u0),length(u0));
for n=1:N
    u_delta{n+1}=u_delta{n}+dt*df(t(1+n),U(n+1,:)); 
end
end

