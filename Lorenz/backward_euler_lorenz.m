function [t, U] = backward_euler_lorenz(sigma, beta, rho, tspan, u0, N)

f = @(t, X) [sigma*(X(2)-X(1)); ...
             X(1)*(rho-X(3))-X(2); ...
             X(1)*X(2)-beta*X(3)];

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
    fn1 = f(tn+1, Xn1);
    J = eye(3) - dt * df(tn+1, Xn1, sigma, rho, beta);
    for k = 1:20
        fn = f(tn+1, Xn1);
        Xn1 = Xn + J \ (dt * fn);
        if norm(fn - fn1) < 1e-10
            break;
        end
        fn1 = fn;
    end
    U(n+1, :) = Xn1';
end

function J = df(t, X, sigma, rho, beta)

J = [-sigma, sigma, 0; ...
    rho-X(3), -1, -X(1); ...   
    X(2), X(1), -beta];
end

end
