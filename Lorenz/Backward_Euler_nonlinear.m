function [t, Y] = Backward_Euler_nonlinear(df,f,tspan, u0, N)
dt = (tspan(2) - tspan(1)) / N;
t = linspace(tspan(1), tspan(2), N+1);

% functions for Euler method
G = @(Y,Yn,n) Y - dt*(f(n*dt+dt,Y)) - Yn;
Gp = @(Y,n) eye(3) - dt*df(n*dt+dt,Y);

% initialization
Y = zeros(N+1,3);
Y(1,:) = u0;
k = zeros(1,N);

for n=1:N
    Yn = Y(n,:);
    Yk = Yn;
    err = 1;
    while (err>1e-10)&&(k(n)<20)
        Ykp1 = Yk - (Gp(Yk,n)\G(Yk,Yn,n)')';
        err = norm(Ykp1-Yk);
        Yk = Ykp1;
        k(n) = k(n)+1;
    end
    if k==20
        disp("!!Maximum nb of iterations reached!!")
    end
    Y(n+1,:) = Ykp1;
end
