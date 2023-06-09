%% script to run advection parareal analysis
clear all; close all; clc
parent_folder = fileparts(pwd);
addpath(parent_folder);
%%  data 
% \partial_t u + a \partial_x u =0    (-inf,+inf)x(0,T]
f=@(x,t) 0;
a=[0.5,0.7,1,2]; 

T=8; N=16; K=16; J=20; % parareal parameters
dx=1/J; x=0:dx:1; % spatial mesh
u0=sin(2*pi*x); % initial condition

MG=1; % MG no of coarse steps
MF=20; % MF no of fine steps

for ii=1:length(a)
    % solve the problem
    G=@(t0,t1,u0) STransportBE(f,a(ii),[t0 t1],[0 1],u0,MG); % G coarse solver
    F=@(t0,t1,u0) STransportBE(f,a(ii),[t0 t1],[0 1],u0,MF); % F fine solver
    U=Parareal(F,G,T,u0,N,K);
    u=TransportBE(f,a(ii),[0 T],[0 1],u0,N*MF); % fine solution
    dt=T/(MF*N);dT=T/N; t=(0:dt:T)'; TT=(0:dT:T)';
    
    for k=1:K 
        
        up=[]; 
        for n=1:N % reconstruct fine
            up((n-1)*MF+1:n*MF+1,:)=TransportBE(f,a(ii),... % solution from
            [(n-1)*dT n*dT],[0 1],U{k}(n,:),MF); % parareal for plotting
        end
        
        up(N*MF+1,:)=U{k}(end,:);

        if (rem(k,3)==0)
            figure
            %subplot(3,1,1)
            mesh(x,t,u)
            axis([0 1 0 T -1 1])
            %title(["real solution a="+a(ii)])
            figure
            %subplot(3,1,2)
            mesh(x,t,up)
            axis([0 1 0 T -1 1])
            %title(["parareal k= "+k])
            figure
            %subplot(3,1,3)
            mesh(x,t,u-up)
            %title('error')
        end
        mesh(x,t,up); xlabel('x'); ylabel('t'); % plot parareal approx.
        axis([0 1 0 T -1 1])
        mesh(x,t,u-up); xlabel('x'); ylabel('t'); % plot parareal error
        err(ii,k)=max(max(abs(u-up)));
    end
    
end

figure
semilogy([1:K],err(1,:)/err(1,1),[1:K],err(2,:)/err(2,1),[1:K],err(3,:)/err(3,1),[1:K],err(4,:)/err(4,1),'LineWidth',2)
legend(["a="+a(1)],["a="+a(2)],["a="+a(3)],["a="+a(4)],'Fontsize',18)
%title(["N="+N+", K="+K])
xlabel('k'); ylabel('error');
set(gca, 'FontSize', 16);