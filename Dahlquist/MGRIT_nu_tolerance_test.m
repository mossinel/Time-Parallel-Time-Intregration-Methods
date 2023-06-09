function [step]=MGRIT_nu_tolerance_test(la,T,MF,MG,N,u0,K,nu,tol)
%%%% Study of the finite step thm
%%%% la: parameter
%%%% T: final time 
%%%% MF,MG: F and G time steps for the backward euler 
%%%% N: divisions
%%%% u0:initial guess
%%%% K: parareal iterations 
%%%% nu: nb of overlapping \Delta T

[t,u]=DahlquistBE(la,[0 T],u0,MF*N); % fine solution
%% MGRIT
dT=T/N; TT=0:dT:T; % coarse time mesh

% Backward Euler
G_=1/(1-la*dT/MG);   
F_=(1-la*dT/MF)^(-MF/MG); 
M_=speye(MG*N+1)-diag(G_*ones(MG*N,1),-1);
A_=speye(MG*N+1)-diag(F_*ones(MG*N,1),-1);
f=zeros(MG*N+1,1);
U{1}(:,1)=zeros(1,length(f)); U{1}(1,1)=u0; % cell 
Y{1}(:,1)=zeros(1,length(f));  % cell 
f(1)=u0;

for k=1:K+1
    
    U{k+1}(1,:)=u0;
    Y{1}=U{k};
    for mu=1:nu
        Y{mu+1}=Y{mu}+f-A_*Y{mu};
    end
    U{k+1}=M_\(M_*Y{end}+f-A_*Y{end});
end
    
%%
for k=1:K+1
    err(k)=max(abs(u(1:MF:end)-U{k+1}(1:MG:N*MG+1,:))); % compute error
end

DT=T/N; R0=abs(1/(1-la*DT)); if R0<1, R0=1; end % compute error bounds

step=find(err<=tol);
if length(step)==0
    step=K+1;
else
    step=step(1);
end
    

    