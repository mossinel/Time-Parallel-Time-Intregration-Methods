function [err]=Dahlquist_Parareal_theta_BE(la,T, MF,MG,N,u0,K,graph)
%%%% Dahlquist equation, Backward Euler 
%%%% la: parameter
%%%% T: final time 
%%%% MF,MG: F and G time steps for the backward euler 
%%%% N: divisions
%%%% u0:initial guess
%%%% K: parareal iterations  
%%%% graph: flag to generate plots  

f=@(t,x) la*x; % Dahlquist rhs

theta=(MF/MG-1)/(2*MF/MG);
F=@(t0,t1,u0) SDahlquistBE(la,t0,t1,u0,MF); % fine solver F
G=@(t0,t1,u0) SDahlquist_theta(la,t0,t1,u0,MG,theta,MF); % coarse solver G

[U]=Parareal(F,G,T,u0,N,K); % apply parareal
[t,u]=DahlquistBE(la,[0 T],u0,MF*N); % fine solution

%% subplots of the solutions
if strcmp(graph,"Y")
figure 
num=0;
if mod(K,2)
    flag=1;
else
    flag=0;
end

for i = 1:K
   
    if i==5||i==(K-flag)||i==K   
        num=num+1;
        subplot(1,3,num);                       
        plot(t,u,'b','linewidth',2);  
        
        hold on
        plot((0:(T/N):T),U{i},'r--','linewidth',2);
        
        %legend("exact","parareal")
        title(["iteration "+i])

    end  
    
end

r=MG/MF;
t_str1=["Solution \partial_t u =("+la+")u, T="+T+", N="+N+", grid ratio "+r]; 
set(gca, 'FontSize', 10);
sgtitle(t_str1)
end
%% error computation
dT=T/N;
for k=1:K+1
    err(k)=max(abs(u(1:MF:end)-U{k})); % compute error
end

DT=T/N; R0=abs((1+theta*la*DT)/(1-(1-theta)*la*DT)); if R0<1, R0=1; end % compute error bounds
errsup(1)=err(1); errlin(1)=err(1);

%% error comparison
if strcmp(graph,"Y")
for k=1:K
    errsup(k+1)=err(1)*abs(exp(la*DT)-(1+theta*la*DT)/(1-(1-theta)*la*DT))^k/factorial(k)...
    *R0^(N-k-1)*prod(N-(1:k));
    errlin(k+1)=err(1)*(abs(exp(la*DT)-(1+theta*la*DT)/(1-(1-theta)*la*DT))/(1-abs((1+theta*la*DT)/(1-(1-theta)*la*DT))))^k;
end
%% plot
err_diff=abs(diff(errsup));
index=find(err_diff<sqrt(eps)); 
if isempty(index)
    index=[K];
end
figure
semilogy(0:K,err,'--',0:K,errsup,'-',0:K,errlin,'-')
xlim([0 index(1)]);
r=MG/MF;
t_str=["Dahlquist equation \partial_t u =("+la+")u, T="+T+", N="+N+", grid ratio "+r];%" grid ratio "+r+
title(t_str)
xlabel('k'); ylabel('error');
set(gca, 'FontSize', 10);
legend('parareal error','superlinear bound','linear bound')

end