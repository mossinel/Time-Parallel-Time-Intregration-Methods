function [err]=MGRIT_FCF(la,T,MF,MG,N,u0,K)
%%%% Dahlquist equation MGRIT, FCF relaxation
%%%% functioning version which follows the AMG procedure
%%%% la: parameter
%%%% T: final time 
%%%% MF,MG: F and G time steps for the backward euler 
%%%% N: divisions
%%%% u0:initial guess
%%%% K: parareal iterations  

[t,u]=DahlquistBE(la,[0 T],u0,MF*N); % fine solution
%% parareal
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
    Y{k}=U{k}+f-A_*U{k};
    U{k+1}(1,:)=u0;
    U{k+1}=M_\(M_*Y{k}+f-A_*Y{k});
end

%%
figure 
num=0;
for i = 1:K

    if i==5||i==25||i==50           
        num=num+1;
        subplot(1,3,num);                       
        plot(t,u,'b','linewidth',2);  
        plot(t,u,'b');  
        hold on
        plot((0:(T/(N*MG)):T),U{i+1},'r--','linewidth',2);
        plot((0:(T/(N*MG)):T),U{i+1},'r--');
        %legend("exact","parareal")
        title(["iteration "+i])

    end  

end

r=MG/MF;
t_str1=["Solution \partial_t u =("+la+")u, T="+T+", N="+N+", grid ratio "+r]; 
sgtitle(t_str1)
%% error computation

for k=1:K+1
    err(k)=max(abs(u(1:MF:end)-U{k+1}(1:MG:N*MG+1,:))); % compute error
end

DT=T/N; R0=abs(1/(1-la*DT)); if R0<1, R0=1; end % compute error bounds
errsup(1)=err(1); errlin(1)=err(1);

%% error comparison
for k=1:K
    errsup(k+1)=err(1)*abs(exp(la*DT)-1/(1-la*DT))^k/factorial(k)...
    *R0^(N-k-1)*prod(N-(1:k));
    errlin(k+1)=err(1)*(abs(exp(la*DT)-1/(1-la*DT))/(1-abs(1/(1-la*DT))))^k;
end
%% plot
err_diff=abs(diff(errsup));
index=find(err_diff<sqrt(eps)); %or (eps)?
if isempty(index)
    index=[K];
end
figure
semilogy(0:K,err,'--',0:K,errsup,'-',0:K,errlin,'-','Linewidth',3)
xlim([0 index(1)]);
r=MG/MF;
t_str=["Dahlquist equation \partial_t u =("+la+")u, T="+T+", N="+N+", grid ratio "+r];%" grid ratio "+r+
%title(t_str)
xlabel('k'); ylabel('error');
set(gca, 'FontSize', 15)
legend('MGRIT error','superlinear bound','linear bound','Fontsize',18)





