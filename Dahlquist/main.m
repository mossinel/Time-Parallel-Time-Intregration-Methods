%% script to run Dahlquist parareal analysis
clear all; close all; clc
parent_folder = fileparts(pwd);
addpath(parent_folder);
%%

N=[30];
l=[5];

period=2*pi/l;
T=N*period;    % chosen to cover exactly one period

MF=20*2;
MG =[1,5,10,20,30];     
u0=1;
K=30;
for la=l
    for t=T
        for mg =MG
            [err]=Dahlquist_Parareal_BE(la,t,MF,mg,N,u0,K,"Y");
            [err]=Dahlquist_Parareal_FE(la,t,MF,mg,N,u0,K,"Y");
        end
    end
end

tol=1e-8;
store=zeros(MF,1);
store2=zeros(MF,1);
for mg = 1:MF
    [err1]=Dahlquist_Parareal_BE(lambda_,T,MF,mg,N,u0,K,"N");
    
    if isempty(find(err1<tol))
        store(mg)=K+1;
    else
        step=find(err1<=tol);
        store(mg)=step(1);
    end
end

for mg = 1:MF
    [err1]=Dahlquist_Parareal_FE(lambda_,T,MF,mg,N,u0,K,"N");
    
    if isempty(find(err1<tol))
        store2(mg)=K+1;
    else
        step=find(err1<=tol);
        store2(mg)=step(1);
    end
end

figure
plot([1:MF],store,'r','linewidth',2)
hold on
plot([1:MF],store2,'b','linewidth',2)
legend('BE','FE','Fontsize',18)
xlabel('MG'); ylabel('iteration');
set(gca, 'FontSize', 15)
%% finite step convergence MGRIT
N=[50];
l=5;
lambda_=[-l*i];  
period=2*pi/l;
T=N*period;    % chosen to cover exactly one period
MF=40;
MG =1;     
u0=1;
K=50;

tol=1e-12;
for nu=0:10
    step1(nu+1)=MGRIT_nu_tolerance_test(lambda_,T,MF,1,N,u0,K,nu,tol);
    step2(nu+1)=ceil(N/(nu+1));
end
figure
plot(0:10,step1,'-o','Linewidth',3)
hold on
plot(0:10,step2,'-x','Linewidth',3)
legend("MGRIT F(CF)^{\nu}",'threshold','Fontsize',18)
xlabel('\nu'); ylabel('iteration');
set(gca, 'FontSize', 15)

%% confront MGRIT_F with theta method and Parareal
N=[30];% 500 dà problemi
l=5;
lambda_=[-l*i]; % higher value to increase oscillations 
period=2*pi/l;
T=N*period;    % chosen to cover exactly one period
MF=20*2;
MG =1;     
u0=1;
K=30;

tol=1e-8;
store=zeros(MF,1);
store2=zeros(MF,1);
for mg = 1:MF
    [err1]=Dahlquist_Parareal_BE(lambda_,T,MF,mg,N,u0,K,"N");
    
    if isempty(find(err1<tol))
        store(mg)=K+1;
    else
        step=find(err1<=tol);
        store(mg)=step(1);
    end
end



for mg = 1:MF
    [err2]=Dahlquist_Parareal_theta_BE(lambda_,T,MF,mg,N,u0,K,"N");
    
    if isempty(find(err2<tol))
        store2(mg)=K+1;
    else
        step=find(err2<=tol);
        store2(mg)=step(1);
    end
end
figure
plot([1:MF],store,'r','linewidth',2)
hold on
plot([1:MF],store2,'b','linewidth',2)
legend('BE','theta')
xlabel('M'); ylabel('iteration');
set(gca, 'FontSize', 15)
store3=zeros(MF,1);
store4=zeros(MF,1);
for mg = 1:MF
    [err3]=Dahlquist_Parareal_FE(lambda_,T,MF,mg,N,u0,K,"N");
    
    if isempty(find(err3<tol))
        store3(mg)=K+1;
    else
        step=find(err3<=tol);
        store3(mg)=step(1);
    end
end



for mg = 1:MF
    [err4]=Dahlquist_Parareal_theta_FE(lambda_,T,MF,mg,N,u0,K,"N");
    
    if isempty(find(err4<tol))
        store4(mg)=K+1;
    else
        step=find(err4<=tol);
        store4(mg)=step(1);
    end
end
figure
plot([1:MF],store3,'r','linewidth',2)
hold on
plot([1:MF],store4,'b','linewidth',2)
legend('FE','theta')
xlabel('M'); ylabel('iteration');
set(gca, 'FontSize', 15)

       
