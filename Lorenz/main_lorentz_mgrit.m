%% script to run Lorentz parareal analysis
clear all; close all; clc

parent_folder = fileparts(pwd);
addpath(parent_folder);
%% data 

sigma=10;r=28;b=8/3;
T_lambda=log(10)/0.9;
T=4*T_lambda;
MF=20;
MG=1;
N=1024;
u0=[20;5;-5];
K=50;
err_inf=zeros(4,K);
err_2=zeros(4,K);
residual_inf=zeros(4,K);

%% comparison F vs FCF
N=[512];
T=[4*T_lambda];
tol=1e-8;

[~,~,err_inf(1,:),err_2(1,:)]=MGRIT_Lorentz_F(sigma,r,b,T,MF,N,u0,K,"FE","Y");
[~,~,err_inf(2,:),err_2(2,:)]=MGRIT_Lorentz_FCF_version(sigma,r ,b, T,MF,MG,N,u0,K);

treshold_inf=zeros(2);
treshold_L2=zeros(2);
for i=linspace(1,2,2)
            if isempty(find(err_inf(i,:)<=tol, 1))
                treshold_inf(i)=K+1;
            else
                temp=find(err_inf(i,:)<=tol);
                treshold_inf(i)=temp(1);
            end
end
        figure
        semilogy([1:K],err_inf(1,:),'-*','Linewidth',2)
        hold on
        semilogy([1:K],err_inf(2,:),'-*','Linewidth',2)
       
        legend("FE, \theta method","BE, \theta method",'Fontsize',15)
        %t_str=["FE T = "+T(ii)+", N ="+N(ii)];
        %title(t_str)
        xlabel('iteration'); ylabel('error L^{\infty}');
        set(gca, 'FontSize', 15);

%% comparison Lyapunov times and N
%BE
N=[512*2,512*4,512*8,512*16];
T=[4*T_lambda,8*T_lambda,16*T_lambda,32*T_lambda];
treshold_inf=zeros(4,4);
treshold_L2=zeros(4,4);
for ii=1:4
tol=1e-8;
        [~,~,err_inf(1,:),err_2(1,:)]=MGRIT_Lorentz_F(sigma,r,b,T(ii),MF,N(ii),u0,K,"BE","N");

        [~,~,err_inf(2,:),err_2(2,:)]=MGRIT_Lorentz_F_Delta2(sigma,r,b,T(ii),MF,N(ii),u0,K,"BE","N","N");

        [~,~,err_inf(3,:),err_2(3,:)]=MGRIT_Lorentz_F_theta(sigma,r ,b, T(ii),MF,N(ii),u0,K,"BE","N");

        [~,~,err_inf(4,:),err_2(4,:),Y]=MGRIT_Lorentz_F_Delta_theta2(sigma,r,b,T(ii),MF,N(ii),u0,K,"BE","N","N");

        figure

        semilogy([1:K],err_inf(1,:),'-*','Linewidth',2)
        hold on
        semilogy([1:K],err_inf(2,:),'-*','Linewidth',2)
        semilogy([1:K],err_inf(3,:),'-*','Linewidth',2)
        semilogy([1:K],err_inf(4,:),'-d','Linewidth',2)
       
        legend("Parareal","Delta","theta","Delta-theta",'Fontsize',12)
        xlabel('iteration'); ylabel('error L^{\infty}');
        set(gca, 'FontSize', 15);
        
        for i=linspace(1,4,4)
            if isempty(find(err_inf(i,:)<=tol, 1))
                treshold_inf(i,ii)=K+1;
            else
                temp=find(err_inf(i,:)<=tol);
                treshold_inf(i,ii)=temp(1);
            end
        end
        
        for i=linspace(1,4,4)
            if isempty(find(err_2(i,:)<=tol, 1))
                treshold_L2(i,ii)=K+1;
            else
                temp=find(err_2(i,:)<=tol);
                treshold_L2(i,ii)=temp(1);
            end
        end
end

% comparison, FE
N=[512*2,512*4,512*8];
%N=[1,2,3,4]
T=[4*T_lambda,8*T_lambda,16*T_lambda,32*T_lambda];
%T=[0.1,0.2,0.3,0.4];
treshold_inf=zeros(4,3);
treshold_L2=zeros(4,3);
for ii=1:3
tol=1e-8;
        [~,~,err_inf(1,:),err_2(1,:)]=MGRIT_Lorentz_F(sigma,r,b,T(ii),MF,N(ii),u0,K,"FE","N");

        [~,~,err_inf(2,:),err_2(2,:)]=MGRIT_Lorentz_F_Delta2(sigma,r,b,T(ii),MF,N(ii),u0,K,"FE","N","N");

        [~,~,err_inf(3,:),err_2(3,:)]=MGRIT_Lorentz_F_theta(sigma,r ,b, T(ii),MF,N(ii),u0,K,"FE","N");

        [~,~,err_inf(4,:),err_2(4,:),Y]=MGRIT_Lorentz_F_Delta_theta2(sigma,r,b,T(ii),MF,N(ii),u0,K,"FE","N","N");

        figure        

        semilogy([1:K],err_inf(1,:),'-*','Linewidth',2)
        hold on
        semilogy([1:K],err_inf(2,:),'-*','Linewidth',2)
        semilogy([1:K],err_inf(3,:),'-*','Linewidth',2)
        semilogy([1:K],err_inf(4,:),'-d','Linewidth',2)
       
        legend("Parareal","Delta","theta","Delta-theta",'Fontsize',12)
        xlabel('iteration'); ylabel('error L^{\infty}');
        set(gca, 'FontSize', 15);
        %t_str=["FE T = "+T(ii)+", N ="+N(ii)];
        %title(t_str)
        
        for i=linspace(1,4,4)
            if isempty(find(err_inf(i,:)<=tol, 1))
                treshold_inf(i,ii)=K+1;
            else
                temp=find(err_inf(i,:)<=tol);
                treshold_inf(i,ii)=temp(1);
            end
        end
        
        for i=linspace(1,4,4)
            if isempty(find(err_2(i,:)<=tol, 1))
                treshold_L2(i,ii)=K+1;
            else
                temp=find(err_2(i,:)<=tol);
                treshold_L2(i,ii)=temp(1);
            end
        end
end

%% comparison N
N=[512*2,512*4,512*8,512*16];
T=[4*T_lambda];
treshold_inf=zeros(4,4);
treshold_L2=zeros(4,4);
for ii=1:4
        tol=1e-8;
        [~,~,err_inf(1,:),err_2(1,:)]=MGRIT_Lorentz_F(sigma,r,b,T,MF,N(ii),u0,K,"BE","N");

        [~,~,err_inf(2,:),err_2(2,:)]=MGRIT_Lorentz_F_Delta2(sigma,r,b,T,MF,N(ii),u0,K,"BE","N","N");

        [~,~,err_inf(3,:),err_2(3,:)]=MGRIT_Lorentz_F_theta(sigma,r ,b, T,MF,N(ii),u0,K,"BE","N");

        [~,~,err_inf(4,:),err_2(4,:),Y]=MGRIT_Lorentz_F_Delta_theta2(sigma,r,b,T,MF,N(ii),u0,K,"BE","N","N");
           
        figure
        xlabel('iteration','FontSize',20)
        ylabel('err_{L^{\infty}}','FontSize',20)
        

        semilogy([1:K],err_inf(1,:),'-*','Linewidth',2)
        hold on
        semilogy([1:K],err_inf(2,:),'-*','Linewidth',2)
        semilogy([1:K],err_inf(3,:),'-*','Linewidth',2)
        semilogy([1:K],err_inf(4,:),'-d','Linewidth',2)
       
        legend("Parareal","Delta","theta","Delta-theta",'Fontsize',12)
        %t_str=["FE T = "+T+", N ="+N(ii)];
        %title(t_str)
        xlabel('iteration'); ylabel('error L^{\infty}');
        set(gca, 'FontSize', 15);
        
        for i=linspace(1,4,4)
            if isempty(find(err_inf(i,:)<=tol, 1))
                treshold_inf(i,ii)=K+1;
            else
                temp=find(err_inf(i,:)<=tol);
                treshold_inf(i,ii)=temp(1);
            end
        end
        
        for i=linspace(1,4,4)
            if isempty(find(err_2(i,:)<=tol, 1))
                treshold_L2(i,ii)=K+1;
            else
                temp=find(err_2(i,:)<=tol);
                treshold_L2(i,ii)=temp(1);
            end
        end
end

