%% script to run Lorentz parareal analysis with Forward euler with truncated input
clear all; close all; clc
parent_folder = fileparts(pwd);
addpath(parent_folder);
%% data 

sigma=10;r=28;b=8/3;

T_lambda=log(10)/0.9;
T=4*T_lambda;
MF=20;
MG=1;
N=512;
u0=[20;5;-5];
K=50;

%% problematic case
Lorenz_time_test(sigma,r,b,T,MF,MG,N,u0,K)
%% bounded

B=50; % L2 con fine e coarse
MF=20;
[U,u,err,err2]=Lorentz_Parareal_bound(sigma, r , b, T,MF,MG,N,u0,K,B,'Y','N');
%%
B=100; % L2 con fine e coarse
MF=20;
[U,u,err,err2]=Lorentz_Parareal_bound(sigma, r , b, T,MF,MG,N,u0,K,B,'Y','Y');
%%
B=100; % L2 con fine e coarse
MF=20;
[U,u,err,err2]=Lorentz_Parareal_bound(sigma, r , b, T,MF,MG,N,u0,K,B,'Y','N');
%%
B=1000; % L2 con fine e coarse
MF=20;
[U,u,err,err2]=Lorentz_Parareal_bound(sigma, r , b, T,MF,MG,N,u0,K,B,'Y','N');
%%
B=30; % L2 con fine e coarse
MF=20;
[U,u,err,err2]=Lorentz_Parareal_bound(sigma, r , b, T,MF,MG,N,u0,K,B,'Y','N');
%%
B=50;
MF=20;
[U,u,err,err2]=Lorentz_Parareal_bound(sigma, r , b, T,MF,MG,N,u0,K,B,'Y','Y');
%%
B=70;
MF=20;
[U,u,err,err2]=Lorentz_Parareal_bound(sigma, r , b, T,MF,MG,N,u0,K,B,'Y','Y');