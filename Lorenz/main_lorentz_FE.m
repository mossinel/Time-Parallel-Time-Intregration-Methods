%% script to run Lorentz parareal analysis with Forward euler
clear all; close all; clc
parent_folder = fileparts(pwd);
addpath(parent_folder);
%% data and test for varying T 

sigma=10;r=28;b=8/3;
T_lambda=log(10)/0.9;

T=[1*T_lambda,2*T_lambda,3*T_lambda,4*T_lambda];

MF=40;
MG=1;
N=512;
u0=[20;5;-5];
K=50;
for t=T
   Lorenz_time_test(sigma,r,b,t,MF,10,N,u0,K) 
end

