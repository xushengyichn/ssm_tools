%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-07-15 11:37:02
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-07-19 16:42:17
%FilePath: \ssm_tools\examples\AKF_LFM_example.m
%Description: 
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initial the problem
% image a spring-mass system
% mx''+cx'+kx=F
% x''+2*zeta*omega0*x'+omega0^2*x=F/m
% X=[x;x']
% X'=AX+Bu
% Z=HX
% A=[0 1;-omega0^2 -2*zeta*omega0]
% B=[0;1]
% H=[1 0]
% u=F/m


% initialization
clc; clear; close all;
addpath(genpath("D:\Users\xushe\Documents\GitHub\ssm_tools"))
addpath(genpath("F:\git\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\Function_shengyi_package"))
addpath(genpath("F:\git\Function_shengyi_package"))
subStreamNumberDefault = 2132;
run('InitScript.m');


%% 0 绘图参数
fig_bool = ON;
num_figs_in_row = 3; %每一行显示几个图
figPos = figPosSmall; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
gap_between_images = [0,0];
figureIdx = 0;

% initial the parameters
f=0.1;
omega0=2*pi*f;
m = 1;
F = 10;
zeta=0.1;
Ac=[0 1;-omega0^2 -2*zeta*omega0];
Bc=[0;1];
Gc=[1,0;0,1];
Jc = [0;0];
% k = m*omega0^2;
% c=2*zeta*omega0*m;
% Gc=[1,0;1,0];
% Jc = [0;0];
% Gc=[1,0;0,1;-k/m,-c/m];
% Jc=[0;0;1/m];
Q=0.001*eye(2);
Qxd=Q;


x0=[0;0];

dt = 0.01;
T = 100;
t = 0:dt:T;

%% LFM model
% Choose the matern kernel
lambda =0.001;
sigma_p = 2;
% [Fc,Lc,Hc,sigma_w]=ssmod_matern(lambda,p_order,sigma_p);
[Fc, Lc, Hc, sigma_w] =ssmod_quasiperiod(lambda, sigma_p, omega0);
[Fac,Bac,Hac,Jac,Fad,Bad,Had,Jad,Qad]=ssmod_lfm_aug(Ac,Bc,Gc,Jc,Fc,Hc,Lc,Qxd,sigma_w,dt);

nx=size(Ac,1);
ns=size(Fc,1);
nmeasure=size(Jc,1);
Rad=zeros(nmeasure,nmeasure);
R = zeros(nmeasure,nmeasure);
RR=0.01;
for k1 =1:nmeasure
    Rad(k1,k1)=RR;
    R(k1,k1)=RR;
end


% S = ones(1,1)*0.01;
% initial the input and output
N=length(t);
[A, B, G, J ,~]=ssmod_c2d(Ac,Bc,Gc,Jc,dt);
Ft=F*sin(2*pi*f*t);
u=Ft;
x=zeros(2,N);
z=zeros(2,N);
w = sqrt(Q)*randn(2,N);
v = sqrt(R)*randn(2,N);
x00= x0;
for k1=1:N
    x(:,k1)=A*x00+B*u(k1)+w(:,k1);
    z(:,k1)=G*x(:,k1)+v(:,k1);
    x00=x(:,k1);
end



P_0_0=eye(nx+ns);
xa0=[x0;zeros(ns,1)];
[x_k_k,x_k_kmin,P_k_k,P_k_kmin]=KalmanFilterNoInput(Fad,Had,Qad,Rad,z,xa0,P_0_0);

Jp=[1];
Gp=[zeros(1,size(Gc,1)),Jp*Hc];
p_k_k=Gp*x_k_k;
%% plot
[figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images);
plot(t, x(1,:),  'Color', 'r','LineWidth', lineWidthThin);
hold on
plot(t, z(1,:),  'Color', 'g','LineWidth', lineWidthThin);
plot(t, x_k_k(1,:),  'Color', 'b','LineWidth', lineWidthThin);
xlabel('time (s)')
ylabel('displacement (m)')
legend('real','measure','filter')

[figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images);
plot(t, Ft,  'Color', 'r','LineWidth', lineWidthThin);
hold on
plot(t,p_k_k,  'Color', 'b','LineWidth', lineWidthThin);

xlabel('time (s)')
ylabel('Force (N)')
legend('real','estimate')



function [figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images)
    figureIdx=figureIdx+1;
    row_idx = floor((figureIdx - 1) / num_figs_in_row);
    col_idx = mod((figureIdx - 1), num_figs_in_row);
    figPos_temp = figPos;
    figPos_temp(1:2) = figPos(1:2) + [col_idx * (figPos(3) + gap_between_images(1)) row_idx * (figPos(4) + gap_between_images(2))];
    hFigure = figure('Position', figPos_temp);
end