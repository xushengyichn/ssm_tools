%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-07-15 11:37:02
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-07-15 11:38:11
%FilePath: \ssm_tools\examples\KalmanFilterNoInput_example.m
%Description: track the state of a spring-mass system with Kalman filter
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initial the problem
% image a spring-mass system
% mx''+cx'+kx=0
% x''+2*zeta*omega0*x'+omega0^2*x=0
% X=[x;x']
% X'=AX+Bu
% Z=HX
% A=[0 1;-omega0^2 -2*zeta*omega0]
% B=[0;omega0^2]
% H=[1 0]



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
f=1;
omega0=2*pi*f;
m = 1;
F = 2;
zeta=0.01;
Ac=[0 1;-omega0^2 -2*zeta*omega0];
Bc=[0;omega0^2];
Hc=[1 0];
Q=0.001*eye(2);
R=2;
P_0_0=eye(2);
x0=[50;0];

dt = 0.01;
T = 100;
t = 0:dt:T;

[A, B, H, ~ ,~]=ssmod_c2d(Ac,Bc,Hc,[],dt);

% initial the input and output
N=length(t);


x=zeros(2,N);
z=zeros(1,N);
w = sqrt(Q)*randn(2,N);
v = sqrt(R)*randn(1,N);
x00= x0;
for k1=1:N
    x(:,k1)=A*x00+w(:,k1);
    z(k1)=H*x(:,k1)+v(k1);
    x00=x(:,k1);
end

[x_k_k,x_k_kmin,P_k_k,P_k_kmin]=KalmanFilterNoInput(A,H,Q,R,z,x0,P_0_0,'steadystate',true);
% [x_k_k,x_k_kmin,P_k_k,P_k_kmin]=KalmanFilterNoInput(A,H,Q,R,z,x0,P_0_0,'steadystate',false);
% S = zeros(2,1);
% y=z;
% [x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss]=KalmanFilter(A,H,Q,R,S,y,x0,P_0_0);

%% plot
[figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images);


plot(t, z(1,:),  'Color', 'g','LineWidth', lineWidthNormal);
hold on
plot(t, x(1,:),  'Color', 'r','LineWidth', lineWidthNormal);
plot(t, x_k_k(1,:),  'Color', 'b','LineWidth', lineWidthNormal,'LineStyle','-.');
xlabel('time (s)')
ylabel('displacement (m)')
legend('measure','real','filter')



function [figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images)
    figureIdx=figureIdx+1;
    row_idx = floor((figureIdx - 1) / num_figs_in_row);
    col_idx = mod((figureIdx - 1), num_figs_in_row);
    figPos_temp = figPos;
    figPos_temp(1:2) = figPos(1:2) + [col_idx * (figPos(3) + gap_between_images(1)) row_idx * (figPos(4) + gap_between_images(2))];
    hFigure = figure('Position', figPos_temp);
end