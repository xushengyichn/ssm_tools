
%% initial the problem
% track a car 
% v2=v1 + w1
% s2 = s1 +v1*t +w2

% X=[s;v]
% X'=AX
% Z=HX
% A=[0 dt;0 1]
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
num_figs_in_row = 4; %每一行显示几个图
figPos = figPosSmall; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
gap_between_images = [0,0];
figureIdx = 0;

% initial the parameters


Q=0.01*eye(2);
R=10;
P_0_0=eye(2);


dt = 0.01;
T = 5;
t = 0:dt:T;
A=[1 dt;0 1];
H = [1 0];

x0=[10;10];
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
% S = zeros(2,1);
% y=z;
% [x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss]=KalmanFilter(A,H,Q,R,S,y,x0,P_0_0);
% [x_k_N,P_k_N]=RTSFixedInterval(A,x_k_k,x_k_kmin,P_k_k,P_k_kmin);
[x_k_N,P_k_N,P_klag_N]=RTSSmoother(A,x_k_k,x_k_kmin,P_k_k(:,:,1),P_k_kmin(:,:,1),'steadystate',true);
%% plot,
[figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images);

plot(t, x(1,:),  'Color', 'r','LineWidth', lineWidthThin);
hold on
plot(t, z(1,:),  'Color', 'g','LineWidth', lineWidthThin);
plot(t, x_k_k(1,:),  'Color', 'b','LineWidth', lineWidthThin);
plot(t, x_k_N(1,:),  'Color', 'cyan','LineWidth', lineWidthNormal,'LineStyle','--');
xlabel('time (s)')
ylabel('displacement (m)')
legend('real','measure','filter','smoothing')
title("估计值")

[figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images);

plot(t, x(1,:),  'Color', 'r','LineWidth', lineWidthThin);
title("真实值")
xlabel('time (s)')
ylabel('displacement (m)')
legend('real')

[figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images);
plot(t, z(1,:),  'Color', 'g','LineWidth', lineWidthThin);
title("测量值")
xlabel('time (s)')
ylabel('displacement (m)')
legend('measure')

[figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images);
plot(t,squeeze(P_k_k(1,1,:)), 'Color', 'b','LineWidth', lineWidthNormal)
hold on
plot(t,squeeze(P_k_N(1,1,:)), 'Color', 'cyan','LineWidth', lineWidthNormal,'LineStyle','-')
legend('filter','smoothing')
xlabel('time (s)')
ylabel('covariance')


function [figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images)
    figureIdx=figureIdx+1;
    row_idx = floor((figureIdx - 1) / num_figs_in_row);
    col_idx = mod((figureIdx - 1), num_figs_in_row);
    figPos_temp = figPos;
    figPos_temp(1:2) = figPos(1:2) + [col_idx * (figPos(3) + gap_between_images(1)) row_idx * (figPos(4) + gap_between_images(2))];
    hFigure = figure('Position', figPos_temp);
end