%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-07-15 11:37:02
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-07-17 18:47:46
%FilePath: \ssm_tools\examples\AKF_example.m
%Description: add the Accelerometer information compared with AKF_example
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
k = m*omega0^2;
c=2*zeta*omega0*m;
Ac=[0 1;-omega0^2 -2*zeta*omega0];
Bc=[0;1];
Hc=[1,0;0,1;-k/m,-c/m];
Jc=[0;0;1/m];
Q=0.001*eye(2);
R=0.02;
P_0_0=eye(2);
x0=[0;0];

dt = 0.01;
T = 100;
t = 0:dt:T;

[A, B, H, J ,~]=ssmod_c2d(Ac,Bc,Hc,Jc,dt);
% J = zeros(2,1);
S = ones(1,1)*0.5;
% initial the input and output
N=length(t);

Ft=F*sin(2*pi*f*t);
u=Ft;
x=zeros(2,N);
z=zeros(3,N);
w = sqrt(Q)*randn(2,N);
v = sqrt(R)*randn(3,N);
x00= x0;
for k1=1:N
    x(:,k1)=A*x00+B*u(k1)+w(:,k1);
    z(:,k1)=H*x(:,k1)+J*u(k1)+v(:,k1);
    x00=x(:,k1);
end

y=z;
p=u;
G=H;
np=size(Ft,1);
Pp_0_0=eye(np)*1;
p0 = zeros(np,1);
% [x_k_k,x_k_kmin,P_k_k,P_k_kmin]=KalmanFilterWithInput_shengyi(A,B,G,J,Q,R,y,p,x0,P_0_0)
[x_k_k,x_k_kmin,P_k_k,P_k_kmin]=AKF(A,B,G,J,Q,R,S,y,x0,p0,P_0_0,Pp_0_0);
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
plot(t, x_k_k(end,:),  'Color', 'b','LineWidth', lineWidthThin);

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