
%% initial the problem
% image a spring-mass system
% mx''+cx'+kx=F
% x''+2*zeta*omega0*x'+omega0^2*x=F/m
% X=[x;x']
% X'=AX+Bu
% Z=HX
% A=[0 1;-omega0^2 -2*zeta*omega0]
% B=[0;omega0^2]
% H=[1 0]
% u=F/m


% initialization
clc; clear; close all;
addpath(genpath("D:\Users\xushe\Documents\GitHub\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\Function_shengyi_package"))
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
zeta=0.1;
Ac=[0 1;-omega0^2 -2*zeta*omega0];
Bc=[0;omega0^2];
Hc=[1 0];
Q=0.1*eye(2);
R=2;
P_0_0=eye(2);
x0=[0;0];

dt = 0.01;
T = 25;
t = 0:dt:T;

[A, B, H, ~ ,~]=ssmod_c2d(Ac,Bc,Hc,[],dt);

% initial the input and output
N=length(t);

Ft=F*sin(2*pi*f*t);
u=Ft/m;
x=zeros(2,N);
z=zeros(1,N);
w = sqrt(Q)*randn(2,N);
v = sqrt(R)*randn(1,N);
x00= x0;
for k1=1:N
    x(:,k1)=A*x00+B*u(k1)+w(:,k1);
    z(k1)=H*x(:,k1)+v(k1);
    x00=x(:,k1);
end

[x_k_k,x_k_kmin,P_k_k,P_k_kmin]=KalmanFilterOneInput(A,B,H,Q,R,z,u,x0,P_0_0);



%% plot
[figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images);
hFigure = figure('Position', figPos_temp);
plot(t, x(1,:),  'Color', 'r','LineWidth', lineWidthThin);
hold on
plot(t, z(1,:),  'Color', 'g','LineWidth', lineWidthThin);
plot(t, x_k_k(1,:),  'Color', 'b','LineWidth', lineWidthThin);
xlabel('time (s)')
ylabel('displacement (m)')
legend('real','measure','filter')




function [figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images)
    figureIdx=figureIdx+1;
    row_idx = floor((figureIdx - 1) / num_figs_in_row);
    col_idx = mod((figureIdx - 1), num_figs_in_row);
    figPos_temp = figPos;
    figPos_temp(1:2) = figPos(1:2) + [col_idx * (figPos(3) + gap_between_images(1)) row_idx * (figPos(4) + gap_between_images(2))];
end