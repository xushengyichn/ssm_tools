%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-07-15 11:37:02
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-07-17 19:09:57
%FilePath: \ssm_tools\examples\AKF_example2.m
%Description: case in Øyvind's slide
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

%% 0 number of modes to consider
nmodes=4;
ns = nmodes*2;
%% 1 former exercise code to calculate the matrix
m1=1.5; m2=1; m3=1; m4=1;
k1=25; k2=20; k3=15; k4=15;
M=[m1 0 0 0; 0 m2 0 0; 0 0 m3 0; 0 0 0 m4];
K=[k1+k2 -k2 0 0; -k2 k2+k3 -k3 0; 0 -k3 k3+k4 -k4; 0 0 -k4 k4];
C= 0.1*M+0.005*K;

[V,D]=eig(K,M);
wn=sqrt(diag(D));
f = wn/2/pi;
zeta = diag(V'*C*V)./2./wn;



w1=wn(1); w2=wn(2); w3=wn(3); w4=wn(4);
xi1=zeta(1); xi2=zeta(2); xi3=zeta(3); xi4=zeta(4);

data = table(wn,zeta);
data.Properties.VariableNames = {'Undamped_Natural_Frequency', 'Damping_Ratio'};
% disp(data);

%% 2 calculate the modal shape matrix
phi=V;
phi1=phi(:,1); phi2=phi(:,2); phi3=phi(:,3); phi4=phi(:,4);

%% 3 establish continuous time matrices
T = phi'*C*phi;
omega2=phi'*K*phi;
T = T(1:nmodes,1:nmodes);
omega2 = omega2(1:nmodes,1:nmodes);
phi = phi(:,1:nmodes);

% S_p
S_p=[1 0; 0 1; 0 0; 0 0];
% S_a
S_a=[1 0 0 0; 0 1 0 0 ; 0 0 1 0; 0 0 0 0; 0 0 0 0];
% S_v
S_v= [0 0 0 0; 0 0 0 0 ; 0 0 0 0; 0 0 0 0; 0 0 0 0];
% S_d
S_d= [0 0 0 0; 0 0 0 0 ; 0 0 0 0; 1 0 0 0; 0 0 1 0];


A_c = [zeros(nmodes,nmodes) eye(nmodes,nmodes); -omega2 -T];
B_c = [zeros(nmodes,2); phi'*S_p];
G_c = [S_d*phi-S_a*phi*omega2 S_v*phi-S_a*phi*T];
J_c = [S_a*phi*phi'*S_p];

%% 4 establish discrete time matrices
Time=20;
dt = 0.05;
t = 0:dt:Time;

A_d = expm(A_c*dt);
B_d = A_c\(A_d-eye(ns))*B_c;
G_d = G_c;
J_d = J_c;


%% 5 calculate the response of the system
A1=10; A2=20; A3=10; A4=20;
f1=0.2; f2=0.3; f3=0.7; f4=1.1;
phi1=0; phi2=pi; phi3=0; phi4=pi;
p1=A1*sin(2*pi*f1*t-phi1)+A2*sin(2*pi*f2*t-phi2);
p2=A3*sin(2*pi*f3*t-phi3)+A4*sin(2*pi*f4*t-phi4);
p = [p1;p2];


A = A_d;
B = B_d;
C = G_d;
D = J_d;

H = zeros(5,ns);
G = eye(ns);


Q = 10^(-4)*eye(ns);
R = 10^(-1)*eye(5);
S = zeros(8,5);


% %% 测试用代码
% Q = 10^(-8)*eye(ns);
% R = 10^(-8)*eye(5);
% %%

N = length(t);
w = sqrt(Q)*randn(ns,N);
v = sqrt(R)*randn(5,N);
C = cov(w.');
xn = zeros(ns,length(t));
xn_true = zeros(ns,length(t));

for t1 = 1:length(t)-1
    xn(:,t1+1) = A_d*xn(:,t1)+B_d*p(:,t1)+w(:,t1);
    yn(:,t1) = G_d*xn(:,t1)+J_d*p(:,t1)+v(:,t1);
    xn_true(:,t1+1) = A_d*xn_true(:,t1)+B_d*p(:,t1);
end
yn(:,end+1) = G_d*xn(:,end)+J_d*p(:,end)+v(:,end);

x = [phi zeros(4,nmodes); zeros(4,nmodes) phi]*xn;
x_true = [phi zeros(4,nmodes); zeros(4,nmodes) phi]*xn_true;
clear x;

% figure(1);
% plot(t,x_true(3,:),'r');
% hold on
% plot(t(1:end),yn(5,:),'b','Linestyle','--');
% 
% legend('u1_{true}','u1_{measured}');
% xlabel('time(s)');
% ylabel('displacement(m)');
% title('displacement response of the system');



%% 6 Augmented Kalman Filter
% initialize the state estimate

% x = zeros(ns,1);
% P = 10^(-4)*eye(ns);



% collect=[];
% collect2=[];
% collect_cov=[];
% NN=N;
% for k1 = 1:NN
%     Omega = G_d*P*G_d'+R;
%     % estimate
%     x = x + P*G_d'*inv(Omega)*(yn(:,k1)-G_d*x-J_d*p(:,k1));
%     P = P - P*G_d'*inv(Omega)*G_d*P;
%     collect(1:8,k1)=x;
%     collect_cov(1:8,k1)=diag(P);
%     % predict
%     K = (A_d*P*G_d'+S)*inv(Omega);
%     x = A_d*x + B_d*p(:,k1) + K*(yn(:,k1)-G_d*x-J_d*p(:,k1));
%     P = A_d*P*A_d'+(A_d*P*G_d'+S)*inv(Omega)*(A_d*P*G_d'+S)'+Q;
%     collect2(1:8,k1)=x;
% %     scatter(k1,x(1))

% end


% x_est = [phi zeros(4,nmodes); zeros(4,nmodes) phi]*collect;
% x_pred = [phi zeros(4,nmodes); zeros(4,nmodes) phi]*collect2;
% x_pred = [zeros(ns,1) x_pred(:,1:end-1)];

% x1_upper_bound = x_est(1,:)+3*sqrt(collect_cov(1,:));
% x1_lower_bound = x_est(1,:)-3*sqrt(collect_cov(1,:));

% t_fill = [t(1:NN) fliplr(t(1:NN))];
% x1_fill = [x1_upper_bound fliplr(x1_lower_bound)];

% Augmented Kalman Filter
% initialize the state estimate

np = 2;
x_ak = zeros(ns+np,1);
P_ak = 10^(-4)*eye(ns+np);

A_a = [A_d B_d; zeros(np,ns) eye(np,np)];
G_a= [G_d  J_d];

SS = 10^(1)*eye(np);
Q_a = [Q zeros(ns,np); zeros(np,ns) SS];

yn_a = [yn]; %外力是不作观测值的，所以不需要加入到yn_a中

NN = N;
xa_history = zeros(ns+np,NN);
pa_history = zeros(ns+np,NN);
for k1 =1:NN
    % measurement update
    Lk=P_ak*G_a'*inv(G_a*P_ak*G_a'+R);
%     yn_a(end-1:end,k1) = x_ak(end-1:end);
    dk = yn_a(:,k1);
    % dk = G_a*yn_a(:,k1)+v(:,k1);
    x_ak = x_ak + Lk*(dk - G_a*x_ak);
    P_ak = P_ak - Lk*G_a*P_ak;

    xa_history(:,k1) = x_ak;
    pa_history(:,k1) = diag(P_ak);
    % time update
    x_ak = A_a*x_ak;
    P_ak = A_a*P_ak*A_a'+Q_a;

    
end

x0 = zeros(ns,1);
p0 = zeros(np,1);
P_0_0 = 10^(-4)*eye(ns);
Pp_0_0 = 10^(-4)*eye(np);
% [x_k_k,x_k_kmin,P_k_k,P_k_kmin]=AKF(A_d,B_d,G_d,J_d,Q,R,SS,yn,x0,p0,P_0_0,Pp_0_0);
[x_k_k,x_k_kmin,P_k_k,P_k_kmin,p_k_k,Pp_k_k]=JIS(A_d,B_d,G_d,J_d,Q,R,yn,x0,P_0_0);
xa_history=x_k_k;

x_history = [phi zeros(4,nmodes); zeros(4,nmodes) phi]*xa_history(1:ns,:);
Px_history = abs([phi zeros(4,nmodes); zeros(4,nmodes) phi])*pa_history(1:ns,:);
p_history = p_k_k;
Pp_history = Pp_k_k;


for k1 =1:4

[figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images);
plot(t,x_true(k1,:),'r');
hold on
plot(t,x_history(k1,:),'b','Linestyle','--');
xlabel('time (s)')
ylabel('displacement (m)')
legend(["u"+num2str(k1)+"_{true}"],"u"+num2str(k1)+"_{predit}")
end


for k1 =1:4

[figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images);
plot(t,x_true(k1+4,:),'r');
hold on
plot(t,x_history(k1+4,:),'b','Linestyle','--');
xlabel('time (s)')
ylabel('velocity (m/s)')
legend(["u"+num2str(k1)+"_{true}"],"u"+num2str(k1)+"_{predit}")
end

for k1 = 1:2

[figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images);
plot(t,p(1,:),'r');
hold on
plot(t,p_history(1,:),'b','Linestyle','--');
xlabel('time(s)');
ylabel('Force(N)');
legend(["p"+num2str(k1)+"_{true}"],"p"+num2str(k1)+"_{predit}")
end

function [figureIdx,figPos_temp] = create_figure(figureIdx, num_figs_in_row,figPos,gap_between_images)
    figureIdx=figureIdx+1;
    row_idx = floor((figureIdx - 1) / num_figs_in_row);
    col_idx = mod((figureIdx - 1), num_figs_in_row);
    figPos_temp = figPos;
    figPos_temp(1:2) = figPos(1:2) + [col_idx * (figPos(3) + gap_between_images(1)) row_idx * (figPos(4) + gap_between_images(2))];
    hFigure = figure('Position', figPos_temp);
end
