%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-07-15 12:44:21
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-07-17 18:37:33
%FilePath: \ssm_tools\algorithm\AKF.m
%Description: based on Lourens, E., et al. "An augmented Kalman filter for force identification in structural dynamics." Mechanical systems and signal processing 27 (2012): 446-460.
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_k_k,x_k_kmin,P_k_k,P_k_kmin]=AKF(A,B,G,J,Q,R,S,y,x0,p0,P_0_0,Pp_0_0)

%% Kalman filter
%
% Model:
% x(k+1)=A*x(k)+B*p(k)+w(k);
% y=G*x(k)+J*p(k)+v(k);
%
% Inputs:
% A: state matrix
% B: input matrix
% G: output matrix
% Q: state noise covariance
% R: output noise covariance
% S: covariance of the force
% y: output vector
% p: input vector
% x0: initial state estimate
% P_0_0: initial state error covariance
% Pp_0_0: initial state error covariance
%
% Outputs:
% x_k_k: filter state estimate
% x_k_kmin: prediction state estimate
% P_k_k: filter error covariance
% P_k_kmin: prediction error covariance

N = size(y,2);
ns=size(A,2);
np=size(B,2);

A_a = [A,B;zeros(np,ns),eye(np)];
G_a = [G,J];
Q_a = [Q,zeros(ns,np);zeros(np,ns),S];

% x=x0;
x=[x0;p0];
P_0_0 = [P_0_0,zeros(ns,np);zeros(np,ns),Pp_0_0];
P_k = P_0_0;
for k1 =1:N
    %% Prediction
    x_=A_a*x;
    P_=A_a*P_k*A_a.'+Q_a;

    x_k_kmin(:,k1)=x_;
    P_k_kmin(:,:,k1)=P_;
    %% Correction
    Kk=P_*G_a.'/(G_a*P_*G_a.'+R);
    x=x_+Kk*(y(:,k1)-G_a*x_);

    x_k_k(:,k1)=x;
    %% Update
    P_k=forcesym((eye(size(A_a))-Kk*G_a)*P_);
    P_k_k(:,:,k1)=P_k;

end

% P_ =P_0_0;
% x_ = x;
% for k1 =1:N
%     %% measurement update
%     Kk=P_*G_a.'/(G_a*P_*G_a.'+R);
%     x=x_+Kk*(y(:,k1)-G_a*x_);
%     P_k=forcesym((eye(size(A_a))-Kk*G_a)*P_);
%     x_k_k(:,k1)=x;    
%     P_k_k(:,:,k1)=P_k;
% 
%     % time update
%     x_=A_a*x;
%     P_=A_a*P_k*A_a.'+Q_a;
% 
%     x_k_kmin(:,k1)=x_;
%     P_k_kmin(:,:,k1)=P_;
% 
% 
% end


end
%% Other functions

function B=forcesym(A)

B=(A+A.')/2;

end
