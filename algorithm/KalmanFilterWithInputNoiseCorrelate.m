%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-07-15 12:44:21
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-07-17 13:53:39
%FilePath: \ssm_tools\algorithm\KalmanFilterWithInputNoiseCorrelate.m
%Description: 
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_k_k,x_k_kmin,P_k_k,P_k_kmin]=KalmanFilterWithInputNoiseCorrelate(A,B,G,J,Q,R,S,y,p,x0,P_0_0)

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
% S: mixed noise covariance (E(wv')=S)
% y: output vector
% p: input vector
% x0: initial state estimate
% P_0_0: initial state error covariance
%
% Outputs:
% x_k_k: filter state estimate
% x_k_kmin: prediction state estimate
% P_k_k: filter error covariance
% P_k_kmin: prediction error covariance

N = size(y,2);
x=x0;
for k1 =1:N
    %% Prediction
    x_=A*x+B*p(:,k1)+S*inv(R)*(y(:,k1)-J*p(:,k1)-G*x);
    P_=A*P_0_0*A.'+Q-S*inv(R)*S.';

    x_k_kmin(:,k1)=x_;
    P_k_kmin(:,:,k1)=P_;
    %% Correction
    Kk=P_*G.'/(G*P_*G.'+R);
    x=x_+Kk*(y(:,k1)-J*p(:,k1)-G*x_);

    x_k_k(:,k1)=x;
    %% Update
    P_k=forcesym((eye(size(A))-Kk*G)*P_);
    P_k_k(:,:,k1)=P_k;

end
end
%% Other functions

function B=forcesym(A)

B=(A+A.')/2;

end
