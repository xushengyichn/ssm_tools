%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-07-15 12:44:21
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-07-18 15:44:50
%FilePath: \ssm_tools\algorithm\JIS.m
%Description: based on The derivative process can be found in
% *Gillijns, Steven, and Bart De Moor. "Unbiased minimum-variance input and state estimation for linear discrete-time systems with direct feedthrough." Automatica 43.5 (2007): 934-937.*
% The case study can be found in
% *Lourens, E., et al. "Joint input-response estimation for structural systems based on reduced-order models and vibration data from a limited number of sensors." Mechanical Systems and Signal Processing 29 (2012): 310-327.*
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_k_k,x_k_kmin,P_k_k,P_k_kmin,p_k_k,Pp_k_k]=JIS(A,B,G,J,Q,R,y,x0,P_0_0)

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
% y: output vector
% x0: initial state estimate
% P_0_0: initial state error covariance
%
% Outputs:
% x_k_k: filter state estimate
% x_k_kmin: prediction state estimate
% P_k_k: filter error covariance
% P_k_kmin: prediction error covariance

N = size(y,2);
ns=size(A,2);
np=size(B,2);


% initial state
x_=x0;
Px_ = P_0_0 ;
x_k_kmin(:,1)=x_;
P_k_kmin(:,:,1)=Px_;
for k1 =1:N
    %% input estimation:
    R_bar = G*Px_*G.'+R;
    M=inv(J'*inv(R_bar)*J)*J'*inv(R_bar);
    p=M*(y(:,k1)-G*x_);
    Pp=inv(J'*inv(R_bar)*J);
    p_k_k(:,k1)=p;
    Pp_k_k(:,:,k1)=Pp;
    %% measurement update:
    L = Px_*G'*inv(R_bar);
    x = x_ + L*(y(:,k1)-G*x_-J*p);
    Px = forcesym(Px_ - L*(R_bar - J*Pp*J')*L');
    Pxp = -L*J*Pp;
    x_k_k(:,k1)=x;
    P_k_k(:,:,k1)=Px;
    %% time update:
    x_ = A*x + B*p;
    Px_ = [A,B]*[Px,Pxp;Pxp.',Pp]*[A,B]'+Q;
    x_k_kmin(:,k1+1)=x_;
    P_k_kmin(:,:,k1+1)=Px_;
end
% for k1 =1:N
%     %% Prediction
%     x_=A_a*x;
%     P_=A_a*P_k*A_a.'+Q_a;

%     x_k_kmin(:,k1)=x_;
%     P_k_kmin(:,:,k1)=P_;
%     %% Correction
%     Kk=P_*G_a.'/(G_a*P_*G_a.'+R);
%     x=x_+Kk*(y(:,k1)-G_a*x_);

%     x_k_k(:,k1)=x;
%     %% Update
%     P_k=forcesym((eye(size(A_a))-Kk*G_a)*P_);
%     P_k_k(:,:,k1)=P_k;

% end

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
