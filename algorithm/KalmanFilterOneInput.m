function [x_k_k,x_k_kmin,P_k_k,P_k_kmin]=KalmanFilterOneInput(A,B,H,Q,R,z,u,x0,P_0_0)

%% Kalman filter
%
% Model:
% x(k+1)=A*x(k)+B*u(k)+w(k);
% z=H*x(k)+v(k);
%
% Inputs:
% A: state matrix
% B: input matrix
% H: output matrix
% Q: state noise covariance
% R: output noise covariance
% S: mixed noise covariance
% z: output vector
% u: input vector
% x0: initial state estimate
% P_0_0: initial state error covariance
%
% Outputs:
% x_k_k: filter state estimate
% x_k_kmin: prediction state estimate
% P_k_k: filter error covariance
% P_k_kmin: prediction error covariance

N = size(z,2);
x=x0;
for k1 =1:N
    %% Prediction
    x_=A*x+B*u(k1);
    P_=A*P_0_0*A.'+Q;

    x_k_kmin(:,k1)=x_;
    P_k_kmin(:,:,k1)=P_;
    %% Correction
    Kk=P_*H.'/(H*P_*H.'+R);
    x=x_+Kk*(z(k1)-H*x_);

    x_k_k(:,k1)=x;
    %% Update
    P_k=forcesym((eye(size(A))-Kk*H)*P_);
    P_k_k(:,:,k1)=P_k;

end
end
%% Other functions

function B=forcesym(A)

B=(A+A.')/2;

end
