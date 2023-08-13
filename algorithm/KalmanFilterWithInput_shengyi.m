%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-07-15 12:44:21
%LastEditors: xushengyichn xushengyichn@outlook.com
%LastEditTime: 2023-08-13 10:52:35
%FilePath: \ssm_tools\algorithm\KalmanFilterWithInput_shengyi.m
%Description: 
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_k_k,x_k_kmin,P_k_k,P_k_kmin]=KalmanFilterWithInput_shengyi(A,B,G,J,Q,R,y,p_det,x0,P_0_0, varargin)

%% Kalman filter
%
% Model:
% x(k+1)=A*x(k)+B*u(k)+w(k);
% y=G*x(k)+J*p(k)+v(k);
%
% Inputs:
% A: state matrix
% B: input matrix
% G: output matrix
% Q: state noise covariance
% R: output noise covariance
% S: mixed noise covariance
% y: output vector
% p_det: known (deterministic) input
% x0: initial state estimate
% P_0_0: initial state error covariance
%
% Outputs:
% x_k_k: filter state estimate
% x_k_kmin: prediction state estimate
% P_k_k: filter error covariance
% P_k_kmin: prediction error covariance

%% Parse inputs

p = inputParser;
p.KeepUnmatched = true;

addParameter(p, 'steadystate', true, @islogical)
addParameter(p, 'showtext', true, @islogical)
addParameter(p, 'noscaling', false, @islogical)

parse(p, varargin{1:end});

steadystate = p.Results.steadystate;
showtext = p.Results.showtext;
noscaling = p.Results.noscaling;


%% Initialization
N = size(y,2);
p = p_det;
x=x0;
P_k = P_0_0;

%% Conventional
if steadystate == false
    t0 = tic;

    for k1 =1:N
        %% Prediction
        x_=A*x+B*p(:,k1);
        P_=A*P_k*A.'+Q;

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
    telapsed = toc(t0);
end




    %% Steady state
    if steadystate == true

        if noscaling == true
            [P_k_kmin_ss, ~, ~, info] = idare(A.', G.', Q, R, [], [], 'noscaling');
        else
            [P_k_kmin_ss, ~, ~, info] = idare(A.', G.', Q, R);
        end

        if info.Report ~= 0

            if info.Report == 1
                disp('***** DARE solution accuracy poor, running with scaling');
                warning('***** DARE solution accuracy poor, running with scaling');
            end

            if info.Report == 2
                warning('***** DARE solution not finite, running with scaling');
            end

            if info.Report == 3
                warning('***** DARE solution not found, running with scaling');
            end

            [P_k_kmin_ss, ~, ~, info] = idare(A.', G.', Q, R);

        end

        if info.Report == 1
            disp('***** DARE solution accuracy poor');
            warning('***** DARE solution accuracy poor');
        elseif info.Report == 2
            error('DARE solution not finite');
        elseif info.Report == 3
            error('DARE solution not found');
        end

        P_k_kmin_ss = forcesym(P_k_kmin_ss);
        Kk_ss = P_k_kmin_ss * G.' / (G * P_k_kmin_ss * G.' + R);
        t0 = tic;


        for k1 =1:N
            %% Prediction
            x_=A*x+B*p(:,k1);
    
            x_k_kmin(:,k1)=x_;
            %% Correction
            x=x_+Kk_ss*(y(:,k1)-J*p(:,k1)-G*x_);
    
            x_k_k(:,k1)=x;
            %% Update
    
        end
        [m, n] = size(P_k_kmin_ss);  % 获取P_k_kmin_ss的尺寸
        P_k_kmin=repmat(P_k_kmin_ss, [1, 1, N]);
        P_k_k=P_k_kmin;
        telapsed = toc(t0);
    end
    if showtext == true
        disp(['Kalman filter calculated in ' sprintf('%2.1f', telapsed) ' seconds, ' sprintf('%2.1f', telapsed * 10 ^ 5 ./ N) ' seconds per 1M steps']);
    end


end
%% Other functions

function B=forcesym(A)

B=(A+A.')/2;

end
