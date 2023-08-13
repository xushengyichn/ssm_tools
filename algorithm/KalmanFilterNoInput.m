%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn xushengyichn@outlook.com
%Date: 2023-07-19 23:19:38
%LastEditors: xushengyichn xushengyichn@outlook.com
%LastEditTime: 2023-08-12 12:07:21
%FilePath: \ssm_tools\algorithm\KalmanFilterNoInput.m
%Description:
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_k_k, x_k_kmin, P_k_k, P_k_kmin] = KalmanFilterNoInput(A, H, Q, R, z, x0, P_0_0, varargin)

    %% Kalman filter
    %
    % Model:
    % x(k+1)=A*x(k)+w(k);
    % z=H*x(k)+v(k);
    %
    % Inputs:
    % A: state matrix
    % H: output matrix
    % Q: state noise covariance
    % R: output noise covariance
    % S: mixed noise covariance
    % z: output vector
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

    N = size(z, 2);
    x = x0;
    P_k = P_0_0;

    %% Conventional
    if steadystate == false
        t0 = tic;

        for k1 = 1:N
            %% Prediction
            x_ = A * x;
            P_ = A * P_k * A.' + Q;

            x_k_kmin(:, k1) = x_;
            P_k_kmin(:, :, k1) = P_;
            %% Correction
            Kk = P_ * H.' / (H * P_ * H.' + R);
            x = x_ + Kk * (z(:, k1) - H * x_);

            x_k_k(:, k1) = x;
            %% Update
            P_k = forcesym((eye(size(A)) - Kk * H) * P_);
            P_k_k(:, :, k1) = P_k;

        end

        telapsed = toc(t0);

    end



    %% Steady state
    if steadystate == true

        if noscaling == true
            [P_k_kmin_ss, ~, ~, info] = idare(A.', H.', Q, R, [], [], 'noscaling');
        else
            [P_k_kmin_ss, ~, ~, info] = idare(A.', H.', Q, R);
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

            [P_k_kmin_ss, ~, ~, info] = idare(A.', H.', Q, R);

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
        Kk_ss = P_k_kmin_ss * H.' / (H * P_k_kmin_ss * H.' + R);
        t0 = tic;

        for k1 = 1:N
            %% Prediction
            x_ = A * x;

            x_k_kmin(:, k1) = x_;
            % P_k_kmin(:, :, k1) = P_k_kmin_ss;
            %% Correction

            x = x_ + Kk_ss * (z(:, k1) - H * x_);

            x_k_k(:, k1) = x;
            %% Update
            % P_k = forcesym((eye(size(A)) - Kk * H) * P_);
            % P_k_k(:, :, k1) = P_k_kmin_ss;

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

function B = forcesym(A)

    B = (A + A.') / 2;

end
