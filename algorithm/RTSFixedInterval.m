%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-07-20 16:55:21
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-07-20 17:11:01
%FilePath: \ssm_tools\algorithm\RTSFixedInterval.m
%Description:
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_k_N, P_k_N] = RTSFixedInterval(A, x_k_k, x_k_kmin, P_k_k, P_k_kmin)

    %% RTS smoother
    %
    % Inputs:
    % A: state matrix
    % x_k_k: filter state estimate from Kalman filter
    % x_k_kmin: prediction state estimate from Kalman filter
    % P_k_k: filter error covariance from Kalman filter
    % P_k_kmin: prediction error covariance from Kalman filter
    %
    % Outputs:
    % x_k_N: smoothed state estimate
    % P_k_N: smoothed state error covariance
    %
    % Note:
    % For S~=0, the system must be transformed so that the process and maesurement noise is uncorrelated. A_star=A-S/R*G;
    % See Niu (2011) and Verify_KF_RTS.m

    N = size(x_k_k, 2);
    x_k_N = zeros(size(x_k_k));
    P_k_N = zeros(size(P_k_k));

    for k = N:-1:1

        if k == N
            x_k_N(:, k) = x_k_k(:, k);
            P_k_N(:, :, k) = P_k_k(:, :, k);
        else
            xkplus1_star = x_k_N(:, k + 1);
            Pkplus1_star = P_k_N(:, :, k + 1);
            xk = x_k_k(:, k);
            xkplus1_ = x_k_kmin(:, k + 1);
            Pkplus1_ = P_k_kmin(:, :, k + 1);
            Pk = P_k_k(:, :, k);
            Lk = Pk * A' / Pkplus1_;
            x_k_N(:, k) = xk + Lk * (xkplus1_star - xkplus1_);
            P_k_N(:, :, k) = Pk + Lk * (Pkplus1_star - Pkplus1_) * Lk';
        end

    end
