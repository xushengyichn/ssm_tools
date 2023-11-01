%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-07-04 12:48:06
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-07-04 15:30:11
%FilePath: \ssm_tools\lfm\ssmod_quasiperiod.m
%Description:have not been verified
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fc, Lc, Hc, sigma_w] = ssmod_quasiperiod(lambda, sigma_p, omega_0)

    %% State space model for Matern
    %
    % ds/dt = F*s(t) + L*w(t)
    % p(t) = H*s(t)
    %
    % Inputs:
    % lambda: hyperparameter, inverse length scale
    % omega_0: the frequency of the periodic component
    % sigma_p: standard deviation of p
    %
    % Outputs:
    % Fc: state matrix (cont)
    % Lc: input matrix (cont)
    % Hc: output matrix (cont)
    % sigma_w: standard deviation (cont)
    %
    %
    %% Assign state space matrices

    Lc = eye(2, 2);
    Hc = [1, 0];
    Fc = [-lambda, -omega_0; omega_0, -lambda];

    %% Variance of p:

    sigma_w = sqrt(sigma_p .^ 2 * (2 * lambda));
