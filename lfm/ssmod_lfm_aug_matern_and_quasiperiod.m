%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-11-01 23:10:34
%LastEditors: xushengyichn xushengyichn@outlook.com
%LastEditTime: 2023-12-04 11:38:21
%FilePath: /ssm_tools_sy/lfm/ssmod_lfm_aug_matern_and_quasiperiod.m
%Description: 
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fac,Bac,Gac,Jac,Fad,Bad,Gad,Jad,Qad]=ssmod_lfm_aug_matern_and_quasiperiod(Ac,Bc,Gc,Jc,Fc,Hc,Lc,Ec,Tc,Kc,Qxd,sigma_w,sigma_z,dt)

%% Augmented model with modal states and latent states
%
% xa(t)=[ x(t) ; s(t) ; gamma(t)]
%
% dx/dt=Ac*x(t)+Bc*p(t) +Bc*eta(t)
% y=Gc*x(t)+Jc*p(t)+Jc*eta(t)
%
% ds/dt=Fc*s(t)+Lc*w(t)
% dgamma/dt=Ec*gamma(t)+Kc*z(t)
% p (t)=Lc*s(t)
% eta(t)=Tc*gamma(t)
%
% dxa/dt = [Ac Bc*Hc Bc*Tc; 0 Fc 0; 0 Ec 0]*xa(t)+[0; Lc*w(t); Kc*z(t)]
% y(t)=[ Gc Jc*Hc Jc*Tc]*xa(t)
%
% Inputs:
% Ac: state matrix (cont)
% Bc: input matrix (cont)
% Gc: output matrix (cont)
% Jc: direct matrix (cont)
% Fc: state matrix for LFM1 (cont)
% Lc: input matrix for LFM1 (cont)
% Hc: output matrix for LFM1 (cont)
% Ec: state matrix for LFM2 (cont)
% Kc: input matrix for LFM2 (cont)
% Tc: output matrix for LFM2 (cont)
% Qxd: additional covariance on state equation
% sigma_w: vector with standard deviations of LFM1 (cont)
% sigma_z: vector with standard deviations of LFM2 (cont)
% dt: time discretization
%
% Outputs:
% Fac: augmented state matrix (cont)
% Bac: 
% Gac: augmented output matrix (cont)
% Jac:
% Fad: augmented state matrix (disc)
% Bad: 
% Gad: augmented output matrix (disc)
% Jad:
% Qad: augmented covariance (disc)
%
%% 

% ns=size(Ac,2);
% ns2=size(Fc,1);

% nw=size(Lc,2);


if isvector(sigma_w)
    Sigma_w_squared=diag(sigma_w).^2;
elseif ismatrix(sigma_w)
    Sigma_w_squared=sigma_w.^2;
else
    Sigma_w_squared=sigma_w.^2;
end

if isvector(sigma_z)
    Sigma_z_squared=diag(sigma_z).^2;
elseif ismatrix(sigma_z)
    Sigma_z_squared=sigma_z.^2;
else
    Sigma_z_squared=sigma_z.^2;
end

Qwc=Lc*Sigma_w_squared*Lc.';
Qzc=Kc*Sigma_z_squared*Kc.';
Qxc=zeros(size(Ac,2));

Qac=blkdiag(Qxc,Qwc,Qzc);

Fac=[Ac Bc*Hc Bc*Tc ; zeros(size(Fc,1),size(Ac,2)),Fc,zeros(size(Fc,1),size(Tc,2)); zeros(size(Ec,1),size(Ac,2)),zeros(size(Ec,1),size(Fc,2)),Ec];

Gac=[Gc Jc*Hc Jc*Tc];

Bac=[];

Jac=[];

Qd=cov_c2d(Fac,Qac,dt);

if isempty(Qxd)
    Qxd=zeros(size(Qd));
end


Qad=Qd+blkdiag(Qxd,zeros(size(Qwc,1)+size(Qzc,1),size(Qwc,2)+size(Qzc,2)));

Qad=(Qad+Qad.')/2;

[Fad,Bad,Gad,Jad]=ssmod_c2d(Fac,Bac,Gac,Jac,dt);

