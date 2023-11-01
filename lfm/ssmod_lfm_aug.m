%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-11-01 23:10:34
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-11-01 23:14:32
%FilePath: \Exercises-for-Techniques-for-estimation-in-dynamics-systemsf:\git\ssm_tools\lfm\ssmod_lfm_aug.m
%Description: 
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fac,Bac,Hac,Jac,Fad,Bad,Had,Jad,Qad]=ssmod_lfm_aug(Ac,Bc,Gc,Jc,Fc,Hc,Lc,Qxd,sigma_w,dt)

%% Augmented model with modal states and latent states
%
% xa(t)=[ x(t) ; s(t) ]
%
% dx/dt=Ac*x(t)+Bc*p(t)
% y=Gc*x(t)+Jc*p(t)
%
% ds/dt=Fc*s(t)+Lc*w(t)
% p(t)=Lc*s(t)
%
% dxa/dt=[ Ac Bc*Hc ; 0 Fc ]*xa(t)+[ 0 ;  Lc*w(t)]
% y(t)=[ Gc Jc*Hc ]*xa(t)
%
% Inputs:
% Ac: state matrix (cont)
% Bc: input matrix (cont)
% Gc: output matrix (cont)
% Jc: direct matrix (cont)
% Fc: state matrix for LFM (cont)
% Hc: input matrix for LFM (cont)
% Lc: output matrix for LFM (cont)
% Qxd: additional covariance on state equation
% sigma_w: vector with standard deviations (cont)
% dt: time discretization
%
% Outputs:
% Fac: augmented state matrix (cont)
% Bac: 
% Hac: augmented output matrix (cont)
% Jac:
% Fad: augmented state matrix (disc)
% Bad: 
% Had: augmented output matrix (disc)
% Jad:
% Qad: augmented covariance (disc)
%
%% 

ns=size(Ac,2);
ns2=size(Fc,1);

nw=size(Lc,2);


if isvector(sigma_w)
    Sigma_w_squared=diag(sigma_w).^2;
elseif ismatrix(sigma_w)
    Sigma_w_squared=sigma_w.^2;
else
    Sigma_w_squared=sigma_w.^2;
end



Qwc=Lc*Sigma_w_squared*Lc.';
Qxc=zeros(ns);

Qac=blkdiag(Qxc,Qwc);

Fac=[Ac Bc*Hc ; zeros(ns2,ns) Fc];

Hac=[Gc Jc*Hc];

Bac=[Bc ; zeros(size(Fc,1),size(Bc,2))];

Jac=Jc;

Qd=cov_c2d(Fac,Qac,dt);

if isempty(Qxd)
    Qxd=zeros(size(Qd));
end


if isempty(Qxd)
    Qxd=zeros(size(Qd));
end


Qad=Qd+blkdiag(Qxd,zeros(size(Qwc)));

Qad=(Qad+Qad.')/2;

[Fad,Bad,Had,Jad]=ssmod_c2d(Fac,Bac,Hac,Jac,dt);

