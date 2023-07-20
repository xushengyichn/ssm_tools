function [Hpy Hx0y Hx1y]=JIS_tf(A,B,G,J,Q,R,S,dt,omega_axis)
%% Transfer function for steady state operation of joint input and state estimator
%
% NOT VERIFIED YET
%
% Model
% x(k+1)=A*x(k)+B*p(k)+w(k);
% y=G*x(k)+J*p(k)+v(k);
%
% Inputs:
% A: state matrix
% B: input matrix
% G: output matrix
% J: direct transmission matrix
% Q: state noise covariance
% R: output noise covariance
% S: mixed noise covariance
% dt: time step
% omega_axis: frequency axis
%
% Outputs:
% Hpy: matrix with TF, output-to-force estimate
% Hx0y: matrix with TF, output-to-filter estimate
% Hx1y: matrix with TF, output-to-prediction estimate
%
%%

ns=size(A,1);
ny=size(G,1);
np=size(B,2);

y_dummy=nan(ny,1);
x0=[];
P0=[];

[~,~,~,~,M_ss,L_ss] = JIS_ss(A,B,G,J,y_dummy,x0,Q,R,S,P0,'trunc',false,'convtol',1e-8);

%%

M2=[M_ss ; L_ss ; zeros(ns,ny) ];

for k=1:length(omega_axis)

    M1=[ eye(np) zeros(np,ns) M_ss*G ;
        L_ss*J eye(ns) (-eye(ns)+L_ss*G);
        B A -exp(1i.*omega_axis(k)*dt)*eye(ns) ];

    M3=eye(size(M1)) / M1 * M2;

    M4(:,:,k)=M3; %M4 is [Hpd;Hx0d;Hx1d]

end

Hpy=M4(1:np,:,:);
Hx0y=M4(np+[1:ns],:,:);
Hx1y=M4(np+ns+[1:ns],:,:);

