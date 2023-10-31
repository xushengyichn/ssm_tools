function Qd=cov_c2d(Fac,Qc,dt,method)

%% Convert covariance from cont to disc 
%
% Model:
% x(k+1)=F*x(k)+w(k);
% y=H*x(k)+v(k);
%
% Inputs:
% Fac: state matrix (cont)
% Qc: covariance matrix (cont)
% dt: time discretization (cont)
% method: 'matrix' or 'integral', should give same result but the former is faster
%
% Outputs:
% Qd: covariance matrix (disc)
%
%%

% Default matrix method
if nargin==3
    method='matrix';
end

ns=size(Fac,1);

%% Matrix method

% From Särka, Recursive Bayesian inference on stochastic processes, eq 2.89

if strcmpi(method,'matrix')

	Matrix_1= expm( [Fac Qc ; zeros(ns) -Fac.']*dt);
	Matrix_2=[ zeros(ns) ; eye(ns) ];

	Matrix_3=Matrix_1*Matrix_2;

	C=Matrix_3(1:ns,:);
	D=Matrix_3(ns+[1:ns],:);

	Qd=C/D;
	
end

%% Old way, numerical integration

if strcmpi(method,'integral')

	Psi_a_func = @(tau) expmnorm(Fac*tau);

	n=300;
	tau_axis=linspace(0,dt,n);

	integrand=zeros(ns,ns,n);
	for k=1:length(tau_axis)
		Psi_a_temp=Psi_a_func(dt-tau_axis(k));
		integrand(:,:,k)=Psi_a_temp*Qc*Psi_a_temp.';
	end

	Qd=trapz(tau_axis,integrand,3);

end

%%

ratio=norm(Qd-Qd.')./norm(Qd);

if ratio>1e-6
	ratio
	warning('Loss of symmetry in covariance, take a look here')
end

Qd=forcesym(Qd);

%% Tests


% ratio=norm(Qd-Qd_old)./norm(Qd_old)
% ratio
% 
% [v,d]=eig(Qd);
% [v_old,d_old]=eig(Qd_old);
% 
% figure();  hold on;
% plot(diag(d_old)./diag(d)-1);


%%

% Fac=magic(50); Fac=Fac+Fac.';
% Fac=diag(-[0.2:0.2:2]); Fac=Fac-(magic(size(Fac))+magic(size(Fac)).')/1000;
% 
% Qc=eye(size(Fac)); a=randn(size(Qc)); a=(a+a.')/10; Qc=Qc+a;
% dt=0.05
% 
% plotcorr(Qc)
