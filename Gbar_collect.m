% The script collects and analyzes the outputs of GbarLMGPU_worker.m.

parameters;
fres = sprintf('Gbar_collect_%02d.mat', komega);

%% Collect the results of Gbar*_workers
clear Tbars TFbars TF2bars VFbars VTbars VTFbars 
Tbars = 0;
TFbars = 0;
TF2bars = 0;
Ntot = 0;

jids = 1:NGPU;
for jid = jids
	fn = sprintf('./GbarLM_%02d_%02d.mat', komega, jid);
	load(fn);
	Ntot = Ntot + NS;				% Total number of ion configurations
	Tbars = Tbars + Tbar;			
	TFbars = TFbars + TFbar;
	TF2bars = TF2bars + TF2bar;
end

Vol = det(C);

%% Determine the averages of quantities
Tbars = Tbars / Ntot;
TFbars = TFbars / Ntot * NBEAD;

% TF2bars_full: the fluctuation of the T-matrix, in the larger k shell
% with 0.5kF < k < 1.5kF.
TF2bars_full = TF2bars / Ntot * NBEAD - abs(TFbars).^2;   

%% Set of momenta in the larger k-shell
idF_full = idF;
gF_full = g(idF_full);
nF_full = length(gF_full);

%% Index to the smaller k-shell, relative to the whole set of momenta
idF = g > 0.9*kF & g < 1.1*kF;

% Index to the smaller k-shell, relative to the larger k-shell
idF1 = gF_full > 0.9*kF & gF_full < 1.1*kF;

% TF2bars: the fluctuation of the T-matrix, in the smaller k shell
% with 0.9kF < k < 1.1kF.
TF2bars = TF2bars_full(:, idF1, :);

% Determine the self-energy.  We assume the self-energy is diagonal.
% This is valid only for liquids.
Gbar = G0 + G0 * Tbars * G0;
DSigma = inv(G0) - inv(Gbar);
Sigma = Sigma + DSigma;
gbar = diag(Gbar);
g0 = 1./ (eF2 + DeF2 - g.^2/2 + 1i * eta);
sigma = 1./g0 - 1./gbar;

save(fres, 'g', 'gbar', 'sigma', 'TF2bars', 'Vol');

% Average the self energy for k-points with same magnitudes.
[gc, gia, gic] = uniquetol(g/kF, 1e-6); 
ngc = length(gc);
sigma_PIMD = [];
for i = 1:ngc
	sigma_PIMD(i) = median(sigma(gic == i)/eF2);
end

save(fres, 'gc', 'sigma_PIMD', '-append');

% Effective interaction
gvF = gv(idF, :); % Wave vectors on the Fermi surface
nF = length(gvF);

%% Determine effective pairing interaction by solving the Bethe-Salpeter
%% equation. We use the direct matrix inverse in the larger truncated k-space.
TF2barstau = real(fft(TF2bars_full(idF_full, :, :), NBEAD, 3));
V2tau_full = [];
for m = 1:NBEAD
	V2tau_full(:,:,m) = (eye(nF_full) + TF2barstau(:,:,m) ...
				* diag(abs(gbar(idF_full)).^2) ) \ TF2barstau(:,:,m);
end
V2_full = ifft(V2tau_full, NBEAD, 3);

%% Scattering amplitudes and effective interaction in the smaller k-shell
%% They are converted from the atomic units to our units.
TF2 = reshape(TF2bars(idF,:,:), nF^2, NBEAD) * 9 * pi^4 * Natoms / eF2^2;
V2 = reshape(V2_full(idF1, idF1, :), nF^2, NBEAD) * 9 * pi^4 * Natoms / eF2^2;
save(fres, 'TF2', 'V2', '-append');

%% Convert the resulting effective interaction in the matrix form
%% to a function of q  -- We assume that the liquid phase is 
%% isotropic and uniform.
gxm = repmat(gv(:,1), 1, ng);
gym = repmat(gv(:,2), 1, ng);
gzm = repmat(gv(:,3), 1, ng);
qxm = gxm - gxm';
qym = gym - gym';
qzm = gzm - gzm';

qm = sqrt(qxm.^2 + qym.^2 + qzm.^2);

qF = reshape(qm(idF, idF)/kF, 1, nF^2);
save(fres, 'qF', '-append');

qvF = [reshape(qxm(idF, idF), nF^2, 1), reshape(qym(idF, idF), nF^2, 1), ...
		reshape(qzm(idF, idF), nF^2, 1)];

%% Unique set of q-vectors
[qvFc, qvFia, qvFic] = uniquetol(qvF, 1e-5, 'ByRows', true);

%% Determine a unique set of q, and determine the pair distribution function
save('./qvFc.mat', 'qvFc');
TCPA2GPU_qvFc_PDF;
load('./TCPA2_qvFc_PDF.mat');

clear chiq
for m = 1:NBEAD
	chiq(:, m) = double(omegaq(qvFic, m) + hq(qvFic, m));
end

%% Average over q-vectors with the same |q|
[qFc, qFia, qFic] = uniquetol(qF, 0.003); % Unique set of qF

clear TF2c V2c DTF2c DV2c M2c DM2c
for m = 1:NBEAD
	TF2c(:,m)   = accumarray(qFic, TF2(:, m),       [], @median);
	V2c(:,m)    = accumarray(qFic, V2(:, m),        [], @median);
	M2c(:, m)   = accumarray(qFic, V2(:, m)./chiq(:,m), [], @median);
	DTF2c(:,m)  = accumarray(qFic, TF2(:, m),       [], @std);
	DV2c(:,m)   = accumarray(qFic, V2(:, m),        [], @std);
	DM2c(:, m)  = accumarray(qFic, V2(:, m)./chiq(:,m), [], @std);
end

%% bare lambda values
lambdaV = [];
for m = 1:NBEAD
	xV2cf = @(x) interp1(qFc/2, V2c(:,m), x, 'pchip') .* x;
	lambdaV(m) = integral(xV2cf, 0, 1) * beta / (6*pi^4);
end

%% bare lambda values + 1 sigma
lambdaVa = [];
for m = 1:NBEAD
	xV2cf = @(x) interp1(qFc/2, V2c(:,m) + DV2c(:,m), x, 'pchip') .* x;
	lambdaVa(m) = integral(xV2cf, 0, 1) * beta / (6*pi^4);
end

%% bare lambda values - 1 sigma
lambdaVb = [];
for m = 1:NBEAD
	xV2cf = @(x) interp1(qFc/2, V2c(:,m) - DV2c(:,m), x, 'pchip') .* x;
	lambdaVb(m) = integral(xV2cf, 0, 1) * beta / (6*pi^4);
end

%% Estimate the standard deviation of lambda
DlambdaV = (lambdaVa - lambdaVb)/2;
save(fres, 'qFc', 'TF2c', 'V2c', 'DV2c', 'lambdaV', 'DlambdaV', '-append');

%% Determine q-dependence of the correlation functions.
omegac = zeros(length(qFc), NBEAD);
hqc = zeros(length(qFc), NBEAD);
dchic = zeros(length(qFc), NBEAD);
dhqc = zeros(length(qFc), NBEAD);
for m = 1:NBEAD
	omegac(:, m) = accumarray(qFic, omegaq(qvFic,m), [], @median);
	hqc(:, m)    = accumarray(qFic, hq(qvFic,m),     [], @median);
	dchic(:, m)  = accumarray(qFic, hq(qvFic,m) + omegaq(qvFic,m), [], @std);
	dhqc(:, m)  = accumarray(qFic, hq(qvFic,m), [], @std);
end
chic = double(omegac + hqc);

%% Model of the interaction matrix element
Mf = @(xv) -8*pi*alpha*rs ./ q2epsUI_et(2*xv, rs);
M2chi = diag(abs(Mf(qFc/2)).^2) * chic;

%% Fitting to the effective pairing interaction
scalsV = ones(5, NBEAD); % The scaling coefficient for the VT matrix elements
nscals = size(scalsV, 1);
qmin = 0.2;
% xscal = linspace(qmin/2, 1, nscals);
xscal = [qmin, 0.5, 0.75, 0.875, 1];
id1 = qF < 2 & qF > qmin;
lb = 0.1*ones(5, 1);
ub = 3  *ones(5, 1);

for m = 1:NBEAD
	fitfun = @(c0, xdata) interp1(xscal, c0, xdata, 'pchip', 'extrap').^2 ...
			 .* Mf(xdata).^2;
	scalsV(:,m) = lsqcurvefit(fitfun, scalsV(:,m), qF(id1)'/2, ...
			V2(id1, m)./chiq(id1,m), lb, ub);
end

scalVq = interp1(xscal, scalsV, qFc/2, 'pchip', 'extrap');
scalM2chi = scalVq.^2 .* M2chi;

%% Lamda values from the fitting model
lambdaVM = [];
for m = 1:NBEAD
	xVMf = @(xv) interp1(qFc/2, scalM2chi(:,m), xv, 'pchip') .* xv;
	lambdaVM(m) = integral(xVMf, 0, 1) * beta / (6*pi^4);
end

save(fres, 'xscal', 'scalsV', 'scalM2chi', 'lambdaVM', '-append');

% Interpolate omegas in tau dependence, so that it gives rise to a correct 
% large q behavior.
NBEAD1 = NBEAD * 16;
[chitau1, omegatau1, hqtau1] = extrapolatePDF(NBEAD, NBEAD1, omegac, hqc);
omega1q = real(ifft(omegatau1, NBEAD1, 2));
omega1f = @(xv) interp1(qFc/2, omega1q, xv);

NG = ceil((NBEAD+1)/2);

hqc1 = [hqc(:, 1:NG), zeros(length(hqc), NBEAD1-NBEAD), hqc(:, (NG+1):NBEAD)];
chi1f = @(xv) omega1f(xv) + interp1(qFc/2, hqc1, xv, 'pchip');
M2chi1 = diag(abs(Mf(qFc/2)).^2) * chi1f(qFc/2);

scalsVq1 = [scalVq(:,1:NG), repmat(scalVq(:, NG), 1, NBEAD1-NBEAD), ...
		scalVq(:, (NG+1):NBEAD)];

scalM2chi1 = abs(scalsVq1).^2 .* M2chi1;
lambdaVM1 = [];
for m = 1:NBEAD
	xVM1f = @(xv) interp1(qFc/2, scalM2chi1(:,m), xv, 'pchip') .* xv;
	lambdaVM1(m) = integral(xVM1f, 0, 1) * beta / (6*pi^4);
end

DlambdaVM1 = DlambdaV;
DlambdaVM1((NG+1):NBEAD) = lambdaVM1((NG+1):NBEAD) * DlambdaV(NG)/ lambdaV(NG); 

save(fres, 'NBEAD1', 'lambdaVM1', 'scalM2chi1', '-append');

% Solve the linearized Eliahsberg equation
Ncut = floor(NBEAD/2);
[rho, omega2] = Gbar_Eliashberg(lambdaVM1, Ncut, mustar);
[rhoa, omega2a] = Gbar_Eliashberg(lambdaVM1+DlambdaVM1, Ncut, mustar);
[rhob, omega2b] = Gbar_Eliashberg(lambdaVM1-DlambdaVM1, Ncut, mustar);
drho = (rhoa-rhob)/2;
domega2 = (omega2a-omega2b)/2;
save(fres, 'rho', 'omega2', 'drho', 'domega2', '-append');