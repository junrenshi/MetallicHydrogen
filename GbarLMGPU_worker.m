%% Calculating the T-matrices of ionic configurations 
%% from a PIMD simulation of metallic Hydrogen. 
%% One needs to provide a jobid and kssolvpath
%% The script will be launched by runGbar.m.

parameters;

setenv('KSSOLVPATH', kssolvpath);
addpath(kssolvpath);

%% Matsubara frequency for which the average Green's function to be evaluated.
%% We adopt the quasi-static approximation, which is valid only for 
%% large Matsubara frequencies. 
omegan = (2*komega+1) *pi / beta1;
eta = omegan;

%% Reading in PIMD simulation results.  Parameters and headers.
fid = fopen('./XDATCAR');
data = textscan(fid, '%f %f %f %f %f', 1, 'headerLines', 1);
lats = [data{2}  data{3}  data{4}];
data = textscan(fid, '%f %f %f', 0, 'headerLines', 4);

C = diag(lats)/aB;
CI = inv(C);
Vol = det(C);

%% Setup k-mesh
mol = Molecule();
mol = set(mol,'supercell',C);
mol = set(mol,'ecut', ecut);  % kinetic energy cut off
gmask = FreqMask(mol);
gkx = get(gmask, 'gkx');
gky = get(gmask, 'gky');
gkz = get(gmask, 'gkz');
gv = [gkx gky gkz];
ng = get(gmask, 'ng');
g = sqrt(gkx.^2 + gky.^2 + gkz.^2);

%% Generate a list of unique q-vectors
qij = zeros(ng, ng);
qvij = zeros(ng, ng, 3);

for j = 1:ng
	for i = (j+1):ng
		qvij(i, j, :) = gv(i,:) - gv(j, :);
		qij(i, j) = norm(gv(i, :) - gv(j, :));
		qij(j, i) = qij(i, j);
		qvij(j, i, :) = -qvij(i, j, :);
	end
end

[qc, qia, qic] = unique(reshape(qij, ng^2, 1));
[qvc, qvia, qvic] = unique(reshape(qvij, ng^2, 3), 'rows'); 

% Setup self-energy, Fermi energy and eigenstates
fn_Self = sprintf('./Gbar_Self_%02d.mat', komega);
load(fn_Self);

idF = g > 0.5*kF & g < 1.5*kF;  % Energy shell cutoff a +/-0.5kF
gF = g(idF);
nF = length(gF);
if exist('DeF2') == 0
	% Hbar = diag(g.^2/2) + (Sigma + Sigma')/2;
	% ebar = eig(Hbar);
	% DeF2 = (ebar(Ne2-1) + ebar(Ne2) + ebar(Ne2+1))/3 - eF2;

	% The shift of the chemical potential is estimated from
	% the diagonal of the self-energy.
	sigma = diag(Sigma);
	DeF2 = median(real(sigma(g > 0.95*kF & g < 1.05*kF)));
end

%% Free Green's function
G0 = (diag(eF2 + DeF2 - (g.^2)/2  + 1i*eta) - Sigma) \ eye(ng);


% Coulomb potential between electrons and ions
vc = -4*pi./q2epsUI_et(qc/kF, rs) / kF^2;
vc(1) = 0;  % q=0 component of the Coulomb potential only shifts the chemical potential

Tbar = zeros(ng);				% Average T-matrix, full 
TFbar = zeros(ng, nF, NBEAD);	% Average T-matrix, partial for WVs close to FS
TF2bar = zeros(ng, nF, NBEAD);	% T-T correlation: <TT>

% Setup GPU
gpuDevice(jobid);
NGPU = gpuDeviceCount;

G0_gpu = gpuArray(single(G0));
Sigma_gpu = gpuArray(single(Sigma));

qic_gpu = gpuArray(qic);
qvic_gpu = gpuArray(qvic);
qvc_gpu = gpuArray(single(qvc));

vc_gpu = gpuArray(single(vc));
vcarray_gpu = vc_gpu(qic_gpu);
eye_gpu = eye(ng, 'single', 'gpuArray');

Tbar_gpu = gpuArray(single(Tbar));
TFbar_gpu = gpuArray(single(TFbar));
TF2bar_gpu = gpuArray(single(TF2bar));

idF_gpu = gpuArray(idF);
TF_gpu = zeros(ng, nF, NBEAD, 'single', 'gpuArray');

%% Reading in PIMD trajectories, and calculating matrix elements
skip = jobid - 1 + NSKIP;
xyzmats = [];
while 1
	%% Reading in a ion configuration
	data = textscan(fid, '%f %f %f', Natoms*NBEAD, 'headerLines', (Natoms+1)*NBEAD*skip);
	% skip the trailing blank line
	textscan(fid, '%f %f %f', 0, 'headerLines', 1);

	skip = NGPU-1;

	if feof(fid)
		break;
	end

	%% Setting up ion positions
	coefs = [data{1}  data{2}  data{3}];
	xyzmat = coefs*C';
	xyzmats = [xyzmats; xyzmat];
end

xyzmats_gpu = gpuArray(single(xyzmats));

NS = 0
tic;
for i = 1:length(xyzmats)/Natoms/NBEAD
	for m = 1:NBEAD
		%% Start sampling the matrix elements
		NS = NS + 1;

		%% Pick a sample
		xyzmatm_gpu = xyzmats_gpu((1:Natoms) + (m-1 + (i-1)*NBEAD)*Natoms, :);

		%% Determine ion potential -- Linear screening approximation
		rhoc_gpu = sum(exp(-1i * qvc_gpu * xyzmatm_gpu'), 2)/Vol;
		V_gpu = reshape(vcarray_gpu .* rhoc_gpu(qvic_gpu), ng, ng);
		
		%% Subtract self energy from V -- The potential in the effective medium
		V_gpu = V_gpu - Sigma_gpu;

		%% Get the T-matrix by matrix inverse
		T_gpu = (eye_gpu - V_gpu * G0_gpu) \ V_gpu;

		%% Pick the elements of the T-matrix for the cutoff energy shell 
		TF_gpu(:,:,m) = T_gpu(:, idF);

		%% Accumulate Tbar
		Tbar_gpu = Tbar_gpu + T_gpu;
	end

	%% Transform the T-matrices from the time-domain to the frequency domain
	TFw_gpu = ifft(TF_gpu, NBEAD, 3);

	%% Accumulate the TFbar, in the frequency domain
	TFbar_gpu = TFbar_gpu + TFw_gpu; 

	%% Accumulate |T|^2, which is related to the scattering amplitude of EPC
	TF2bar_gpu = TF2bar_gpu + abs(TFw_gpu).^2;
end

%% Collect results from GPU
Tbar = gather(Tbar_gpu);
TFbar = gather(TFbar_gpu);
TF2bar = gather(TF2bar_gpu);

%% File to save results
fn = sprintf('GbarLM_%02d_%02d.mat', komega, jobid);
save(fn, 'Tbar', 'TFbar', 'TF2bar', 'ng', ...
		'gv', 'g', 'G0', 'C', 'DeF2', 'NS', 'idF', ...
		'Sigma', 'eta');

fclose(fid);
gpuDevice([]);
clear