%% Determine the density correlation functions from PIMD simulation
%% utilizing GPUs

%% Parameters
parameters;
setenv('KSSOLVPATH', kssolvpath);
addpath(kssolvpath);


%% Reading in PIMD simulation results.  Parameters and headers.
fid = fopen('XDATCAR');
data = textscan(fid, '%f %f %f %f %f', 1, 'headerLines', 1);
lats = [data{2}  data{3}  data{4}];
data = textscan(fid, '%f %f %f', 0, 'headerLines', 4);

C = diag(lats)/aB;
CI = inv(C);
Omega = det(C);

%% Load the set of q-vectors for which the PDF is evaluated
load('qvFc.mat');

gv1 = qvFc;
ng1 = length(gv1);
g1 = sqrt(gv1(:,1).^2 + gv1(:,2).^2 + gv1(:,3).^2)/kF;

%% Reading in PIMD trajectories, and calculating matrix elements
NS = 0;
rho2self = zeros(ng1, NBEAD);
rho2 = zeros(ng1, NBEAD);
rhobar = zeros(ng1, NBEAD);

gv1_gpu = gpuArray(single(gv1));
rho2self_gpu = gpuArray(single(rho2self));
rho2_gpu = gpuArray(single(rho2));
rhobar_gpu = gpuArray(single(rhobar));

%% Skipping 
if NSKIP > 0
	textscan(fid, '%f %f %f', 0, 'headerLines', (Natoms+1)*NBEAD*NSKIP);
end
tic;
NS = 0;
while 1
	%% Reading in an ion configuration
	data = textscan(fid, '%f %f %f', Natoms*NBEAD);

	if feof(fid)
		break;
	end

	%% Setting up ion positions
	coefs = [data{1}  data{2}  data{3}];
	xyzmat = coefs*C';

	if size(xyzmat,1) < Natoms*NBEAD
		break;
	end 

	NS = NS + 1;
	xyzmat_gpu = gpuArray(single(xyzmat));

	%% Determine rho(q, t)
	rhoqit_gpu = reshape(exp(-1i* gv1_gpu * xyzmat_gpu'), ng1, Natoms, NBEAD);
	%% Transform from the time domain to the frequency domain
	rhoqiw_gpu = ifft(rhoqit_gpu, NBEAD, 3);
	
	%% Accumulate the self-correlation function
	rho2self_gpu = rho2self_gpu + squeeze(sum(abs(rhoqiw_gpu).^2, 2));

	%% Accumulate for hq
	rhoqw_gpu = squeeze(sum(rhoqiw_gpu, 2));
	rho2_gpu = rho2_gpu + abs(rhoqw_gpu).^2;
	rhobar_gpu = rhobar_gpu + rhoqw_gpu;
end

rho2 = gather(rho2_gpu);
rhobar = gather(rhobar_gpu);
rho2self = gather(rho2self_gpu);

rho2 = rho2 / NS;
rhobar = rhobar / NS;
rho2self = rho2self / NS;

hq = (rho2  -  abs(rhobar).^2) - rho2self;
hq = hq / Natoms;

omegaq = rho2self / Natoms;

save('TCPA2_qvFc_PDF.mat', 'hq', 'omegaq');

gpuDevice([]);