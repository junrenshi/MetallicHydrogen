% The script dispatches GbarLMGPU_worker, and iterates until converged.
% It determinies the average of the Green's function.

parameters;

%% kssolv is used to contruct the k-mesh.
addpath(kssolvpath);

%% Reading in PIMD parameters
fid = fopen('XDATCAR');		% Read in PIMD result

% Read the configuration of the simulation supercell
data = textscan(fid, '%f %f %f %f %f', 1, 'headerLines', 1);
lats = [data{2}  data{3}  data{4}];
fclose(fid);

% Invoke kssolve to construct the k-mesh
C = diag(lats)/aB;
mol = Molecule();
mol = set(mol,'supercell',C);
mol = set(mol,'ecut', ecut);  % kinetic energy cut off

%% Obtain k-mesh
gmask = FreqMask(mol);
ng = get(gmask, 'ng');	% Number of k-points

% Set the output filename
fn_Self = sprintf('./Gbar_Self_%02d.mat', komega);

err = 1;
jbatch = {};
while err > 1e-3
	if ~exist(fn_Self, 'file')  % For a new set of data
		Sigma = zeros(ng);
		save(fn_Self, 'Sigma');		
	end
	
	prefix = 'GbarLM';
	execfn = 'GbarLMGPU_worker';

	% Launch workers to GPUs
	for jobid = 1:NGPU
		workspace = struct('jobid', jobid, 'AdditionalPaths', '../../');
		jbatch{jobid} = batch(execfn, 'Workspace', workspace);
	end

	for jobid = 1:NGPU
		wait(jbatch{jobid});
	end

	%% Collect the results of Gbar*_workers
	Tbars = zeros(ng);
	Ntot = 0;

	for jobid = 1:NGPU
		fn = sprintf('%s_%02d_%02d.mat', prefix, komega, jobid);
		load(fn);
		Ntot = Ntot + NS;
		Tbars = Tbars + Tbar;
	end

	%% Determine the average Green's function, T matrix, and Self energy
	Tbars = Tbars / Ntot;			% Average T-matrix
	Gbar = G0 + G0 * Tbars * G0;	% Average Green's functioon 
	DSigma = inv(G0) - inv(Gbar);	% Change of the self-energy
	Sigma = Sigma + DSigma;			% The total self-energy

	err = norm(DSigma)

	save(fn_Self, 'Sigma');
end
