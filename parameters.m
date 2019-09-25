%% Parameters
rs = 1.17;  		% Density parameter
Temp = 150; 		% Temperature in K

%% PIMD parameters
NBEAD = 24;			% Number of beads discretizing the imaginary time
Natoms = 200;   	% Number of atoms
Ne2 = Natoms/2;		% Number of electrons per spin species.
NSKIP = 1000;   	% Number of samples to be skipped

%% Local running environment
NGPU = 8;			% Number of GPUs 
%% kssolv is used to contruct the k-mesh.
kssolvpath = '/home/shi/Documents/MATLAB/Add-Ons/kssolv';	% The path of kssolv

%% Parameters to evaluate Tc
komega = 16;		% The frequency index for the quasi-static approximation
mustar = 0.089;		% Coulomb pseudopotential parameter in determining Tc
ecut = 30;			% Energy cutoff in Ry to set up the k-meash


%% Constants and derived parameters
alpha = (4/(9*pi))^(1/3);
kF = 1/(alpha*rs); 	% Fermi wave number in 1/aB
eF = kF^2;			% Fermi energy in Ry
eF2 = kF^2/2;		% in Hartree = 2*Ry
aB = 5.2917721067e-11;	% Bohr radius in m
Ry = 13.605693;		% Ry in eV
beta = eF*Ry * 1.16045221e4 / Temp;  % in unit of eF
beta1 = 2*Ry * 1.16045221e4 / Temp;  % in unit of 2Ry