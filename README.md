# MetallicHydrogen
A MATLAB program for analyzing outputs of PIMD simulations to determine instability toward forming a superconductivitng state.

## Additional software
The program utilizes KSSOLV--a MATLAB toolbox for solving the Kohn-Sham equation. The description of the toolbox can be found in https://dl.acm.org/citation.cfm?id=1499096.1499099. It could be downloaded in http://www.crd.lbl.gov/~chao/KSSOLV/.

The toolbox should be installed somewhere. The variable 'kssolvpath' in 'parameters.m' should be modifty to reflect its installed path.

Our program only uses the toolbox to generate the mesh of wave-vectors.

## Running environment
The program is tested and run under MATLAB R2015b on a GPU workstation which equips with 4 Nvidia Tesla K80 (8 GPUs).

The MATLAB has the "Distributed Computing Toolbox" installed.  The toolbox is used in runGbar.m to dispatch sun-tasks.

The GPUs are heavily used. The number of GPUs should be properly set in 'parameters.m'.

## Setting parameters
Various running parameters should be set in 'parameters.m'.  See the file for the explanations of the parameters.

## Workflow
Under the MATLAB prompt:

>> runGbar
>> Gbar_collect

The calculation results can then be checked by typing in variable names. Some of them are:

* rho -- the maximal eigenvalue the Eliashberg equations

* drho -- the standard deviation of rho

* omega2 -- the average phonon frenqiency, in $k_B*T$ 

* domega2 -- the standard deviation of omega2

* lambdaVM -- the values of $\lambda(n)$

* lambdaVM1 -- the values of $\lambda(n)$ after eliminating the discretization error of PIMD by oversampling.

