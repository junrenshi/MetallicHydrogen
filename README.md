# MetallicHydrogen
A MATLAB program for analyzing outputs of PIMD simulations to determine instability toward forming a superconducting state.

## Additional software
The program uses KSSOLV--a MATLAB toolbox for solving the Kohn-Sham equation. The description of the toolbox can be found in https://dl.acm.org/citation.cfm?id=1499096.1499099. It could be downloaded in http://www.crd.lbl.gov/~chao/KSSOLV/.

The toolbox should be installed somewhere loadable by MATLAB. The variable 'kssolvpath' in 'parameters.m' should be modify to reflect its installation path.

Our program only uses the toolbox to generate the mesh of wave-vectors.

## Running environment
The program is tested and run under MATLAB R2015b on a GPU workstation equipping with 4 Nvidia Tesla K80 GPU cards (8 GPUs).

The MATLAB should have the "Distributed Computing Toolbox" installed.  The toolbox is used in runGbar.m to dispatch sub-tasks.

The GPUs are heavily utilized to accelerate calculations. The number of GPUs may be set in 'parameters.m'.

## Setting parameters
Various parameters can be be set in 'parameters.m'.  See comments in the file for explanations of the parameters.

## Work flow
Under the MATLAB prompt:

```
>> runGbar
>> Gbar_collect
```

Calculation results can then be checked by typing in variable names in MATLAB prompt. Some of them are:

* rho--the maximal eigenvalue of the Eliashberg equations
* drho--the standard deviation of rho
* omega2--the average phonon frequency $\omega_2$, in $k_B T$ 
* domega2--the standard deviation of $\omega_2$
* lambdaVM--the values of $\lambda(n)$
* lambdaVM1--the values of $\lambda(n)$ after oversampling to eliminate the discretization error of PIMD.

