# DIRHBZ_BLV

Authors:
A. Ravlic[^1], E. Yuksel[^2], T. Niksic[^1], N. Paar[^1]

Correspondance to: ```aravlic@phy.hr```

[^1]: University of Zagreb, Faculty of Science, Physics Department, Croatia
[^2]: Department of Physics, University of Surrey, United Kingdom 
-------------------------------------------

The FORTRAN77 relativistic finite-temperature Hartree-Bogoliubov (FT-RHB) solver in harmonic oscillator basis with an explicit treatment of the particle continuum. The code assumes axial symmetry, conserved parity, and time-reversal symmetry. It is suitable for the description of hot even-even nuclei, including pairing and deformation effects. The code is based on zero-temperature RHB code DIRHBZ described in Ref. [^3]. The extension to finite temperature RHB is derived in Ref. [^4], while the Bonche-Levit-Vautherin subtraction procedure is detailed in Ref. [^5]. The DIRHBZ_BLV code implements the latter two frameworks.

[^3]: T. Niksic, N. Paar, D. Vretenar, P. Ring, Computer Physics Communications 185, 6, 1808–1821 (2014).
[^4]: A. L. Goodman, Nuclear Physics A 352, 1, 30–44 (1981).
[^5]: P. Bonche, S. Levit, D. Vautherin, Nuclear Physics A 436, 2, 265–293 (1985).

The code uses the Intel compilers, MPI for parallel calculations, and MKL library for matrix operations. We recommend downloading the oneAPI tools from:

https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html

Make sure that ```mpiifort``` is installed. To compile the calculations use the provided Makefile by typing ```make```.

There are four main files containing the code:
* ```dirhbz.f```- the main file containing the FT-RHB solver routines
* ```dirhbz_vapor.f```- the FT-RHB solver for the vapor states
* ```dirhb.dat```- input file
* ```dirhb.par```- parameter file

The structure of the input ```dirhb.dat```file is as follows:
* ```n0f      =   16  16```
Number of harmonic oscillator shells $N_{osc}$ for fermion and boson states
* ```def-bas  =    -0.000```
The quadrupole deformation of the axial harmonic oscillator basis
* ```b0       =    -1.500```
The harmonic oscillator length parameter, if negative, it is calculated from the prescription $\hbar \omega = 41 A^{-1/3}$ where $A$ is the nucleon number and $\omega$ oscillator frequency.
* ```def-ini  =    -0.000```
Initial quadrupole deformation of the Woods-Saxon basis used to initialize nucleus and vapor potentials.
* ```inin     =    1    1```
Initialization of mean-field and pairing potential (1 - yes)
* ```Sm 180 ```
Nucleus under consideration
* ```Init.Gap =    1.000     1.000```
Initial value of the pairing gap $\Delta$.
* ```Force    =  DD-PC1```
The relativistic functional under consideration. Can choose among: DD-PC1, DD-PCX, DD-PCJ30-36, and DD-ME2
* ```icstr    =  1```
Constrained calculation with quadrupole constraint. It is released after the first 20 calculations.
* ```cquad    =  0.005 ```
Strength of the quadrupole constraint.
* ```temp     =  1.01 ```
Temperature in MeV.
* ```matel    =  0```
Calculate pairing matrix element (0), or read from the ```wnn.del```file (1).

-------------------------------------------

The code used MPI to split the calculation of potential energy surface on different CPUs. The number of CPUs is hardcoded to 11, however, this can be changed with minimal modification of the code. To run the code, after successful compilation, use:
```mpirun -np 11 ./run```

-------------------------------------------

The results of the calculations are stored in the ```results.out```file. It has the following column structure, where each column corresponds to 1 out of 11 CPU cores with its own initial quadrupole deformation $\beta_2$:

* neutron number
* proton number
* temperature $T$ (MeV)
* subtracted binding energy $E$ (MeV)
* subtracted entropy $S$
* isovector quadrupole deformation $\beta_2$
* neutron root-mean-square radius $r_n$
* proton root-mean-square radius $r_p$
* neutron chemical potential $\lambda_n$ (MeV)
* proton chemical potential $\lambda_p$ (MeV)
* neutron pairing gap $\Delta_n$ (MeV)
* proton pairing gap $\Delta_p$ (MeV)
* number of iterations required for convergence
* density of vapor states
* neutron emission lifetime $T_n$

Other details regarding the iteration-by-iteration convergence for each CPU are stored in ```dirhb_5xx.out```files, where xx ranges from 01 to 10.

The vector densities are plotted in an $(r,z)$ plane for the Nucleus+Vapor system ```dirhb.plo``` and Vapor only system ```dirhb_vapor.plo``` for quadrupole deformation $\beta_2$ which minimizes the free energy $F = E - TS$.

-------------------------------------------

To benchmark the compilation we have included the ```results.out```file which contains an example run with the provided ```dirhb.dat``` file.




  









