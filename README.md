# Themis

Themis is a statistichal mechanics software designed to obtain the association
thermodynamics of two structures (ions, molecules, crystals, nanoparticles, etc). 
It generates a configurational partition function by systematically sampling 
the phase space using discrete grids to perform translations and rotations of 
one structure around another. Interaction energy for each microstate can be 
obtained by one of the potentials implemented or by using external softwares.  

Themis is a free software written in Fortran 2003 language, being available at
https://github.com/colombarifm/themis under the GPLv3+ License. It runs under 
Linux environment with gfortran/gcc 5.4+ compilers.  

Since it was written in modules, new potential functions and analysis routines 
can be easily implemented.

# Install Themis

See [installation instructions](./INSTALL.md)  

# Links for useful programs

* [COM](https://github.com/colombarifm/com) program to calculate the center of mass of molecular structures
* [SAS_GRID](https://github.com/colombarifm/sas_grid) program to generate the solvent accessible surface translation grid for molecular structures

# Directory organization

* [`src`](./src): The source code
* [`utils`](./utils): Useful scripts that could be used with Themis
* [`tests`](./tests): Input files for test (water dimer using spherical grids) 
* [`manual`](./manual): Themis user guide 
* [`doc`](./doc): Themis documentation

