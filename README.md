# Themis

Themis is a statistichal mechanics software designed to obtain the association
thermodynamics of two structures (ions, molecules, crystals, nanoparticles, etc). 
It generates a configurational partition function by systematically sampling 
the phase space using discrete grids to perform translations and rotations of 
one structure around another. Interaction energy for each microstate can be 
obtained by one of the potentials implemented or by using external softwares.

Themis is a free software written in Fortran 2003 language, being available at
http://www.lqt.dq.ufscar.br/lqt/lqt_software-pt.html under the GPLv3+ License.
It runs under Linux environment with gfortran/gcc 5.4+ compilers.

Since it was written in modules, new potential functions and analysis routines 
can be easily implemented.

# Downloading Themis source code

The source code of Themis, some utilities and documentation is available on 
http://www.lqt.dq.ufscar.br/lqt/lqt_software-pt.html  

# Install Themis

See [installation instructions](./INSTALL.md)

# Links

* [LQT webpage](http://www.lqt.dq.ufscar.br) for information about the UFSCAR Laboratory of Theoretical Chemistry
etc...

# Directory organization

* [`src`](./src): The source code
* [`utils`](./utils): Useful programs and scripts that could be used with Themis
* [`tests`](./tests): Input files for test (water dimer using spherical grids) 
* [`manual`](./manual): Themis user guide 
* [`doc`](./manual): Themis documentation

