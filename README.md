# 2D-Ising-Model
A set of c++ and python codes which simulate and analyze the 2D Ising Model

To use:

The file 2DIsingSimulation.cpp is used to generate data files to analyze the behavior of the 2D Ising model. To determine the side length of the lattice (for a total size of L*L for the lattice) change the variable int L. To increase accuracy (and, obviously, runtime) of the simulation, change int nData in the function MonteCarlo(). The lattice is sweeped over nData + nTherm times. nTherm allows the lattice to approach equilibrium before taking data.

The .cpp file will create 5 files, one each for Energy, Energy^2, Magnetization, Mag^2, and Mag^4, which are used to analyze the model. If you change the L variable, the file name will change as well. 

The .ipynb file is used to generate relevant physical quantities such as magnetization, magnetic susceptibility, specific heat, and binder cumulant using the data from the .cpp file. It also plots these quantities, as well as computes and plots the finite-size scaling using critical exponents.
