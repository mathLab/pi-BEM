#Parallel BEM Solver for deal.II

[![Build Status](https://travis-ci.org/mathLab/pi-BEM.svg?branch=master)](https://travis-ci.org/mathLab/pi-BEM)

The library represents a parallel solver for the Laplace equation through Boundary Element Methods. We have developed the software on top of the deal.II library. 

## Provided features

We provide the following capabilities

- Read of the grid through an external file in one of the following formats (.prm, .msh, .vtk), the files specifies the kind of boundary condition to be applied on the nodes.
- Possibility of solving mixed Dirichlet Neumann boundary value problem.
- Automatic treatment of sharp edges via the double nodes technique.
- Usage of Lagrangian Finite Elements of arbitrary order. We also provide interfaces with discontinuous elements.
- Distributed memory (MPI) parallelisation of the standard collocation BEM for the Laplace equation.
- Coupling with a FMM to get a performance improvement.
- Hybrid Distributed (MPI) - Shared (TBB) memory parallelisation for the BEM-FMM code

## Code Structure

## Install Procedure
In order to successfully compile the code you need 

- to install the [deal.II](https://github.com/dealii/dealii) library allowing both for multiprocessors and multithreaded environment.
- to install the Trilinos and Metis wrappers of the library, see the official [instructions](https://www.dealii.org/developer/readme.html) 
- to install the pi-BEM library.

After you have compiled your application, you can run 

	make test

or
	
	ctest 

to start the testsuite.

Take a look at
https://www.dealii.org/developer/developers/testsuite.html for more
information on how to create tests and add categories of tests.

#Notice to developers

Before making a pull request, please make sure you run the script

    ./scripts/indent

from the base directory of this project, to ensure that no random 
white space changes are inserted in the repository.

The script requires Artistic Style Version 2.04 (astyle) to work 
properly.

#Licence

Please see the file [./LICENSE](https://github.com/mathLab/pi-BEM/blob/master/LICENSE) for details



