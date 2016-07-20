#Parallel BEM Solver for deal.II

[![Build Status](https://travis-ci.org/mathLab/pi-BEM.svg?branch=master)](https://travis-ci.org/mathLab/pi-BEM)

The library represents a parallel solver for the Laplace equation through Boundary Element Methods. We have developed the software in C++ on top of the [deal.II](https://github.com/dealii/dealii) library. 

## Provided features

We provide the following capabilities

- Read of the grid through an external file in one of the following formats (.prm, .msh, .vtk), the files specifies the kind of boundary condition to be applied on the nodes.
- Possibility of solving mixed Dirichlet Neumann boundary value problem.
- Automatic treatment of sharp edges via the double nodes technique.
- Usage of Lagrangian Finite Elements of arbitrary order. We also provide interfaces with discontinuous elements.
- Distributed memory (MPI) parallelisation of the standard collocation BEM for the Laplace equation.
- Coupling with a Fast Multiple Method (FMM) to get a performance improvement.
- Hybrid Distributed (MPI) - Shared (Intel Threaded Building Block) memory parallelisation for the BEM-FMM code
- Recovery of both primal (potential) and dual (potential normal derivative) unknowns.
- L2 projection of the full 3D potential gradient for post processing.
- Extensive tuning via parameter file using the [deal2lkit](https://github.com/mathLab/deal2lkit) library 


## Code Structure
We have subdivided the code in main classes to handle the many different aspects of a complete BEM simulation.

- Driver. This class is in charge of organising the overall BEM simulation. It has interfaces with all the other classes in order to perform a complete simulation.
- ComputationalDomain. This class handles, and provides to the other classes, ONLY the geometry of the problem. In particular
 - it handles the domain decomposition using a graph partitioning tool (METIS);
 - it reads the domain from an external file.
- BoundaryCondition. The class handles the boundary conditions. In particular
	- it reads the boundary conditions for the potential and its normal derivative; 
	- given the peculiarities of the BEM, the boundary conditions represent the actual unknowns, thus it creates the vectors containing the two variables and fills them with the proper data;
	- it performs the error analysis on both unknowns.

- BEMProblem. This class is the core of the BEM simulation
	- it receives the variables vector filled in with the proper boundary condition;
	- it creates the codimension 1 functional space setting up the BEM;
	- it solves the system using a preconditioned parallel GMRES solver;
	- it eventually interacts with the FMM accelerator.
- BEMFMA. This class handles the accelerator, in particular
	- it sets up an hierarchical 3D space subdivision (octree);
	- it receives two distributed vectors representing the unknowns and performs a full FMM matrix vector product.

## Install Procedure
In order to successfully compile the code you need 

- to install the Trilinos and Metis wrappers of the library, see the official [instructions](https://www.dealii.org/developer/readme.html) 
- to install the [deal.II](https://github.com/dealii/dealii) library allowing both for multiprocessors and multithreaded environment.
- to install the [deal2lkit](https://github.com/mathLab/deal2lkit) library allowing both for multiprocessors and multithreaded environment.


Then you can clone the repository and compile it

	git clone https://github.com/mathLab/pi-BEM.git
	cd pi-BEM
	mkdir build
	cd build
	cmake ../
	make -j4

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



