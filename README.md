# $\pi$-BEM: Parallel BEM Solver 

[![Build Status](https://travis-ci.org/mathLab/pi-BEM.svg?branch=master)](https://travis-ci.org/mathLab/pi-BEM)

The library represents a parallel solver for the Laplace equation through Boundary Element Methods. We have developed the software in C++ on top of many high performance libraries, the [deal.II](https://github.com/dealii/dealii) library for Finite Element Handling, the [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/download) project and [Trilinos](https://trilinos.org/) library for automatic Workload balance, [OpenCASCADE](https://www.opencascade.com/) for CAD integration, and [deal2lkit](https://github.com/luca-heltai/deal2lkit) for parameter handling. 

## Provided features

We provide the following capabilities

- Read of the grid through an external file in one of the following formats (.prm, .msh, .vtk), the files specifies the kind of boundary condition to be applied on the nodes.
- Possibility of solving mixed Dirichlet Neumann boundary value problem.
- Automatic treatment of sharp edges via the double nodes technique.
- Integration of complex geometry descriptors (CAD files through OpenCASCADE).
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

## Install Procedure using CANDI
To install from scratch all the needed library you can look to the automatic installation procedure provided by [CANDI](https://github.com/koecher/candi) developed Uwe Köcher.

## Install Procedure using spack
Just follow the [instructions](https://github.com/dealii/dealii/wiki/deal.II-in-Spack) to install dealii@develop using spack.

## Install Procedure using Docker
We provide the possibility of using Docker as a tool to provide a fully operational environment for our library. To use such tool you need to install Docker following the [instructions](https://docs.docker.com/engine/installation/) provided by its authors. Then you can execute the following command line instruction 

	docker run -v `pwd`:/pwd_to_your_own_directory/ -i -t mathlab/deal2lkit:latest bash

to retrieve the environment. In such a shell you can easily compile the $\pi$-BEM library following its own instructions.

## Install Procedure from scratch
In order to successfully compile the code you need 

- to install the Trilinos and Metis wrappers of the library, see the official [instructions](https://www.dealii.org/developer/readme.html) 
- to install the [deal.II](https://github.com/dealii/dealii) library allowing both for multiprocessors and multithreaded environment.
- to install the [deal2lkit](https://github.com/mathLab/deal2lkit) library allowing both for multiprocessors and multithreaded environment.



### deal.II Installation procedure
Follow the detailed [instruction](https://www.dealii.org/developer/readme.html) to install deal with METIS and Trilinos wrappers. We highlight that in order to fully exploit $\pi$-BEM you need to properly install the following additional packages: MPI, TBB, METIS, TRILINOS ans OPENCASCADE. For more detailed instruction you can look to the the deal.ii install procedures. In the following we provide an example of the installation of all the proper packages.


### OpenCASCADE Installation procedure
- Download the latest version at [OpenCASCADE](https://github.com/tpaviot/oce) 
- Follow the instructions


### METIS-PARMETIS Installation procedure
- Download the latest version at [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/download)
- Follow the [instructions](http://glaros.dtc.umn.edu/gkhome/metis/metis/download) for the correct installation of the package


### Trilinos Installation procedure
- Download the latest version at [Trilinos](https://trilinos.org/download/)
- This is a possible configuration file

		cmake -D Trilinos_ENABLE_OPTIONAL_PACKAGES:BOOL=ON \
		-D Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON \
		-D CMAKE_CXX_FLAGS:STRING="-O3" \
		-D CMAKE_C_FLAGS:STRING="-O3" \
		-D CMAKE_VERBOSE_MAKEFILE:BOOL=FALSE \
		-D Trilinos_VERBOSE_CONFIGURE:BOOL=FALSE \
		-D TPL_ENABLE_MPI:BOOL=ON \
		-D CMAKE_BUILD_TYPE:STRING=RELEASE \
		-D Trilinos_ENABLE_Fortran:BOOL=ON \
		-D BLAS_LIBRARY_NAMES:STRING="blas" \
		-D BLAS_LIBRARY_DIRS:PATH=/usr/lib/libblas/ \
		-D TPL_BLAS_LIBRARIES:PATH=/usr/lib/libblas/ \
		-D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="" \
		-D CMAKE_INSTALL_PREFIX:PATH=/usr/local/ \
		-D BUILD_SHARED_LIBS:BOOL=ON ../
		
		make install

- if some packages create conflicts they can be disabled as seen on the Trilinos webpage





### deal2lkit Installation procedure
Follow the detailed [instruction](https://https://github.com/mathLab/deal2lkit) to install deal2lkit.


### $\pi$-BEM Installation procedure

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

If you want you can run a preliminary execution in the build library typing
	
	mpirun -np 1 bem_fma_2d
	
this will automatically generate the parameter file for the bi-dimensional run while 

	mpirun -np 1 bem_fma_3d
	
will create a proper parameter file for a 3 dimensional simulation.

#Notice to developers

Before making a pull request, please make sure you run the script

    ./scripts/indent

from the base directory of this project, to ensure that no random 
white space changes are inserted in the repository.

The script requires Artistic Style Version 2.04 (astyle) to work 
properly.

#Licence

Please see the file [./LICENSE](https://github.com/mathLab/pi-BEM/blob/master/LICENSE) for details



