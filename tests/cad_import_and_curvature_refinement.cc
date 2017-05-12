//-----------------------------------------------------------
//
//    Copyright (C) 2014 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

// Read a file in iges format, and write it out again in the same
// format.

#include "tests.h"
#include "computational_domain.h"

int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  //initlog();
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);
  ComputationalDomain<3> computational_domain(mpi_communicator);
  ParameterAcceptor::initialize(SOURCE_DIR "/parameter_cad_import_and_curvature_refinement.prm","used.prm");
  computational_domain.read_domain();
  computational_domain.refine_and_resize(computational_domain.n_cycles);
}
