// Verify that the unsupported 2D hypersingular gradient path fails before
// accessing mesh/vector data, in both Debug and Release.
#include "bem_problem.h"
#include "tests.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  initlog();
  ComputationalDomain<2>        domain;
  BEMProblem<2>                 problem(domain);
  TrilinosWrappers::MPI::Vector phi, normal_derivative;
  try
    {
      problem.compute_gradients_hypersingular(phi, normal_derivative);
    }
  catch (const ExcNotImplemented &)
    {
      deallog << "2D hypersingular gradients rejected safely" << std::endl;
      return 0;
    }
  AssertThrow(false, ExcMessage("Expected ExcNotImplemented"));
}
