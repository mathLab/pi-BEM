#ifndef driver_h
#define driver_h

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>

#include "bem_problem.h"
#include "bem_fma.h"
#include "boundary_conditions.h"
#include "computational_domain.h"

#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/utilities.h>
using namespace dealii;

/**
 * This class is in charge of organising the overall BEM simulation,
 * while timing the various phases. It has interfaces with all the
 * other classes in order to have a complete simulation.
 */
template <int dim>
class Driver : public ParameterAcceptor
{
public:
  /**
   * Initialize variables and setup parameters.
   */
  Driver();

  /**
   * Sets up BEMProblem and ComputationalDomain, then solve BEM.
   */
  void run();

private:

  /// Monitor output only on first processor.
  ConditionalOStream pcout;

  /// Mpi communicator object.
  MPI_Comm mpi_communicator;

  /// Domain related object.
  ComputationalDomain<dim> computational_domain;

  /// The actual bem solver.
  BEMProblem<dim> bem_problem;

  /// Boundary conditions.
  BoundaryConditions<dim> boundary_conditions;

  /// Flag to use global refinement.
  bool global_refinement;

  /// Total number of processors
  const unsigned int n_mpi_processes;

  /// The id of this processor.
  const unsigned int this_mpi_process;

  /// Time monitor object.
  TimeMonitor timer;
};

#endif
