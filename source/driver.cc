#include "../include/driver.h"

using namespace std;

template <int dim>
Driver<dim>::Driver() :
  pcout (std::cout),
  mpi_communicator (MPI_COMM_WORLD),
  computational_domain(mpi_communicator),
  bem_problem(computational_domain,mpi_communicator),
  boundary_conditions(computational_domain, bem_problem),
  global_refinement(true),
  n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
  this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator))
{
  pcout.set_condition(this_mpi_process == 0);
  add_parameter(global_refinement, "Set Global Refinement");
}

template <int dim>
void Driver<dim>::run()
{
  {
    auto t_ = timer.scoped_timer("Total Time");
    unsigned int local_refinement_cycles = 0;
    {
      auto t_ = timer.scoped_timer("Mesh Time");
      computational_domain.read_domain();
      if (global_refinement)
        {
          computational_domain.refine_and_resize(computational_domain.n_cycles);
        }
      else
        {
          computational_domain.refine_and_resize(1);
          local_refinement_cycles=computational_domain.n_cycles;
        }
      //computational_domain.generate_octree_blocking();
    }
    computational_domain.update_triangulation();
    for (unsigned int i = 0; i<=local_refinement_cycles; ++i)
      {
        {
          auto t_ = timer.scoped_timer("Solve Time");
          bem_problem.reinit();
          boundary_conditions.solve_problem();
        }
        if (!global_refinement && i<local_refinement_cycles)
          {
            // Compute error estimator and local refinement strategy
            bem_problem.adaptive_refinement(boundary_conditions.get_phi());
            computational_domain.update_triangulation();
          }

      }

    std::string filename = ( boundary_conditions.output_file_name);
    boundary_conditions.compute_errors();
    {
      auto t_ = timer.scoped_timer("Output Time");
      boundary_conditions.output_results(filename);
    }
  }
}

template class Driver<3>;
