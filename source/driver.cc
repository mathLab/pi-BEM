#include "../include/driver.h"
#include <sys/time.h>


#include "Teuchos_TimeMonitor.hpp"

using Teuchos::Time;
using Teuchos::TimeMonitor;
using Teuchos::RCP;

RCP<Time> TotalTime = TimeMonitor::getNewTimer("Total Time");
RCP<Time> MeshTime = TimeMonitor::getNewTimer("Mesh Time");
RCP<Time> OutputTime = TimeMonitor::getNewTimer("Output Time");
RCP<Time> SolveTime = TimeMonitor::getNewTimer("Solve Time");

using namespace std;

template <int dim>
Driver<dim>::Driver(int argc, char **argv) :
  pcout (std::cout),
  computational_domain(1,1),
  fma(computational_domain),
  bem_problem(computational_domain, fma),
  boundary_conditions(computational_domain, bem_problem),
  prm(),
  mpi_communicator (MPI_COMM_WORLD),
  n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
  this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator))
{
  pcout.set_condition(this_mpi_process == 0);
  //PathSearch search_prm("PARAMETER");

  // Declare the parameter entries..
  DeclareParameters();

  deallog.depth_console(this_mpi_process == 0 ? 3 : 0);

  vector<string> args;
  for (int i=0; i<argc; ++i)
    args.push_back (argv[i]);

  string default_prm;
  // The default parameter file is the name of the application plus prm
  default_prm = args[0] + ".prm";

  prm.read_input(default_prm, false, true);

  for (int i=1; i<argc; ++i)
    prm.read_input(args[i], true);

  // Now that we have the final version of the parameters, parse them.
  ParseParameters();

  // And write the used ones.
  default_prm = args.front() + "_used.prm";
  ofstream outprm(default_prm.c_str());
  prm.print_parameters(outprm, ParameterHandler::ShortText);
}

template <int dim>
Driver<dim>::~Driver()
{}


template <int dim>
void Driver<dim>::run()
{
  {
    TimeMonitor LocalTimer(*TotalTime);
    {
      TimeMonitor LocalTimer(*MeshTime);
      //computational_domain.create_initial_mesh();
      computational_domain.read_domain();
      computational_domain.refine_and_resize(computational_domain.n_cycles);
      //computational_domain.generate_octree_blocking();
    }

    {
      TimeMonitor LocalTimer(*SolveTime);
      bem_problem.reinit();
      boundary_conditions.solve_problem();
    }

    // {
    //   TimeMonitor LocalTimer(*OutputTime);
    std::string filename = ( boundary_conditions.output_file_name + ".vtk" );
    boundary_conditions.compute_errors();
    boundary_conditions.output_results(filename);
    // }
  }
  // Write a summary of all timers
  TimeMonitor::summarize();
}

template <int dim>
void Driver<dim>::DeclareParameters()
{
  computational_domain.declare_parameters(prm);
  boundary_conditions.declare_parameters(prm);
  bem_problem.declare_parameters(prm);
  fma.declare_parameters(prm);
}

template <int dim>
void Driver<dim>::ParseParameters()
{
  computational_domain.parse_parameters(prm);
  boundary_conditions.parse_parameters(prm);
  bem_problem.parse_parameters(prm);
  fma.parse_parameters(prm);
}

//template class Driver<3>;
