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
Driver<dim>::Driver() :
  pcout (std::cout),
  mpi_communicator (MPI_COMM_WORLD),
  computational_domain(mpi_communicator),
  bem_problem(computational_domain,1,mpi_communicator),
  boundary_conditions(computational_domain, bem_problem),
  prm(),
  n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
  this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator))
{
  pcout.set_condition(this_mpi_process == 0);
  //PathSearch search_prm("PARAMETER");

  // Declare the parameter entries..
  // DeclareParameters();
  //
  // deallog.depth_console(this_mpi_process == 0 ? 3 : 0);
  //
  // vector<string> args;
  // for (int i=0; i<argc; ++i)
  //   args.push_back (argv[i]);
  //
  // string default_prm;
  // // The default parameter file is the name of the application plus prm
  // default_prm = args[0] + ".prm";
  //
  // prm.read_input(default_prm, false, true);
  //
  // for (int i=1; i<argc; ++i)
  //   prm.read_input(args[i], true);
  //
  // // Now that we have the final version of the parameters, parse them.
  // ParseParameters();
  //
  // // And write the used ones.
  // default_prm = args.front() + "_used.prm";
  // ofstream outprm(default_prm.c_str());
  // prm.print_parameters(outprm, ParameterHandler::ShortText);
}

template <int dim>
Driver<dim>::~Driver()
{}


template <int dim>
void Driver<dim>::declare_parameters (ParameterHandler &prm)
{
  prm.declare_entry("Set Global Refinement", "true",
                    Patterns::Bool());


}

template <int dim>
void Driver<dim>::parse_parameters (ParameterHandler &prm)
{
  global_refinement = prm.get_bool("Set Global Refinement");


}


template <int dim>
void Driver<dim>::run()
{
  {
    TimeMonitor LocalTimer(*TotalTime);
    {
      TimeMonitor LocalTimer(*MeshTime);
      //computational_domain.create_initial_mesh();
      computational_domain.read_domain();
      if (global_refinement)
        computational_domain.refine_and_resize(computational_domain.n_cycles);
      else
        computational_domain.conditional_refine_and_resize(computational_domain.n_cycles);
      //computational_domain.generate_octree_blocking();
    }

    {
      TimeMonitor LocalTimer(*SolveTime);
      bem_problem.reinit();
      boundary_conditions.solve_problem();
    }

    std::string filename = ( boundary_conditions.output_file_name);
    boundary_conditions.compute_errors();
    boundary_conditions.output_results(filename);
    // }
  }
  // Write a summary of all timers
  TimeMonitor::summarize();
}


template class Driver<3>;
