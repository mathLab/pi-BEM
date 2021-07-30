#include <sys/time.h>

#include "../include/driver.h"
#include "Teuchos_TimeMonitor.hpp"

using Teuchos::RCP;
using Teuchos::Time;
using Teuchos::TimeMonitor;

RCP<Time> TotalTime  = Teuchos::TimeMonitor::getNewTimer("Total Time");
RCP<Time> MeshTime   = Teuchos::TimeMonitor::getNewTimer("Mesh Time");
RCP<Time> OutputTime = Teuchos::TimeMonitor::getNewTimer("Output Time");
RCP<Time> SolveTime  = Teuchos::TimeMonitor::getNewTimer("Solve Time");

using namespace std;

template <int dim>
Driver<dim>::Driver(unsigned int n_components)
  : pcout(std::cout)
  , mpi_communicator(MPI_COMM_WORLD)
  , computational_domain(mpi_communicator)
  , bem_problem(computational_domain, mpi_communicator, n_components)
  , boundary_conditions(computational_domain,
                        bem_problem,
                        mpi_communicator,
                        n_components)
  , prm()
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , n_components(n_components)
{
  pcout.set_condition(this_mpi_process == 0);
  // PathSearch search_prm("PARAMETER");

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
void
Driver<dim>::declare_parameters(ParameterHandler &prm)
{
  pcout << "Driver::declare_parameters" << std::endl;
  prm.declare_entry("Set Global Refinement", "true", Patterns::Bool());
  prm.declare_entry("Potential components", "1", Patterns::Integer());
  prm.declare_entry(
    "Complex components",
    "",
    Patterns::List(Patterns::Integer(0)),
    "Lists the indices of the real components of complex-valued problems; empty means no complex problem is present");
}

template <int dim>
void
Driver<dim>::parse_parameters(ParameterHandler &prm)
{
  global_refinement = prm.get_bool("Set Global Refinement");
  n_components      = prm.get_integer("Potential components");

  // what's the order of parse_parameter calls?
  boundary_conditions.set_n_phi_components(n_components);
  bem_problem.set_n_phi_components(n_components);

  complex_problems.clear();
  std::vector<std::string> complex_problems_list =
    Utilities::split_string_list(prm.get("Complex components"));
  if (complex_problems_list.size())
    {
      for (unsigned int i = 0; i < complex_problems_list.size(); ++i)
        {
          std::istringstream key_reader(complex_problems_list[i]);
          types::manifold_id id;
          key_reader >> id;
          if (id)
            {
              Assert(
                complex_problems.count(id) == 0 &&
                  complex_problems.count(id - 1) == 0,
                ExcMessage(
                  "Found soperposition of real and imaginary components of complex problems"));
            }
          complex_problems.insert(id);
        }
    }

  for (auto i : complex_problems)
    {
      pcout << "Found complex problem " << i << ", " << (i + 1) << std::endl;
    }
}

template <int dim>
void
Driver<dim>::run()
{
  {
    Teuchos::TimeMonitor LocalTimer(*TotalTime);
    unsigned int         local_refinement_cycles = 0;
    {
      Teuchos::TimeMonitor LocalTimer(*MeshTime);
      // computational_domain.create_initial_mesh();
      computational_domain.read_domain();
      if (global_refinement)
        {
          computational_domain.refine_and_resize(computational_domain.n_cycles);
        }
      else
        {
          // computational_domain.conditional_refine_and_resize(1);
          computational_domain.refine_and_resize(
            computational_domain.pre_global_refinements);
          local_refinement_cycles = computational_domain.n_cycles;
        }

      // computational_domain.generate_octree_blocking();
      computational_domain.update_triangulation();
    }

    for (unsigned int i = 0; i <= local_refinement_cycles; ++i)
      {
        {
          Teuchos::TimeMonitor LocalTimer(*SolveTime);
          bem_problem.reinit();
          MPI_Barrier(MPI_COMM_WORLD);

          for (unsigned int i = 0; i < boundary_conditions.n_phi_components();
               ++i)
            {
              boundary_conditions.set_current_phi_component(i);
              if (!complex_problems.count(i))
                {
                  pcout << "solving for component " << i << std::endl;
                  // this is a purely real problem
                  boundary_conditions.solve_problem(i == 0);
                }
              else
                {
                  pcout << "solving for components " << i << " real and "
                        << i + 1 << " imaginary" << std::endl;
                  // this is a complex-valued problem
                  boundary_conditions.solve_complex_problem(i == 0);

                  ++i;
                }
              MPI_Barrier(MPI_COMM_WORLD);
            }

          boundary_conditions.set_current_phi_component(0);
        }

        if (!global_refinement && i < local_refinement_cycles)
          {
            // Compute error estimator and local refinement strategy
            MPI_Barrier(MPI_COMM_WORLD);
            if (boundary_conditions.n_phi_components() == 1)
              {
                bem_problem.adaptive_refinement(boundary_conditions.get_phi());
              }
            else
              {
                // drive refinement using phi's norm
                bem_problem.adaptive_refinement(
                  boundary_conditions.get_phi_components_norm());
              }

            computational_domain.update_triangulation();
          }
      }

    for (unsigned int i = 0; i < boundary_conditions.n_phi_components(); ++i)
      {
        pcout << "output component " << i + 1 << std::endl;

        std::string filename = boundary_conditions.output_file_name;
        if (i)
          {
            filename += "_" + Utilities::int_to_string(i + 1);
          }

        boundary_conditions.set_current_phi_component(i);
        if (!complex_problems.count(i))
          {
            pcout << "error for component " << i + 1 << std::endl;
            boundary_conditions.compute_errors(false, true);
          }
        else
          {
            pcout << "error for components " << i + 1 << " real and " << i + 2
                  << " imaginary" << std::endl;
            boundary_conditions.compute_errors(true, true);
            boundary_conditions.output_results(filename);

            ++i;
            filename += "_" + Utilities::int_to_string(i + 1);
            boundary_conditions.set_current_phi_component(i);
            boundary_conditions.compute_errors(true, false);
          }
        boundary_conditions.output_results(filename);
      }
    boundary_conditions.set_current_phi_component(0);

    /* sample output using the norm vectors
    if (boundary_conditions.n_phi_components() > 1)
      {
        const auto phi_comp_normed =
          Vector<double>(boundary_conditions.get_phi_components_norm());
        const auto dphi_dn =
          Vector<double>(boundary_conditions.get_dphi_dn_components_norm());

        if (!this_mpi_process)
          {
            // do something with it
            DataOut<dim - 1, DoFHandler<dim - 1, dim>> dataout_scalar;

            dataout_scalar.attach_dof_handler(bem_problem.dh);

            dataout_scalar.add_data_vector(
              phi_comp_normed,
              "phi",
              DataOut<dim - 1, DoFHandler<dim - 1, dim>>::type_dof_data);

            dataout_scalar.add_data_vector(
              dphi_dn,
              "dphi_dn",
              DataOut<dim - 1, DoFHandler<dim - 1, dim>>::type_dof_data);

            //dataout_scalar.add_data_vector(
            //  alpha,
            //  "alpha",
            //  DataOut<dim - 1, DoFHandler<dim - 1, dim>>::type_dof_data);

            dataout_scalar.build_patches(
              *bem_problem.mapping,
              bem_problem.mapping_degree,
              DataOut<dim - 1, DoFHandler<dim - 1, dim>>::curved_inner_cells);

            std::ofstream file_scalar("result_scalar_components_norm.vtu");

            dataout_scalar.write_vtu(file_scalar);
          }
        }
    */
    // }
  }
  // Write a summary of all timers
  Teuchos::TimeMonitor::summarize();
}

template class Driver<2>;
template class Driver<3>;
