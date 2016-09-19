
#ifndef boundary_conditions_h
#define boundary_conditions_h
#include<deal.II/base/smartpointer.h>
#include<deal.II/base/convergence_table.h>
#include<deal.II/base/quadrature_lib.h>
#include<deal.II/base/quadrature_selector.h>
#include<deal.II/base/parsed_function.h>
#include<deal.II/base/utilities.h>

#include<deal.II/lac/full_matrix.h>
#include<deal.II/lac/sparse_matrix.h>
#include<deal.II/lac/matrix_lib.h>
#include<deal.II/lac/vector.h>
#include<deal.II/lac/solver_control.h>
#include<deal.II/lac/solver_gmres.h>
#include<deal.II/lac/precondition.h>

#include<deal.II/grid/tria.h>
#include<deal.II/grid/tria_iterator.h>
#include<deal.II/grid/tria_accessor.h>
#include<deal.II/grid/grid_generator.h>
#include<deal.II/grid/grid_in.h>
#include<deal.II/grid/grid_out.h>
#include<deal.II/grid/tria_boundary_lib.h>

#include<deal.II/dofs/dof_handler.h>
#include<deal.II/dofs/dof_accessor.h>
#include<deal.II/dofs/dof_tools.h>
#include<deal.II/dofs/dof_renumbering.h>

#include<deal.II/fe/fe_q.h>
#include<deal.II/fe/fe_values.h>
#include<deal.II/fe/fe_system.h>
#include<deal.II/fe/mapping_q1_eulerian.h>
#include<deal.II/fe/mapping_q1.h>

#include<deal.II/numerics/data_out.h>
#include<deal.II/numerics/vector_tools.h>
#include<deal.II/numerics/solution_transfer.h>
#include <deal.II/base/types.h>

// And here are a few C++ standard header
// files that we will need:
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <set>

#include "../include/bem_problem.h"
#include "../include/computational_domain.h"

#include <deal2lkit/parsed_data_out.h>
#include <deal2lkit/utilities.h>

/**
* - BoundaryCondition. The class handles the boundary conditions. In particular
*   - it reads the boundary conditions for the potential and its normal derivative;
*   - given the peculiarities of the BEM, the boundary conditions represent the actual unknowns, thus it creates the vectors containing the variables and fills them with the proper data;
*   - it performs the error analysis on both unknowns.
*/
template <int dim>
class BoundaryConditions : public ParameterAcceptor
{
public:
  BoundaryConditions(ComputationalDomain<dim> &comp_dom, BEMProblem<dim> &bem, const MPI_Comm comm = MPI_COMM_WORLD) :
    wind(dim), comp_dom(comp_dom), bem(bem),
    mpi_communicator (comm),
    n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
    this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator)),
    pcout(std::cout,
          (this_mpi_process
           == 0)),
    data_out_scalar("Scalar data out", "vtu"),
    data_out_vector("Vector data out", "vtu")
  {
    dofs_number = 0,
    output_frequency = 1;
  }


  typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;

  virtual void declare_parameters(ParameterHandler &prm);

  virtual void parse_parameters(ParameterHandler &prm);

  void prepare_bem_vectors();

  void solve_problem();

  void output_results(const std::string);

  void compute_errors();

  const TrilinosWrappers::MPI::Vector &get_phi();

  const TrilinosWrappers::MPI::Vector &get_dphi_dn();


  std::string output_file_name;

private:

  Functions::ParsedFunction<dim> wind;

  Functions::ParsedFunction<dim> potential;

  std::string node_displacement_type;

  SolverControl solver_control;

  ComputationalDomain<dim> &comp_dom;

  BEMProblem<dim> &bem;

  types::global_dof_index dofs_number;

  unsigned int output_frequency;

  TrilinosWrappers::MPI::Vector        tmp_rhs;
  TrilinosWrappers::MPI::Vector        phi;
  TrilinosWrappers::MPI::Vector        dphi_dn;

  MPI_Comm mpi_communicator;

  unsigned int n_mpi_processes;

  unsigned int this_mpi_process;

  bool have_dirichlet_bc;

  IndexSet this_cpu_set;

  ConditionalOStream pcout;

  ParsedDataOut<dim-1, dim> data_out_scalar;
  ParsedDataOut<dim-1, dim> data_out_vector;

};

#endif
