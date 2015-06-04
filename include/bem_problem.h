//----------------------------  step-34.cc  ---------------------------
//    $Id: step-34.cc 18734 2009-04-25 13:36:48Z heltai $
//    Version: $Name$
//
//    Copyright (C) 2009, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//    Authors: Luca Heltai, Cataldo Manigrasso
//
//----------------------------  step-34.cc  ---------------------------

#ifndef bem_problem_h
#define bem_problem_h
// @sect3{Include files}

// The program starts with including a bunch
// of include files that we will use in the
// various parts of the program.


#include<deal.II/base/smartpointer.h>
#include<deal.II/base/convergence_table.h>
#include<deal.II/base/quadrature_lib.h>
#include<deal.II/base/quadrature_selector.h>
#include<deal.II/base/parsed_function.h>
#include<deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>

#include<deal.II/lac/full_matrix.h>
#include<deal.II/lac/sparse_matrix.h>
#include<deal.II/lac/constraint_matrix.h>
#include<deal.II/lac/matrix_lib.h>
#include<deal.II/lac/vector.h>
#include<deal.II/lac/solver_control.h>
#include<deal.II/lac/solver_gmres.h>
#include<deal.II/lac/precondition.h>
#include<deal.II/lac/compressed_sparsity_pattern.h>
#include<deal.II/lac/sparse_direct.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

//#include <deal.II/lac/petsc_vector.h>
//#include <deal.II/lac/petsc_parallel_vector.h>
//#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
//#include <deal.II/lac/petsc_solver.h>
//#include <deal.II/lac/petsc_precondition.h>

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



// And here are a few C++ standard header
// files that we will need:
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>

#include "../include/octree_block.h"
#include "../include/local_expansion.h"
#include "../include/multipole_expansion.h"
#include "../include/ass_leg_function.h"
#include "../include/computational_domain.h"
#include "../include/bem_fma.h"
#include "../include/constrained_matrix.h"
#include "parameter_acceptor.h"


#include <mpi.h>

using namespace dealii;
//using namespace TrilinosWrappers;
//using namespace TrilinosWrappers::MPI;

template <int dim>
class BEMProblem : public ParameterAcceptor
{
public:

  typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;

  BEMProblem(ComputationalDomain<dim> &comp_dom,
             const unsigned int fe_degree,
             const MPI_Comm comm = MPI_COMM_WORLD);

  void solve(TrilinosWrappers::MPI::Vector &phi, TrilinosWrappers::MPI::Vector &dphi_dn,
             const TrilinosWrappers::MPI::Vector &tmp_rhs);

  void reinit();

  void compute_constraints(ConstraintMatrix &constraints, const TrilinosWrappers::MPI::Vector &tmp_rhs);

  //  private:

  virtual void declare_parameters(ParameterHandler &prm);

  virtual void parse_parameters(ParameterHandler &prm);

  // To be commented

  void compute_alpha();

  void assemble_system();

  // The next three methods are
  // needed by the GMRES solver:
  // the first provides result of
  // the product of the system
  // matrix (a combination of Neumann
  // and Dirichlet matrices) by the
  // vector src. The result is stored
  // in the vector dst.

  void vmult(TrilinosWrappers::MPI::Vector &dst, const TrilinosWrappers::MPI::Vector &src) const;

  // The second method computes the
  // right hand side vector of the
  // system.

  void compute_rhs(TrilinosWrappers::MPI::Vector &dst, const TrilinosWrappers::MPI::Vector &src) const;

  // The third method computes the
  // product between the solution vector
  // and the (fully populated) sytstem
  // matrix.

  //TODO CHECK
  void assemble_preconditioner();

  void solve_system(TrilinosWrappers::MPI::Vector &phi, TrilinosWrappers::MPI::Vector &dphi_dn,
                    const TrilinosWrappers::MPI::Vector &tmp_rhs);


  void output_results(const std::string);

  //TODO PARALLELIZE
  void compute_surface_gradients(const TrilinosWrappers::MPI::Vector &tmp_rhs);
  //TODO PARALLELIZE
  void compute_gradients(const TrilinosWrappers::MPI::Vector &phi, const TrilinosWrappers::MPI::Vector &dphi_dn);
  //TODO PARALLELIZE
  void compute_normals();

  // this method is needed to
  // separate Dirichlet dofs from
  // Neumann nodes.

  void compute_dirichlet_and_neumann_dofs_vectors();


  // in the imported mesh, the nodes on the
  // domain edges are doubled: this routine
  // creates a std::vector of std::set which
  // allows to relate each node to their
  // double(s)

  void compute_double_nodes_set();




  ConditionalOStream pcout;
  ComputationalDomain<dim> &comp_dom;

  FE_Q<dim-1,dim>                   fe;
  DoFHandler<dim-1,dim>             dh;
  FESystem<dim-1,dim>      gradient_fe;
  DoFHandler<dim-1,dim>    gradient_dh;



  // An Eulerian Mapping is created to deal
  // with the free surface and boat mesh
  // deformation

  MappingQ<dim-1, dim>      mapping;
  Vector<double> map_points;


  // these are the std::vectors of std::sets
  // containing informations on multiple
  // nodes on the edges: one vector is
  // created for the points associated with
  // the degrees of freedom of the potential
  // function, and one is created for the
  // points associated with the degrees of
  // freedom of its gradient (a vector field)

  std::vector <std::set<unsigned int> >   double_nodes_set;
  std::vector <std::set<unsigned int> >   gradient_double_nodes_set;



  // the following vectors are needed to
  // treat Dirichlet and Neumann nodes
  // differently. Each component of the
  // first one is null if it corresponds
  // to a Dirichlet node, and zero if
  // it corresponds to a Neumann node.
  // The second vector has instead null
  // entries for Dirichlet nodes, and ones
  // for Neumann nodes

  // the number of standard quadrature points
  // and singular kernel quadrature to be
  // used

  std_cxx1x::shared_ptr<Quadrature<dim-1> > quadrature;
  unsigned int singular_quadrature_order;


  TrilinosWrappers::SparsityPattern full_sparsity_pattern;
  TrilinosWrappers::SparseMatrix neumann_matrix;
  TrilinosWrappers::SparseMatrix dirichlet_matrix;
  //FullMatrix<double>    neumann_matrix;
  //FullMatrix<double>    dirichlet_matrix;

  TrilinosWrappers::MPI::Vector        system_rhs;

  TrilinosWrappers::MPI::Vector              sol;
  TrilinosWrappers::MPI::Vector              alpha;

  mutable TrilinosWrappers::MPI::Vector              serv_phi;
  mutable  TrilinosWrappers::MPI::Vector              serv_dphi_dn;
  TrilinosWrappers::MPI::Vector              serv_tmp_rhs;

  ConstraintMatrix     constraints;

  std::string solution_method;

  SolverControl solver_control;

  TrilinosWrappers::PreconditionILU preconditioner;

  TrilinosWrappers::SparsityPattern preconditioner_sparsity_pattern;

  TrilinosWrappers::SparseMatrix band_system;

  int preconditioner_band;

  bool is_preconditioner_initialized;

  MPI_Comm mpi_communicator;

  unsigned int n_mpi_processes;

  unsigned int this_mpi_process;

  ///TODO all the story of the double nodes imposes all the procs to know almost everything.
  TrilinosWrappers::MPI::Vector dirichlet_nodes;
  TrilinosWrappers::MPI::Vector neumann_nodes;
  // Vector<double> dirichlet_nodes;
  // Vector<double> neumann_nodes;

  IndexSet this_cpu_set;

  IndexSet ghosted_set;

  IndexSet vector_this_cpu_set;

  // std::vector<Point<dim> > node_surface_gradients;

  // std::vector<Point<dim> > node_gradients;
  //
  // std::vector<Point<dim> > node_normals;

  TrilinosWrappers::MPI::Vector vector_gradients_solution;

  TrilinosWrappers::MPI::Vector vector_surface_gradients_solution;

  TrilinosWrappers::MPI::Vector vector_normals_solution;

  std::vector<types::global_dof_index> start_per_process;

  std::vector<types::global_dof_index> vector_start_per_process;

  TrilinosWrappers::SparsityPattern vector_sparsity_pattern;// TrilinosWrappers::SparsityPattern

  ConstraintMatrix  vector_constraints;

  BEMFMA<dim> fma;
};

#endif
