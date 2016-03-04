#ifndef preconditioner_bem_h
#define preconditioner_bem_h

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

#include <deal.II/base/types.h>

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
#include <deal.II/fe/mapping_fe_field.h>

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



#include <mpi.h>

#include <bem_problem_access.h>

#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_finite_element.h>
#include <deal2lkit/parsed_grid_refinement.h>
#include <deal2lkit/utilities.h>




using namespace dealii;
using namespace deal2lkit;

template <int dim>//, Type V>
class PreconditionerBEM : public BEMProblemAccess<dim>
{
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * This function is only present to provide the interface of a
   * preconditioner to be handed to a smoother.  This does nothing.
   */
  struct AdditionalData
  {
    /**
     * Constructor.
     */
    AdditionalData () {}
  };

  /**
   * Constructor, sets the domain and range sizes to their defaults.
   */
  // PreconditionerBEM(){};

  // void initialize(const unsigned int quad_order, const unsigned int sing_quad_order);
  /**
   * The matrix argument is ignored and here just for compatibility with more
   * complex preconditioners.
   */
public:
  void build_coarse_inverse(unsigned int level_mg=0);

  void build_projector();//const DoFHandler<dim-1, dim> &dh_fine_in


  /**
   * Apply preconditioner.
   */
  void vmult (TrilinosWrappers::MPI::Vector &out, const TrilinosWrappers::MPI::Vector &in) const;

  // /**
  //  * Apply transpose preconditioner. to be implemented
  //  */
  // void Tvmult (TrilinosWrappers::MPI::Vector &out, const TrilinosWrappers::MPI::Vector &in) const;
  //
  // /**
  //  * Apply preconditioner, adding to the previous value.
  //  */
  // void vmult_add (TrilinosWrappers::MPI::Vector &out, const TrilinosWrappers::MPI::Vector &in) const;
  //
  // /**
  //  * Apply transpose preconditioner, adding. to be implemented
  //  */
  // void Tvmult_add (TrilinosWrappers::MPI::Vector &out, const TrilinosWrappers::MPI::Vector &in) const;

  /**
   * This function is only present to provide the interface of a
   * preconditioner to be handed to a smoother.  This does nothing.
   */
  void clear () {}

  /**
   * Return the dimension of the codomain (or range) space. To remember: the
   * matrix is of dimension $m \times n$.
   *
   * @note This function should only be called if the preconditioner has been
   * initialized.
   */
  size_type m () const;

  /**
   * Return the dimension of the domain space. To remember: the matrix is of
   * dimension $m \times n$.
   *
   * @note This function should only be called if the preconditioner has been
   * initialized.
   */
  size_type n () const;

  private:
  /**
   * The dimension of the range space.
   */
  size_type n_rows;

  /**
   * The dimension of the domain space.
   */
  size_type n_columns;

  SparseMatrix<double> prec_matrix;
  SparsityPattern sp_coarse_full;
  SparseDirectUMFPACK prec_solver;
  // TrilinosWrappers::SparseMatrix projector;
  // TrilinosWrappers::SparsityPattern projector_sp;
  ConstraintMatrix constraints_coarse;
  DoFHandler<dim-1, dim> *dh_coarse;
  size_type n_coarse;

};

#endif
