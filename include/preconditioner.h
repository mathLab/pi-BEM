#ifndef preconditioner_bem_h
#define preconditioner_bem_h
// @sect3{Include files}

// The program starts with including a bunch
// of include files that we will use in the
// various parts of the program.


#include <bem_problem_access.h>

#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_finite_element.h>
#include <deal2lkit/parsed_grid_refinement.h>
#include <deal2lkit/utilities.h>



#include <mpi.h>

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
  TrilinosWrappers::SparseMatrix projector;
  TrilinosWrappers::SparsityPattern projector_sp;
  ConstraintMatrix constraints_coarse;
  DoFHandler<dim-1, dim> dh_coarse;
  size_type n_coarse;

};

#endif
