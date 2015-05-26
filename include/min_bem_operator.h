#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include "bem_fma.h"

namespace Operator
{
  using namespace dealii;

  template <int dim>
  class MinBEMOperator
  {
  public:
    template<int fdim>
    friend void test(MinBEMOperator<fdim> &);

    MinBEMOperator(const BEMFMA<dim> &fma_in, const MPI_Comm comm_in, const IndexSet &cpu_set_in, const unsigned int mpi_process_in);
    void set_alpha();
    void vmult(TrilinosWrappers::MPI::Vector &dst, const TrilinosWrappers::MPI::Vector &src)const;
    void compute_rhs(TrilinosWrappers::MPI::Vector &dst, const TrilinosWrappers::MPI::Vector &src_dir, const TrilinosWrappers::MPI::Vector &src_neum);
    void increase_fma_order();
  private:
    const BEMFMA<dim> &op_fma;
    TrilinosWrappers::MPI::Vector alpha;
    MPI_Comm mpi_communicator;
    unsigned int this_mpi_process;

    IndexSet this_cpu_set;
    unsigned int n_dofs;
  };
}
