// @sect3{Include files}

// The program starts with including a bunch of include files that we will use
// in the various parts of the program. Most of them have been discussed in
// previous tutorials already:
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/base/types.h>

// And here are a few C++ standard header files that we will need:
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <mpi.h>
#include "bem_fma.h"
#include "min_bem_operator.h"
#include <deal2lkit/parsed_grid_generator.h>
#include <deal2lkit/parsed_finite_element.h>
#include <deal2lkit/parsed_function.h>
#include <deal2lkit/parsed_data_out.h>
#include <deal2lkit/error_handler.h>
#include <deal2lkit/utilities.h>
// The last part of this preamble is to import everything in the dealii
// namespace into the one into which everything in this program will go:
namespace MinFmm
{
  using namespace dealii;


  /// @sect3{Single and double layer operator kernels}

  // First, let us define a bit of the boundary integral equation machinery.

  // The following two functions are the actual calculations of the single and
  // double layer potential kernels, that is $G$ and $\nabla G$. They are well
  // defined only if the vector $R = \mathbf{y}-\mathbf{x}$ is different from
  // zero.
  namespace LaplaceKernel
  {
    template <int dim>
    double single_layer(const Tensor<1,dim> &R)
    {
      switch (dim)
        {
        case 2:
          return (-std::log(R.norm()) / (2*numbers::PI) );

        case 3:
          return (1./( R.norm()*4*numbers::PI ) );

        default:
          Assert(false, ExcInternalError());
          return 0.;
        }
    }



    template <int dim>
    Tensor<1,dim> double_layer(const Tensor<1,dim> &R)
    {
      switch (dim)
        {
        case 2:
          return R / ( -2*numbers::PI * R.norm_square());
        case 3:
          return R / ( -4*numbers::PI * R.norm_square() * R.norm() );

        default:
          Assert(false, ExcInternalError());
          return Tensor<1,dim>();
        }
    }
  }

// @sect3{The StepFMA class}

// The structure of a boundary element method code is very similar to the
// structure of a finite element code, and so the member functions of this
// class are like those of most of the other tutorial programs. In
// particular, by now you should be familiar with reading parameters from an
// external file, and with the splitting of the different tasks into
// different modules. The same applies to boundary element methods, and we
// won't comment too much on them, except on the differences.
  template <int dim>
  class StepFMA : public ParameterAcceptor
  {
  public:
    StepFMA(const unsigned int fe_degree = 1, bool fmm_method = true, const MPI_Comm comm = MPI_COMM_WORLD);

    void run();

    void run_for_octree();

  private:

    virtual void declare_parameters (ParameterHandler &prm);

    virtual void parse_parameters (ParameterHandler &prm);

    void read_domain();

    void refine_and_resize(const bool ref=true);

    void compute_boundary_condition();

    void compute_double_nodes_set();

    void assemble_direct_system();

    void solve_system();

    void compute_errors(const unsigned int cycle);

    void output_results(const unsigned int cycle);

    void save_direct_solution(const unsigned int cycle);

    void read_direct_solution(const unsigned int cycle, Vector<double> &phi_direct, Vector<double> &alpha_dir);

    const Quadrature<dim-1> & get_singular_quadrature(const unsigned int index) const;

    MPI_Comm mpi_communicator;

    // shared_ptr<Triangulation<dim> >       tria;
    // shared_ptr<FiniteElement<dim,dim> >   fe;
    // shared_ptr<DoFHandler<dim> >          dh;

    Triangulation<dim-1, dim>   tria;
    FE_Q<dim-1,dim>             fe;
    DoFHandler<dim-1,dim>       dh;
    MappingQ<dim-1, dim>      mapping;


    TrilinosWrappers::MPI::Vector              phi;
    TrilinosWrappers::MPI::Vector              dphi_dn;


    // The following variables are the ones that we fill through a parameter
    // file.  The new objects that we use in this example are the
    // Functions::ParsedFunction object and the QuadratureSelector object.
    //
    // The Functions::ParsedFunction class allows us to easily and quickly
    // define new function objects via parameter files, with custom
    // definitions which can be very complex (see the documentation of that
    // class for all the available options).
    //
    // We will allocate the quadrature object using the QuadratureSelector
    // class that allows us to generate quadrature formulas based on an
    // identifying string and on the possible degree of the formula itself. We
    // used this to allow custom selection of the quadrature formulas for the
    // standard integration, and to define the order of the singular
    // quadrature rule.
    //
    // We also define a couple of parameters which are used in case we wanted
    // to extend the solution to the entire domain.

    Functions::ParsedFunction<dim> exact_dphi_dn_solution;
    Functions::ParsedFunction<dim> exact_phi_solution;

    unsigned int singular_quadrature_order;
    std_cxx11::shared_ptr<Quadrature<dim-1> > quadrature;

    SolverControl solver_control;

    unsigned int initial_ref;
    unsigned int n_cycles;
    unsigned int n_cycles_octree;
    unsigned int external_refinement;

    bool from_theory;
    bool run_in_this_dimension;
    bool extend_solution;
    IndexSet this_cpu_set;
    TrilinosWrappers::MPI::Vector dirichlet_nodes;
    TrilinosWrappers::MPI::Vector              dirichlet_values;
    TrilinosWrappers::MPI::Vector              neumann_values;
    std::vector <std::set<types::global_dof_index> >   double_nodes_set;
    bool fmm_sol;
    ErrorHandler<2> eh;
    unsigned int n_mpi_processes;
    unsigned int this_mpi_process;
    ConditionalOStream pcout;

    TrilinosWrappers::SparseMatrix system_matrix;//(const size_type m, const size_type n, const unsigned int n_max_entries_per_row);
    TrilinosWrappers::SparsityPattern tril_sp;
    ConstraintMatrix constraints;
    TrilinosWrappers::MPI::Vector              system_rhs;
    TrilinosWrappers::MPI::Vector              system_alpha;
    BEMFMA<dim> fma;

  };
}
