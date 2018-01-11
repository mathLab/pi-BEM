
/// This class contains all the methods
/// of the Fast Multipole Algorithm
/// applied to the Boundary Element
/// Method




#ifndef bem_fma_h
#define bem_fma_h

#include "../include/octree_block.h"
#include "../include/local_expansion.h"
#include "../include/multipole_expansion.h"
#include "../include/ass_leg_function.h"
#include "../include/computational_domain.h"
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/types.h>
#include <deal.II/base/std_cxx11/bind.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>

#include <mpi.h>
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/utilities.h>

namespace Operator
{
  template <int dim> class MinBEMOperator;
}

using namespace dealii;
using namespace deal2lkit;

/**
* A class for the handling of the Fast Multiple Method coupled with the Bundary Element Method.
* It is derived from ParameterAcceptor to have a common interface with parameter files. In
* particular this class performs the matrix vector products approximating long range interactions
* via Multipole and Local Expansions that are contained in dedicated classes.

*/
template <int dim>//, Type V>
class BEMFMA : public ParameterAcceptor
{
public:

  /// Function to be used in the tests, it can access everything inside bemfma
  // template<int fdim>
  // friend void test(Operator::MinBEMOperator<fdim> &);
  /// Class that implements a minimal operator
  friend class Operator::MinBEMOperator<dim>;

  template<int gdim>
  friend void test(BEMFMA<gdim> &);

  /// Just renaming the cell iterator type

  typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;

  /// Contructor: needs only a MPI communicator
  /// to be consistent with the other classes of
  /// the problem.
  BEMFMA(MPI_Comm mpi_commy=MPI_COMM_WORLD);

  /// Destructor: does not need anything. It frees all
  /// the block pointer in the octree lists.
  ~BEMFMA();

  /// Initialization function. It sets up the DoFHandler, FiniteElement,
  /// Mapping in the BEMFMA class. It also sets up some useful vector for,
  /// mixed boundary conditions and double nodes handling.
  void init_fma(const DoFHandler<dim-1,dim> &input_dh,
                const std::vector<std::set<types::global_dof_index> > &db_in,
                const TrilinosWrappers::MPI::Vector &input_sn,
                const Mapping<dim-1,dim> &input_mapping = StaticMappingQ1<dim-1, dim>::mapping,
                unsigned int quad_order=4, unsigned int sing_quad_order=5);

  /// Parameters declaration: we take the number of octree level and
  /// the level of approximation of the Kernels
  virtual void declare_parameters(ParameterHandler &prm);

  /// Parameter parsing from input file
  virtual void parse_parameters(ParameterHandler &prm);

  /// Method computing the parts of the
  /// BEM system matrices in which the
  /// integrals have to be performed
  /// directly. It also sets up the preconditioner. We
  /// don't consider any constraints in this function.
  /// We have set up the function to work with two levels
  /// of parallelism, TBB and MPI. We build up the direct
  /// contribution for all those block that doesn't satisfy
  /// the accuracy bounds for multipole (and consequently)
  /// local expansions.
  void direct_integrals();

  /// Method computing the multipole
  /// expansion containing the integrals
  /// values for each bottom level block.
  ///  It is called once for each
  /// GMRES solved. This function set up
  /// the structure to build the actual
  /// multipoles associated with quadrature
  /// points. For this reason it has to be called just once.
  void multipole_integrals();

  /// [TODO] TO BE MOVED INSIDE multipole_matr_vect_products and made private
  /// Ascending phase of the FMA method.
  /// Multipole expansions are genarated
  ///  at the bottom level blocks, and then
  /// translated to their parent blocks up
  /// to the highest level. It is
  /// called once per matrix-vector multiplication. We need to
  /// ensure that the proper structure as been set through the
  /// function multipole_integral(). We have chosen to use only
  /// one level of parallelism inside this function: multithreading with TBB.
  void generate_multipole_expansions(const TrilinosWrappers::MPI::Vector &phi_values, const TrilinosWrappers::MPI::Vector &dphi_dn_values) const;

  /// Descending phase of the FMA method. Local
  /// Expansions are obtained for each block
  /// starting from the top level and down
  /// to the bottom ones, where they are used
  /// to approximete the values of the
  /// integrals, i.e. the BEM matrix-vector
  /// product values.
  /// This function is called once per every matrix-vector product.
  /// We need to call this function after having generating the multipole expansions.
  /// Since this function takes a considerable amount of time we have chosen to parallelise
  /// it using multithreaded TBB on a single node and MPI to allow even a greater level of
  /// parallelism. In this case we don't need any communication because we have made sure that
  /// the multipole expansions are replicated on each node. Thus we can safely split the
  /// descending phase.
  void multipole_matr_vect_products(const TrilinosWrappers::MPI::Vector &phi_values, const TrilinosWrappers::MPI::Vector &dphi_dn_values,
                                    TrilinosWrappers::MPI::Vector &matrVectProdN,    TrilinosWrappers::MPI::Vector &matrVectProdD) const;


  // void compute_m2l_flags();

  /// this methods creates the adaptive
  /// octree partitioning of the domain,
  /// needed by the FMA algorithm.
  /// Since it takes a very small relative amount of time, and it is
  /// called just once in our program we have chosen not to parallelise it.
  void generate_octree_blocking();


  // In this function we have grouped some geometrical computation that
  // are useful for setting up the octree blocking. For example the maps
  // containing the surrounding elements for each element, or the component
  // associated with each vectorial degree of freedom.
  void compute_geometry_cache();

  /// Method for the assembling of the
  /// sparse preconitioning matrix for FMA.
  /// In this function we actually build up the final preconditioner taking great care of all the
  /// constraints, through a constraint matrix, we need to impose on our system. We need the vector alpha to properly build up the
  /// band of the BEM system we are going to invert (with ILU in this case).
  /// This function takes a lot of time to properly compute the preconditioner sparsity pattern and fill it.
  /// For this very reason we have taken great care in providing a double parallelisation strategy. Firstly we split
  /// among different nodes the computations with MPI and then we use TBB to ensure a parallelisation between different
  /// threads on multicore architectures.
  TrilinosWrappers::PreconditionILU &FMA_preconditioner(const TrilinosWrappers::MPI::Vector &alpha, ConstraintMatrix &c);

private:

  /// Three pointers to the problem parameters to be set equal to
  /// the ones in the calling problem through the init function
  SmartPointer<const DoFHandler<dim-1,dim> >         fma_dh;
  SmartPointer<const Mapping<dim-1,dim> >  fma_mapping;


  /// Truncation order for the multipole
  /// and local expansion series: it is
  /// read from the parameters input file.

  unsigned int trunc_order;

  /// Sparsity pattern for the
  /// initial preconditioning matrix
  /// (no constraints applied)

  TrilinosWrappers::SparsityPattern init_prec_sparsity_pattern;

  /// Sparsity pattern for the
  /// final preconditioning matrix
  /// (constraints applied)

  TrilinosWrappers::SparsityPattern final_prec_sparsity_pattern;

  /// Sparse Neumann matrix containing
  /// only direct integrals contributions
  TrilinosWrappers::SparseMatrix prec_neumann_matrix;
  //SparseMatrix<double> prec_neumann_matrix;

  /// Sparse Dirichlet matrix containing
  /// only direct integrals contributions

  TrilinosWrappers::SparseMatrix prec_dirichlet_matrix;
  //SparseMatrix<double> prec_dirichlet_matrix;

  /// Initial sparse preconditioning matrix (without constraints)

  TrilinosWrappers::SparseMatrix init_preconditioner;

  /// Final sparse preconditioning matrix (with constraints)

  TrilinosWrappers::SparseMatrix final_preconditioner;

  /// Structures where the Dirichlet
  /// matrix multipole
  /// integrals are stored: for each cell
  /// there are as many multipole
  /// expansions as the lowest level
  /// blocks in which element's quad
  /// points lie.

  mutable std::map <types::global_dof_index, std::map <cell_it, std::vector <MultipoleExpansion > > > elemMultipoleExpansionsKer1;

  /// Structures where the Neumann
  /// matrix multipole
  /// integrals are stored: for each cell
  /// there are as many multipole
  /// expansions as the lowest level
  /// blocks in which element's quad
  /// points lie.


  mutable std::map <types::global_dof_index, std::map <cell_it, std::vector <MultipoleExpansion > > > elemMultipoleExpansionsKer2;

  /// Vector storing the Dirichlet
  /// integrals multipole expansions
  /// for each block

  mutable std::vector <MultipoleExpansion > blockMultipoleExpansionsKer1;

  /// Vector storing the Neumann
  /// integrals multipole expansions
  /// for each block

  mutable std::vector <MultipoleExpansion > blockMultipoleExpansionsKer2;

  /// Vector storing the Dirichlet
  /// integrals local expansions
  /// for each block

  mutable std::vector <LocalExpansion > blockLocalExpansionsKer1;

  /// Vector storing the Neumann
  /// integrals local expansions
  /// for each block

  mutable std::vector <LocalExpansion > blockLocalExpansionsKer2;

  /// Associated Legendre functions class

  AssLegFunction assLegFunction;

  /// the preconditioner to be passed to bem_problem
  /// contributi diretti del multipolo, assembla la parte
  /// diretta in una vera matrice.
  TrilinosWrappers::PreconditionILU preconditioner;

  unsigned int quadrature_order;

  unsigned int singular_quadrature_order;

  /// mpi related variables

  MPI_Comm mpi_communicator;

  unsigned int n_mpi_processes;

  unsigned int this_mpi_process;

  /// maximum number of collocation points
  /// that can be contained by a childless block:
  /// if a block contains less nodes than this number
  /// it is not further refined

  unsigned int max_num_nodes_per_block;

  /// number of levels of the octree
  /// partitioning

  unsigned int num_octree_levels;

  /// Granularity for TBB cycles. Simple parallel_for at the moment.
  unsigned int tbb_granularity;

  /// here are declared dome structures which
  /// will be created in the framework of the
  /// octree partitioning of the mesh, and
  /// will be used in the FMA
  /// a map associating each DoF with the cells
  /// it belongs to

  std::map<types::global_dof_index, std::vector<cell_it> > dof_to_elems;

  /// a map associating each gradient DoF
  /// with the cells it belongs to

  std::map<types::global_dof_index, std::vector<cell_it> > gradient_dof_to_elems;

  /// a vector associating each gradient DoF
  /// with the component it represents

  std::vector<unsigned int > gradient_dof_components;

  /// a map associating each DoF to the
  /// block it belongs to
  /// for each level

  std::map<types::global_dof_index, std::vector<types::global_dof_index> > dof_to_block;

  /// a map associating each quad point to the
  /// block it belongs to for
  /// each level

  std::map<cell_it, std::vector<std::vector <types::global_dof_index> > > quad_point_to_block;

  /// a map associating each cell with a std::set
  /// containing the surrounding
  /// cells

  std::map <cell_it, std::set <cell_it> > elem_to_surr_elems;

  /// a vector to store all OctreeBlocks
  /// in which the geometry is divided

  mutable std::vector<OctreeBlock<dim> *> blocks;

  /// the total blocks number

  types::global_dof_index num_blocks;

  /// the indices in the blocks vector, at which
  /// each of the levels start or end

  std::vector <types::global_dof_index> endLevel;
  std::vector <types::global_dof_index> startLevel;

  /// a list of the indices of all the childless
  /// blocks

  std::vector <types::global_dof_index> childlessList;

  /// a list of the number of parent blocks
  /// for each level
  std::vector <types::global_dof_index> numParent;

  /// a std::vector containing the list of
  /// parent blocks for each level

  std::vector <std::vector<types::global_dof_index> > parentList;

  /// a std::map of std::vectors containing the
  /// list of quadrature points

  std::map <cell_it, std::vector <Point <dim> > > quadPoints;

  /// a std::map of std::vectors containing the
  /// list of normals at quadrature points

  // std::map <cell_it, std::vector <Point <dim> > > quadNormals;
  std::map <cell_it, std::vector <Tensor <1, dim> > > quadNormals;

  /// a std::map of std::vectors containing the
  /// list of shape function values at
  /// quadrature points

  std::map <cell_it, std::vector <std::vector<double> > > quadShapeFunValues;

  /// a std::map of std::vectors containing the
  /// list of JxW values at
  /// quadrature points

  std::map <cell_it, std::vector <double > > quadJxW;

  /// a std::vector containing std::vectors with
  /// the IDs of blocks with at least one dof,
  /// for each level

  std::vector< std::vector<types::global_dof_index> > dofs_filled_blocks;

  /// a std::vector containing std::vectors with
  /// the IDs of blocks with at least one
  /// quad point, for each level

  std::vector< std::vector<types::global_dof_index> > quad_points_filled_blocks;

  ConditionalOStream pcout;


  /// TODO parsed quadrature?
  shared_ptr<Quadrature<dim-1> > quadrature;
  SmartPointer<const Vector<double> > dirichlet_nodes;
  /// This should be erased by the usage of the constraint matrix.
  const std::vector <std::set<types::global_dof_index> >   *double_nodes_set;

  std::vector<std::vector<types::global_dof_index> > m2l_flags;
  IndexSet this_cpu_set;

};

#endif
