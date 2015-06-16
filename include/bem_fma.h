
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
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>

#include <mpi.h>
#include "parameter_acceptor.h"
#include "utilities.h"

namespace Operator
{
  template <int dim> class MinBEMOperator;
}

using namespace dealii;

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
                const std::vector<std::set<unsigned int> > &db_in,
                const TrilinosWrappers::MPI::Vector &input_sn,
                const Mapping<dim-1,dim> &input_mapping = StaticMappingQ1<dim-1, dim>::mapping);

  /// Parameters declaration: we take the number of octree level and
  /// the level of approximation of the Kernels
  virtual void declare_parameters(ParameterHandler &prm);

  /// Parameter parsing from input file
  virtual void parse_parameters(ParameterHandler &prm);

  /// Method computing the parts of the
  /// BEM system matrices in which the
  /// integrals have to be performed
  /// directly. It also sets up the preconditioner.
  void direct_integrals();

  /// Method computing the multipole
  /// expansion containing the integrals
  /// values for each bottom level block.
  ///  It is called once for each
  /// GMRES solved.

  void multipole_integrals();

  /// [TODO] TO BE MOVED INSIDE multipole_matr_vect_products and made private
  /// Ascending phase of the FMA method.
  /// Multipole expansions are genarated
  ///  at the bottom level blocks, and then
  /// translated to their parent blocks up
  /// to the highest level. It is
  /// called once per GMRES iteration.

  void generate_multipole_expansions(const TrilinosWrappers::MPI::Vector &phi_values, const TrilinosWrappers::MPI::Vector &dphi_dn_values) const;

  /// Descending phase of the FMA method. Local
  /// Expansions are obtained for each block
  /// starting from the top level and down
  /// to the bottom ones, where they are used
  /// to approximete the values of the
  /// integrals, i.e. the BEM matrix-vector
  /// product values

  void multipole_matr_vect_products(const TrilinosWrappers::MPI::Vector &phi_values, const TrilinosWrappers::MPI::Vector &dphi_dn_values,
                                    TrilinosWrappers::MPI::Vector &matrVectProdN,    TrilinosWrappers::MPI::Vector &matrVectProdD) const;


  void compute_m2l_flags();

  /// this methods creates the adaptive
  /// octree partitioning of the domain,
  /// needed by the FMA algorithm

  void generate_octree_blocking();

  void compute_geometry_cache();

  /// Method for the assembling of the
  /// sparse preconitioning matrix for FMA

  TrilinosWrappers::PreconditionILU &FMA_preconditioner(const TrilinosWrappers::MPI::Vector &alpha, ConstraintMatrix &c);

private:

  /// Three pointers to the problem parameters to be set equal to
  /// the ones in the calling problem through the init function
  SmartPointer<const DoFHandler<dim-1,dim> >         fma_dh;
  SmartPointer<const FiniteElement<dim-1,dim> >  fma_fe;
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

  mutable std::map <unsigned int, std::map <cell_it, std::vector <MultipoleExpansion > > > elemMultipoleExpansionsKer1;

  /// Structures where the Neumann
  /// matrix multipole
  /// integrals are stored: for each cell
  /// there are as many multipole
  /// expansions as the lowest level
  /// blocks in which element's quad
  /// points lie.


  mutable std::map <unsigned int, std::map <cell_it, std::vector <MultipoleExpansion > > > elemMultipoleExpansionsKer2;

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

  unsigned int singular_quadrature_order;

  /// mpi related variables

  MPI_Comm mpi_communicator;

  unsigned int n_mpi_processes;

  unsigned int this_mpi_process;


  /// number of levels of the octree
  /// partitioning

  unsigned int num_octree_levels;

  /// here are declared dome structures which
  /// will be created in the framework of the
  /// octree partitioning of the mesh, and
  /// will be used in the FMA
  /// a map associating each DoF with the cells
  /// it belongs to

  std::map<unsigned int, std::vector<cell_it> > dof_to_elems;

  /// a map associating each gradient DoF
  /// with the cells it belongs to

  std::map<unsigned int, std::vector<cell_it> > gradient_dof_to_elems;

  /// a vector associating each gradient DoF
  /// with the component it represents

  std::vector<unsigned int > gradient_dof_components;

  /// a map associating each DoF to the
  /// block it belongs to
  /// for each level

  std::map<unsigned int, std::vector<unsigned int> > dof_to_block;

  /// a map associating each quad point to the
  /// block it belongs to for
  /// each level

  std::map<cell_it, std::vector<std::vector <unsigned int> > > quad_point_to_block;

  /// a map associating each cell with a std::set
  /// containing the surrounding
  /// cells

  std::map <cell_it, std::set <cell_it> > elem_to_surr_elems;

  /// a vector to store all OctreeBlocks
  /// in which the geometry is divided

  mutable std::vector<OctreeBlock<dim> *> blocks;

  /// the total blocks number

  unsigned int num_blocks;

  /// the indices in the blocks vector, at which
  /// each of the levels start or end

  std::vector <unsigned int> endLevel;
  std::vector <unsigned int> startLevel;

  /// a list of the indices of all the childless
  /// blocks

  std::vector <unsigned int> childlessList;

  /// a list of the number of parent blocks
  /// for each level
  std::vector <unsigned int> numParent;

  /// a std::vector containing the list of
  /// parent blocks for each level

  std::vector <std::vector<unsigned int> > parentList;

  /// a std::map of std::vectors containing the
  /// list of quadrature points

  std::map <cell_it, std::vector <Point <dim> > > quadPoints;

  /// a std::map of std::vectors containing the
  /// list of normals at quadrature points

  std::map <cell_it, std::vector <Point <dim> > > quadNormals;

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

  std::vector< std::vector<unsigned int> > dofs_filled_blocks;

  /// a std::vector containing std::vectors with
  /// the IDs of blocks with at least one
  /// quad point, for each level

  std::vector< std::vector<unsigned int> > quad_points_filled_blocks;

  ConditionalOStream pcout;


  /// TODO parsed quadrature?
  shared_ptr<Quadrature<dim-1> > quadrature;
  SmartPointer<const TrilinosWrappers::MPI::Vector > dirichlet_nodes;
  /// This should be erased by the usage of the constraint matrix.
  const std::vector <std::set<unsigned int> >   *double_nodes_set;

  std::vector<std::vector<unsigned int> > m2l_flags;

};

#endif
