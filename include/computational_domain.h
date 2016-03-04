
// We start with including a bunch
// of include files: they might be more than
// needed, we might want to check, some time
#ifndef computational_domain_h
#define computational_domain_h

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/matrix_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q_eulerian.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/base/types.h>

#include <mpi.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>


#include "../include/local_expansion.h"
#include "../include/multipole_expansion.h"
#include "../include/ass_leg_function.h"

#include <deal2lkit/parameter_acceptor.h>
#include <mpi.h>

using namespace dealii;
using namespace deal2lkit;

template <int dim>
class ComputationalDomain : public ParameterAcceptor
{
public:

  /// constructor: since this is the
  /// class containing all the geometry and
  /// the base instruments needed by all the
  /// other classes, it is created first and
  /// the constructor does not need
  /// arguments.
  /// For the same reason, most of the class
  /// attributes are public: we can leter
  /// make them public end introduce suitable
  /// Get and Set methods, if needed

  ComputationalDomain(MPI_Comm comm = MPI_COMM_WORLD);


  ~ComputationalDomain();

  /// method to declare the parameters
  /// to be read from the parameters file

  virtual void declare_parameters(ParameterHandler &prm);

  /// method to parse the needed parameters
  /// from the parameters file

  virtual void parse_parameters(ParameterHandler &prm);

  /// method to create initial mesh

  void create_initial_mesh();
  /// alternative method to read initial mesh
  /// from file

  void read_domain();

  /// method to refine the imported mesh
  /// according to the level requested in
  /// the parameters file

  void refine_and_resize(const unsigned int refinement_level);

  void conditional_refine_and_resize(const unsigned int refinement_level);

  void update_triangulation();
  /// Here are the members of the class:
  /// they are all public, as the upper level
  /// classes (bem_problem, bem_fma,
  /// free_surface) will all need to perform
  /// operations based on the greometry (and
  /// the tools to handle it) contained in
  /// this class

  /// here are some basic classes needed by
  /// the program: a triangulation, and the
  /// FiniteElement and DoFHandler classes.
  /// A second DoF handler and FiniteElement
  /// must be created in order to compute
  /// the solution gradients, which are
  /// vectorial functions


  void make_edges_conformal(const bool with_double_nodes = true);

  void compute_double_vertex_cache();
  //const unsigned int fe_degree;
  //const unsigned int mapping_degree;

  Triangulation<dim-1, dim>             tria;

  /// here we are just renaming the cell
  /// iterator



  // values to be imported from the
  // parameters file:

  /// number of refining cycles

  unsigned int n_cycles, mg_cycles;

  /// Strings identifying the input grid name and format
  std::string input_grid_name;
  std::string input_grid_format;

  /// the material ID numbers in the mesh
  /// input file, for the dirichlet_nodes
  std::vector<unsigned int> dirichlet_boundary_ids;
  /// the material ID numbers in the mesh
  /// input file, for the neumann_nodes
  std::vector<unsigned int> neumann_boundary_ids;

  MPI_Comm mpi_communicator;

  unsigned int n_mpi_processes;

  unsigned int this_mpi_process;

  // to deal with conformity on edges with double nodes
  std::vector<bool> vertex_on_boundary;
  std::vector<std::vector<unsigned int> > double_vertex_vector;
  std::map<unsigned int, std::vector<typename Triangulation<dim-1,dim>::active_cell_iterator> > vert_to_elems;
  std::set<typename Triangulation<dim-1,dim>::active_cell_iterator> edge_cells;
  Manifold<dim-1, dim> *manifold;
  ConditionalOStream pcout;


};

#endif
