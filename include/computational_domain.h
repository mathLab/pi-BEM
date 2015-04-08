//----------------------------  step-34.cc  ---------------------------
//    $Id: step-34.cc 18734 2009-04-25 13:36:48Z heltai $
//    Version: $Name$ 
//
//    Copyright (C) 2009, 2010 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//    Authors: Luca Heltai, Cataldo Manigrasso
//
//----------------------------  step-34.cc  ---------------------------

				 // We start with including a bunch
				 // of include files: they might be more than
				 // needed, we might want to check, some time
#ifndef computational_domain_h
#define computational_domain_h

#include<deal.II/base/smartpointer.h>
#include<deal.II/base/convergence_table.h>
#include<deal.II/base/quadrature_lib.h>
#include<deal.II/base/quadrature_selector.h>
#include<deal.II/base/parsed_function.h>
#include<deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>

#include<deal.II/lac/full_matrix.h>
#include<deal.II/lac/sparse_matrix.h>
#include<deal.II/lac/matrix_lib.h>
#include<deal.II/lac/vector.h>
#include<deal.II/lac/solver_control.h>
#include<deal.II/lac/solver_gmres.h>
#include<deal.II/lac/precondition.h>
#include<deal.II/lac/compressed_sparsity_pattern.h>
#include<deal.II/lac/sparse_direct.h>

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
#include<deal.II/fe/mapping_q_eulerian.h>
#include<deal.II/fe/mapping_q.h>

#include<deal.II/numerics/data_out.h>
#include<deal.II/numerics/vector_tools.h>
#include<deal.II/numerics/solution_transfer.h>

#include <mpi.h>

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

#include <mpi.h>

using namespace dealii;

template <int dim>
class ComputationalDomain 
{
  public:
  
                                      // constructor: since this is the
				      // class containing all the geometry and
				      // the base instruments needed by all the
				      // other classes, it is created first and
				      // the constructor does not need
				      // arguments.
				      // For the same reason, most of the class
				      // attributes are public: we can leter
				      // make them public end introduce suitable
				      // Get and Set methods, if needed
    
    ComputationalDomain(const unsigned int fe_degree = 1,
	                const unsigned int mapping_degree = 1);
    

    ~ComputationalDomain();

                                      // method to declare the parameters
				      // to be read from the parameters file
    
    void declare_parameters(ParameterHandler &prm);
    
                                      // method to parse the needed parameters
				      // from the parameters file
  
    void parse_parameters(ParameterHandler &prm);
    
                                      // method to create initial mesh

    void create_initial_mesh();    
                                      // alternative method to read initial mesh
                                      // from file

    void read_domain();    

                                      // method to refine the imported mesh
				      // according to the level requested in
				      // the parameters file

    void refine_and_resize(const unsigned int refinement_level);
    
                                      // in the imported mesh, the nodes on the
				      // domain edges are doubled: this routine
				      // creates a std::vector of std::set which
				      // allows to relate each node to their
				      // double(s)
    
    void generate_double_nodes_set();

                                      // this method is needed to
				      // separate Dirichlet nodes from
				      // Neumann nodes.
 
    void compute_phi_nodes();
                                      // computation of normals at collocation points (dofs)
    void compute_normals();
    
                                      // this methods creates the adaptive 
				      // octree partitioning of the domain,
				      // needed by the FMA algorithm 
    
    void generate_octree_blocking(); 
				     
                                      // Here are the members of the class:
				      // they are all public, as the upper level
				      // classes (bem_problem, bem_fma,
				      // free_surface) will all need to perform
				      // operations based on the greometry (and
				      // the tools to handle it) contained in
				      // this class
				      
                                     // here are some basic classes needed by
				     // the program: a triangulation, and the
				     // FiniteElement and DoFHandler classes. 
                                     // A second DoF handler and FiniteElement
                                     // must be created in order to compute
                                     // the solution gradients, which are
                                     // vectorial functions

    //const unsigned int fe_degree;
    //const unsigned int mapping_degree;
				     
    Triangulation<dim-1, dim>             tria;
    FE_Q<dim-1,dim>                       fe;
    DoFHandler<dim-1,dim>                 dh;
    FESystem<dim-1,dim>                   gradient_fe;
    DoFHandler<dim-1,dim>                 gradient_dh;

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

                                     // An Eulerian Mapping is created to deal
				     // with the free surface and boat mesh
				     // deformation

    Vector<double> map_points;
    //MappingQ<dim-1, dim>	mapping;
    MappingQEulerian<dim-1, Vector<double>, dim> * mapping;

                                     // here we are just renaming the cell
				     // iterator
				     
    typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;

                                     // the following vectors are needed to
				     // treat Dirichlet and Neumann nodes
				     // differently. Each component of the
				     // first one is null if it corresponds
				     // to a Dirichlet node, and zero if
				     // it corresponds to a Neumann node.
				     // The second vector has instead null
				     // entries for Dirichlet nodes, and ones
				     // for Neumann nodes 
    
    Vector<double> surface_nodes;
    Vector<double> other_nodes;

                                     // values to be imported from the
				     // parameters file:
				     
				     // number of refining cycles
				      
    unsigned int n_cycles;


 				     // the number of standard quadrature points
				     // and singular kernel quadrature to be
				     // used 				      
        
    std_cxx1x::shared_ptr<Quadrature<dim-1> > quadrature;
    unsigned int singular_quadrature_order;

 				     // the material ID numbers in the mesh
				     // input file, for the free surface cells
				     // and wall boundary (boat) cells 				      

    unsigned int free_sur_ID1;
    unsigned int free_sur_ID2;
    unsigned int free_sur_ID3;
    unsigned int wall_sur_ID1;
    unsigned int wall_sur_ID2;
    unsigned int wall_sur_ID3;
				     
				     // number of levels of the octree
				     // partitioning
				      
    unsigned int num_octree_levels;
    
                                     // here are declared dome structures which
				     // will be created in the framework of the
				     // octree partitioning of the mesh, and
				     // will be used in the FMA
				     
				     // a map associating each DoF with the cells
				     // it belongs to

    std::map<unsigned int, std::vector<cell_it> > dof_to_elems;

				     // a map associating each gradient DoF
				     // with the cells it belongs to

    std::map<unsigned int, std::vector<cell_it> > gradient_dof_to_elems;

				     // a vector associating each gradient DoF
				     // with the component it represents

    std::vector<unsigned int > gradient_dof_components;
    
				     // a map associating each DoF to the
				     // block it belongs to
				     // for each level

    std::map<unsigned int, std::vector<unsigned int> > dof_to_block;

				     // a map associating each quad point to the
				     // block it belongs to for
				     // each level

    std::map<cell_it, std::vector<std::vector <unsigned int> > > quad_point_to_block;
    
				     // a map associating each cell with a std::set
				     // containing the surrounding
				     // cells
				     
    std::map <cell_it, std::set <cell_it> > elem_to_surr_elems; 

                                     // a vector to store all OctreeBlocks
				     // in which the geometry is divided 
				     
    mutable std::vector<OctreeBlock<dim> *> blocks;
    
                                     // the total blocks number
    
    unsigned int num_blocks;
    
                                     // the indices in the blocks vector, at which
				     // each of the levels start or end 
    
    std::vector <unsigned int> endLevel;
    std::vector <unsigned int> startLevel;
    
                                     // a list of the indices of all the childless
				     // blocks
    
    std::vector <unsigned int> childlessList;
    
                                     // a list of the number of parent blocks
				     // for each level
    std::vector <unsigned int> numParent;	
                                     
				     // a std::vector containing the list of
				     // parent blocks for each level

    std::vector <std::vector<unsigned int> > parentList;
    
				     // a std::map of std::vectors containing the
				     // list of quadrature points    

    std::map <cell_it, std::vector <Point <dim> > > quadPoints;
    
				     // a std::map of std::vectors containing the
				     // list of normals at quadrature points    

    std::map <cell_it, std::vector <Point <dim> > > quadNormals;
    
				     // a std::map of std::vectors containing the
				     // list of shape function values at
				     // quadrature points    

    std::map <cell_it, std::vector <std::vector<double> > > quadShapeFunValues;
    
				     // a std::map of std::vectors containing the
				     // list of JxW values at
				     // quadrature points    

    std::map <cell_it, std::vector <double > > quadJxW;
    
                                     // a std::vector containing std::vectors with
				     // the IDs of blocks with at least one dof,
				     // for each level				     

    std::vector< std::vector<unsigned int> > dofs_filled_blocks;
 
                                     // a std::vector containing std::vectors with
				     // the IDs of blocks with at least one
				     // quad point, for each level

    std::vector< std::vector<unsigned int> > quad_points_filled_blocks;
    

                                     // a std::set containing nodes of
				     // free surface with doubles on the boat
				     // needed to calculate their displacement

    std::set<unsigned int> free_surf_and_boat_nodes;

                                     // a std::set containing nodes of
				     // the boat with doubles on the boat
				     // needed to calculate the boat
				     // surface smoothing

     std::set<unsigned int> boat_keel_nodes;

     MPI_Comm mpi_communicator;

     unsigned int n_mpi_processes;

     unsigned int this_mpi_process;

     ConditionalOStream pcout;

				     // vector of Point<dim> containing the node normals 
    
    std::vector<Point<dim> > node_normals;

    //______________________
    //just for check: to be removed
    
    std::map <unsigned int, std::map<cell_it,unsigned int> > integralCheck;
     
};

#endif
