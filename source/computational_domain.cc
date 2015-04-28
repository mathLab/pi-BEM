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

#define TOLL 0.001
#define MAXELEMENTSPERBLOCK 1

#include "../include/computational_domain.h"
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_reordering.h>

				 // @sect4{ComputationalDomain::ComputationalDomain and
				 // ComputationalDomain::read_parameters}
				 // The constructor initializes the
				 // variuous object in much the same
				 // way as done in the finite element
				 // programs such as step-4 or
				 // step-6. The only new ingredient
				 // here is the ParsedFunction object,
				 // which needs, at construction time,
				 // the specification of the number of
				 // components.
				 //
				 // For the exact solution the number
				 // of vector components is one, and
				 // no action is required since one is
				 // the default value for a
				 // ParsedFunction object. The wind,
				 // however, requires dim components
				 // to be specified. Notice that when
				 // declaring entries in a parameter
				 // file for the expression of the
				 // Functions::ParsedFunction, we need
				 // to specify the number of
				 // components explicitly, since the
				 // function
				 // Functions::ParsedFunction::declare_parameters
				 // is static, and has no knowledge of
				 // the number of components.
template <int dim>
ComputationalDomain<dim>::ComputationalDomain(const unsigned int fe_degree,
	                                      const unsigned int mapping_degree)
		:
		fe(fe_degree),
		dh(tria),
                gradient_fe(FE_Q<dim-1,dim>(fe_degree), dim),
                gradient_dh(tria),
		mapping(NULL),
                mpi_communicator (MPI_COMM_WORLD),
                n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
                this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator)),
                pcout(std::cout)

{
// Only output on first processor.
pcout.set_condition(this_mpi_process == 0);
}

template <int dim>
ComputationalDomain<dim>::~ComputationalDomain()
{
if (blocks.size() > 0)
   {
   for (unsigned int ii = 0; ii < num_blocks;  ii++)
        delete blocks[ii];
   }
 if(mapping != NULL)
   {
     delete mapping;
     mapping = NULL;
   }
}



template <int dim> 
void ComputationalDomain<dim>::declare_parameters (ParameterHandler &prm)
{
  prm.declare_entry("Number of cycles", "4",
		    Patterns::Integer()); 
    
  prm.enter_subsection("Quadrature rules");
  {
    prm.declare_entry("Quadrature type", "gauss", 
		      Patterns::Selection(QuadratureSelector<(dim-1)>::get_quadrature_names()));
    prm.declare_entry("Quadrature order", "4", Patterns::Integer());
    prm.declare_entry("Singular quadrature order", "5", Patterns::Integer());
  }
  prm.leave_subsection();
    

  
  prm.enter_subsection("Boundary Conditions ID Numbers");
  {    
   prm.declare_entry("Free Surface 1 ID", "110", Patterns::Integer());
   prm.declare_entry("Free Surface 2 ID", "110", Patterns::Integer());
   prm.declare_entry("Free Surface 3 ID", "110", Patterns::Integer());
   prm.declare_entry("Wall Surface 1 ID", "110", Patterns::Integer());
   prm.declare_entry("Wall Surface 2 ID", "110", Patterns::Integer());
   prm.declare_entry("Wall Surface 3 ID", "110", Patterns::Integer());
  }
  prm.leave_subsection();
  
  prm.enter_subsection("Octree Params");
  {    
   prm.declare_entry("Number of Octree Levels", "1", Patterns::Integer());
  }
  prm.leave_subsection();// to be moved

}

template <int dim> 
void ComputationalDomain<dim>::parse_parameters (ParameterHandler &prm)
{
  n_cycles = prm.get_integer("Number of cycles");                   
    
  prm.enter_subsection("Quadrature rules");
  {
    quadrature =
      std_cxx1x::shared_ptr<Quadrature<dim-1> >
      (new QuadratureSelector<dim-1> (prm.get("Quadrature type"),
				      prm.get_integer("Quadrature order")));
    singular_quadrature_order = prm.get_integer("Singular quadrature order");
  }
  prm.leave_subsection();
    
  prm.enter_subsection("Boundary Conditions ID Numbers");
  {    
   free_sur_ID1 = prm.get_integer("Free Surface 1 ID");
   free_sur_ID2 = prm.get_integer("Free Surface 2 ID");
   free_sur_ID3 = prm.get_integer("Free Surface 3 ID");
   wall_sur_ID1 = prm.get_integer("Wall Surface 1 ID");
   wall_sur_ID2 = prm.get_integer("Wall Surface 2 ID");
   wall_sur_ID3 = prm.get_integer("Wall Surface 3 ID");
  }
  prm.leave_subsection();
  
  prm.enter_subsection("Octree Params");
  {
   num_octree_levels = prm.get_integer("Number of Octree Levels");
  }
  prm.leave_subsection();

}


				 // @sect4{ComputationalDomain::read_domain}
    
				 // A boundary element method
				 // triangulation is basically the
				 // same as a (dim-1) dimensional
				 // triangulation, with the difference
				 // that the vertices belong to a
				 // (dim) dimensional space.
				 //
				 // Some of the mesh formats supported
				 // in deal.II use by default three
				 // dimensional points to describe
				 // meshes. These are the formats
				 // which are compatible with the
				 // boundary element method
				 // capabilities of deal.II. In
				 // particular we can use either UCD
				 // or GMSH formats. In both cases, we
				 // have to be particularly careful
				 // with the orientation of the mesh,
				 // because, unlike in the standard
				 // finite element case, no reordering
				 // or compatibility check is
				 // performed here.  All meshes are
				 // considered as oriented, because
				 // they are embedded in a higher
				 // dimensional space. (See the
				 // documentation of the GridIn and of
				 // the Triangulation for further
				 // details on orientation of cells in
				 // a triangulation.) In our case, the
				 // normals to the mesh are external
				 // to both the circle in 2d or the
				 // sphere in 3d.
				 //
				 // The other detail that is required
				 // for appropriate refinement of the
				 // boundary element mesh, is an
				 // accurate description of the
				 // manifold that the mesh is
				 // approximating. We already saw this
				 // several times for the boundary of
				 // standard finite element meshes
				 // (for example in step-5 and
				 // step-6), and here the principle
				 // and usage is the same, except that
				 // the HyperBallBoundary class takes
				 // an additional template parameter
				 // that specifies the embedding space
				 // dimension. The function object
				 // still has to be static to live at
				 // least as long as the triangulation
				 // object to which it is attached.
        
 template <int dim>
  void ComputationalDomain<dim>::read_domain()
  {


    std::ifstream in;
    switch (dim)
      {
      case 2:
        in.open ("coarse_circle.inp");
        break;

      case 3:
        in.open ("coarse_cube_double_nodes.inp");
        break;

      default:
        Assert (false, ExcNotImplemented());
      }

    GridIn<dim-1, dim> gi;
    gi.attach_triangulation (tria);
    gi.read_ucd (in);
    dh.distribute_dofs(fe);


  }


template <int dim>
void ComputationalDomain<dim>::create_initial_mesh()
{

	  std::vector<Point<dim> > vertices;
	  std::vector<CellData<dim-1> > cells;
	  SubCellData subcelldata;

	  vertices.resize(34);

	  vertices[0](0) = 0;  vertices[0](1) = 0;  vertices[0](2) = 1;
	  vertices[1](0) = 0.57735;  vertices[1](1) = 0.57735;  vertices[1](2) = 0.57735;
	  vertices[2](0) = -0.57735;  vertices[2](1) = 0.57735;  vertices[2](2) = 0.57735;
	  vertices[3](0) = -0.57735;  vertices[3](1) = -0.57735;  vertices[3](2) = 0.57735;
	  vertices[4](0) = 0.57735;  vertices[4](1) = -0.57735;  vertices[4](2) = 0.57735;
	  vertices[5](0) = 0.70711;  vertices[5](1) = 0;  vertices[5](2) = 0.70711;
	  vertices[6](0) = 0;  vertices[6](1) = 0.70711;  vertices[6](2) = 0.70711;
	  vertices[7](0) = -0.70711;  vertices[7](1) = 0;  vertices[7](2) = 0.70711;
	  vertices[8](0) = 0;  vertices[8](1) = -0.70711;  vertices[8](2) = 0.70711;
	  vertices[9](0) = 1;  vertices[9](1) = 0;  vertices[9](2) = 0;
	  vertices[10](0) = 0; vertices[10](1) = 1; vertices[10](2) = 0;
	  vertices[11](0) = -1; vertices[11](1) = 0; vertices[11](2) = 0;
	  vertices[12](0) = 0; vertices[12](1) = -1; vertices[12](2) = 0;
	  vertices[13](0) = 0.70711; vertices[13](1) = 0.70711; vertices[13](2) = 0;
	  vertices[14](0) = -0.70711; vertices[14](1) = 0.70711; vertices[14](2) = 0;
	  vertices[15](0) = -0.70711; vertices[15](1) = -0.70711; vertices[15](2) = 0;
	  vertices[16](0) = 0.70711; vertices[16](1) = -0.70711; vertices[16](2) = 0;
	  vertices[17](0) = 0;  vertices[17](1) = 0;  vertices[17](2) = -1;
	  vertices[18](0) = 0.57735;  vertices[18](1) = 0.57735;  vertices[18](2) = -0.57735;
	  vertices[19](0) = -0.57735;  vertices[19](1) = 0.57735;  vertices[19](2) = -0.57735;
	  vertices[20](0) = -0.57735;  vertices[20](1) = -0.57735;  vertices[20](2) = -0.57735;
	  vertices[21](0) = 0.57735;  vertices[21](1) = -0.57735;  vertices[21](2) = -0.57735;
	  vertices[22](0) = 0.70711;  vertices[22](1) = 0;  vertices[22](2) = -0.70711;
	  vertices[23](0) = 0;  vertices[23](1) = 0.70711;  vertices[23](2) = -0.70711;
	  vertices[24](0) = -0.70711;  vertices[24](1) = 0;  vertices[24](2) = -0.70711;
	  vertices[25](0) = 0;  vertices[25](1) = -0.70711;  vertices[25](2) = -0.70711;
	  vertices[26](0) = 1;  vertices[26](1) = 0;  vertices[26](2) = 0;
	  vertices[27](0) = 0; vertices[27](1) = 1; vertices[27](2) = 0;
	  vertices[28](0) = -1; vertices[28](1) = 0; vertices[28](2) = 0;
	  vertices[29](0) = 0; vertices[29](1) = -1; vertices[29](2) = 0;
	  vertices[30](0) = 0.70711; vertices[30](1) = 0.70711; vertices[30](2) = 0;
	  vertices[31](0) = -0.70711; vertices[31](1) = 0.70711; vertices[31](2) = 0;
	  vertices[32](0) = -0.70711; vertices[32](1) = -0.70711; vertices[32](2) = 0;
	  vertices[33](0) = 0.70711; vertices[33](1) = -0.70711; vertices[33](2) = 0;


	  cells.resize(24);



	  cells[0].vertices[0]=0; cells[0].vertices[1]=5; cells[0].vertices[2]=1; cells[0].vertices[3]=6;
	  cells[1].vertices[0]=7; cells[1].vertices[1]=0; cells[1].vertices[2]=6; cells[1].vertices[3]=2;
	  cells[2].vertices[0]=7; cells[2].vertices[1]=3; cells[2].vertices[2]=8; cells[2].vertices[3]=0;
	  cells[3].vertices[0]=4; cells[3].vertices[1]=5; cells[3].vertices[2]=0; cells[3].vertices[3]=8;
	  cells[4].vertices[0]=9; cells[4].vertices[1]=13; cells[4].vertices[2]=1; cells[4].vertices[3]=5;
	  cells[5].vertices[0]=13; cells[5].vertices[1]=10; cells[5].vertices[2]=6; cells[5].vertices[3]=1;
	  cells[6].vertices[0]=10; cells[6].vertices[1]=14; cells[6].vertices[2]=2; cells[6].vertices[3]=6;
	  cells[7].vertices[0]=14; cells[7].vertices[1]=11; cells[7].vertices[2]=7; cells[7].vertices[3]=2;
	  cells[8].vertices[0]=11; cells[8].vertices[1]=15; cells[8].vertices[2]=3; cells[8].vertices[3]=7;
	  cells[9].vertices[0]=15; cells[9].vertices[1]=12; cells[9].vertices[2]=8; cells[9].vertices[3]=3;
	  cells[10].vertices[0]=12; cells[10].vertices[1]=16; cells[10].vertices[2]=4; cells[10].vertices[3]=8;
	  cells[11].vertices[0]=16; cells[11].vertices[1]=9; cells[11].vertices[2]=5; cells[11].vertices[3]=4;

	  cells[12].vertices[0]=29; cells[12].vertices[1]=32; cells[12].vertices[2]=20; cells[12].vertices[3]=25;
	  cells[13].vertices[0]=33; cells[13].vertices[1]=29; cells[13].vertices[2]=25; cells[13].vertices[3]=21;
	  cells[14].vertices[0]=32; cells[14].vertices[1]=28; cells[14].vertices[2]=24; cells[14].vertices[3]=20;
	  cells[15].vertices[0]=28; cells[15].vertices[1]=31; cells[15].vertices[2]=19; cells[15].vertices[3]=24;
	  cells[16].vertices[0]=31; cells[16].vertices[1]=27; cells[16].vertices[2]=23; cells[16].vertices[3]=19;
	  cells[17].vertices[0]=27; cells[17].vertices[1]=30; cells[17].vertices[2]=18; cells[17].vertices[3]=23;
	  cells[18].vertices[0]=30; cells[18].vertices[1]=26; cells[18].vertices[2]=22; cells[18].vertices[3]=18;
	  cells[19].vertices[0]=26; cells[19].vertices[1]=33; cells[19].vertices[2]=21; cells[19].vertices[3]=22;
	  cells[20].vertices[0]=25; cells[20].vertices[1]=20; cells[20].vertices[2]=24; cells[20].vertices[3]=17;
	  cells[21].vertices[0]=22; cells[21].vertices[1]=21; cells[21].vertices[2]=25; cells[21].vertices[3]=17;
	  cells[22].vertices[0]=24; cells[22].vertices[1]=19; cells[22].vertices[2]=23; cells[22].vertices[3]=17;
	  cells[23].vertices[0]=23; cells[23].vertices[1]=18; cells[23].vertices[2]=22; cells[23].vertices[3]=17;

	  
	  cells[0].material_id=1;
	  cells[1].material_id=1;
	  cells[2].material_id=1;
	  cells[3].material_id=1;
	  cells[4].material_id=1;
	  cells[5].material_id=1;
	  cells[6].material_id=1;
	  cells[7].material_id=1;
	  cells[8].material_id=1;
	  cells[9].material_id=1;
	  cells[10].material_id=1;
	  cells[11].material_id=1;
	  cells[12].material_id=0;
	  cells[13].material_id=0;
	  cells[14].material_id=0;
	  cells[15].material_id=0;
	  cells[16].material_id=0;
	  cells[17].material_id=0;
	  cells[18].material_id=0;
	  cells[19].material_id=0;
	  cells[20].material_id=0;
	  cells[21].material_id=0;
	  cells[22].material_id=0;
	  cells[23].material_id=0;

	  
	  GridTools::delete_unused_vertices(vertices,cells,subcelldata);
	  GridReordering<dim-1,dim>::reorder_cells(cells);
          tria.create_triangulation_compatibility(vertices,cells,subcelldata);

          static const Point<dim> center = Point<dim>();
          static const HyperBallBoundary<dim-1, dim> boundary(center,1.);
          tria.set_boundary(1, boundary);
          tria.set_boundary(0, boundary);


}

				 // @sect4{ComputationalDomain::refine_and_resize}

				 // This function globally refines the
				 // mesh, distributes degrees of
				 // freedom, and resizes matrices and
				 // vectors.

template <int dim>
void ComputationalDomain<dim>::refine_and_resize(const unsigned int refinement_level)
{
  pcout<<"Refining and resizing mesh as required"<<std::endl;    

  tria.refine_global(refinement_level);
  GridTools::partition_triangulation(n_mpi_processes, tria);
  std::string filename0 = ( "meshResult.inp" );
  std::ofstream logfile0(filename0.c_str());
  GridOut grid_out0;
  grid_out0.write_ucd(tria, logfile0);   


  dh.distribute_dofs(fe);
  DoFRenumbering::subdomain_wise(dh);
  gradient_dh.distribute_dofs(gradient_fe);
  DoFRenumbering::subdomain_wise(dh);
  map_points.reinit(gradient_dh.n_dofs(),false);

  const types::global_dof_index n_local_dofs = DoFTools::count_dofs_with_subdomain_association (dh, this_mpi_process);
  pcout<<"Proc: "<<this_mpi_process<<"  Dofs: "<<n_local_dofs<<std::endl;

  if(mapping == NULL)
    mapping = new MappingQEulerian<dim-1, Vector<double>, dim>
	      (gradient_fe.degree, map_points, gradient_dh);

  //Vector<double> map_points_old(map_points);
  //soltrans.refine_interpolate(map_points_old, map_points);
  
  generate_double_nodes_set();
    

  pcout<<"...done refining and resizing mesh"<<std::endl;  
}    


template <int dim>
void ComputationalDomain<dim>::generate_double_nodes_set()
{ 
pcout<<"Generating double nodes set..."<<std::endl;

				 // @sect5{ComputationalDomain::generate_double_nodes_set}

				 // The following is the function
				 // which creates a set containing
				 // the double nodes.
				 

  double tol = 1e-8;
  double_nodes_set.clear();
  double_nodes_set.resize(dh.n_dofs());
  std::vector<Point<dim> > support_points(dh.n_dofs());

  DoFTools::map_dofs_to_support_points<dim-1, dim>( *mapping,
						    dh, support_points);
 
  for (unsigned int i=0; i<dh.n_dofs(); ++i)  
      {
      for (unsigned int j=0; j<dh.n_dofs(); ++j) 
	  {
	  //pcout<<"i "<<i<<" ("<<support_points[i]<<")  j "<<j<<" ("<<support_points[j]<<")  distance "<<support_points[i].distance(support_points[j])<<std::endl;
	  if (support_points[i].distance(support_points[j]) < tol)
	     {
	     double_nodes_set[i].insert(j);
	     //pcout<<"i "<<i<<" double "<<double_nodes_set[i].count(j)<<std::endl;
	     }
	  }
	
      }
                                 // for the gradient dofs finding coupled
                                 // dofs is a little bit difficult, as the
                                 // gradient is a vectorial function: usually
                                 // in such case the dofs are numbered
                                 // so that for each support point dim
                                 // consecutive dofs represent each component
                                 // of the vector field: in this case
                                 // (and only in this case) the following
                                 // piece of code works
   
  gradient_double_nodes_set.clear();
  gradient_double_nodes_set.resize(gradient_dh.n_dofs());
  std::vector<Point<dim> > gradient_support_points(gradient_dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *mapping,
						    gradient_dh, gradient_support_points);
  for (unsigned int k = 0; k < dim; k++)
      {
      for(unsigned int i=k; i<gradient_dh.n_dofs(); i=i+dim)  
	 {
	 for(unsigned int j=k; j<gradient_dh.n_dofs(); j=j+dim) 
	    {
             //pcout<<"i "<<i<<" ("<<support_points[i]<<")  j "<<j<<" ("<<support_points[j]<<")  distance "<<support_points[i].distance(support_points[j])<<std::endl;
            if (gradient_support_points[i].distance(gradient_support_points[j]) < tol )
               {
		gradient_double_nodes_set[i].insert(j);
		//pcout<<"i "<<i<<" double "<<j<<std::endl;
		}
	    }
	 }
      }
      
  cell_it
  gradient_cell = gradient_dh.begin_active(),
  gradient_endc = gradient_dh.end();
  
  cell_it
  cell = dh.begin_active(),
  endc = dh.end();
     
  std::vector<unsigned int> dofs(fe.dofs_per_cell);
  std::vector<unsigned int> gradient_dofs(gradient_fe.dofs_per_cell);
    // mappa che associa ad ogni dof le celle cui esso appartiene
  dof_to_elems.clear();
  
  // mappa che associa ad ogni gradient dof le celle cui esso appartiene
  gradient_dof_to_elems.clear();
  
  // vettore che associa ad ogni gradient dof la sua componente
  gradient_dof_components.clear();
  gradient_dof_components.resize(gradient_dh.n_dofs());
  
  // mappa che associa ad ogni cella un set contenente le celle circostanti
  elem_to_surr_elems.clear();
  
  // set che raccoglie i nodi della free surface che stanno sulla barca
  free_surf_and_boat_nodes.clear();
  
  // mappa raccoglie i nodi degli edges della barca
  boat_keel_nodes.clear();  
  
   
   for (; cell!=endc; ++cell,++gradient_cell)
    {
    Assert(cell->index() == gradient_cell->index(), ExcInternalError());
    
    cell->get_dof_indices(dofs);
    for(unsigned int j=0; j<fe.dofs_per_cell; ++j)
        {
	dof_to_elems[dofs[j]].push_back(cell);
	}
    gradient_cell->get_dof_indices(gradient_dofs);
    for(unsigned int j=0; j<gradient_fe.dofs_per_cell; ++j)
        {
	gradient_dof_to_elems[gradient_dofs[j]].push_back(gradient_cell);
	gradient_dof_components[gradient_dofs[j]] = gradient_fe.system_to_component_index(j).first;
	}
    }
  
  // qui viene creata la mappa dei elmenti che circondano ciascun elemento    
  for (cell = dh.begin_active(); cell != endc; ++cell)    
    {
    cell->get_dof_indices(dofs);
    for(unsigned int j=0; j<fe.dofs_per_cell; ++j)
        {
	std::set <unsigned int> duplicates = double_nodes_set[dofs[j]];
	for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
	    {
	    std::vector<cell_it>
	    dof_cell_list = dof_to_elems[*pos];
	    for (unsigned int k=0; k<dof_cell_list.size(); ++k)
                elem_to_surr_elems[cell].insert(dof_cell_list[k]);
	    }
	}
    }
  
  gradient_cell = gradient_dh.begin_active();
  cell = dh.begin_active();
  
  for (; cell!=endc,gradient_cell!=gradient_endc; ++cell,++gradient_cell)
      {
      Assert(cell->index() == gradient_cell->index(), ExcInternalError());
  
      if (cell->material_id() == wall_sur_ID1 ||
          cell->material_id() == wall_sur_ID2 ||
          cell->material_id() == wall_sur_ID3)
	 {
         // This is a free surface node.
         gradient_cell->get_dof_indices(gradient_dofs);
         for (unsigned int j=0; j<gradient_fe.dofs_per_cell; ++j)
             {
	     std::set <unsigned int> duplicates = gradient_double_nodes_set[gradient_dofs[j]];
	     for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
	         {
	         cell_it duplicate_cell = gradient_dof_to_elems[*pos][0];
		 if (duplicate_cell->material_id() == free_sur_ID1 ||
                     duplicate_cell->material_id() == free_sur_ID2 ||
                     duplicate_cell->material_id() == free_sur_ID3)
		    {
		    free_surf_and_boat_nodes.insert(gradient_dofs[j]);
		    }
	         
	         }
	     }
             //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"  surfNodes: "<<surface_nodes(dofs[i])<<"  otherNodes: "<<other_nodes(dofs[i])<<std::endl;
         }
      }    

  gradient_cell = gradient_dh.begin_active();
  cell = dh.begin_active();
  
  for (; cell!=endc,gradient_cell!=gradient_endc; ++cell,++gradient_cell)
      {
      Assert(cell->index() == gradient_cell->index(), ExcInternalError());
  
      if (gradient_cell->material_id() == wall_sur_ID1)
	 {
         gradient_cell->get_dof_indices(gradient_dofs);
         for (unsigned int j=0; j<gradient_fe.dofs_per_cell; ++j)
             {
	     std::set <unsigned int> duplicates = gradient_double_nodes_set[gradient_dofs[j]];
	     for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
	         {
	         cell_it duplicate_cell = gradient_dof_to_elems[*pos][0];
		 if (
                     duplicate_cell->material_id() == wall_sur_ID2 ||
                     duplicate_cell->material_id() == wall_sur_ID3)
		    {
		    boat_keel_nodes.insert(gradient_dofs[j]);
		    }
	         
	         }
	     }
             //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"  surfNodes: "<<surface_nodes(dofs[i])<<"  otherNodes: "<<other_nodes(dofs[i])<<std::endl;
         }
      }
  gradient_cell = gradient_dh.begin_active();
  cell = dh.begin_active();
  
  for (; cell!=endc,gradient_cell!=gradient_endc; ++cell,++gradient_cell)
      {
      Assert(cell->index() == gradient_cell->index(), ExcInternalError());

      if (gradient_cell->material_id() == wall_sur_ID2)
	 {
         // This is a free surface node.
         gradient_cell->get_dof_indices(gradient_dofs);
         for (unsigned int j=0; j<gradient_fe.dofs_per_cell; ++j)
             {
	     std::set <unsigned int> duplicates = gradient_double_nodes_set[gradient_dofs[j]];
	     for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
	         {
	         cell_it duplicate_cell = gradient_dof_to_elems[*pos][0];
		 if (duplicate_cell->material_id() == wall_sur_ID3)
		    {
		    boat_keel_nodes.insert(gradient_dofs[j]);
		    }
	         
	         }
	     }
             //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"  surfNodes: "<<surface_nodes(dofs[i])<<"  otherNodes: "<<other_nodes(dofs[i])<<std::endl;
         }
      }

      
pcout<<"...done"<<std::endl;
}


template <int dim>
void ComputationalDomain<dim>::compute_phi_nodes() {

  cell_it
  gradient_cell = gradient_dh.begin_active(),
  gradient_endc = gradient_dh.end();

  cell_it
  cell = dh.begin_active(),
  endc = dh.end();
     
  surface_nodes.reinit(dh.n_dofs());
  other_nodes.reinit(dh.n_dofs());
  other_nodes.add(1);
  std::vector<unsigned int> dofs(fe.dofs_per_cell);
  std::vector<unsigned int> gradient_dofs(gradient_fe.dofs_per_cell); 
  
  for(; cell != endc; ++cell) {
    if(cell->material_id() == free_sur_ID1 ||
       cell->material_id() == free_sur_ID2 ||
       cell->material_id() == free_sur_ID3) {
      // This is a free surface node.
      cell->get_dof_indices(dofs);
      for(unsigned int i=0; i<fe.dofs_per_cell; ++i) {
	surface_nodes(dofs[i]) = 1;
	other_nodes(dofs[i]) = 0;
        //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"  surfNodes: "<<surface_nodes(dofs[i])<<"  otherNodes: "<<other_nodes(dofs[i])<<std::endl;
      }
    }
    else
    {
      for(unsigned int i=0; i<fe.dofs_per_cell; ++i) {
        cell->get_dof_indices(dofs);
        //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"  surfNodes: "<<surface_nodes(dofs[i])<<"  otherNodes: "<<other_nodes(dofs[i])<<std::endl;
      }
    
    }
    
  }
    
  //for (unsigned int i=0; i<dh.n_dofs(); ++i)
  //    if (this_mpi_process == 1)
  //       pcout<<i<<" "<<surface_nodes(i)<<" "<<other_nodes(i)<<std::endl;

}

template <int dim>
void ComputationalDomain<dim>::generate_octree_blocking()
{

pcout<<"Generating octree blocking... "<<std::endl;

				 // @sect5{BEMProblem::generate_double_nodes_set}

				 // The following is the function
				 // which creates the octree blocking
				 // for the fast multipole algorithm


std::vector<Point<dim> > support_points(dh.n_dofs());
DoFTools::map_dofs_to_support_points<dim-1, dim>( *mapping,
                                                  dh, support_points);

FEValues<dim-1,dim> fe_v(*mapping,fe, *quadrature,
			 update_values |
			 update_cell_normal_vectors |
			 update_quadrature_points |
			 update_JxW_values);

double max_coor_value = 0;

for (unsigned int i=0; i < dh.n_dofs(); i++)
    {
    //for printout
    //pcout<<"Node "<<i<< "["<<support_points[i]<<"] "<<std::endl;
    for (unsigned int j=0; j < dim; j++)
        {
	max_coor_value = std::max(max_coor_value,std::abs(support_points[i](j)));
	}
    }
    
if (blocks.size() > 0)
   {
   for (unsigned int ii = 0; ii < num_blocks;  ii++)
        delete blocks[ii];
   }
    
unsigned int maxNumBlocks = num_octree_levels*tria.n_active_cells()*fe_v.n_quadrature_points;
//unsigned int maxNumBlocks = 0;
//for (unsigned int ii = 0; ii < num_octree_levels + 1;  ii++)
//	{
//	 maxNumBlocks += int(pow(8.,double(ii)));
//	}

blocks.clear();
blocks.reserve(maxNumBlocks);
blocks.resize(maxNumBlocks);

unsigned int blocksCount = 0;
startLevel.resize(num_octree_levels+1);
endLevel.resize(num_octree_levels+1);

//for (unsigned int j=0; j < num_octree_levels + 1; j++)
//     parentList[j].clear();
parentList.clear();
parentList.resize(num_octree_levels+1);
parentList[0].push_back(0);


childlessList.clear();
unsigned int numChildless = 0;
numParent.resize(num_octree_levels+1);

//qui di seguito vengono reinizializzate strutture utili al multipolo

// mappa che associa ad ogni dof un vettore con i blocchi cui essa appartiene per ogni livello
dof_to_block.clear();

// mappa che associa ad ogni quad point un vettore con i blocchi cui essa appartiene per ogni livello
quad_point_to_block.clear();

// vettore di vettori contenente per ogni livello, gli ids dei blocchi
// contenenti almeno un dof
dofs_filled_blocks.clear();

// vettore di vettori contenente per ogni livello, gli ids dei blocchi
// contenenti almeno un quad point
quad_points_filled_blocks.clear();

quadPoints.clear();
quadNormals.clear();
quadShapeFunValues.clear();
quadJxW.clear();

dofs_filled_blocks.resize(num_octree_levels+1);

quad_points_filled_blocks.resize(num_octree_levels+1);



for (unsigned int ii = 0; ii < num_octree_levels + 1 ;  ii++)
	{
	 numParent[ii] = 0;
	}



Point<dim> pMin;
for (int i=0; i<dim; i++)
   pMin(i) = -1.1*max_coor_value;

				 // delta e' il lato del kazzo di kubo...
double delta = 2.2*max_coor_value;

OctreeBlock<dim>* block = new OctreeBlock<dim>(0, 0, pMin, delta);

std::vector<unsigned int> local_dof_indices(fe.dofs_per_cell);
cell_it
cell = dh.begin_active(),
endc = dh.end();
for (cell = dh.begin_active(); cell != endc; ++cell)    
    {
    fe_v.reinit(cell);
    const unsigned int n_q_points = fe_v.n_quadrature_points;
    quadPoints[cell] = fe_v.get_quadrature_points();
    quadNormals[cell] = fe_v.get_normal_vectors();
    quadJxW[cell].resize(n_q_points);
    quadShapeFunValues[cell].resize(n_q_points);
    for(unsigned int q=0; q<n_q_points; ++q)
       {
       quadJxW[cell][q] = fe_v.JxW(q);
       for(unsigned int j=0; j<fe.dofs_per_cell; ++j)
          quadShapeFunValues[cell][q].push_back(fe_v.shape_value(j,q));
       }
    
    quad_point_to_block[cell].resize(n_q_points);		     
    for (unsigned int j=0; j<n_q_points; ++j)
        {
	block->AddQuadPoint(cell,j);
	quad_point_to_block[cell][j].push_back(0);
	}
    
    cell->get_dof_indices(local_dof_indices);
    for(unsigned int j=0; j<fe.dofs_per_cell; ++j)
        {
	dof_to_elems[local_dof_indices[j]].push_back(cell);
	}
    }
    
for (unsigned int ii = 0; ii < dh.n_dofs(); ii++)
    {
    block->AddNode(ii);
    dof_to_block[ii].push_back(0);
    }    


// just for output    
/*for (cell = dh.begin_active(); cell != endc; ++cell)    
    {
    std::set<cell_it> surr_elems = elem_to_surr_elems[cell];
    pcout<<std::endl<<"cell "<<cell<<"  surrounded by: "; 
    for (typename std::set<cell_it>::iterator pos = surr_elems.begin(); pos !=surr_elems.end(); pos++)
         pcout<<" "<<*pos;
    }*/
    
blocks[0] = block;
numParent[0] = 1;

//pcout<<"blocks[0].GetBlockChildrenNum() "<<blocks[0].GetBlockChildrenNum()<<std::endl;

/*pcout<<std::endl;
pcout<<blocks[0].GetPMin()<<std::endl;
pcout<<pMin<<std::endl;
pcout<<block.GetDelta()<<std::endl;
pcout<<block.GetBlockNodeList()[0]<<std::endl;
pcout<<block.GetBlockElementsList()[1]<<std::endl;
pcout<<delta<<std::endl;
pcout<<std::endl;//*/

unsigned int quadPointsInChildless = 0;
unsigned int nodesInChildless = 0;

for (unsigned int level = 1; level < num_octree_levels + 1;  level++)

    {
    unsigned int quadPointsCheck = quadPointsInChildless;
    unsigned int nodesCheck = nodesInChildless;
    delta /= 2.;

    for (unsigned int kk = 0; kk < numParent[level-1];  kk++)

        {
	unsigned int jj = parentList[level-1][kk];
        //pcout<<" level "<<level<<"     block "<<jj<<std::endl;
	OctreeBlock<dim> *parent = blocks[jj];
        //pcout<<"parent.GetBlockChildrenNum() "<<parent.GetBlockChildrenNum()<<std::endl;
        //pcout<<" Pmin "<<parent.GetPMin()(0)<<", "<<parent.GetPMin()(1)<<", "<<parent.GetPMin()(2)<<" "<<std::endl;
        //pcout<<" delta "<<parent.GetDelta()<<" "<<std::endl;

        pMin = parent->GetPMin();
	unsigned int num_children_per_block = int(pow((double)2,(double)dim));
        std::vector<OctreeBlock<dim> *> children(num_children_per_block);
	
	if (dim == 3)
	   {
           children[0] = new OctreeBlock<dim>(level, jj, pMin                              , delta);
           children[1] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta   ,0.,   0.), delta);
           children[2] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta,delta,   0.), delta);
           children[3] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(0.   ,delta,   0.), delta);
           children[4] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(0.   ,   0.,delta), delta);
           children[5] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta,   0.,delta), delta);
           children[6] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta,delta,delta), delta);
           children[7] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(   0.,delta,delta), delta);
           }

	if (dim == 2)
	   {
           children[0] = new OctreeBlock<dim>(level, jj, pMin                        , delta);
           children[1] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta   ,0.), delta);
           children[2] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta,delta), delta);
           children[3] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(0.   ,delta), delta);
           
           }

           std::map <cell_it, std::vector <unsigned int> > blockQuadPointsList =
	                                   parent->GetBlockQuadPointsList();
					   
           std::vector <unsigned int> blockNodeList = parent->GetBlockNodeList();
        
	if (dim == 3) 
	{
        for (unsigned int i = 0; i < blockNodeList.size(); i++)    
            {
	    Point <dim> node = support_points[blockNodeList[i]];
	    // assegnamento nodi del blocco padre ai blocchi figli
	    
	    if (node(2) <= parent->GetPMin()(2)+delta)
	       {
	       if (node(1) <= parent->GetPMin()(1)+delta)
		  {
		  if (node(0) <= parent->GetPMin()(0)+delta)
		     {
		     //pcout<<" Sono in 1 "<<std::endl;
		     children[0]->AddNode(blockNodeList[i]);
		     }
		  else
		     {
		     //pcout<<" Sono in 2 "<<std::endl;
		     children[1]->AddNode(blockNodeList[i]);
		     }
		  }
	       else
		  {
                  if (node(0) <= parent->GetPMin()(0)+delta)
		     {
		     //pcout<<" Sono in 4 "<<std::endl;
		     children[3]->AddNode(blockNodeList[i]);
		     }
		  else
		     {
		     //pcout<<" Sono in 3 "<<std::endl;
		     children[2]->AddNode(blockNodeList[i]);
		     }
		  }
	       }
            else
	       {
	       if (node(1) <= parent->GetPMin()(1)+delta)
		  {
	          if (node(0) <= parent->GetPMin()(0)+delta)
		     {
		     //pcout<<" Sono in 5 "<<std::endl;
		     children[4]->AddNode(blockNodeList[i]);
		     }
		  else
		     {
		     //pcout<<" Sono in 6 "<<std::endl;
		     children[5]->AddNode(blockNodeList[i]);
		     }
		  }
               else
		  {
		  if (node(0) <= parent->GetPMin()(0)+delta)
		     {
		     //pcout<<" Sono in 8 "<<std::endl;
		     children[7]->AddNode(blockNodeList[i]);
		     }
	          else
	             {
		     //pcout<<" Sono in 7 "<<std::endl;
		     children[6]->AddNode(blockNodeList[i]);
		     }
		  }
	       } //fine assegnazione nodi del padre ai blocchi figli
	       
	    } //fine loop nodi del blocco 
	
	typename std::map <cell_it, std::vector<unsigned int> >::iterator it;    
        for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)    
            {
	    for (unsigned int pp = 0; pp < (*it).second.size(); pp++)
	        {
	        Point<dim> quadPoint = quadPoints[(*it).first][(*it).second[pp]];
	        // assegnamento punti quadratura del blocco padre ai blocchi figli
	        if (quadPoint(2) <= parent->GetPMin()(2)+delta)
	           {
	           if (quadPoint(1) <= parent->GetPMin()(1)+delta)
		      {
		      if (quadPoint(0) <= parent->GetPMin()(0)+delta)
		         {
		         //pcout<<" Sono in 1 "<<std::endl;
		         children[0]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
		      else
		         {
		         //pcout<<" Sono in 2 "<<std::endl;
		         children[1]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
		      }
	           else
		      {
                      if (quadPoint(0) <= parent->GetPMin()(0)+delta)
		         {
		         //pcout<<" Sono in 4 "<<std::endl;
		         children[3]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
		      else
		         {
		         //pcout<<" Sono in 3 "<<std::endl;
		         children[2]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
		      }
	           }
                else
	           {
	           if (quadPoint(1) <= parent->GetPMin()(1)+delta)
		      {
	              if (quadPoint(0) <= parent->GetPMin()(0)+delta)
		         {
		         //pcout<<" Sono in 5 "<<std::endl;
		         children[4]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
		      else
		         {
		         //pcout<<" Sono in 6 "<<std::endl;
		         children[5]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
		      }
                   else
		      {
		      if (quadPoint(0) <= parent->GetPMin()(0)+delta)
		         {
		         //pcout<<" Sono in 8 "<<std::endl;
		         children[7]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
	              else
	                 {
		         //pcout<<" Sono in 7 "<<std::endl;
		         children[6]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
		      }
	           } //fine assegnazione punti quadratura del padre ai blocchi figli
	        }
	    } //fine loop punti quadratura del blocco

        for (unsigned int j=0; j < num_children_per_block; j++ )
	    {   
	    if (children[j]->GetBlockNodeList().size() + 
	        children[j]->GetBlockQuadPointsList().size()> 0)
	       {
	       blocksCount += 1;
	       blocks[blocksCount] = new OctreeBlock<dim>();
	       blocks[blocksCount]->CopyContent(children[j]);
	       delete children[j];
	       
	       parent->AddChild(blocksCount);
               std::map <cell_it, std::vector<unsigned int> >
	       blockQuadPointsList = blocks[blocksCount]->GetBlockQuadPointsList();
	       typename std::map <cell_it, std::vector<unsigned int> >::iterator it;    
               for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
	           {
		   cell_it cell = (*it).first;
		   for(unsigned int kk = 0; kk < (*it).second.size(); kk++)
		      {
		      quad_point_to_block[(*it).first][(*it).second[kk]].push_back(blocksCount);
		      }
		   }
               std::vector<unsigned int> blockNodesList = blocks[jj]->GetBlockNodeList();
	       for(unsigned int k = 0; k < blockNodesList.size(); k++)
	          dof_to_block[blockNodesList[k]].push_back(jj);
		   
	       }
	    else
	       {
	       delete children[j];
	       }

	    } // fine loop sui blocchi figlio appena creati
	
	} //fine ramo dim = 3 dell'if
	else
	{
	for (unsigned int i = 0; i < blockNodeList.size(); i++)    
            {
	    
	    // assegnamento nodi del blocco padre ai blocchi figli
            Point <dim> node = support_points[blockNodeList[i]];
	    
	       if (node(1) <= parent->GetPMin()(1)+delta)
		  {
		  if (node(0) <= parent->GetPMin()(0)+delta)
		     {
		     //pcout<<" Sono in 1 "<<std::endl;
		     children[0]->AddNode(blockNodeList[i]);
		     }
		  else
		     {
		     //pcout<<" Sono in 2 "<<std::endl;
		     children[1]->AddNode(blockNodeList[i]);
		     }
		  }
	       else
		  {
                  if (node(0) <= parent->GetPMin()(0)+delta)
		     {
		     //pcout<<" Sono in 4 "<<std::endl;
		     children[3]->AddNode(blockNodeList[i]);
		     }
		  else
		     {
		     //pcout<<" Sono in 3 "<<std::endl;
		     children[2]->AddNode(blockNodeList[i]);
		     }
		  }//fine assegnazione blocchi del padre ai blocchi figli
	       
	    } //fine loop nodi del blocco	
	
	typename std::map <cell_it, std::vector<unsigned int> >::iterator it;    
        for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)    
            {
	    for (unsigned int pp = 0; pp < (*it).second.size(); pp++)
	        {
                // assegnamento quad points del blocco padre ai blocchi figli
               Point<dim> quadPoint = quadPoints[(*it).first][(*it).second[pp]];
	       if (quadPoint(1) <= parent->GetPMin()(1)+delta)
		  {
		  if (quadPoint(0) <= parent->GetPMin()(0)+delta)
		     {
		     //pcout<<" Sono in 1 "<<std::endl;
		     children[0]->AddQuadPoint((*it).first,(*it).second[pp]);
		     }
		  else
		     {
		     //pcout<<" Sono in 2 "<<std::endl;
		     children[1]->AddQuadPoint((*it).first,(*it).second[pp]);
		     }
		  }
	       else
		  {
                  if (quadPoint(0) <= parent->GetPMin()(0)+delta)
		     {
		     //pcout<<" Sono in 4 "<<std::endl;
		     children[3]->AddQuadPoint((*it).first,(*it).second[pp]);
		     }
		  else
		     {
		     //pcout<<" Sono in 3 "<<std::endl;
		     children[2]->AddQuadPoint((*it).first,(*it).second[pp]);
		     }
		  }//fine assegnazione blocchi del padre ai blocchi figli

		}    
  	    }
	    
        for (unsigned int j=0; j < num_children_per_block; j++ )
	    {   
	    if (children[j]->GetBlockNodeList().size() + 
	        children[j]->GetBlockQuadPointsList().size()> 0)
	       {
	       blocksCount += 1;
	       blocks[blocksCount] = new OctreeBlock<dim>();
	       blocks[blocksCount]->CopyContent(children[j]);
	       delete children[j];
	       
	       parent->AddChild(blocksCount);
               std::map <cell_it, std::vector<unsigned int> >
	       blockQuadPointsList = blocks[blocksCount]->GetBlockQuadPointsList();
	       typename std::map <cell_it, std::vector<unsigned int> >::iterator it;    
               for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
	           {
		   cell_it cell = (*it).first;
		   for(unsigned int kk = 0; kk < (*it).second.size(); kk++)
		      {
		      quad_point_to_block[(*it).first][(*it).second[kk]].push_back(blocksCount);
		      }
		   }
               std::vector<unsigned int> blockNodesList = blocks[jj]->GetBlockNodeList();
	       for(unsigned int k = 0; k < blockNodesList.size(); k++)
	          dof_to_block[blockNodesList[k]].push_back(jj);
		   
	       }
	    else
	       {
	       delete children[j];
	       }

	    } // fine loop sui blocchi figlio appena creati
	    
	
	} // fine ramo dim == 2 dell'if
	
	} //fine loop blocchi livello precedente	
    
	
    //double elemCheck = numChildless;

    startLevel[level] = endLevel[level-1] + 1;
    endLevel[level] = blocksCount;
	
	// here we loop over the blocks of the newly created level and
	// we decide if each block is to be split again in the next level:
	// if it contains more
	// than a node or quad point, it will be placed in the parent list.
	// Instead, if it only contains a node or quad point, it goes in the
	// childless list, and not be refined any more. It is important to
	// account for the presence of double nodes: if not, blocks will be
	// always refined    
    for (unsigned int jj = startLevel[level]; jj < endLevel[level]+1;  jj++)
	{

	// here we get the number of nodes in the block 
	std::vector<unsigned int> nodesId = blocks[jj]->GetBlockNodeList();
	int blockNumNodes = (int) nodesId.size();
	
	// now we compute the number of the nodes that are double of others
	int numDoubleNodes = 0;  
	for (unsigned int kk = 0; kk < nodesId.size();  kk++)
	    {
	    int a = (int) double_nodes_set[nodesId[kk]].size();
	    numDoubleNodes += a - 1;
	    }   
	
	// here we compute the number of quad points in the block
	int blockNumQuadPoints = 0;
	std::map <cell_it, std::vector<unsigned int> >
	blockQuadPointsList = blocks[jj]->GetBlockQuadPointsList();
	typename std::map <cell_it, std::vector<unsigned int> >::iterator it;    
        for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
	     blockNumQuadPoints += (int) (*it).second.size();
	//pcout<<"Level "<<level<<" Block "<<jj<<"  nodes "<<blockNumNodes<<" + quadPoints "<<blockNumQuadPoints<<std::endl;     
	
	quadPointsCheck += blockNumQuadPoints;
	nodesCheck += blockNumNodes;
	// here we decide if a block is to be placed in the parent
	// or childless list
	//if (blockNumNodes + blockNumQuadPoints - numDoubleNodes < 2)
	if (blockNumNodes - numDoubleNodes < 2)
	   {
	   numChildless += 1;
	   childlessList.push_back(jj);
	   quadPointsInChildless += blockNumQuadPoints; 
	   nodesInChildless += blockNumNodes;
	   
	   
	   // if a block is childless, we must assign now the nodes and quad points
	   // that belong to it for all the next levels
	    
	for (unsigned int kk = 0; kk < nodesId.size();  kk++)   
	    for (unsigned int j = level+1; j < num_octree_levels+1; j++)
	        dof_to_block[nodesId[kk]].push_back(jj);
	    
	for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
	    for (unsigned int i = 0; i < (*it).second.size(); i++)    
	        for (unsigned int j = level+1; j < num_octree_levels+1; j++)
	             quad_point_to_block[(*it).first][(*it).second[i]].push_back(jj);      
	   
	   }
	else
	   {
	   numParent[level] += 1;
	   parentList[level].push_back(jj);
	   }
	
	   // let's update the list of node filled block
	if (blockNumNodes > 0)
	   dofs_filled_blocks[level].push_back(jj);

	   // let's update the list of quad point filled block	   
	if (blockNumQuadPoints > 0)
	   quad_points_filled_blocks[level].push_back(jj);   
	//elemCheck += blockNumNodes + blockNumQuadPoints;  
        }


     pcout<<" Total nodes at level "<<level<<" of "<<num_octree_levels<<" are "<<nodesCheck<<std::endl;
     pcout<<" Total quad points at level "<<level<<" of "<<num_octree_levels<<" are "<<quadPointsCheck<<std::endl;
     pcout<<" Blocks at level "<<level<<" of "<<num_octree_levels<<" are "<<endLevel[level]-endLevel[level-1]<<std::endl;
     pcout<<" Total blocks at level "<<level<<" of "<<num_octree_levels<<" are "<<endLevel[level] + 1<<std::endl;
     pcout<<std::endl;	   
	
     } //fine loop livelli

childlessList.resize(childlessList.size()+parentList[num_octree_levels].size());

for (unsigned int jj = 0; jj < parentList[num_octree_levels].size();  jj++)
	{
	 childlessList[numChildless + jj] = parentList[num_octree_levels][jj];
	}




num_blocks = blocksCount+1;

	
pcout<<"...done generating octree blocking"<<std::endl;

pcout<<"Computing proximity lists for blocks"<<std::endl;

//just for output
/*for (cell = dh.begin_active(); cell != endc; ++cell)    
    {
    unsigned int levelCheck = elem_to_blocks[cell].size();
    pcout<<std::endl<<"Elem "<<cell<<" is in the "<<levelCheck<<" blocks: ";
    for (unsigned int zz = 0; zz < levelCheck; zz++)
         pcout<<elem_to_blocks[cell][zz]<<" ";
    }*/

/*for (cell_it cell = dh.begin_active(); cell != endc; ++cell)
    for (unsigned int j=0; j < quadPoints[cell].size(); j++)
        pcout<<"Cell "<<cell<<"  QP "<<j<<"  of "<<quadPoints[cell].size()<<": "<<quadPoints[cell][j]<<std::endl;//*/

/*for (cell_it cell = dh.begin_active(); cell != endc; ++cell)
    for (unsigned int j=0; j < quad_point_to_block[cell].size(); j++)
        {
	pcout<<"Cell "<<cell<<"  QP "<<j<<"  of "<<quad_point_to_block[cell].size()<<": ";
        for (unsigned int i=0; i < quad_point_to_block[cell][j].size(); i++)
            pcout<<quad_point_to_block[cell][j][i]<<" ";
	pcout<<std::endl;    
	}  //*/  


/*for (unsigned int i=0; i < dh.n_dofs(); i++)
    {
    pcout<<"Node "<<i<<"  doubles: ";
    std::set <unsigned int> doubleNodes = double_nodes_set[i]; 
    for (std::set<unsigned int>::iterator pos = doubleNodes.begin(); pos != doubleNodes.end(); pos++)
        {
	pcout<<*pos<<"( ";
	for (unsigned int j=0; j < dof_to_elems[*pos].size(); j++)
	    pcout<<" "<<dof_to_elems[*pos][j];
	pcout<<") ";
        }
    pcout<<std::endl;
    } //*/



// ricerca blocchi nearest neighbors

for (unsigned int ii = startLevel[1]; ii < endLevel[1] + 1;  ii++)
	{
	for (unsigned int jj = startLevel[1]; jj < endLevel[1] + 1;  jj++)
		{
		blocks[ii]->AddNearNeigh(0,jj);
		}
	}


for (unsigned int level = 2; level < num_octree_levels + 1;  level++)

	{
	for (unsigned int kk = startLevel[level]; kk < endLevel[level]+1;  kk++)

		{
		OctreeBlock<dim> *block1 = blocks[kk];
		block1->AddNearNeigh(0,kk); // a block is NearNeigh of itself

		double delta1 = block1->GetDelta();
		Point<dim> PMin1 = block1->GetPMin();
		Point<dim> Center1;
		for (unsigned int iii = 0; iii < dim; iii++)
		     Center1(iii) = delta1/2.;
		Point<dim> PMax1 = 2.*Center1;     
		PMax1 += PMin1;
		Center1 += PMin1;
		unsigned int parentId = block1->GetParentId();
		std::set <unsigned int> parentNNeighs = blocks[parentId]->GetNearNeighs(0);

                // the nearest neighbors are searched among the father's nearest neighbors children
		for (std::set <unsigned int>::iterator pos = parentNNeighs.begin(); pos != parentNNeighs.end();  pos++)

			{
		        if (blocks[*pos]->GetBlockChildrenNum() == 0) // if a parent's near neigh is childless, he can be a near neigh: let's check
				{
				unsigned int block2Id = *pos;
				OctreeBlock<dim> *block2 = blocks[block2Id];
				double delta2 = block2->GetDelta();
		                Point<dim> PMin2 = block2->GetPMin();
		                Point<dim> Center2;
		                for (unsigned int iii = 0; iii < dim; iii++)
		                    Center2(iii) = delta2/2.;
		                Point<dim> PMax2 = 2.*Center2;     
		                PMax2 += PMin2;
		                Center2 += PMin2;
				
				if (dim == 3)
				{
				if 	((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
					{
					if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
						{
						if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
							{
							block1->AddNearNeigh(0,block2Id);
							//pcout<<" *"<<block2Id;
							}
						}
					}

				if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
					{
					if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
						{
						if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
							{
							block1->AddNearNeigh(0,block2Id);
							//pcout<<" *"<<block2Id;
							}
						}
					}

				if ((fabs(PMin1(2) - PMax2(2)) <= TOLL) || (fabs(PMax1(2) - PMin2(2)) <= TOLL))
					{
					if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
						{
						if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
							{
							block1->AddNearNeigh(0,block2Id);
							//pcout<<" *"<<block2Id;
							}
						}
					}
				} //fine caso dim ==3
						
				else if (dim == 2)
				{
				if 	((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
					{
					if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
						{
						block1->AddNearNeigh(0,block2Id);
						//pcout<<block2Id<<" ";
						}
					}

				if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
					{
					if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
						{
						block1->AddNearNeigh(0,block2Id);
						//pcout<<block2Id<<" ";
						}
					}

                                } // fine caso dim == 2 
				
			        }

			for (unsigned int ii = 0; ii < blocks[*pos]->GetBlockChildrenNum();  ii++)
				{
				unsigned int block2Id = blocks[*pos]->GetChildId(ii);
				OctreeBlock<dim> *block2 = blocks[block2Id];
				double delta2 = block2->GetDelta();
		                Point<dim> PMin2 = block2->GetPMin();
		                Point<dim> Center2;
		                for (unsigned int iii = 0; iii < dim; iii++)
		                    Center2(iii) = delta2/2.;
		                Point<dim> PMax2 = 2.*Center2;     
		                PMax2 += PMin2;
		                Center2 += PMin2;
				
				if (dim == 3)
				{
				if 	((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
					{
					if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
						{
						if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
							{
							block1->AddNearNeigh(0,block2Id);
							//pcout<<" "<<block2Id;
							}
						}
					}

				if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
					{
					if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
						{
						if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
							{
							block1->AddNearNeigh(0,block2Id);
							//pcout<<" "<<block2Id;
							}
						}
					}

				if ((fabs(PMin1(2) - PMax2(2)) <= TOLL) || (fabs(PMax1(2) - PMin2(2)) <= TOLL))
					{
					if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
						{
						if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
							{
							block1->AddNearNeigh(0,block2Id);
							//pcout<<" "<<block2Id;
							}
						}
					} 
				} //fine caso dim ==3
						
				else if (dim == 2)
				{
				if 	((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
					{
					if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
						{
						block1->AddNearNeigh(0,block2Id);
						//pcout<<block2Id<<" ";
						}
					}

				if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
					{
					if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
						{
						block1->AddNearNeigh(0,block2Id);
						//pcout<<block2Id<<" ";
						}
					}

                                } // fine caso dim == 2 
				
				} // fine loop sui figli di un nearest neighbor del padre


			} // fine loop sui nearest neighbors del padre

                        

		if ((block1->GetBlockChildrenNum() == 0))  // if the block is childless we must compute now its nearneigh at all residual levels
			{
			block1->SetNearNeighSize(num_octree_levels-level+1);
			block1->SetIntListSize(num_octree_levels-level+1); // intList is a vector of sets with the same number of members of nearNeigh
			block1->SetNonIntListSize(num_octree_levels-level+1); // nonIntList is a vector of sets with the same number of members of nearNeigh
			
			for (unsigned int subLevel = 1; subLevel < num_octree_levels - level + 1; subLevel++)

				{
				
				std::set <unsigned int> upperLevelNNeighs = block1->GetNearNeighs(subLevel-1);
				for (std::set <unsigned int>::iterator pos = upperLevelNNeighs.begin(); pos != upperLevelNNeighs.end();  pos++)

					{
					if (blocks[*pos]->GetBlockChildrenNum() == 0) // if nearneigh is childless, it will stay a near neigh
						block1->AddNearNeigh(subLevel,*pos);

					for (unsigned int ii = 0; ii < blocks[*pos]->GetBlockChildrenNum();  ii++)
					    {
						unsigned int block2Id = blocks[*pos]->GetChildId(ii);
						OctreeBlock<dim>* block2 = blocks[block2Id];
						double delta2 = block2->GetDelta();
		                		Point<dim> PMin2 = block2->GetPMin();
		                		Point<dim> Center2;
		                		for (unsigned int iii = 0; iii < dim; iii++)
		                    		    Center2(iii) = delta2/2.;
		                		Point<dim> PMax2 = 2.*Center2;     
		                		PMax2 += PMin2;
		                		Center2 += PMin2;
						
						if (dim == 3)
						{
						if 	((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
							{
							if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
								{
								if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
									{
									block1->AddNearNeigh(subLevel,block2Id);
									//pcout<<block2Id<<" ";
									}
								}
							}

						if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
							{
							if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
								{
								if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
									{
									block1->AddNearNeigh(subLevel,block2Id);
									//pcout<<block2Id<<" ";
									}
								}
							}

						if ((fabs(PMin1(2) - PMax2(2)) <= TOLL) || (fabs(PMax1(2) - PMin2(2)) <= TOLL))
							{
							if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
								{
								if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
									{
									block1->AddNearNeigh(subLevel,block2Id);
									//pcout<<block2Id<<" ";
									}
								}
							} 
                                                } //fine caso dim ==3
						
						else if (dim == 2)
						{
						if 	((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
							{
							if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
								{
								block1->AddNearNeigh(subLevel,block2Id);
								//pcout<<block2Id<<" ";
								}
							}

						if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
							{
							if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
								{
								block1->AddNearNeigh(subLevel,block2Id);
								//pcout<<block2Id<<" ";
								}
							}

                                                } // fine caso dim == 2 

					    } // fine loop sui figli di ciascun nearest neighbor del blocco childless


					} // fine loop sui nearest neighbors del blocco childless


				} // fine loop sui subLevels (da quello del blocco childless all'ultimo)


			} // fine if (il blocco e' childless?)



	} // fine loop sui blocchi di un livello

} // fine loop sui livelli

//for printout                      

/*pcout<<std::endl;
pcout<<"-------------------------------- "<<std::endl;
pcout<<"-------------------------------- "<<std::endl;
pcout<<std::endl;
		
//for printout
for (cell=dh.begin_active();cell!=endc; cell++)
    {
    pcout<<std::endl;
    pcout<<"-------------------------------- "<<std::endl;
    pcout<<"Cell "<<cell<<"  elementPlot(";
    cell->get_dof_indices(local_dof_indices);
    for (unsigned int j = 0; j<local_dof_indices.size(); j++)
        pcout<<"["<<support_points[local_dof_indices[j]]<<"],";
        pcout<<"'r')";
    }
			
pcout<<std::endl;
pcout<<"-------------------------------- "<<std::endl;
pcout<<"-------------------------------- "<<std::endl;
pcout<<std::endl;
pcout<<std::endl;//*/


// search for interaction list blocks (NearNeigh + NearNeighOfNearNeigh)
// and for non interaction list blocks (nonIntList for a block B is composed by blocks that are children
// of blocks being in intList of B's parent, but are not in intList of B)

for (unsigned int ii = startLevel[1]; ii < endLevel[1] + 1;  ii++) // at level 1, add all blocks to intList
	{
	for (unsigned int jj = startLevel[1]; jj < endLevel[1] + 1;  jj++)
		{
		blocks[ii]->AddBlockToIntList(0,jj);
		}
	}


for (unsigned int level = 2; level < num_octree_levels + 1;  level++) // loop over levels

	{
	for (unsigned int jj = startLevel[level]; jj < endLevel[level] + 1;  jj++) // loop over blocks of each level
		{
		OctreeBlock<dim>* block1 = blocks[jj];
		//pcout<<"??Out "<<jj<<" "<<block1.GetNearNeighSize()<<std::endl;
		for (unsigned int subLevel = 0; subLevel < block1->NumNearNeighLevels(); subLevel++)
                    {
		    std::set <unsigned int> NNList = block1->GetNearNeighs(subLevel);
		    
		    for (std::set <unsigned int>::iterator pos1 = NNList.begin(); pos1 != NNList.end();  pos1++) //loop over blocks in NN list and get their NNs
		        {
			block1->AddBlockToIntList(subLevel,*pos1);
			}
		    //pcout<<std::endl<<"Sublevel "<<subLevel<<" elem("<<block1.GetBlockElementsList()[0]<<") NearNeighs: ";
		    std::vector <unsigned int> nodeIds = block1->GetBlockNodeList();
		    //pcout<<std::endl<<"Level "<<level<<"  Block1: "<<kk<<"  NumElements: "<<block1.GetBlockElementsNum()<<" "<<elemIds.size()<<std::endl;
		    //pcout<<"Nearest Neighbors Found:"<<std::endl;
		    for (unsigned int pp = 0; pp < nodeIds.size();  pp++)
		        {
			//pcout<<"Node "<<nodeIds[pp]<<std::endl;
			std::set <unsigned int> doubleNodes = double_nodes_set[nodeIds[pp]];
			for (std::set<unsigned int>::iterator pos = doubleNodes.begin();
			     pos != doubleNodes.end(); pos++)
			     {
			     //pcout<<"Node Double"<<*pos<<std::endl;
			     //std::vector<cell_it > surrCellIds = dof_to_elems[*pos];
			     for (unsigned int k=0; k < dof_to_elems[*pos].size(); k++)
			         {
				 cell_it cell = dof_to_elems[*pos][k];
				 //pcout<<cell<<std::endl;
				 for (unsigned int j=0; j < quadPoints[cell].size(); j++)
				     {
				     block1->AddBlockToIntList(subLevel,quad_point_to_block[cell][j][level+subLevel]);
				     }
				 }  
			     }
		        }
                    }
		for (unsigned int subLevel = 0; subLevel < block1->GetNearNeighSize();  subLevel++) // for each block, loop over all sublevels in his NN list (to account for childless blocks)
			{
			// now use intList to compute nonIntList
			std::set <unsigned int> intList = block1->GetIntList(subLevel);
			std::set <unsigned int> parentIntList; // intList at the  previous level
			if (subLevel == 0) // if a block is childless we get its intList at the previous level, otherwise we get its parent's intList
				parentIntList = blocks[block1->GetParentId()]->GetIntList(0);
			else
				parentIntList = block1->GetIntList(subLevel-1);

			for (std::set <unsigned int>::iterator pos1 = parentIntList.begin(); pos1 != parentIntList.end();  pos1++) // loop over blocks in parentIntList
				{
			 	OctreeBlock<dim>* block2 = blocks[*pos1];
			 	if (block2->GetBlockChildrenNum() == 0) // if blocks in parentIntList are childless, don't look for their children, but see if they are in nonIntList
			 		{
			 		if (intList.count(*pos1) == 0) // if these blocks are not in intList
			 			block1->AddBlockToNonIntList(subLevel,*pos1); // then they go in nonIntList
			 		}
			 	else // if blocks in parentIntList are not childless, do the same test on all their children
			 		{
			 		for (unsigned int kk = 0; kk < block2->GetBlockChildrenNum();  kk++) // loop over children of blocks in parentIntList
			 			{
						//pcout<<"Sublevel "<<subLevel<<" Block1 "<<jj<<"  Block2 "<<*pos1<<" child(kk) "<<block2.GetChildId(kk)<<std::endl;
			 			if (intList.count(block2->GetChildId(kk)) == 0)   // if these blocks are not in intList
			 				block1->AddBlockToNonIntList(subLevel,block2->GetChildId(kk)); // then they go in nonIntList
			 			} // end loop over children of blocks in parentIntList
			 		}
				}	// loop over blocks in parentIntList

			}	// end loop over subLevels of each block's intList

//			for printout
                        
/*			
			pcout<<std::endl;
                        pcout<<"-------------------------------- "<<std::endl;
			pcout<<"Block "<<jj<<": Parent "<<block1->GetParentId();
			pcout<<"  cubePlot(["<<block1->GetPMin()<<"],"<<block1->GetDelta()<<",'m')"<<std::endl;
			pcout<<"Quad points of block "<<jj<<": "<<std::endl;
			std::map <cell_it,std::vector<unsigned int> > quadPointsList = block1->GetBlockQuadPointsList();
			typename std::map <cell_it, std::vector<unsigned int> >::iterator it;    
                        for (it = quadPointsList.begin(); it != quadPointsList.end(); it++)
	                    {
			    pcout<<(*it).first<<"( ";
			    for (unsigned int zz = 0; zz < (*it).second.size();  zz++)
			        pcout<<(*it).second[zz]<<" ";
			    pcout<<")  ";
			    }
			pcout<<std::endl;

			pcout<<"Nodes of block "<<jj<<" :"<<std::endl;
			std::vector <unsigned int> nodeList = block1->GetBlockNodeList();
			for (unsigned int zz = 0; zz < block1->GetBlockNodesNum();  zz++)
				{
				pcout<<nodeList.at(zz)<<" ";
				}
			pcout<<std::endl;

			for (unsigned int subLevel = 0; subLevel < block1->GetNearNeighSize();  subLevel++)
				{
				std::set <unsigned int> NNList = block1->GetNearNeighs(subLevel);
				pcout<<"NearNeigh for block "<<jj<<" at level "<<level+subLevel<<":"<<std::endl;
				for (std::set <unsigned int>::iterator pos1 = NNList.begin(); pos1 != NNList.end();  pos1++)
					{
					pcout<<*pos1<<" ";
					}
				pcout<<std::endl;

				std::set <unsigned int> intList = block1->GetIntList(subLevel);
				pcout<<"IntList for block "<<jj<<" at level "<<level+subLevel<<":"<<std::endl;
				for (std::set <unsigned int>::iterator pos1 = intList.begin(); pos1 != intList.end();  pos1++)
					{
					pcout<<*pos1<<" ";
					}
				pcout<<std::endl;

				std::set <unsigned int> nonIntList = block1->GetNonIntList(subLevel);
				pcout<<"NonIntList for block "<<jj<<" at level "<<level+subLevel<<":"<<std::endl;
				for (std::set <unsigned int>::iterator pos1 = nonIntList.begin(); pos1 != nonIntList.end();  pos1++)
					{
					pcout<<*pos1<<" ";
					}
				pcout<<std::endl;
				}
//*/



		} // end loop over blocks of a level

	} // end loop over levels

//if (blocks.size() > 0)
//   {
//   for (unsigned int ii = 0; ii < num_blocks;  ii++)
//        delete blocks[ii];
//   }

integralCheck.clear(); 
for (unsigned int i = 0; i < dh.n_dofs(); i++)
    {
    for (cell_it cell = dh.begin_active(); cell != endc; ++cell)
        {
	//pcout<<i<<" "<<cell<<" "<<integralCheck[i][cell]<<std::endl;
	integralCheck[i][cell] = 0;
	}
    }



pcout<<"Done computing proximity lists for blocks"<<std::endl;
} //end method for octree blocking generation






  template <int dim>
  void ComputationalDomain<dim>::compute_normals()
  {
   typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;

   SparsityPattern      normals_sparsity_pattern;
   normals_sparsity_pattern.reinit(gradient_dh.n_dofs(),
                                   gradient_dh.n_dofs(),
                                   gradient_dh.max_couplings_between_dofs());
   ConstraintMatrix  vector_constraints;
   vector_constraints.clear();
   DoFTools::make_hanging_node_constraints (gradient_dh,vector_constraints);
   vector_constraints.close();
   DoFTools::make_sparsity_pattern (gradient_dh, normals_sparsity_pattern, vector_constraints);
   normals_sparsity_pattern.compress();
   Vector<double> vector_normals_solution(gradient_dh.n_dofs());
                                   
   SparseMatrix<double> vector_normals_matrix;
   Vector<double> vector_normals_rhs;

   vector_normals_matrix.reinit (normals_sparsity_pattern);
   vector_normals_rhs.reinit(gradient_dh.n_dofs());
   vector_normals_solution.reinit(gradient_dh.n_dofs());


   FEValues<dim-1,dim> gradient_fe_v(*mapping, gradient_fe, *quadrature,
			     update_values | update_cell_normal_vectors |  
			     update_JxW_values);

   const unsigned int vector_n_q_points = gradient_fe_v.n_quadrature_points;
   const unsigned int   vector_dofs_per_cell   = gradient_fe.dofs_per_cell;
   std::vector<unsigned int> vector_local_dof_indices (vector_dofs_per_cell);

   FullMatrix<double>   local_normals_matrix (vector_dofs_per_cell, vector_dofs_per_cell);
   Vector<double>       local_normals_rhs (vector_dofs_per_cell);

   cell_it
   vector_cell = gradient_dh.begin_active(),
   vector_endc = gradient_dh.end();
            
   for (; vector_cell!=vector_endc; ++vector_cell)
     {
       gradient_fe_v.reinit (vector_cell);
       local_normals_matrix = 0;
       local_normals_rhs = 0;
       const std::vector<Point<dim> > &vector_node_normals = gradient_fe_v.get_normal_vectors();
       unsigned int comp_i, comp_j;
       
       for (unsigned int q=0; q<vector_n_q_points; ++q)
	 for (unsigned int i=0; i<vector_dofs_per_cell; ++i)
	   {
	     comp_i = gradient_fe.system_to_component_index(i).first;
	     for (unsigned int j=0; j<vector_dofs_per_cell; ++j)
	       {
		 comp_j = gradient_fe.system_to_component_index(j).first;
		 if (comp_i == comp_j) 
		   {
		     local_normals_matrix(i,j) += gradient_fe_v.shape_value(i,q)*
						  gradient_fe_v.shape_value(j,q)*
						  gradient_fe_v.JxW(q);
		   }
	       }
	   local_normals_rhs(i) += (gradient_fe_v.shape_value(i, q)) *
                                    vector_node_normals[q](comp_i) * gradient_fe_v.JxW(q);
	   }
        
       vector_cell->get_dof_indices (vector_local_dof_indices);
       
       vector_constraints.distribute_local_to_global
       (local_normals_matrix,
	local_normals_rhs,
	vector_local_dof_indices,
	vector_normals_matrix,
	vector_normals_rhs);
     }



   SparseDirectUMFPACK normals_inverse;
   normals_inverse.initialize(vector_normals_matrix);
   normals_inverse.vmult(vector_normals_solution, vector_normals_rhs);
   vector_constraints.distribute(vector_normals_solution);

   node_normals.resize(dh.n_dofs());
 
   for (unsigned int i=0; i<gradient_dh.n_dofs()/dim; ++i)
       { 
       for (unsigned int d=0; d<dim; d++)
           node_normals[i](d) = vector_normals_solution(3*i+d);
       node_normals[i]/= node_normals[i].distance(Point<dim>(0.0,0.0,0.0));
       //cout<<i<<" Gradient: "<<node_normals[i]<<endl;
       for (unsigned int d=0; d<dim; d++)
           vector_normals_solution(3*i+d) = node_normals[i](d);
       }


  }


template class ComputationalDomain<3>;
