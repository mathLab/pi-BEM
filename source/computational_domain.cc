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
  if (mapping != NULL)
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

  vertices[0](0) = 0;
  vertices[0](1) = 0;
  vertices[0](2) = 1;
  vertices[1](0) = 0.57735;
  vertices[1](1) = 0.57735;
  vertices[1](2) = 0.57735;
  vertices[2](0) = -0.57735;
  vertices[2](1) = 0.57735;
  vertices[2](2) = 0.57735;
  vertices[3](0) = -0.57735;
  vertices[3](1) = -0.57735;
  vertices[3](2) = 0.57735;
  vertices[4](0) = 0.57735;
  vertices[4](1) = -0.57735;
  vertices[4](2) = 0.57735;
  vertices[5](0) = 0.70711;
  vertices[5](1) = 0;
  vertices[5](2) = 0.70711;
  vertices[6](0) = 0;
  vertices[6](1) = 0.70711;
  vertices[6](2) = 0.70711;
  vertices[7](0) = -0.70711;
  vertices[7](1) = 0;
  vertices[7](2) = 0.70711;
  vertices[8](0) = 0;
  vertices[8](1) = -0.70711;
  vertices[8](2) = 0.70711;
  vertices[9](0) = 1;
  vertices[9](1) = 0;
  vertices[9](2) = 0;
  vertices[10](0) = 0;
  vertices[10](1) = 1;
  vertices[10](2) = 0;
  vertices[11](0) = -1;
  vertices[11](1) = 0;
  vertices[11](2) = 0;
  vertices[12](0) = 0;
  vertices[12](1) = -1;
  vertices[12](2) = 0;
  vertices[13](0) = 0.70711;
  vertices[13](1) = 0.70711;
  vertices[13](2) = 0;
  vertices[14](0) = -0.70711;
  vertices[14](1) = 0.70711;
  vertices[14](2) = 0;
  vertices[15](0) = -0.70711;
  vertices[15](1) = -0.70711;
  vertices[15](2) = 0;
  vertices[16](0) = 0.70711;
  vertices[16](1) = -0.70711;
  vertices[16](2) = 0;
  vertices[17](0) = 0;
  vertices[17](1) = 0;
  vertices[17](2) = -1;
  vertices[18](0) = 0.57735;
  vertices[18](1) = 0.57735;
  vertices[18](2) = -0.57735;
  vertices[19](0) = -0.57735;
  vertices[19](1) = 0.57735;
  vertices[19](2) = -0.57735;
  vertices[20](0) = -0.57735;
  vertices[20](1) = -0.57735;
  vertices[20](2) = -0.57735;
  vertices[21](0) = 0.57735;
  vertices[21](1) = -0.57735;
  vertices[21](2) = -0.57735;
  vertices[22](0) = 0.70711;
  vertices[22](1) = 0;
  vertices[22](2) = -0.70711;
  vertices[23](0) = 0;
  vertices[23](1) = 0.70711;
  vertices[23](2) = -0.70711;
  vertices[24](0) = -0.70711;
  vertices[24](1) = 0;
  vertices[24](2) = -0.70711;
  vertices[25](0) = 0;
  vertices[25](1) = -0.70711;
  vertices[25](2) = -0.70711;
  vertices[26](0) = 1;
  vertices[26](1) = 0;
  vertices[26](2) = 0;
  vertices[27](0) = 0;
  vertices[27](1) = 1;
  vertices[27](2) = 0;
  vertices[28](0) = -1;
  vertices[28](1) = 0;
  vertices[28](2) = 0;
  vertices[29](0) = 0;
  vertices[29](1) = -1;
  vertices[29](2) = 0;
  vertices[30](0) = 0.70711;
  vertices[30](1) = 0.70711;
  vertices[30](2) = 0;
  vertices[31](0) = -0.70711;
  vertices[31](1) = 0.70711;
  vertices[31](2) = 0;
  vertices[32](0) = -0.70711;
  vertices[32](1) = -0.70711;
  vertices[32](2) = 0;
  vertices[33](0) = 0.70711;
  vertices[33](1) = -0.70711;
  vertices[33](2) = 0;


  cells.resize(24);



  cells[0].vertices[0]=0;
  cells[0].vertices[1]=5;
  cells[0].vertices[2]=1;
  cells[0].vertices[3]=6;
  cells[1].vertices[0]=7;
  cells[1].vertices[1]=0;
  cells[1].vertices[2]=6;
  cells[1].vertices[3]=2;
  cells[2].vertices[0]=7;
  cells[2].vertices[1]=3;
  cells[2].vertices[2]=8;
  cells[2].vertices[3]=0;
  cells[3].vertices[0]=4;
  cells[3].vertices[1]=5;
  cells[3].vertices[2]=0;
  cells[3].vertices[3]=8;
  cells[4].vertices[0]=9;
  cells[4].vertices[1]=13;
  cells[4].vertices[2]=1;
  cells[4].vertices[3]=5;
  cells[5].vertices[0]=13;
  cells[5].vertices[1]=10;
  cells[5].vertices[2]=6;
  cells[5].vertices[3]=1;
  cells[6].vertices[0]=10;
  cells[6].vertices[1]=14;
  cells[6].vertices[2]=2;
  cells[6].vertices[3]=6;
  cells[7].vertices[0]=14;
  cells[7].vertices[1]=11;
  cells[7].vertices[2]=7;
  cells[7].vertices[3]=2;
  cells[8].vertices[0]=11;
  cells[8].vertices[1]=15;
  cells[8].vertices[2]=3;
  cells[8].vertices[3]=7;
  cells[9].vertices[0]=15;
  cells[9].vertices[1]=12;
  cells[9].vertices[2]=8;
  cells[9].vertices[3]=3;
  cells[10].vertices[0]=12;
  cells[10].vertices[1]=16;
  cells[10].vertices[2]=4;
  cells[10].vertices[3]=8;
  cells[11].vertices[0]=16;
  cells[11].vertices[1]=9;
  cells[11].vertices[2]=5;
  cells[11].vertices[3]=4;

  cells[12].vertices[0]=29;
  cells[12].vertices[1]=32;
  cells[12].vertices[2]=20;
  cells[12].vertices[3]=25;
  cells[13].vertices[0]=33;
  cells[13].vertices[1]=29;
  cells[13].vertices[2]=25;
  cells[13].vertices[3]=21;
  cells[14].vertices[0]=32;
  cells[14].vertices[1]=28;
  cells[14].vertices[2]=24;
  cells[14].vertices[3]=20;
  cells[15].vertices[0]=28;
  cells[15].vertices[1]=31;
  cells[15].vertices[2]=19;
  cells[15].vertices[3]=24;
  cells[16].vertices[0]=31;
  cells[16].vertices[1]=27;
  cells[16].vertices[2]=23;
  cells[16].vertices[3]=19;
  cells[17].vertices[0]=27;
  cells[17].vertices[1]=30;
  cells[17].vertices[2]=18;
  cells[17].vertices[3]=23;
  cells[18].vertices[0]=30;
  cells[18].vertices[1]=26;
  cells[18].vertices[2]=22;
  cells[18].vertices[3]=18;
  cells[19].vertices[0]=26;
  cells[19].vertices[1]=33;
  cells[19].vertices[2]=21;
  cells[19].vertices[3]=22;
  cells[20].vertices[0]=25;
  cells[20].vertices[1]=20;
  cells[20].vertices[2]=24;
  cells[20].vertices[3]=17;
  cells[21].vertices[0]=22;
  cells[21].vertices[1]=21;
  cells[21].vertices[2]=25;
  cells[21].vertices[3]=17;
  cells[22].vertices[0]=24;
  cells[22].vertices[1]=19;
  cells[22].vertices[2]=23;
  cells[22].vertices[3]=17;
  cells[23].vertices[0]=23;
  cells[23].vertices[1]=18;
  cells[23].vertices[2]=22;
  cells[23].vertices[3]=17;


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

  if (mapping == NULL)
    mapping = new MappingQEulerian<dim-1, Vector<double>, dim>
    (gradient_fe.degree, gradient_dh, map_points);

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
      for (unsigned int i=k; i<gradient_dh.n_dofs(); i=i+dim)
        {
          for (unsigned int j=k; j<gradient_dh.n_dofs(); j=j+dim)
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
  //fma.dof_to_elems.clear();

  // mappa che associa ad ogni gradient dof le celle cui esso appartiene
  //fma.gradient_dof_to_elems.clear();

  // vettore che associa ad ogni gradient dof la sua componente
  //fma.gradient_dof_components.clear();
  //fma.gradient_dof_components.resize(gradient_dh.n_dofs());

  // mappa che associa ad ogni cella un set contenente le celle circostanti
  //fma.elem_to_surr_elems.clear();

  // set che raccoglie i nodi della free surface che stanno sulla barca
  free_surf_and_boat_nodes.clear();

  // mappa raccoglie i nodi degli edges della barca
  boat_keel_nodes.clear();


  for (; cell!=endc; ++cell,++gradient_cell)
    {
      Assert(cell->index() == gradient_cell->index(), ExcInternalError());

      cell->get_dof_indices(dofs);
      for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
        {
          //fma.dof_to_elems[dofs[j]].push_back(cell);
        }
      gradient_cell->get_dof_indices(gradient_dofs);
      for (unsigned int j=0; j<gradient_fe.dofs_per_cell; ++j)
        {
          //fma.gradient_dof_to_elems[gradient_dofs[j]].push_back(gradient_cell);
          //fma.gradient_dof_components[gradient_dofs[j]] = gradient_fe.system_to_component_index(j).first;
        }
    }

  // qui viene creata la mappa dei elmenti che circondano ciascun elemento
  for (cell = dh.begin_active(); cell != endc; ++cell)
    {
      cell->get_dof_indices(dofs);
      for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
        {
          std::set <unsigned int> duplicates = double_nodes_set[dofs[j]];
          // TO BE FIXED !!!
          // for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
          //     {
          //     std::vector<cell_it>
          //     dof_cell_list = //fma.dof_to_elems[*pos];
          //     //for (unsigned int k=0; k<dof_cell_list.size(); ++k)
          //               //fma.elem_to_surr_elems[cell].insert(dof_cell_list[k]);
          //     }
        }
    }

  gradient_cell = gradient_dh.begin_active();
  cell = dh.begin_active();

  for (; cell!=endc; ++cell,++gradient_cell)
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
              // TO BE FIXED !!!
              //    for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
              //        {
              //        //cell_it duplicate_cell = //fma.gradient_dof_to_elems[*pos][0];
              //  if (duplicate_cell->material_id() == free_sur_ID1 ||
              //                  duplicate_cell->material_id() == free_sur_ID2 ||
              //                  duplicate_cell->material_id() == free_sur_ID3)
              //     {
              //     free_surf_and_boat_nodes.insert(gradient_dofs[j]);
              //     }
              //
              //        }
            }
          //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"  surfNodes: "<<surface_nodes(dofs[i])<<"  otherNodes: "<<other_nodes(dofs[i])<<std::endl;
        }
    }

  gradient_cell = gradient_dh.begin_active();
  cell = dh.begin_active();

  for (; cell!=endc; ++cell,++gradient_cell)
    {
      Assert(cell->index() == gradient_cell->index(), ExcInternalError());

      if (gradient_cell->material_id() == wall_sur_ID1)
        {
          gradient_cell->get_dof_indices(gradient_dofs);
          for (unsigned int j=0; j<gradient_fe.dofs_per_cell; ++j)
            {
              std::set <unsigned int> duplicates = gradient_double_nodes_set[gradient_dofs[j]];
              // TO BE FIXED !!!
              //    for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
              //        {
              //        //cell_it duplicate_cell = //fma.gradient_dof_to_elems[*pos][0];
              //  if (
              //                  duplicate_cell->material_id() == wall_sur_ID2 ||
              //                  duplicate_cell->material_id() == wall_sur_ID3)
              //     {
              //     boat_keel_nodes.insert(gradient_dofs[j]);
              //     }
              //
              //        }
            }
          //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"  surfNodes: "<<surface_nodes(dofs[i])<<"  otherNodes: "<<other_nodes(dofs[i])<<std::endl;
        }
    }
  gradient_cell = gradient_dh.begin_active();
  cell = dh.begin_active();

  for (; cell!=endc; ++cell,++gradient_cell)
    {
      Assert(cell->index() == gradient_cell->index(), ExcInternalError());

      if (gradient_cell->material_id() == wall_sur_ID2)
        {
          // This is a free surface node.
          gradient_cell->get_dof_indices(gradient_dofs);
          for (unsigned int j=0; j<gradient_fe.dofs_per_cell; ++j)
            {
              std::set <unsigned int> duplicates = gradient_double_nodes_set[gradient_dofs[j]];
              // TO BE FIXED !!!
              //    for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
              //        {
              //        //cell_it duplicate_cell = //fma.gradient_dof_to_elems[*pos][0];
              // //  if (duplicate_cell->material_id() == wall_sur_ID3)
              // //     {
              // //     boat_keel_nodes.insert(gradient_dofs[j]);
              // //     }
              //
              //        }
            }
          //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"  surfNodes: "<<surface_nodes(dofs[i])<<"  otherNodes: "<<other_nodes(dofs[i])<<std::endl;
        }
    }


  pcout<<"...done"<<std::endl;
}


template <int dim>
void ComputationalDomain<dim>::compute_phi_nodes()
{

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

  for (; cell != endc; ++cell)
    {
      if (cell->material_id() == free_sur_ID1 ||
          cell->material_id() == free_sur_ID2 ||
          cell->material_id() == free_sur_ID3)
        {
          // This is a free surface node.
          cell->get_dof_indices(dofs);
          for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
            {
              surface_nodes(dofs[i]) = 1;
              other_nodes(dofs[i]) = 0;
              //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"  surfNodes: "<<surface_nodes(dofs[i])<<"  otherNodes: "<<other_nodes(dofs[i])<<std::endl;
            }
        }
      else
        {
          for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
            {
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
