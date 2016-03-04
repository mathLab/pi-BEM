

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
ComputationalDomain<dim>::ComputationalDomain(MPI_Comm comm)
  :
  mpi_communicator (comm),
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
  tria.set_all_manifold_ids(0);
  tria.set_manifold(0);
  //
  // if(manifold != NULL)
  //   delete manifold;

}



template <int dim>
void ComputationalDomain<dim>::declare_parameters (ParameterHandler &prm)
{

  prm.declare_entry("Input grid name", "../utilities/coarse_cube_double_nodes",
                    Patterns::Anything());

  prm.declare_entry("Input grid format", "inp",
                    Patterns::Anything());

  prm.declare_entry("Number of cycles", "2",
                    Patterns::Integer());

  prm.declare_entry("Number of mg cycles", "1",
                    Patterns::Integer());

  prm.enter_subsection("Boundary Conditions ID Numbers");
  {
    prm.declare_entry("Dirichlet boundary ids", "1,110,110", Patterns::List(Patterns::Integer(0)));
    prm.declare_entry("Neumann boundary ids", "0,110,110", Patterns::List(Patterns::Integer(0)));
    // prm.declare_entry("Dirichlet Surface 1 ID", "1", Patterns::Integer());
    // prm.declare_entry("Dirichlet Surface 2 ID", "110", Patterns::Integer());
    // prm.declare_entry("Dirichlet Surface 3 ID", "110", Patterns::Integer());
    // prm.declare_entry("Neumann Surface 1 ID", "0", Patterns::Integer());
    // prm.declare_entry("Neumann Surface 2 ID", "110", Patterns::Integer());
    // prm.declare_entry("Neumann Surface 3 ID", "110", Patterns::Integer());
  }
  prm.leave_subsection();


}

template <int dim>
void ComputationalDomain<dim>::parse_parameters (ParameterHandler &prm)
{

  input_grid_name = prm.get("Input grid name");
  input_grid_format = prm.get("Input grid format");
  n_cycles = prm.get_integer("Number of cycles");
  mg_cycles = prm.get_integer("Number of mg cycles");




  prm.enter_subsection("Boundary Conditions ID Numbers");
  {
    std::vector<std::string> dirichlet_string_list = Utilities::split_string_list(prm.get("Dirichlet boundary ids"));
    dirichlet_boundary_ids.resize(dirichlet_string_list.size());
    for (unsigned int i=0; i<dirichlet_string_list.size(); ++i)
      {
        std::istringstream reader(dirichlet_string_list[i]);
        reader >> dirichlet_boundary_ids[i];
      }

    std::vector<std::string> neumann_string_list = Utilities::split_string_list(prm.get("Neumann boundary ids"));
    neumann_boundary_ids.resize(neumann_string_list.size());
    for (unsigned int i=0; i<neumann_string_list.size(); ++i)
      {
        std::istringstream reader(neumann_string_list[i]);
        reader >> neumann_boundary_ids[i];
      }

    // dirichlet_sur_ID1 = prm.get_integer("Dirichlet Surface 1 ID");
    // dirichlet_sur_ID2 = prm.get_integer("Dirichlet Surface 2 ID");
    // dirichlet_sur_ID3 = prm.get_integer("Dirichlet Surface 3 ID");
    // neumann_sur_ID1 = prm.get_integer("Neumann Surface 1 ID");
    // neumann_sur_ID2 = prm.get_integer("Neumann Surface 2 ID");
    // neumann_sur_ID3 = prm.get_integer("Neumann Surface 3 ID");
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
  in.open (input_grid_name + "." + input_grid_format);
  GridIn<dim-1, dim> gi;
  gi.attach_triangulation (tria);
  if (input_grid_format=="vtk")
    gi.read_vtk (in);
  else if (input_grid_format=="msh")
    gi.read_msh (in);
  else if (input_grid_format=="inp")
    gi.read_ucd (in);
  else
    Assert (false, ExcNotImplemented());

  //
  // std::ifstream in;
  // switch (dim)
  //   {
  //   case 2:
  //     in.open ("coarse_circle.inp");
  //     break;
  //
  //   case 3:
  //     in.open ("coarse_cube_double_nodes.inp");
  //     break;
  //
  //   default:
  //     Assert (false, ExcNotImplemented());
  //   }
  //
  // GridIn<dim-1, dim> gi;
  // gi.attach_triangulation (tria);
  // gi.read_ucd (in);
  //
  // manifold = new SphericalManifold<dim-1, dim>;

  if (input_grid_name == "../utilities/coarse_sphere" || input_grid_name == "../utilities/coarse_sphere_double_nodes" )
    {
      manifold = new SphericalManifold<dim-1, dim>;
      tria.set_all_manifold_ids(0);
      tria.set_manifold(0, *manifold);

    }

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
  pcout<<"We have a tria of "<<tria.n_active_cells()<<" cells."<<std::endl;
  GridTools::partition_triangulation(n_mpi_processes, tria);
  std::string filename0 = ( "meshResult.inp" );
  std::ofstream logfile0(filename0.c_str());
  GridOut grid_out0;
  grid_out0.write_ucd(tria, logfile0);


  pcout<<"...done refining and resizing mesh"<<std::endl;
}


template <int dim>
void ComputationalDomain<dim>::conditional_refine_and_resize(const unsigned int refinement_level)
{
  pcout<<"Conditionally refining and resizing mesh as required"<<std::endl;

  const Point<dim> center (0,0,0);
  compute_double_vertex_cache();
  make_edges_conformal();

  for (unsigned int step=0; step < refinement_level; ++step)
    {
      auto cell = tria.begin_active();
      auto endc = tria.end();
      for (; cell!=endc; ++cell)
        {
          for (unsigned int v=0;
               v < GeometryInfo<dim-1>::vertices_per_cell;
               ++v)
            {
              const double distance_from_center
                = center.distance (cell->vertex(v));
              if (std::fabs(distance_from_center) < 1.)
                {
                  cell->set_refine_flag ();
                  break;
                }
            }
        }
      tria.prepare_coarsening_and_refinement ();
      tria.execute_coarsening_and_refinement ();
      // compute_double_vertex_cache();
      make_edges_conformal();
    }
  update_triangulation();
}

template<int dim>
void ComputationalDomain<dim>::update_triangulation()
{
  // compute_double_vertex_cache();
  make_edges_conformal();
  // tria.execute_coarsening_and_refinement ();

  pcout<<"We have a tria of "<<tria.n_active_cells()<<" cells."<<std::endl;
  GridTools::partition_triangulation(n_mpi_processes, tria);
  std::string filename0 = ( "meshResult.inp" );
  std::ofstream logfile0(filename0.c_str());
  GridOut grid_out0;
  grid_out0.write_ucd(tria, logfile0);

  std::ostringstream filename;
  filename << "mesh.vtu";
  std::ofstream output (filename.str().c_str());

  FE_Q<dim-1, dim> fe_dummy(1);

  DoFHandler<dim-1,dim> dof_handler(tria);
  dof_handler.distribute_dofs (fe_dummy);
  DataOut<dim-1,DoFHandler<dim-1, dim>> data_out;
  data_out.attach_dof_handler (dof_handler);
  std::vector<unsigned int> partition_int (tria.n_active_cells());
  GridTools::get_subdomain_association (tria, partition_int);
  const Vector<double> partitioning(partition_int.begin(),
                                    partition_int.end());
  data_out.add_data_vector (partitioning, "partitioning");
  data_out.build_patches ();
  data_out.write_vtu (output);


  pcout<<"...done refining and resizing mesh"<<std::endl;


}


template<int dim>
void ComputationalDomain<dim>::make_edges_conformal(const bool with_double_nodes)
{
  if (with_double_nodes==false)
    {
      auto cell = tria.begin_active();
      auto endc = tria.end();

      for (cell=tria.begin_active(); cell!=endc; ++cell)
        {
          for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
            if (cell->face(f)->at_boundary()) //material_id()!=numbers::invalid_material_id)//dovrei essere su un buondary
              {
                //TriaIterator<CellAccessor<dim-1,dim> > cell_neigh = cell->neighbor(f);
                if (cell->neighbor_is_coarser(f))
                  {
                    TriaIterator<CellAccessor<dim-1,dim> > cell_neigh = cell->neighbor(f);
                    cell_neigh->set_refine_flag(RefinementCase<dim-1>::isotropic_refinement);
                    //std::cout<<"mammina..."<<std::endl;

                  }
              }

        }
      tria.prepare_coarsening_and_refinement();
      tria.execute_coarsening_and_refinement();
    }
  else
    {

      pcout<<"Restoring mesh conformity on edges..."<<std::endl;
      pcout<<"cells before : "<<tria.n_active_cells()<<std::endl;
      //pcout<<"dofs before: "<<dhh.n_dofs()<<std::endl;
      bool to_restore = true;

      while (to_restore)
        {
          compute_double_vertex_cache();
          types::global_dof_index n_vertex=tria.n_vertices();
          auto all_vertices=tria.get_vertices();

          double tol=1e-7;

          to_restore=false;
          pcout<<to_restore<<std::endl;
          for (types::global_dof_index i=0; i<n_vertex; ++i)
            {
              if (vertex_on_boundary[i]==true && double_vertex_vector[i].size()==1)
                {
                  std::vector<Point<dim> > nodes(GeometryInfo<dim-1>::vertices_per_face);
                  for (types::global_dof_index kk=0; kk<vert_to_elems[i].size(); ++kk) //ogni faccia ha due estremi
                    {
                      auto cell = vert_to_elems[i][kk];//mi riconduco alla cella con il nodo non conforme
                      for (unsigned int f=0; f<GeometryInfo<dim-1>::faces_per_cell; ++f)
                        {
                          if (cell->face(f)->at_boundary())// ritrovo la faccia con l'edge non doppio
                            {
                              //std::cout<<cell->face(f)->vertex(1)<<"  "<<cell->face(f)->vertex(0)<<std::endl;
                              if (all_vertices[i].distance(cell->face(f)->vertex(0)) <tol)
                                nodes[kk] = cell->face(f)->vertex(1);
                              else if (all_vertices[i].distance(cell->face(f)->vertex(1)) <tol)
                                nodes[kk] = cell->face(f)->vertex(0);
                            }
                        }
                      //std::cout<<std::endl;
                    }
                  //std::cout<<nodes[0]<<"    "<<ref_points[i]<<"   "<<nodes[1]<<std::endl;
                  // we can now compute the center of the parent cell face
                  Point<3> parent_face_center = 0.5*(nodes[0]+nodes[1]);
                  for (auto jt=edge_cells.begin(); jt != edge_cells.end(); ++jt)
                    for (unsigned int d=0; d<GeometryInfo<2>::faces_per_cell; ++d)
                      if ((*jt)->face(d)->at_boundary())
                        {
                          //cout<<parent_face_center.distance((*jt)->face(d)->center())<<" "<<tol<<endl;
                          if ( parent_face_center.distance(((*jt)->face(d)->vertex(0)+(*jt)->face(d)->vertex(1))/2) < tol)
                            {
                              (*jt)->set_refine_flag();
                              to_restore=true;
                            }
                        }
                }

            }

          if (to_restore)
            {
              tria.prepare_coarsening_and_refinement();
              tria.execute_coarsening_and_refinement();
              pcout<<"found non conformity, new cell number : "<<tria.n_active_cells()<<std::endl;

            }
          // pcout<<"pippo"<<std::endl;

        }
      pcout<<"cells after : "<<tria.n_active_cells()<<std::endl;
      pcout<<"...Done restoring mesh conformity"<<std::endl;
    }

}

template <int dim>
void ComputationalDomain<dim>::compute_double_vertex_cache()
{
  pcout<<"Computing infos for the double_vertex"<<std::endl;
  double toll=1e-7;
  double_vertex_vector.clear();
  types::global_dof_index n_vertex=tria.n_vertices();
  double_vertex_vector.resize(n_vertex);
  vertex_on_boundary.resize(n_vertex);
  std::fill(vertex_on_boundary.begin(), vertex_on_boundary.end(),false);

  auto all_vertices=tria.get_vertices();

  for (types::global_dof_index i = 0; i<n_vertex; ++i)
    {
      for (types::global_dof_index j = 0; j<n_vertex; ++j)
        {
          if (all_vertices[i].distance(all_vertices[j])<=toll)
            {
              double_vertex_vector[i].push_back(j);
            }
        }
    }



  auto cell = tria.begin_active();
  auto endc = tria.end();
  vert_to_elems.clear();
  edge_cells.clear();

  for (cell=tria.begin_active(); cell!=endc; ++cell)
    {
      std::vector<Point<dim> > cell_vertices(GeometryInfo<dim-1>::vertices_per_cell);

      for (unsigned int v=0;
           v < GeometryInfo<dim-1>::vertices_per_cell;
           ++v)
        {
          vert_to_elems[cell->vertex_index(v)].push_back(cell);
          cell_vertices[v]=cell->vertex(v);
        }


      if (cell->at_boundary())
        {
          edge_cells.insert(cell);

          for (unsigned int f=0;
               f < GeometryInfo<dim-1>::faces_per_cell;
               ++f)
            {
              if (cell->face(f)->at_boundary())
                {
                  for (unsigned int v=0;
                       v < GeometryInfo<dim-1>::vertices_per_cell;
                       ++v)
                    {
                      if (cell->face(f)->vertex(0)==cell_vertices[v])
                        vertex_on_boundary[cell->vertex_index(v)]=true;
                      else if (cell->face(f)->vertex(1)==cell_vertices[v])
                        vertex_on_boundary[cell->vertex_index(v)]=true;
                    }

                }


            }



        }
    }
  pcout<<"done double_vertex cache"<<std::endl;
}



template class ComputationalDomain<3>;
