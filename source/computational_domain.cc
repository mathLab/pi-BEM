

#include "../include/computational_domain.h"
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal2lkit/utilities.h>

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

  for (unsigned int i=0; i<cad_surfaces.size(); ++i)
    {
      tria.set_manifold(1+i);
    }

  for (unsigned int i=0; i<cad_curves.size(); ++i)
    {
      tria.set_manifold(11+i);
    }

  tria.set_manifold(0);


}



template <int dim>
void ComputationalDomain<dim>::declare_parameters (ParameterHandler &prm)
{
  if (dim == 3)
    prm.declare_entry("Input grid name", "../grids/coarse_cube_double_nodes",
                      Patterns::Anything());
  else
    prm.declare_entry("Input grid name", "../grids/circle",
                      Patterns::Anything());


  prm.declare_entry("Input grid format", "inp",
                    Patterns::Anything());

  prm.declare_entry("Input path to CAD files", "./",
                    Patterns::Anything());

  prm.declare_entry("Number of cycles", "2",
                    Patterns::Integer());

  prm.declare_entry("Max aspect ratio", "3.5",
                    Patterns::Double());

  prm.declare_entry("Use iges surfaces and curves", "false",
                    Patterns::Bool());

  prm.declare_entry("Cad tolerance to projectors tolerance ratio", "100",
                    Patterns::Double());

  prm.declare_entry("Surface curvature adaptive refinement", "false",
                    Patterns::Bool());

  prm.declare_entry("Cells per circle", "12",
                    Patterns::Double());

  prm.declare_entry("Maximum number of curvature adaptive refinement cycles", "5",
                    Patterns::Integer());

  prm.declare_entry("Number of global refinement to be executed before local refinement cycle", "0",
                    Patterns::Integer());

  prm.enter_subsection("Boundary Conditions ID Numbers");
  {
    prm.declare_entry("Dirichlet boundary ids", "1,110,110", Patterns::List(Patterns::Integer(0)));
    prm.declare_entry("Neumann boundary ids", "0,110,110", Patterns::List(Patterns::Integer(0)));
  }
  prm.leave_subsection();


}

template <int dim>
void ComputationalDomain<dim>::parse_parameters (ParameterHandler &prm)
{

  input_grid_name = prm.get("Input grid name");
  input_grid_format = prm.get("Input grid format");
  input_cad_path = prm.get("Input path to CAD files");
  n_cycles = prm.get_integer("Number of cycles");
  max_element_aspect_ratio = prm.get_double("Max aspect ratio");
  use_cad_surface_and_curves = prm.get_bool("Use iges surfaces and curves");
  surface_curvature_refinement = prm.get_bool("Surface curvature adaptive refinement");
  cells_per_circle = prm.get_double("Cells per circle");
  pre_global_refinements = prm.get_integer("Number of global refinement to be executed before local refinement cycle");
  max_curvature_ref_cycles = prm.get_integer("Maximum number of curvature adaptive refinement cycles");
  cad_to_projectors_tolerance_ratio = prm.get_double("Cad tolerance to projectors tolerance ratio");

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

  if (input_grid_name == "../grids/coarse_sphere" || input_grid_name == "../grids/coarse_sphere_double_nodes" || input_grid_name == "../grids/circle" )
    {
      manifold = new SphericalManifold<dim-1, dim>;
      tria.set_all_manifold_ids(0);
      tria.set_manifold(0, *manifold);

    }

}

template <>
void ComputationalDomain<2>::create_initial_mesh()
{
  AssertThrow(true,
              ExcMessage("Create initial mesh only works in 3D"));

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

template <>
void ComputationalDomain<2>::refine_and_resize(const unsigned int refinement_level)
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
void ComputationalDomain<dim>::refine_and_resize(const unsigned int refinement_level)
{
  pcout<<"Refining and resizing mesh as required"<<std::endl;


  double max_tol=0;
  if (use_cad_surface_and_curves)
    {
      pcout<<"Color Files"<<endl;
      unsigned int ii=1;
      bool go_on = true;
      while (go_on == true)
        {
          std::string color_filename = ( input_cad_path + "Color_" +
                                         Utilities::int_to_string(ii) +
                                         ".iges" );
          ifstream f(color_filename);
          if (f.good())
            {
              pcout<<ii<<"-th file exists"<<endl;
              TopoDS_Shape surface = OpenCASCADE::read_IGES(color_filename, 1e-3);
              cad_surfaces.push_back(surface);
            }
          else
            go_on = false;
          ii++;
        }

      pcout<<"Edge Files"<<endl;
      ii=1;
      go_on = true;
      while (go_on == true)
        {
          std::string edge_filename = ( input_cad_path + "Curve_" +
                                        Utilities::int_to_string(ii) +
                                        ".iges" );
          ifstream f(edge_filename);
          if (f.good())
            {
              pcout<<ii<<"-th file exists"<<endl;
              TopoDS_Shape curve = OpenCASCADE::read_IGES(edge_filename, 1e-3);
              cad_curves.push_back(curve);
            }
          else
            go_on = false;
          ii++;
        }


      for (unsigned int i = 0; i<cad_surfaces.size(); ++i)
        {
          pcout<<i<<endl;
          max_tol = fmax(max_tol,OpenCASCADE::get_shape_tolerance(cad_surfaces[i]));
          pcout<<max_tol<<endl;
        }
      for (unsigned int i = 0; i<cad_curves.size(); ++i)
        {
          pcout<<i+cad_surfaces.size()<<endl;
          max_tol = fmax(max_tol,OpenCASCADE::get_shape_tolerance(cad_curves[i]));
          pcout<<max_tol<<endl;
        }

      const double tolerance = cad_to_projectors_tolerance_ratio*max_tol;




      pcout<<"Used tolerance is: "<<tolerance<<endl;
      for (unsigned int i=0; i<cad_surfaces.size(); ++i)
        {
          normal_to_mesh_projectors.push_back(SP(new OpenCASCADE::NormalToMeshProjectionBoundary<2,3>(cad_surfaces[i], tolerance)));
        }
      //static OpenCASCADE::DirectionalProjectionBoundary<2,3>
      //        directional_projector_lat(cad_surfaces[0], Point<3>(0.0,1.0,0.0), tolerance);
      //static OpenCASCADE::NormalProjectionBoundary<2,3>
      //        normal_projector_lat(cad_surfaces[0], tolerance);

      for (unsigned int i=0; i<cad_curves.size(); ++i)
        {
          line_projectors.push_back(SP(new OpenCASCADE::ArclengthProjectionLineManifold<2,3>(cad_curves[i], tolerance)));
        }

      for (unsigned int i=0; i<cad_surfaces.size(); ++i)
        {
          tria.set_manifold(1+i,*normal_to_mesh_projectors[i]);
        }

      for (unsigned int i=0; i<cad_curves.size(); ++i)
        {
          tria.set_manifold(11+i,*line_projectors[i]);
        }

    }


  unsigned int refinedCellCounter = 1;
  unsigned int cycles_counter = 0;
  // we repeat the aspect ratio refinement cycle until no cell has been
  // flagged for refinement, or until we reach a maximum of 10 cycles
  while ( (refinedCellCounter) && (cycles_counter < 10) )
    {
      // the refined cells counter is zeroed at the start of each cycle
      refinedCellCounter = 0;
      // we loop on the all the triangulation active cells
      Triangulation<2,3>::active_cell_iterator cell = tria.begin_active();
      Triangulation<2,3>::active_cell_iterator endc = tria.end();
      for ( ; cell!= endc; ++cell)
        {
          // the following lines determine if the cell is more elongated
          // in its 0 or 1 direction
          unsigned int max_extent_dim = 0;
          unsigned int min_extent_dim = 1;
          if (cell->extent_in_direction(0) < cell->extent_in_direction(1))
            {
              max_extent_dim = 1;
              min_extent_dim = 0;
            }
          // we compute the extent of the cell in its maximum and minimum elongation
          // direction respectively
          double min_extent = cell->extent_in_direction(min_extent_dim);
          double max_extent = cell->extent_in_direction(max_extent_dim);
          // if the aspect ratio exceeds the prescribed maximum value, the cell is refined
          if (max_extent > max_element_aspect_ratio*min_extent)
            {
              cell->set_refine_flag(RefinementCase<2>::cut_axis(max_extent_dim));
              refinedCellCounter++;
            }
        }
      // the number of cells refined in this cycle is reported before
      // proceeding with the next one
      pcout<<"Aspect Ratio Reduction Cycle: "<<cycles_counter<<" ("<<refinedCellCounter<<")"<<endl;
      tria.execute_coarsening_and_refinement();

      // the following commented lines are here for debug puroposes: if
      // something fails during the aspect ratio reduction cycles, they
      // should be uncommented, so that a mesh file per cycle can be
      // produced to document the evolution of the mesh through the
      // refinements. If the make_edges_conformal() function
      // is suspect in creating some error, the lines can also be
      // moved after the make_edges_conformal() function is called

      //std::string filename = ( "meshIntermediateResult_" +
      //         Utilities::int_to_string(int(round(cycles_counter))) +
      //         ".inp" );
      //std::ofstream logfile(filename.c_str());
      //GridOut grid_out;
      //grid_out.write_ucd(tria, logfile);

      make_edges_conformal();
      cycles_counter++;

    }

  // the following refinement cycle is based upon the original CAD
  // geometry curvature. For this reason it can be activated not only
  // when the user requires it with the surface_curvature_refinement
  // option in the input file. Of course, this is possible only if
  // also the use_cad_surface_and_curves flag is set to thrue through
  // the input file. Only in this way in fact, CAD surfaces and curves
  // are prescribed for the triangulation refinements on some of
  // its manifold ids.
  if (use_cad_surface_and_curves && surface_curvature_refinement)
    {
      const double tolerance = cad_to_projectors_tolerance_ratio*max_tol;
      refinedCellCounter = 1;
      cycles_counter = 0;
      // the refinement procedure is recursively repeated until no more cells are flagged
      // for refinement, or until the user specified maximum number of curvature
      // refinement cycles is reached
      while ( (refinedCellCounter) && (cycles_counter < max_curvature_ref_cycles) )
        {
          // the refined cells counter is zeroed at the start of each cycle
          refinedCellCounter = 0;
          // we loop on the all the triangulation active cells
          Triangulation<2,3>::active_cell_iterator cell = tria.begin_active();
          Triangulation<2,3>::active_cell_iterator endc = tria.end();
          for ( ; cell!= endc; ++cell)
            {
              // In the following lines, we try to come up with an estimation
              // of the cell normal. It is obtained from the average of the normal
              // to the 4 triangles in which the cell can be split using the vertices
              // and the center. The commented lines can be used for checks
              // in case something goes wrong.

              //cout<<"center: "<<cell->center()<<endl;
              //cout<<"v0: "<<cell->vertex(0)<<endl;
              //cout<<"v1: "<<cell->vertex(1)<<endl;
              //cout<<"v2: "<<cell->vertex(2)<<endl;
              //cout<<"v3: "<<cell->vertex(3)<<endl;
              Point<3> t0 = cell->vertex(0)+(-1.0)*cell->center();
              Point<3> t1 = cell->vertex(1)+(-1.0)*cell->center();
              Point<3> t2 = cell->vertex(2)+(-1.0)*cell->center();
              Point<3> t3 = cell->vertex(3)+(-1.0)*cell->center();
              //cout<<"t0: "<<t0<<endl;
              //cout<<"t1: "<<t1<<endl;
              //cout<<"t2: "<<t2<<endl;
              //cout<<"t3: "<<t3<<endl;

              Point<3> nn0(t0(1)*t1(2)-t0(2)*t1(1),
                           t0(2)*t1(0)-t0(0)*t1(2),
                           t0(0)*t1(1)-t0(1)*t1(0));
              nn0/=nn0.norm();
              Point<3> nn1(t1(1)*t3(2)-t1(2)*t3(1),
                           t1(2)*t3(0)-t1(0)*t3(2),
                           t1(0)*t3(1)-t1(1)*t3(0));
              nn1/=nn1.norm();
              Point<3> nn2(t3(1)*t2(2)-t3(2)*t2(1),
                           t3(2)*t2(0)-t3(0)*t2(2),
                           t3(0)*t2(1)-t3(1)*t2(0));
              nn2/=nn2.norm();
              Point<3> nn3(t2(1)*t0(2)-t2(2)*t0(1),
                           t2(2)*t0(0)-t2(0)*t0(2),
                           t2(0)*t0(1)-t2(1)*t0(0));
              nn3/=nn3.norm();
              Point<3> n = (nn0+nn1+nn2+nn3)/4.0;
              n/=n.norm();
              //cout<<cell<<endl;
              //cout<<nn0<<endl;
              //cout<<nn1<<endl;
              //cout<<nn2<<endl;
              //cout<<nn3<<endl;
              //cout<<n<<endl;
              //cout<<cell<<"  material id: "<<int(cell->material_id())<<endl;

              // once the cell normal has beed created, we want to use it as the
              // direction of the projection onto the CAD surface
              // first though, let's check that we are using a CAD surface for
              // the refinement of the manifold_id associated with the present
              // cell
              double cell_size;
              if (int(cell->material_id())-1 < cad_surfaces.size())
                {
                  // if so, the cad_surface associated with the present manifold_id is identified...
                  TopoDS_Shape neededShape = cad_surfaces[int(cell->material_id())-1];
                  // ...and used to set up a line intersection to project the cell center
                  // on the CAD surface along the direction specified by the previously computed
                  // cell normal
                  Point<3> projection = OpenCASCADE::line_intersection(neededShape,
                                                                       cell->center(),
                                                                       n,
                                                                       tolerance);
                  // in correspondence with the projected point, we ask all the surface differential forms
                  std_cxx11::tuple<Point<3>,Tensor<1,3>,double,double> tup = OpenCASCADE::closest_point_and_differential_forms(neededShape,
                                                                             projection,
                                                                             tolerance);
                  // among the differential point, we select the maximum absolute curvature
                  double max_abs_curv = fmax(fabs(std_cxx11::get<2>(tup)),fabs(std_cxx11::get<3>(tup)));
                  // this commented line is just for debug purposes
                  // cout<<"Point: "<<std_cxx11::get<0>(tup)<<"  Kmin: "<<std_cxx11::get<2>(tup)<<"  Kmax: "<<std_cxx11::get<3>(tup)<<endl;
                  // the minimum curvature radius is computed from the maximum absolute curvatur
                  double curvature_radius = 1.0/fmax(max_abs_curv,tolerance);
                  // the target cell size is selected so that it corresponds to a cells_per_circle
                  // fraction of the circumference corresponding to the minimum curvature radius
                  cell_size = 2*dealii::numbers::PI/cells_per_circle*curvature_radius;
                }
              else
                {
                  // if the cell manifold_id is not associated to a CAD surface, the
                  // target cell_size is set to and extremely high value, so that the cell is
                  // never refined
                  cell_size = 2*dealii::numbers::PI/cells_per_circle/tolerance;
                }

              // the following line si for debug puropses and should be uncommented if
              // something is not working with the refinement
              //cout<<"Cell Diam: "<<cell->diameter()<<"  Target Cell Size: "<<cell_size<<endl;


              // if the cell diameter is higher than the target cell size, the refinement flag is set
              // (unless the cell is already very small ---which for us means 10xtolerance)
              if ( (cell->diameter() > cell_size)  &&
                   (cell->diameter() > 10*tolerance) )
                {
                  cell->set_refine_flag();
                  refinedCellCounter++;
                }
            }
          // the number of cells to be refined in this cycle is reported, the refinement is carried out
          // and the make_edges_conformal function is called to check no edge presents non comformities
          pcout<<"Curvature Based Local Refinement Cycle: "<<cycles_counter<<" ("<<refinedCellCounter<<")"<<endl;
          tria.execute_coarsening_and_refinement();
          make_edges_conformal();
          cycles_counter++;

          //std::string filename = ( "DTMB_II_meshResult_max_curv" +
          //                         Utilities::int_to_string(int(round(cycles_counter))) +
          //                         ".vtk" );
          //std::ofstream logfile(filename.c_str());
          //GridOut grid_out;
          //grid_out.write_vtk(tria, logfile);
          //std::string stl_filename = ( "DTMB_II_meshResult_max_curv" +
          //                         Utilities::int_to_string(int(round(cycles_counter))) +
          //                         ".stl" );
          //SaveSTL(tria,stl_filename);



        }
    }
//*/


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
  data_out.add_data_vector (partitioning, "partitioning", DataOut<dim-1, DoFHandler<dim-1, dim> >::type_dof_data);
  data_out.build_patches ();
  data_out.write_vtu (output);


  pcout<<"...done refining and resizing mesh"<<std::endl;


}

template<>
void ComputationalDomain<2>::make_edges_conformal(const bool with_double_nodes,
                                                  const bool isotropic_ref_on_opposite_side)
{
}

template<int dim>
void ComputationalDomain<dim>::make_edges_conformal(const bool with_double_nodes,
                                                    const bool isotropic_ref_on_opposite_side)
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
                              if ( isotropic_ref_on_opposite_side )
                                {
                                  (*jt)->set_refine_flag();
                                  to_restore=true;
                                }
                              // otherwise, use anisotropic refinement to make edge mesh conformal
                              else
                                {
                                  if ((d==0) || (d==1))
                                    (*jt)->set_refine_flag(RefinementCase<2>::cut_axis(1));
                                  else
                                    (*jt)->set_refine_flag(RefinementCase<2>::cut_axis(0));
                                  to_restore=true;
                                }
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



template class ComputationalDomain<2>;
template class ComputationalDomain<3>;
