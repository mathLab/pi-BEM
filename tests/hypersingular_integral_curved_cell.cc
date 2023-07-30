//-----------------------------------------------------------
//
//    Copyright (C) 2014 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>
#include "singular_kernel_integral.h"
#include "tests.h"

int
main()
{
  initlog();

  //we start selecting a mapping degree
  double mapping_degree = 2;
  //we also import a CAD geometry of a portion of a cylinder to integrate on curved cell
  std::string cad_filename = std::string(SOURCE_DIR)+"/testData/Revolution_1.iges";
  TopoDS_Shape surface = OpenCASCADE::read_IGES(cad_filename, 1e-3);
  OpenCASCADE::NormalToMeshProjectionManifold<2,3> manifold(surface, 1e-7);
  // the next lines are in case we want to run the test with NURBSPatchManifold
  // and --- possibly --- with MappingManifold. But there is some work to do to make it work
//  std::vector<TopoDS_Face> occ_faces;
//  std::vector<TopoDS_Edge> occ_edges;
//  std::vector<TopoDS_Vertex> occ_vertices;
//  OpenCASCADE::extract_geometrical_shapes(surface, occ_faces, occ_edges, occ_vertices);
//  OpenCASCADE::NURBSPatchManifold<2,3> manifold(occ_faces[0], 1e-7);
 
  
  // and now we generate a grid and use it to distribute dofs
  Triangulation<2,3> triangulation;

  std::vector<Point<3> > vertices;
  std::vector<CellData<2> > cells;
  SubCellData subcelldata;

// Flat

  vertices.resize(4);
  vertices[0](0)=1.0;
  vertices[0](1)=0.0;
  vertices[0](2)=0.0;
  vertices[1](0)=1.0;
  vertices[1](1)=2.0;
  vertices[1](2)=0.0;
  vertices[2](0)=0.0;
  vertices[2](1)=0.0;
  vertices[2](2)=1.0;
  vertices[3](0)=0.0;
  vertices[3](1)=2.0;
  vertices[3](2)=1.0;


  cells.resize(1);

  cells[0].vertices[0]=0;
  cells[0].vertices[1]=1;
  cells[0].vertices[2]=2;
  cells[0].vertices[3]=3;
  
  GridTools::delete_unused_vertices (vertices, cells, subcelldata);
  GridTools::consistently_order_cells(cells);

  // here's the triangulation set up
  triangulation.create_triangulation(vertices, cells, subcelldata );

  triangulation.set_all_manifold_ids(0);
  triangulation.set_manifold(0, manifold);



  // we will need a codimension one finite element and
  // the corresponding dof handler
  // we use zeroth order fe_dgq because in the integral the hypersingular kernel must not be multiplied by the shape function. Thus,
  // we multiply it by a constant shape function
  FE_DGQ<2,3> fe(0);
  DoFHandler<2,3> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  MappingQ<2,3> mapping(mapping_degree);

  // The test cases are taken from
  //M. Guiggiani, G. Krishnasamy, T. J. Rudolphi & F. J. Rizzo,A general algorithm for the numerical solution of hypersingular boundary integral equations, JOURNAL OF APPLIED MECHANICS-TRANSACTIONS OF THE ASME,vol. 59,pp 604-614,tot.pag 11,1992

  // Case 4.2a
  Point<2> singularity_location_in_parametric_plane_1(0.5,0.5);
  // Case 4.2b
  Point<2> singularity_location_in_parametric_plane_2(0.66*.5+.5,0.0*.5+.5);
  // Case 4.2c
  Point<2> singularity_location_in_parametric_plane_3(0.66*.5+.5,0.66*.5+.5);
  

  std::vector<Tensor<1,3> > I;

  for (const auto &cell : dof_handler.active_cell_iterators())
      {
      deallog<<"Case: 4.1a"<<std::endl;
      SingularKernelIntegral sing_kernel_integrator_1(cell, fe, mapping, singularity_location_in_parametric_plane_1);
      I = sing_kernel_integrator_1.evaluate_VkNj_integrals();
      //std::vector<Tensor<1,3> > I = sing_kernel_integrator.evaluate_WkNj_integrals();
      deallog<<std::setprecision(15);
      deallog<<"Result : "<<I[0][2]*4*dealii::numbers::PI<<"    Exact: "<<-0.343807*4*dealii::numbers::PI<<std::endl;
      deallog<<"Abs Err: "<<-0.343807*4*dealii::numbers::PI-I[0][2]*4*dealii::numbers::PI<<std::endl;
      deallog<<"Rel Err: "<<fabs(-0.343807*4*dealii::numbers::PI-I[0][2]*4*dealii::numbers::PI)/fabs(-0.343807*4*dealii::numbers::PI)<<std::endl;
      
      
      deallog<<"Case: 4.1b"<<std::endl;
      SingularKernelIntegral sing_kernel_integrator_2(cell, fe, mapping, singularity_location_in_parametric_plane_2);
      I = sing_kernel_integrator_2.evaluate_VkNj_integrals();
      //std::vector<Tensor<1,3> > I = sing_kernel_integrator.evaluate_WkNj_integrals();
      deallog<<std::setprecision(15);
      deallog<<"Result : "<<I[0][2]*4*dealii::numbers::PI<<"    Exact: "<<-0.497099*4*dealii::numbers::PI<<std::endl;
      deallog<<"Abs Err: "<<-0.497099*4*dealii::numbers::PI-I[0][2]*4*dealii::numbers::PI<<std::endl;
      deallog<<"Rel Err: "<<fabs(-0.497099*4*dealii::numbers::PI-I[0][2]*4*dealii::numbers::PI)/fabs(-0.497099*4*dealii::numbers::PI)<<std::endl;


      deallog<<"Case: 4.1c"<<std::endl;
      SingularKernelIntegral sing_kernel_integrator_3(cell, fe, mapping, singularity_location_in_parametric_plane_3);
      I = sing_kernel_integrator_3.evaluate_VkNj_integrals();
      //std::vector<Tensor<1,3> > I = sing_kernel_integrator.evaluate_WkNj_integrals();
      deallog<<std::setprecision(15);
      deallog<<"Result : "<<I[0][2]*4*dealii::numbers::PI<<"    Exact: "<<-0.877214*4*dealii::numbers::PI<<std::endl;
      deallog<<"Abs Err: "<<-0.877214*4*dealii::numbers::PI-I[0][2]*4*dealii::numbers::PI<<std::endl;
      deallog<<"Rel Err: "<<fabs(-0.877214*4*dealii::numbers::PI-I[0][2]*4*dealii::numbers::PI)/fabs(-0.877214*4*dealii::numbers::PI)<<std::endl;



//    strong integrals case a
//    deallog<<"Abs Err: "<<-2.114175-I[0][0]*4*dealii::numbers::PI<<std::endl;
//    deallog<<"Rel Err: "<<fabs(-2.114175-I[0][0]*4*dealii::numbers::PI)/fabs(-2.114175)<<std::endl;


//    strong integrals case b
//    deallog<<"Abs Err: "<<-1.935711-I[0][0]*4*dealii::numbers::PI<<std::endl;
//    deallog<<"Rel Err: "<<fabs(-1.935711-I[0][0]*4*dealii::numbers::PI)/fabs(-1.935711)<<std::endl;

//    strong integrals case c
//    deallog<<"Abs Err: "<<.8790179-I[0][0]*4*dealii::numbers::PI<<std::endl;
//    deallog<<"Rel Err: "<<fabs(.8790179-I[0][0]*4*dealii::numbers::PI)/fabs(.8790179)<<std::endl;




      }
  


}
