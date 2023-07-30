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


#include "singular_kernel_integral.h"
#include "tests.h"

int
main()
{
  initlog();

  // we start selecting a mapping degree
  double mapping_degree = 1;

  // and now we generate a grid and use it to distribute dofs
  Triangulation<2, 3> triangulation;

  std::vector<Point<3>>    vertices;
  std::vector<CellData<2>> cells;
  SubCellData              subcelldata;

  // Flat

  vertices.resize(4);
  vertices[0](0) = -1.0;
  vertices[0](1) = -1.0;
  vertices[0](2) = 0.0;
  vertices[1](0) = 1.5;
  vertices[1](1) = -1.0;
  vertices[1](2) = 0.0;
  vertices[2](0) = -1.0;
  vertices[2](1) = 1.0;
  vertices[2](2) = 0.0;
  vertices[3](0) = 0.5;
  vertices[3](1) = 1.0;
  vertices[3](2) = 0.0;


  cells.resize(1);

  cells[0].vertices[0] = 0;
  cells[0].vertices[1] = 1;
  cells[0].vertices[2] = 2;
  cells[0].vertices[3] = 3;

  GridTools::delete_unused_vertices(vertices, cells, subcelldata);
  GridTools::consistently_order_cells(cells);

  // here's the triangulation set up
  triangulation.create_triangulation(vertices, cells, subcelldata);


  // we will need a codimension one finite element and
  // the corresponding dof handler
  // we use zeroth order fe_dgq because in the integral the hypersingular kernel
  // must not be multiplied by the shape function. Thus, we multiply it by a
  // constant shape function
  FE_DGQ<2, 3>     fe(0);
  DoFHandler<2, 3> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  MappingQ<2, 3> mapping(mapping_degree);

  // The test cases are taken from
  // M. Guiggiani, G. Krishnasamy, T. J. Rudolphi & F. J. Rizzo,A general
  // algorithm for the numerical solution of hypersingular boundary integral
  // equations, JOURNAL OF APPLIED MECHANICS-TRANSACTIONS OF THE ASME,vol. 59,pp
  // 604-614,tot.pag 11,1992

  // Case 4.1a
  Point<2> singularity_location_in_parametric_plane_1(0.5, 0.5);
  // 4.1b
  Point<2> singularity_location_in_parametric_plane_2(0.66 * .5 + .5, 0.5);
  // 4.1c
  Point<2> singularity_location_in_parametric_plane_3(0.885764071856287,
                                                      0.66 * .5 + .5);

  std::vector<Tensor<1, 3>> I;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      deallog << "Case: 4.1a" << std::endl;
      SingularKernelIntegral sing_kernel_integrator_1(
        cell, fe, mapping, singularity_location_in_parametric_plane_1);
      I = sing_kernel_integrator_1.evaluate_VkNj_integrals();
      // std::vector<Tensor<1,3> > I =
      // sing_kernel_integrator.evaluate_WkNj_integrals();
      deallog << std::setprecision(15);
      deallog << "Result : " << I[0][2] * 4 * dealii::numbers::PI
              << "    Exact: " << -5.749236751228080 << std::endl;
      deallog << "Abs Err: "
              << -5.749236751228080 - I[0][2] * 4 * dealii::numbers::PI
              << std::endl;
      deallog << "Rel Err: "
              << fabs(-5.749236751228080 - I[0][2] * 4 * dealii::numbers::PI) /
                   fabs(-5.749236751228080)
              << std::endl;


      deallog << "Case: 4.1b" << std::endl;
      SingularKernelIntegral sing_kernel_integrator_2(
        cell, fe, mapping, singularity_location_in_parametric_plane_2);
      I = sing_kernel_integrator_2.evaluate_VkNj_integrals();
      // std::vector<Tensor<1,3> > I =
      // sing_kernel_integrator.evaluate_WkNj_integrals();
      deallog << std::setprecision(15);
      deallog << "Result : " << I[0][2] * 4 * dealii::numbers::PI
              << "    Exact: " << -9.154585469918885 << std::endl;
      deallog << "Abs Err: "
              << -9.154585469918885 - I[0][2] * 4 * dealii::numbers::PI
              << std::endl;
      deallog << "Rel Err: "
              << fabs(-9.154585469918885 - I[0][2] * 4 * dealii::numbers::PI) /
                   fabs(-9.154585469918885)
              << std::endl;


      deallog << "Case: 4.1c" << std::endl;
      SingularKernelIntegral sing_kernel_integrator_3(
        cell, fe, mapping, singularity_location_in_parametric_plane_3);
      I = sing_kernel_integrator_3.evaluate_VkNj_integrals();
      // std::vector<Tensor<1,3> > I =
      // sing_kernel_integrator.evaluate_WkNj_integrals();
      deallog << std::setprecision(15);
      deallog << "Result : " << I[0][2] * 4 * dealii::numbers::PI
              << "    Exact: " << -15.32849545090306 << std::endl;
      deallog << "Abs Err: "
              << -15.32849545090306 - I[0][2] * 4 * dealii::numbers::PI
              << std::endl;
      deallog << "Rel Err: "
              << fabs(-15.32849545090306 - I[0][2] * 4 * dealii::numbers::PI) /
                   fabs(-15.32849545090306)
              << std::endl;
    }
}
