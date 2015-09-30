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


// #include "tests.h"
// #include <deal.II/grid/tria.h>
// #include <deal.II/grid/grid_in.h>
//
// #include <deal.II/dofs/dof_handler.h>
// // #include <deal.II/dofs/dof_accessor.h>
//
// #include <deal.II/fe/fe_q.h>
// // #include <deal.II/fe/fe_values.h>
// #include <deal.II/fe/fe_system.h>
// #include <cmath>
// #include <iostream>
// #include <fstream>
// #include <string>
// #include <set>
//
// int main (int argc, char *argv[])
// {
//   initlog();
//   deallog<<"Check on the DoFHandler with codim 1 FESystem and 64 bits indices"<<std::endl;
//
//   Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);
#include "tests.h"
#include "computational_domain.h"
#include <deal.II/grid/grid_tools.h>

#include "bem_problem.h"
int main (int argc, char *argv[])
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);

  const unsigned int dim = 3;
  Triangulation<dim-1, dim> tria;
  std::ifstream in;
  in.open (SOURCE_DIR "/../utilities/coarse_sphere.inp");
  GridIn<dim-1, dim> gi;
  gi.attach_triangulation (tria);
  gi.read_ucd (in);

  FE_Q<dim-1, dim> fe_scalar(1);
  FESystem<dim-1, dim> fe_vector(FE_Q<dim-1, dim> (1), dim);

  DoFHandler<dim-1, dim> dh_scalar(tria);
  DoFHandler<dim-1, dim> dh_vector(tria);
  deallog<<"Distributing scalar dofs"<<std::endl;
  dh_scalar.distribute_dofs(fe_scalar);
  deallog<<"Distributing vector dofs"<<std::endl;
  dh_vector.distribute_dofs(fe_vector);
  deallog<<"Distributed dofs"<<std::endl;



}
