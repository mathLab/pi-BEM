// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#ifndef dealii_occ_my_utilities_h
#define dealii_occ_my_utilities_h

#  include <deal.II/opencascade/utilities.h>
// opencascade needs "HAVE_CONFIG_H" to be exported...
#  define HAVE_CONFIG_H
#  include <Adaptor3d_Curve.hxx>
#  include <Adaptor3d_HCurve.hxx>
#  include <BRepAdaptor_Curve.hxx>
#  include <Bnd_Box.hxx>
#  undef HAVE_CONFIG_H


DEAL_II_NAMESPACE_OPEN

using namespace OpenCASCADE;





  template <int dim>
  std::tuple<Point<dim>, TopoDS_Shape, double, double>
  my_project_point_and_pull_back(const TopoDS_Shape &in_shape,
                              const Point<dim> &  origin,
                              const double        tolerance);



  std::tuple<Point<3>, Tensor<1, 3>, double, double>
  my_closest_point_and_differential_forms(const TopoDS_Shape &in_shape,
                                       const Point<3> &    origin,
                                       const double        tolerance);


  template <int dim>
  Point<dim>
  my_line_intersection(const TopoDS_Shape &  in_shape,
                       const Point<dim> &    origin,
                       const Tensor<1, dim> &direction,
                       const double          tolerance);

DEAL_II_NAMESPACE_CLOSE


#endif // dealii_occ_my_utilities_h
