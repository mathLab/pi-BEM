// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2019 by the pi-BEM authors.
//
// This file is part of the pi-BEM library.
//
// The BEMStokes is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License version 2.1 as published by the Free Software Foundation.
// The full text of the license can be found in the file LICENSE at
// the top level of the BEMStokes distribution.
//
// Authors: Nicola Giuliani, Andrea Mola, Luca Heltai
//
// ---------------------------------------------------------------------

#include <deal.II/base/point.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/utilities.h>

// And here are a few C++ standard header
// files that we will need:
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <string>


namespace LaplaceKernel
{
  template <int dim>
  double
  single_layer(const Point<dim> &R)
  {
    switch (dim)
      {
        case 2:
          return (-std::log(R.norm()) / (2 * numbers::PI));

        case 3:
          return (1. / (R.norm() * 4 * numbers::PI));

        default:
          Assert(false, ExcInternalError());
          return 0.;
      }
  }



  template <int dim>
  Point<dim>
  double_layer(const Point<dim> &R)
  {
    switch (dim)
      {
        case 2:
          return R / (-2 * numbers::PI * R.square());
        case 3:
          return R / (-4 * numbers::PI * R.square() * R.norm());

        default:
          Assert(false, ExcInternalError());
          return Point<dim>();
      }
  }

  template <int dim>
  void
  kernels(const Tensor<1, dim> &R, Tensor<1, dim> &D, double &d)
  {
    double r  = R.norm();
    double r2 = r * r;
    switch (dim)
      {
        case 2:
          d = -std::log(r) / (2 * numbers::PI);
          D = R / (-2 * numbers::PI * r2);
          break;
        case 3:
          d = (1. / (r * 4 * numbers::PI));
          D = R / (-4 * numbers::PI * r2 * r);
          break;
        default:
          Assert(false, ExcInternalError());
      }
  }

} // namespace LaplaceKernel
