// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the pi-BEM authors.
//
// This file is part of the pi-BEM library.
//
// The pi-BEM is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License version 2.1 as published by the Free Software Foundation.
// The full text of the license can be found in the file LICENSE at
// the top level of the pi-BEM distribution.
//
// Authors: Nicola Giuliani, Andrea Mola, Luca Heltai

#ifndef MULTIPOLE_EXPANSION_H_
#define MULTIPOLE_EXPANSION_H_

#include <deal.II/base/point.h>

#include <math.h>

#include <complex>
#include <string>
#include <vector>

#include "ass_leg_function.h"

using namespace dealii;

class MultipoleExpansion
{
public:
  static FullMatrix<double> A_n_m;

  mutable bool is_zero;

private:
  mutable unsigned int p;

  mutable dealii::Point<3> center;

  mutable const AssLegFunction *assLegFunction;

  mutable std::vector<std::complex<double>> _M_n_m;

public:
  MultipoleExpansion();

  MultipoleExpansion(const unsigned int      order,
                     const dealii::Point<3> &center,
                     const AssLegFunction   *assLegFunction);

  void
  Add(const MultipoleExpansion &multipole, const double sol);

  /* unused
  void
  Add(const MultipoleExpansion &         multipole,
      const double                       sol,
      std::vector<std::complex<double>> &cache);
  */

  void
  Add(const double strength, const dealii::Point<3> &point);

  void
  Add(const double                       strength,
      const dealii::Point<3> &           point,
      std::vector<std::complex<double>> &cache);

  void
  Add(const MultipoleExpansion &child);

  void
  Add(const MultipoleExpansion &         child,
      std::vector<std::complex<double>> &cache);

  void
  AddNormDer(const double                strength,
             const dealii::Point<3>     &point,
             const dealii::Tensor<1, 3> &normal);

  void
  AddNormDer(const double                       strength,
             const dealii::Point<3> &           point,
             const dealii::Tensor<1, 3> &       normal,
             std::vector<std::complex<double>> &cache);

  double
  Evaluate(const dealii::Point<3> &evalPoint);

  double
  Evaluate(const dealii::Point<3> &           evalPoint,
           std::vector<std::complex<double>> &cache);

  inline dealii::Point<3>
  GetCenter() const
  {
    return this->center;
  }

  inline void
  SetCenter(const dealii::Point<3> &new_center)
  {
    this->center = new_center;
  }

  inline FullMatrix<double> &
  GetA_n_m() const
  {
    return this->A_n_m;
  }

  inline std::complex<double> &
  GetCoeff(unsigned int n, unsigned int m) const
  {
    return this->_M_n_m[(n) * (n + 1) / 2 + m];
  }

  inline void
  SetCoeff(unsigned int n, unsigned int m, std::complex<double> &value) const
  {
    this->_M_n_m[(n) * (n + 1) / 2 + m] = value;
  }

  inline void
  AddToCoeff(unsigned int n, unsigned int m, std::complex<double> &value) const
  {
    this->_M_n_m[(n) * (n + 1) / 2 + m] += value;
  }

  static void
  spherical_coords(const dealii::Point<3> &center,
                   const dealii::Point<3> &other,
                   dealii::Point<3> &      blockRelPos,
                   double &                rho,
                   double &                cos_alpha,
                   double &                beta)
  {
    blockRelPos = other - center;
    rho         = blockRelPos.norm();
    cos_alpha   = blockRelPos(2) / rho;
    beta        = atan2(blockRelPos(1), blockRelPos(0));
  }

  static FullMatrix<double>
  A_n_m_Matrix(unsigned int dim)
  {
    FullMatrix<double> A_n_m(dim + 1, dim + 1);
    for (unsigned int n = 0; n < dim + 1; n++)
      {
        for (unsigned int m = 0; m < n + 1; m++)
          {
            double f1 = 1.;
            double f2 = 1.;
            /*
            //TODO: validate
            //slightly optimized implementation: less multiplications
            for (unsigned int ii = n - m; ii > 0; ii--)
              {
                f1 *= ii;
              }

            for (unsigned int ii = n + m; ii > n - m; ii--)
              {
                f2 *= (ii);
              }
            A_n_m(n, m) = (n % 2 ? -1. : 1.) / (sqrt(f2) * f1);
            */

            for (int ii = n - m; ii > 0; ii--)
              {
                f1 *= ii;
              }

            for (int ii = n + m; ii > 0; ii--)
              {
                f2 *= (ii);
              }

            A_n_m(n, m) = pow(-1., double(n)) / sqrt(f1 * f2);
          }
      }

    return A_n_m;
  }
};
#endif /*MULTIPOLE_EXPANSION_H_*/
