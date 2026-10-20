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

  /// Accumulate the supplied source or expansion. Cache-taking overloads
  /// overwrite and resize caller-owned scratch storage; use a separate cache
  /// per concurrent operation. Source expansions must have compatible order.
  /// The non-owning AssLegFunction passed at construction must remain alive.
  void
  Add(const MultipoleExpansion &multipole, const double sol);

  /// Accumulate the supplied source or expansion. Cache-taking overloads
  /// overwrite and resize caller-owned scratch storage; use a separate cache
  /// per concurrent operation. Source expansions must have compatible order.
  /// The non-owning AssLegFunction passed at construction must remain alive.
  void
  Add(const double strength, const dealii::Point<3> &point);

  /// Accumulate the supplied source or expansion. Cache-taking overloads
  /// overwrite and resize caller-owned scratch storage; use a separate cache
  /// per concurrent operation. Source expansions must have compatible order.
  /// The non-owning AssLegFunction passed at construction must remain alive.
  void
  Add(const double                       strength,
      const dealii::Point<3>            &point,
      std::vector<std::complex<double>> &cache);

  /// Accumulate the supplied source or expansion. Cache-taking overloads
  /// overwrite and resize caller-owned scratch storage; use a separate cache
  /// per concurrent operation. Source expansions must have compatible order.
  /// The non-owning AssLegFunction passed at construction must remain alive.
  void
  Add(const MultipoleExpansion &child);

  /// Accumulate the supplied source or expansion. Cache-taking overloads
  /// overwrite and resize caller-owned scratch storage; use a separate cache
  /// per concurrent operation. Source expansions must have compatible order.
  /// The non-owning AssLegFunction passed at construction must remain alive.
  void
  Add(const MultipoleExpansion          &child,
      std::vector<std::complex<double>> &cache);

  /// Accumulate a normal-derivative source. The cache overload overwrites and
  /// resizes caller-owned scratch storage; provide one cache per concurrent
  /// call.
  void
  AddNormDer(const double                strength,
             const dealii::Point<3>     &point,
             const dealii::Tensor<1, 3> &normal);

  /// Accumulate a normal-derivative source. The cache overload overwrites and
  /// resizes caller-owned scratch storage; provide one cache per concurrent
  /// call.
  void
  AddNormDer(const double                       strength,
             const dealii::Point<3>            &point,
             const dealii::Tensor<1, 3>        &normal,
             std::vector<std::complex<double>> &cache);

  /// Evaluate the expansion at evalPoint. The cache overload reuses mutable
  /// caller-owned scratch storage and may overwrite/resize it. Do not share
  /// a cache between concurrent calls. Observe the expansion convergence
  /// region.
  double
  Evaluate(const dealii::Point<3> &evalPoint);

  /// Evaluate the expansion at evalPoint. The cache overload reuses mutable
  /// caller-owned scratch storage and may overwrite/resize it. Do not share
  /// a cache between concurrent calls. Observe the expansion convergence
  /// region.
  double
  Evaluate(const dealii::Point<3>            &evalPoint,
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

  /// Return a read-only reference to the shared normalization matrix.
  inline const FullMatrix<double> &
  GetA_n_m() const
  {
    return this->A_n_m;
  }

  /// Return a read-only reference to owned triangular coefficient storage,
  /// indexed by n*(n+1)/2+m. Copy the vector to retain an independent snapshot.
  /// Use SetCoeff()/AddToCoeff() for changes. Assignment can invalidate
  /// references to individual coefficients; the vector object remains owned
  /// here.
  inline const std::vector<std::complex<double>> &
  GetCoeffs() const
  {
    return this->_M_n_m;
  }

  /// Return a read-only coefficient reference for 0 <= m <= n <= order.
  inline const std::complex<double> &
  GetCoeff(unsigned int n, unsigned int m) const
  {
    return this->_M_n_m[(n) * (n + 1) / 2 + m];
  }

  /// Replace coefficient (n,m), with 0 <= m <= n <= order. The value is copied.
  inline void
  SetCoeff(unsigned int                n,
           unsigned int                m,
           const std::complex<double> &value) const
  {
    this->_M_n_m[(n) * (n + 1) / 2 + m] = value;
  }

  /// Accumulate a copied value at coefficient (n,m), with 0 <= m <= n <= order.
  inline void
  AddToCoeff(unsigned int                n,
             unsigned int                m,
             const std::complex<double> &value) const
  {
    this->_M_n_m[(n) * (n + 1) / 2 + m] += value;
  }

  /// Compute other-center, its length rho, cos(polar angle) and azimuth beta.
  /// The points must differ, since cos_alpha divides by rho.
  static void
  spherical_coords(const dealii::Point<3> &center,
                   const dealii::Point<3> &other,
                   dealii::Point<3>       &blockRelPos,
                   double                 &rho,
                   double                 &cos_alpha,
                   double                 &beta)
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
#endif
