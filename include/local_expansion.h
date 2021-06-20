#ifndef LOCALEXPANSION_H_
#define LOCALEXPANSION_H_

#include <deal.II/lac/full_matrix.h>

#include <math.h>

#include <string>
#include <vector>

#include "ass_leg_function.h"
#include "local_expansion_coeff.h"
#include "multipole_expansion.h"

class LocalExpansion
{
public:
  static FullMatrix<double> A_n_m;

  static std::vector<std::vector<std::map<int, std::map<int, double>>>>
    lExp_to_lExp_Coeff;

  static LocalExpansionCoeff mExp_to_lExp_Coeff;

  static bool
  equal(const LocalExpansion &lhs,
        const LocalExpansion &rhs,
        double                tolerance = 1e-10)
  {
    if (lhs.p != rhs.p)
      {
        return false;
      }

    if (lhs.center.distance_square(rhs.center) > tolerance * tolerance)
      {
        return false;
      }

    for (unsigned int i = 0; i < (lhs.p + 1) * (lhs.p + 2) / 2; ++i)
      {
        if (std::abs(lhs._L_n_m[i] - rhs._L_n_m[i]) > tolerance)
          {
            return false;
          }
      }

    return true;
  }

  mutable bool is_zero;

private:
  mutable unsigned int p;

  mutable dealii::Point<3> center;

  mutable const AssLegFunction *assLegFunction;

  mutable std::vector<std::complex<double>> _L_n_m;

public:
  LocalExpansion();

  LocalExpansion(const unsigned int      order,
                 const dealii::Point<3> &center,
                 const AssLegFunction   *assLegFunction);

  void
  Add(const std::vector<double> &real, const std::vector<double> &imag);

  void
  Add(const LocalExpansion &parent);

  void
  Add(const MultipoleExpansion &multipole);

  double
  Evaluate(const dealii::Point<3> &evalPoint);

  inline dealii::Point<3> &
  GetCenter() const
  {
    return this->center;
  }

  inline void
  SetCenter(const dealii::Point<3> &new_center)
  {
    this->center = new_center;
  }

  inline FullMatrix<double>
  GetA_n_m() const
  {
    return this->A_n_m;
  }

  inline unsigned int
  GetOrder() const
  {
    return this->p;
  }

  inline std::complex<double> &
  GetCoeff(unsigned int n, unsigned int m) const
  {
    return this->_L_n_m[(n) * (n + 1) / 2 + m];
  }

  void
  SetCoeff(unsigned int n, unsigned int m, std::complex<double> &value) const
  {
    this->_L_n_m[(n) * (n + 1) / 2 + m] = value;
  }

  void
  AddToCoeff(unsigned int n, unsigned int m, std::complex<double> &value) const
  {
    this->_L_n_m[(n) * (n + 1) / 2 + m] += value;
  }

  static FullMatrix<double>
  A_n_m_Matrix(unsigned int dimension)
  {
    FullMatrix<double> A_n_m(dimension + 1, dimension + 1);
    for (unsigned int n = 0; n < dimension + 1; n++)
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

  static LocalExpansionCoeff
  mExp_to_lExp_Coeff_Build(FullMatrix<double> A_n_m, unsigned int p)
  {
    LocalExpansionCoeff loc_exp_coeff(p);
    for (int n = 0; n < int(p) + 1; n++)
      {
        for (int m = 0; m < n + 1; m++)
          {
            for (int nn = 0; nn < int(p) + 1; nn++)
              {
                for (int mm = -1 * nn; mm < nn + 1; mm++)
                  {
                    double realFact = A_n_m(nn, abs(mm)) /
                                      A_n_m(n + nn, abs(m - mm)) *
                                      A_n_m(n, abs(m));

                    /*
                    //TODO: validate
                    unsigned int steps   = (abs(m - mm) - abs(m) - abs(mm)) % 4;
                    int          rotated = 0;
                    if (steps == 0)
                      {
                        rotated = 1;
                      }
                    else if (steps == 2)
                      {
                        rotated = -1;
                      }
                    rotated *= ((nn % 2) ? -1 : 1);

                    realFact *= rotated;
                    */
                    // reference implementation
                    auto imUnit = std::complex<double>(0, 1);
                    realFact *=
                      (pow(imUnit, double(abs(m - mm) - abs(m) - abs(mm))))
                        .real() /
                      pow(-1., nn);

                    loc_exp_coeff.set(n, m, nn, mm, realFact);
                  }
              }
          }
      }

    return loc_exp_coeff;
  }

  static std::vector<std::vector<std::map<int, std::map<int, double>>>>
  lExp_to_lExp_Coeff_Build(FullMatrix<double> A_n_m, unsigned int p)
  {
    std::vector<std::vector<std::map<int, std::map<int, double>>>> realCoeff;
    realCoeff.resize(p + 1);
    for (int n = 0; n < int(p) + 1; n++)
      {
        realCoeff[n].resize(n + 1);
        for (int m = 0; m < n + 1; m++)
          {
            for (int nn = n; nn < int(p) + 1; nn++)
              {
                for (int mm = -1 * nn; mm < nn + 1; mm++)
                  {
                    double realFact = A_n_m(nn - n, abs(mm - m)) /
                                      A_n_m(nn, abs(mm)) * A_n_m(n, abs(m));
                    /*
                    //TODO: validate
                                        unsigned int steps   = (abs(mm) - abs(mm
                       - m) - abs(m)) % 4; int          rotated = 0; if (steps
                       == 0)
                                          {
                                            rotated = 1;
                                          }
                                        else if (steps == 2)
                                          {
                                            rotated = -1;
                                          }
                                        rotated *= (((nn + n) % 2) ? -1 : 1);

                                        realFact *= rotated;
                                        */
                    // reference implementation
                    auto imUnit = std::complex<double>(0, 1);
                    realFact *=
                      (pow(imUnit, double(abs(mm) - abs(mm - m) - abs(m))))
                        .real() *
                      pow(-1., nn + n);

                    realCoeff[n][m][nn][mm] = realFact;
                  }
              }
          }
      }

    return realCoeff;
  }
};
#endif /*LOCALEXPANSION_H_*/
