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

    if (lhs.center.distance_square(rhs.center) > (tolerance * tolerance))
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
  Add(const LocalExpansion &parent, std::vector<std::complex<double>> &cache);

  void
  Add(const MultipoleExpansion &multipole);

  void
  Add(const MultipoleExpansion          &multipole,
      std::vector<std::complex<double>> &cache);

  double
  Evaluate(const dealii::Point<3> &evalPoint);

  double
  Evaluate(const dealii::Point<3>            &evalPoint,
           std::vector<std::complex<double>> &cache);

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

  inline const std::vector<std::complex<double>> &
  GetCoeffs() const
  {
    return this->_L_n_m;
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

            for (int ii = n - m; ii > 0; ii--)
              {
                f1 *= ii;
              }

            for (int ii = n + m; ii > 0; ii--)
              {
                f2 *= (ii);
              }

            A_n_m(n, m) = std::pow(-1., double(n)) / std::sqrt(f1 * f2);
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
                    double realFact = A_n_m(nn, std::abs(mm)) /
                                      A_n_m(n + nn, std::abs(m - mm)) *
                                      A_n_m(n, std::abs(m));

                    // reference implementation
                    std::complex<double> imUnit(0, 1);
                    int steps = std::abs(m - mm) - std::abs(m) - std::abs(mm);
                    realFact *= std::pow(imUnit, double(steps)).real() /
                                std::pow(-1., nn);

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
                    double realFact = A_n_m(nn - n, std::abs(mm - m)) /
                                      A_n_m(nn, std::abs(mm)) *
                                      A_n_m(n, std::abs(m));

                    // reference implementation
                    std::complex<double> imUnit(0, 1);
                    int steps = std::abs(mm) - std::abs(mm - m) - std::abs(m);
                    realFact *= std::pow(imUnit, double(steps)).real() *
                                std::pow(-1., nn + n);

                    realCoeff[n][m][nn][mm] = realFact;
                  }
              }
          }
      }

    return realCoeff;
  }
};
#endif
