#include <math.h>

#include <iostream>

#include "local_expansion.h"

#define GSL_SIGN(x) (x < 0 ? -1 : (x > 0 ? 1 : 0))

FullMatrix<double> LocalExpansion::A_n_m = LocalExpansion::A_n_m_Matrix(20);

LocalExpansionCoeff LocalExpansion::mExp_to_lExp_Coeff =
  LocalExpansion::mExp_to_lExp_Coeff_Build(LocalExpansion::A_n_m, 10);

std::vector<std::vector<std::map<int, std::map<int, double>>>>
  LocalExpansion::lExp_to_lExp_Coeff =
    LocalExpansion::lExp_to_lExp_Coeff_Build(LocalExpansion::A_n_m, 10);

LocalExpansion::LocalExpansion()
  : is_zero(true)
  , p(0)
  , center(0, 0, 0)
  , assLegFunction(nullptr)
  , _L_n_m()
{}

LocalExpansion::LocalExpansion(const unsigned int      order,
                               const dealii::Point<3> &center,
                               const AssLegFunction *  assLegFunction)
  : is_zero(true)
  , p(order)
  , center(center)
  , assLegFunction(assLegFunction)
  , _L_n_m((order + 1) * (order + 2) / 2, std::complex<double>(0.0, 0.0))
{}

void
LocalExpansion::Add(const std::vector<double> &real,
                    const std::vector<double> &imag)
{
  unsigned int count = 0;
  double       sum   = 0.0;
  for (unsigned int m = 0; m < this->p + 1; m++)
    {
      for (unsigned int n = m; n < this->p + 1; n++)
        {
          Assert(count < real.size() && count < imag.size(),
                 ExcInternalError());
          std::complex<double> a(real[count], imag[count]);
          sum += std::norm(a);
          this->AddToCoeff(n, n + m, a);
          a = std::conj(a);
          this->AddToCoeff(n, n - m, a);
          count++;
        }
    }

  if (sum > 1e-20)
    {
      this->is_zero = false;
    }
}

void
LocalExpansion::Add(
  const LocalExpansion &             other,
  std::vector<std::complex<double>> &cache) // translation of local expansion
{
  if (other.is_zero)
    {
    }
  else
    {
      unsigned int p = this->p;
      if (other.center.distance(this->center) > 1e-7)
        {
          dealii::Point<3> blockRelPos;
          double           rho, cos_alpha, beta;
          MultipoleExpansion::spherical_coords(
            center, other.GetCenter(), blockRelPos, rho, cos_alpha, beta);

          // cache rotations by beta; could we use simple powers, it would be
          // blazingly fast
          cache.reserve(2 * p + 1);
          cache.clear();
          cache.emplace_back(1);
          for (unsigned int i = 1; i < 2 * p + 1; ++i)
            {
              cache.emplace_back(std::cos(i * beta), std::sin(i * beta));
            }

          double P_nn_mm;
          for (int n = 0; n < int(p) + 1; n++)
            {
              for (int m = 0; m < n + 1; m++)
                {
                  std::complex<double> z = std::complex<double>(0., 0.);
                  for (int nn = n; nn < int(p) + 1; nn++)
                    {
                      double rhoFact = pow(rho, double(nn - n));
                      for (int mm = -1 * nn; mm < nn + 1; mm++)
                        {
                          if (abs(mm - m) > nn - n)
                            {
                            }
                          else
                            {
                              std::complex<double> a = std::complex<double>(
                                other.GetCoeff(abs(nn), abs(mm)).real(),
                                GSL_SIGN(mm) *
                                  other.GetCoeff(abs(nn), abs(mm)).imag());

                              P_nn_mm = this->assLegFunction->GetAssLegFunSph(
                                nn - n, abs(mm - m), cos_alpha);
                              double realFact =
                                P_nn_mm * rhoFact *
                                lExp_to_lExp_Coeff[n][m][nn][mm];

                              auto absm    = std::abs(mm - m);
                              auto rotated = ((mm - m) != absm) ?
                                               std::conj(cache[absm]) :
                                               cache[absm];

                              z += a * rotated * realFact;
                            }
                        }
                    }

                  this->AddToCoeff(n, m, z);
                }
            }
        }
      else
        {
          for (int n = 0; n < int(this->p) + 1; n++)
            {
              for (int m = 0; m < n + 1; m++)
                {
                  this->AddToCoeff(n, m, other.GetCoeff(n, m));
                }
            }
        }

      this->is_zero = false;
    }
}

void
LocalExpansion::Add(
  const LocalExpansion &other) // translation of local expansion
{
  if (!other.is_zero)
    {
      std::vector<std::complex<double>> cache;
      this->Add(other, cache);
    }
}

void
LocalExpansion::Add(const MultipoleExpansion &multipole,
                    std::vector<std::complex<double>>
                      &cache) // multipole conversion into local
                              // expansion, and addition to the rest
{
  if (multipole.is_zero)
    {
    }
  else
    {
      dealii::Point<3> blockRelPos;
      double           rho, cos_alpha, beta;
      MultipoleExpansion::spherical_coords(
        center, multipole.GetCenter(), blockRelPos, rho, cos_alpha, beta);

      cache.reserve(2 * p + 1);
      cache.clear();
      cache.emplace_back(1);
      for (unsigned int i = 1; i < 2 * p + 1; ++i)
        {
          cache.emplace_back(std::cos(i * beta), std::sin(i * beta));
        }

      double P_nn_mm;
      for (int n = 0; n < int(this->p) + 1; n++)
        {
          for (int m = 0; m < n + 1; m++)
            {
              std::complex<double> z = std::complex<double>(0., 0.);
              for (int nn = 0; nn < int(this->p) + 1; nn++)
                {
                  double rhoFact = pow(rho, double(-n - nn - 1));
                  for (int mm = -1 * nn; mm < 0; mm++)
                    {
                      std::complex<double> a =
                        std::conj(multipole.GetCoeff(nn, abs(mm)));
                      P_nn_mm =
                        this->assLegFunction->GetAssLegFunSph(nn + n,
                                                              abs(mm - m),
                                                              cos_alpha);
                      double realFact = P_nn_mm * rhoFact *
                                        mExp_to_lExp_Coeff.get(n, m, nn, mm);

                      auto absm    = std::abs(mm - m);
                      auto rotated = ((mm - m) != absm) ?
                                       std::conj(cache[absm]) :
                                       cache[absm];

                      z += a * rotated * realFact;
                    }

                  for (int mm = 0; mm < nn + 1; mm++)
                    {
                      std::complex<double> a = multipole.GetCoeff(nn, abs(mm));
                      P_nn_mm =
                        this->assLegFunction->GetAssLegFunSph(nn + n,
                                                              abs(mm - m),
                                                              cos_alpha);
                      double realFact = P_nn_mm * rhoFact *
                                        mExp_to_lExp_Coeff.get(n, m, nn, mm);

                      auto absm    = std::abs(mm - m);
                      auto rotated = ((mm - m) != absm) ?
                                       std::conj(cache[absm]) :
                                       cache[absm];

                      z += a * rotated * realFact;
                    }
                }

              this->AddToCoeff(n, m, z);
            }
        }

      this->is_zero = false;
    }
}

void
LocalExpansion::Add(
  const MultipoleExpansion &multipole) // multipole conversion into local
                                       // expansion, and addition to the rest
{
  if (!multipole.is_zero)
    {
      std::vector<std::complex<double>> cache;
      Add(multipole, cache);
    }
}

double
LocalExpansion::Evaluate(const dealii::Point<3> &           evalPoint,
                         std::vector<std::complex<double>> &cache)
{
  std::complex<double> fieldValue = std::complex<double>(0., 0.);
  if (this->is_zero)
    {
    }
  else
    {
      dealii::Point<3> blockRelPos;
      double           rho, cos_alpha, beta;
      MultipoleExpansion::spherical_coords(
        center, evalPoint, blockRelPos, rho, cos_alpha, beta);

      cache.reserve(p + 1);
      cache.clear();
      cache.emplace_back(1);
      for (unsigned int i = 1; i < p + 1; ++i)
        {
          // TODO: try to use exp
          cache.emplace_back(std::cos(i * beta), std::sin(i * beta));
        }

      double P_n_m;
      for (int n = 0; n < int(p) + 1; n++)
        {
          P_n_m        = this->assLegFunction->GetAssLegFunSph(n, 0, cos_alpha);
          double rho2n = pow(rho, double(n));
          double realFact = P_n_m * rho2n;

          fieldValue += this->GetCoeff(n, 0) * realFact;
          for (int m = 1; m < n + 1; m++)
            {
              P_n_m = this->assLegFunction->GetAssLegFunSph(n, m, cos_alpha);
              double realFact = P_n_m * rho2n;

              std::complex<double> complexFact = cache[m] * 2. * realFact;

              fieldValue += this->GetCoeff(n, m) * complexFact;
            }
        }
    }

  return fieldValue.real();
}

double
LocalExpansion::Evaluate(const dealii::Point<3> &evalPoint)
{
  std::complex<double> fieldValue = std::complex<double>(0., 0.);
  if (!this->is_zero)
    {
      std::vector<std::complex<double>> cache;
      return Evaluate(evalPoint, cache);
    }

  return fieldValue.real();
}
