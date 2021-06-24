#include <deal.II/base/point.h>

#include <iostream>

#include "multipole_expansion.h"

#define GSL_SIGN(x) (x < 0 ? -1 : (x > 0 ? 1 : 0))

FullMatrix<double> MultipoleExpansion::A_n_m =
  MultipoleExpansion::A_n_m_Matrix(20);

MultipoleExpansion::MultipoleExpansion()
  : is_zero(true)
  , p(0)
  , center(0, 0, 0)
  , assLegFunction(nullptr)
  , _M_n_m()
{}

MultipoleExpansion::MultipoleExpansion(const unsigned int      order,
                                       const dealii::Point<3> &center,
                                       const AssLegFunction *  assLegFunction)
  : is_zero(true)
  , p(order)
  , center(center)
  , assLegFunction(assLegFunction)
  , _M_n_m((order + 1) * (order + 2) / 2, std::complex<double>(0.0, 0.0))
{}

void
MultipoleExpansion::Add(const MultipoleExpansion &multipole, const double sol)
{
  // TODO: sol testing should use a tolerance
  if (!multipole.is_zero && sol != 0)
    {
      this->is_zero = false;
      for (int n = 0; n < int(this->p) + 1; n++)
        for (int m = 0; m < n + 1; m++)
          {
            std::complex<double> a = multipole.GetCoeff(n, m) * sol;
            this->AddToCoeff(n, m, a);
          }
    }
}

void
MultipoleExpansion::Add(const double                       strength,
                        const dealii::Point<3> &           point,
                        std::vector<std::complex<double>> &cache)
{
  // TODO: strength testing should use a tolerance
  if (strength != 0)
    {
      this->is_zero = false;

      dealii::Point<3> pointRelPos;
      double           rho, cos_alpha, beta;
      MultipoleExpansion::spherical_coords(
        center, point, pointRelPos, rho, cos_alpha, beta);

      cache.reserve(p + 1);
      cache.clear();
      cache.emplace_back(1);
      for (unsigned int i = 1; i < p + 1; ++i)
        {
          cache.emplace_back(std::exp(std::complex<double>(0., -(i * beta))));
        }

      double P_n_m;
      for (int n = 0; n < int(this->p) + 1; n++)
        {
          const double rho2n = std::pow(rho, double(n));
          for (int m = 0; m < n + 1; m++)
            {
              P_n_m = this->assLegFunction->GetAssLegFunSph(n, m, cos_alpha);
              double realFact = P_n_m * rho2n * strength;

              std::complex<double> a = cache[m] * realFact;

              this->AddToCoeff(n, m, a);
            }
        }
    }
}
void
MultipoleExpansion::Add(const double strength, const dealii::Point<3> &point)
{
  std::vector<std::complex<double>> cache;
  Add(strength, point, cache);
}

void
MultipoleExpansion::AddNormDer(const double                       strength,
                               const dealii::Point<3> &           point,
                               const dealii::Tensor<1, 3> &       normal,
                               std::vector<std::complex<double>> &cache)
{
  // TODO: strength testing should use a tolerance
  if (strength != 0)
    {
      this->is_zero = false;

      dealii::Point<3> pointRelPos;
      double           rho, cos_alpha, beta;
      MultipoleExpansion::spherical_coords(
        center, point, pointRelPos, rho, cos_alpha, beta);

      const auto normVersor = normal / normal.norm();
      double     dRhodN     = (pointRelPos / rho) * normVersor;
      double     dBetadN =
        (dealii::Point<3>(-pointRelPos(1), pointRelPos(0), 0.) /
         (std::pow(pointRelPos(0), 2.) + std::pow(pointRelPos(1), 2.))) *
        normVersor;

      double sin_alpha = std::sqrt(1. - std::pow(cos_alpha, 2.));
      double dAlphadN  = (dealii::Point<3>(cos_alpha * std::cos(beta),
                                          cos_alpha * std::sin(beta),
                                          -sin_alpha) /
                         rho) *
                        normVersor;

      cache.reserve(p + 1);
      cache.clear();
      cache.emplace_back(1);
      for (unsigned int i = 1; i < p + 1; ++i)
        {
          cache.emplace_back(std::exp(std::complex<double>(0., -(i * beta))));
        }

      double P_n_m;
      double dP_n_m_sin;
      for (int n = 0; n < int(this->p) + 1; n++)
        {
          const double rho2n  = std::pow(rho, double(n));
          const double rho2n1 = std::pow(rho, double(n) - 1.);

          for (int m = 0; m < n + 1; m++)
            {
              P_n_m = this->assLegFunction->GetAssLegFunSph(n, m, cos_alpha);
              dP_n_m_sin =
                this->assLegFunction->GetAssLegFunSphDeriv(n, m, cos_alpha) *
                sin_alpha;

              std::complex<double> z = cache[m];
              z *= std::complex<double>(double(n) * rho2n1 * P_n_m * dRhodN -
                                          rho2n * dP_n_m_sin * dAlphadN,
                                        -double(m) * rho2n * P_n_m * dBetadN);
              z *= strength;

              this->AddToCoeff(n, m, z);
            }
        }
    }
}

void
MultipoleExpansion::AddNormDer(const double                strength,
                               const dealii::Point<3> &    point,
                               const dealii::Tensor<1, 3> &normal)
{
  std::vector<std::complex<double>> cache;
  AddNormDer(strength, point, normal, cache);
}

void
MultipoleExpansion::Add(
  const MultipoleExpansion &other,
  std::vector<std::complex<double>>
    &cache) // translation of a multipole to its parent center
{
  const double tolerance = 1e-7;
  if (!other.is_zero)
    {
      this->is_zero             = false;
      FullMatrix<double> &A_n_m = this->GetA_n_m();
      if (other.center.distance_square(this->center) > (tolerance * tolerance))
        {
          dealii::Point<3> blockRelPos;
          double           rho, cos_alpha, beta;
          MultipoleExpansion::spherical_coords(
            center, other.center, blockRelPos, rho, cos_alpha, beta);

          cache.reserve(p + 1);
          cache.clear();
          cache.emplace_back(1);
          for (unsigned int i = 1; i < p + 1; ++i)
            {
              cache.emplace_back(std::exp(std::complex<double>(0., i * beta)));
            }

          // const double imUnitPow[4] = {1, 0, -1, 0};
          double P_nn_mm;
          for (int n = 0; n < int(this->p) + 1; n++)
            {
              for (int m = 0; m < n + 1; m++)
                {
                  std::complex<double> z(0., 0.);
                  for (int nn = 0; nn < n + 1; nn++)
                    {
                      const double rho2nn = std::pow(rho, double(nn));
                      for (int mm = -1 * nn; mm < nn + 1; mm++)
                        {
                          if (std::abs(m - mm) <= n - nn)
                            {
                              std::complex<double> a(
                                other.GetCoeff(n - nn, std::abs(m - mm)).real(),
                                GSL_SIGN(m - mm) *
                                  (other.GetCoeff(n - nn, std::abs(m - mm)))
                                    .imag());

                              P_nn_mm = this->assLegFunction->GetAssLegFunSph(
                                nn, std::abs(mm), cos_alpha);

                              double realFact =
                                P_nn_mm * rho2nn * A_n_m(nn, std::abs(mm)) *
                                A_n_m(n - nn, std::abs(m - mm)) / A_n_m(n, m);

                              std::complex<double> imUnit(0, 1);
                              int steps = m - abs(mm) - abs(m - mm);
                              realFact *= std::pow(imUnit, steps).real();
                              auto rotated =
                                (mm > 0) ? std::conj(cache[mm]) : cache[-mm];
                              z += realFact * (a * rotated);
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
    }
}

void
MultipoleExpansion::Add(
  const MultipoleExpansion
    &other) // translation of a multipole to its parent center
{
  std::vector<std::complex<double>> cache;
  Add(other, cache);
}

double
MultipoleExpansion::Evaluate(const dealii::Point<3> &           evalPoint,
                             std::vector<std::complex<double>> &cache)
{
  std::complex<double> fieldValue(0., 0.);
  if (!this->is_zero)
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
          cache.emplace_back(std::exp(std::complex<double>(0., i * beta)));
        }

      double P_n_m;

      for (int n = 0; n < int(this->p) + 1; n++)
        {
          P_n_m = this->assLegFunction->GetAssLegFunSph(n, 0, cos_alpha);
          const double rho2n1   = std::pow(rho, double(-n - 1));
          double       realFact = P_n_m * rho2n1;

          fieldValue += this->GetCoeff(n, 0) * realFact;
          for (int m = 1; m < n + 1; m++)
            {
              P_n_m    = this->assLegFunction->GetAssLegFunSph(n, m, cos_alpha);
              realFact = P_n_m * rho2n1;

              fieldValue += this->GetCoeff(n, m) * cache[m] * 2. * realFact;
            }
        }
    }

  return fieldValue.real();
}

double
MultipoleExpansion::Evaluate(const dealii::Point<3> &evalPoint)
{
  std::vector<std::complex<double>> cache;
  return Evaluate(evalPoint, cache);
}
