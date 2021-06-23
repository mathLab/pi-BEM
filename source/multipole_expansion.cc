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
  if (multipole.is_zero || sol == 0)
    {
    }
  else
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
  if (strength == 0)
    {
    }
  else
    {
      this->is_zero = false;

      dealii::Point<3> pointRelPos;
      double           rho, cos_alpha, beta;
      MultipoleExpansion::spherical_coords(
        center, point, pointRelPos, rho, cos_alpha, beta);

      // cache rotations by beta; could we use simple powers, it would be
      // blazingly fast
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
          for (int m = 0; m < n + 1; m++)
            {
              P_n_m =
                this->assLegFunction->GetAssLegFunSph(n, abs(m), cos_alpha);
              double realFact = P_n_m * pow(rho, double(n)) * strength;

              // another option was flipping the nesting of for_n and for_m
              // better due to skipping memory allocation for the cache
              /* reference
              std::complex<double> a =
                exp(std::complex<double>(0., -m * beta)) * realFact;
              */
              std::complex<double> a = cache[m] * realFact;

              this->AddToCoeff(n, m, a);
            }
        }
    }
}
void
MultipoleExpansion::Add(const double strength, const dealii::Point<3> &point)
{
  if (strength != 0)
    {
      std::vector<std::complex<double>> cache;
      Add(strength, point, cache);
    }
}

void
MultipoleExpansion::AddNormDer(const double                       strength,
                               const dealii::Point<3> &           point,
                               const dealii::Tensor<1, 3> &       normal,
                               std::vector<std::complex<double>> &cache)
{
  if (strength == 0)
    {
    }
  else
    {
      this->is_zero = false;

      dealii::Point<3> pointRelPos;
      double           rho, cos_alpha, beta;
      MultipoleExpansion::spherical_coords(
        center, point, pointRelPos, rho, cos_alpha, beta);

      const auto normVersor = normal / normal.norm();
      double     dRhodN     = (pointRelPos / rho) * normVersor;
      double dBetadN = (dealii::Point<3>(-pointRelPos(1), pointRelPos(0), 0.) /
                        (pow(pointRelPos(0), 2.) + pow(pointRelPos(1), 2.))) *
                       normVersor;

      // below, was using cos_alpha_ notice the trailing underscore
      double sin_alpha = std::sqrt(1. - pow(cos_alpha, 2.));
      double dAlphadN  = (dealii::Point<3>(cos_alpha * cos(beta),
                                          cos_alpha * sin(beta),
                                          -sin_alpha) /
                         rho) *
                        normVersor;

      // cache rotations by beta; could we use simple powers, it would be
      // blazingly fast
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
          for (int m = 0; m < n + 1; m++)
            {
              // below, was using cos_alpha_ notice the trailing underscore
              P_n_m =
                this->assLegFunction->GetAssLegFunSph(n, abs(m), cos_alpha);
              dP_n_m_sin =
                this->assLegFunction->GetAssLegFunSphDeriv(n,
                                                           abs(m),
                                                           cos_alpha) *
                sin_alpha;

              // another option was flipping the nesting of for_n and for_m
              // better due to skipping memory allocation for the cache
              /*
              std::complex<double> z =
                exp(std::complex<double>(0., -double(m) * beta));
              */
              std::complex<double> z = cache[m];

              z *= std::complex<double>(
                double(n) * pow(rho, double(n) - 1.) * P_n_m * dRhodN -
                  pow(rho, double(n)) * dP_n_m_sin * dAlphadN,
                -double(m) * pow(rho, double(n)) * P_n_m * dBetadN);
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
  if (strength != 0)
    {
      std::vector<std::complex<double>> cache;
      AddNormDer(strength, point, normal, cache);
    }
}

void
MultipoleExpansion::Add(
  const MultipoleExpansion &other,
  std::vector<std::complex<double>>
    &cache) // translation of a multipole to its parent center
{
  if (other.is_zero)
    {
      this->is_zero = this->is_zero & other.is_zero;
    }
  else
    {
      this->is_zero             = false;
      FullMatrix<double> &A_n_m = this->GetA_n_m();
      if (other.center.distance(this->center) > 1e-7)
        {
          dealii::Point<3> blockRelPos;
          double           rho, cos_alpha, beta;
          MultipoleExpansion::spherical_coords(
            center, other.center, blockRelPos, rho, cos_alpha, beta);

          // cache rotations by beta; could we use simple powers, it would be
          // blazingly fast
          cache.reserve(p + 1);
          cache.clear();
          cache.emplace_back(1);
          for (unsigned int i = 1; i < p + 1; ++i)
            {
              cache.emplace_back(std::exp(std::complex<double>(0., i * beta)));
            }

          double P_nn_mm;
          for (int n = 0; n < int(this->p) + 1; n++)
            {
              for (int m = 0; m < n + 1; m++)
                {
                  std::complex<double> z(0., 0.);
                  for (int nn = 0; nn < n + 1; nn++)
                    {
                      for (int mm = -1 * nn; mm < nn + 1; mm++)
                        {
                          if (abs(m - mm) > n - nn)
                            {
                            }
                          else
                            {
                              std::complex<double> a = std::complex<double>(
                                (other.GetCoeff(abs(n - nn), abs(m - mm)))
                                  .real(),
                                GSL_SIGN(m - mm) *
                                  (other.GetCoeff(abs(n - nn), abs(m - mm)))
                                    .imag());

                              P_nn_mm = this->assLegFunction->GetAssLegFunSph(
                                nn, abs(mm), cos_alpha);

                              double realFact =
                                P_nn_mm * pow(rho, double(nn)) *
                                A_n_m(abs(nn), abs(mm)) *
                                A_n_m(abs(n - nn), abs(m - mm)) /
                                A_n_m(abs(n), abs(m));

                              /*
                              //TODO: validate
                              unsigned int steps =
                                abs(m) - abs(mm) - abs(m - mm);
                              double rotated = 0;
                              if (steps % 4 == 0)
                                {
                                  rotated = 1;
                                }
                              else if (steps % 4 == 2)
                                {
                                  rotated = -1;
                                }

                              realFact *= rotated;
                              */

                              // reference implementation
                              auto imUnit = std::complex<double>(0, 1);
                              realFact *=
                                (pow(imUnit,
                                     double(abs(m) - abs(mm) - abs(m - mm))))
                                  .real();

                              /*
                              z +=
                                realFact *
                                (a * exp(std::complex<double>(0., -mm * beta)));
                              */
                              // TODO: validate
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
          cache.emplace_back(std::exp(std::complex<double>(0., i * beta)));
        }

      double P_n_m;

      for (int n = 0; n < int(this->p) + 1; n++)
        {
          P_n_m = this->assLegFunction->GetAssLegFunSph(n, 0, cos_alpha);
          const double rho_n1 = pow(rho, double(-n - 1));

          double realFact = P_n_m * rho_n1;
          fieldValue += this->GetCoeff(n, 0) * realFact;
          for (int m = 1; m < n + 1; m++)
            {
              P_n_m =
                this->assLegFunction->GetAssLegFunSph(n, abs(m), cos_alpha);
              double realFact = P_n_m * rho_n1;

              /*
              fieldValue += this->GetCoeff(n, m) *
                            exp(std::complex<double>(0., m * beta)) * 2. *
                            realFact;
              */
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
