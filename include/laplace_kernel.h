
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

  template <int dim> // mio//
  Tensor<2, dim>
  hypersingular(const Point<dim> &R)
  {
    Tensor<2, dim> Hyper;
    switch (dim)
      {
        case 2:
          Hyper[0][0] = -R[0] * R[0] + R[1] * R[1];
          Hyper[1][1] = R[0] * R[0] - R[1] * R[1];
          Hyper[0][1] = -2 * R[0] * R[1];
          Hyper[1][0] = -2 * R[0] * R[1];
          Hyper       = Hyper / (-2 * numbers::PI * R.square() * R.square());
          return Hyper;
        case 3:
          {
            for (unsigned int i = 0; i < dim; ++i)
              for (unsigned int j = 0; j < dim; ++j)
                if (i == j)
                  {
                    Hyper[j][j] = 2 * (R[j] * R[j]);
                    for (unsigned int k = 0; k < dim && k != i; ++k)
                      Hyper[j][j] -= (R[k] * R[k]);
                  }
                else
                  Hyper[i][j] = 3 * R[i] * R[j];
            return Hyper /
                   (-4 * numbers::PI * R.square() * R.square() * R.norm());
          }


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

  template <int dim> // mio//
  void
  kernels(const Tensor<1, dim> &R,
          Tensor<2, dim> &      H,
          Tensor<1, dim> &      D,
          double &              d)
  {
    double r  = R.norm();
    double r2 = r * r;
    switch (dim)
      {
        case 2:
          d       = -std::log(r) / (2 * numbers::PI);
          D       = R / (-2 * numbers::PI * r2);
          H       = 0;
          H[0][0] = -R[0] * R[0] + R[1] * R[1];
          H[1][1] = R[0] * R[0] - R[1] * R[1];
          H[0][1] = -2 * R[0] * R[1];
          H[1][0] = -2 * R[0] * R[1];
          H       = H / (-2 * numbers::PI * r2 * r2);
          break;
        case 3:
          d = (1. / (r * 4 * numbers::PI));
          D = R / (4 * numbers::PI * r2 * r);
          H = 0;
          /*	for(unsigned int i=0;i<dim;++i)
                for(unsigned int j=0;j<dim;++j)
                  if(i==j)
              {
              H[j][j]=-2*(R[j]*R[j]);
              for(unsigned int k=0;k<dim;++k)
                            if(k!=i)
                H[j][j]+=(R[k]*R[k]);
              }
                  else H[i][j]=-3*R[i]*R[j];
            H=H/( 4*numbers::PI * r2 *r2*r );
          */
          for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
              H[i][j] = -3 * R[i] * R[j] / r2 / r2 / r;
          for (unsigned int i = 0; i < dim; ++i)
            H[i][i] += 1 / r2 / r;
          H = H / (4 * numbers::PI);

          break;
        default:
          Assert(false, ExcInternalError());
      }
  }

} // namespace LaplaceKernel
