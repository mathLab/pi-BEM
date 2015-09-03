#ifndef MULTIPOLE_EXPANSION_H_
#define MULTIPOLE_EXPANSION_H_

#include <math.h>
#include <string>
#include <vector>
#include <complex>

#include <deal.II/base/point.h>

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

  mutable std::complex <double> *_M_n_m;



public:

  MultipoleExpansion();

  MultipoleExpansion(const unsigned int order,const  dealii::Point<3> &center,const AssLegFunction *assLegFunction);

  MultipoleExpansion(const MultipoleExpansion &other);

  ~MultipoleExpansion();

  void Add(const MultipoleExpansion &multipole,const double sol);

  void Add(const double strength, const dealii::Point<3> &point);

  void Add(const MultipoleExpansion &child);

  void AddNormDer(const double strength, const dealii::Point<3> &point, const dealii::Tensor<1, 3> &normal);

  double Evaluate(const dealii::Point<3> &evalPoint);

  inline dealii::Point<3> GetCenter() const
  {
    return this->center;
  }

  inline void SetCenter(const dealii::Point<3> &new_center)
  {
    this->center = new_center;
  }

  inline FullMatrix<double> &GetA_n_m() const
  {
    return this->A_n_m;
  }

  inline std::complex <double> *GetCoeffs() const
  {
    return this->_M_n_m;
  }

  inline std::complex <double> &GetCoeff(unsigned int n, unsigned int m) const
  {
    return this->_M_n_m[(n)*(n+1)/2+m];
  }

  inline void SetCoeff(unsigned int n, unsigned int m, std::complex <double> &value) const
  {
    this->_M_n_m[(n)*(n+1)/2+m] = value;
  }

  inline void AddToCoeff(unsigned int n, unsigned int m, std::complex <double> &value) const
  {
    this->_M_n_m[(n)*(n+1)/2+m] += value;
  }

  MultipoleExpansion &operator=( const MultipoleExpansion &other );


  static FullMatrix<double> A_n_m_Matrix(unsigned int dim)
  {
    FullMatrix<double> A_n_m(dim+1,dim+1);
    for (unsigned int n = 0; n < dim+1 ; n++)

      {
        for (unsigned int m = 0; m < n+1 ; m++)

          {
            double f1 = 1.;
            double f2 = 1.;

            for (int ii = n-m; ii > 0; ii-- )
              f1 *= ii;

            for (int ii = n+m; ii > 0; ii-- )
              f2 *= (ii);

            A_n_m(n,m) = pow(-1.,double(n))/sqrt(f1*f2);
          }
      }

    return A_n_m;
  }


};



#endif /*MULTIPOLE_EXPANSION_H_*/
