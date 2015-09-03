#include <iostream>

#include "multipole_expansion.h"
#include <deal.II/base/point.h>

#define GSL_SIGN(x) (x<0 ? -1: (x>0 ? 1: 0))

FullMatrix<double> MultipoleExpansion::A_n_m = MultipoleExpansion::A_n_m_Matrix(20);

MultipoleExpansion::MultipoleExpansion()

{

  this->p = 0;
  this->center = Point<3>(0,0,0);
  this->assLegFunction = NULL;
  this->_M_n_m = NULL;
  this->is_zero = true;
}


MultipoleExpansion::MultipoleExpansion(const unsigned int order, const dealii::Point<3> &center, const AssLegFunction *assLegFunction)

{

  this->p = order;
  this->center = center;
  this->assLegFunction = assLegFunction;
  this->_M_n_m= new std::complex <double>[(this->p+1)*(this->p+2)/2];
  for (unsigned int i=0; i<(this->p+1)*(this->p+2)/2; ++i)
    this->_M_n_m[i] = std::complex <double>(0.0,0.0);
  this->is_zero = true;

}

MultipoleExpansion::MultipoleExpansion(const MultipoleExpansion &other)
{
  this->p=other.p;
  this->assLegFunction = other.assLegFunction;
  this->center = other.center;
  this->_M_n_m = new std::complex <double>[(this->p+1)*(this->p+2)/2];
  memcpy(this->_M_n_m, other.GetCoeffs(), sizeof(std::complex <double>)*(this->p+1)*(this->p+2)/2);
  this->is_zero = other.is_zero;
}



MultipoleExpansion &MultipoleExpansion::operator=( const MultipoleExpansion &other )
{
  this->p=other.p;
  this->assLegFunction = other.assLegFunction;
  this->center = other.center;
  if (_M_n_m != NULL)
    delete [] _M_n_m;
  this->_M_n_m = new std::complex <double>[(this->p+1)*(this->p+2)/2];
  memcpy(this->_M_n_m, other.GetCoeffs(), sizeof(std::complex <double>)*(this->p+1)*(this->p+2)/2);
  this->is_zero = other.is_zero;
  return *this;
}


MultipoleExpansion::~MultipoleExpansion()
{
  if (_M_n_m != NULL)
    delete [] _M_n_m;
}

void MultipoleExpansion::Add(const MultipoleExpansion &multipole,const double sol)

{
  if (multipole.is_zero || sol==0)
    {
    }
  else
    {
      this->is_zero = false;
      for (int n = 0; n < int(this->p)+1 ; n++)
        for (int m = 0; m < n + 1 ; m++)
          {
            std::complex <double> a = multipole.GetCoeff(n,m)*sol;
            this->AddToCoeff(n,m,a);
          }
    }

}

void MultipoleExpansion::Add(const double strength, const dealii::Point<3> &point)

{

  if (strength == 0)
    {
    }
  else
    {
      this->is_zero = false;

      dealii::Point<3> pointRelPos = point + (-1.0*this->center);
      double rho = sqrt(pointRelPos.square());
      double cos_alpha_ = pointRelPos(2)/rho;
      double beta = atan2(pointRelPos(1),pointRelPos(0));
      double P_n_m;
      for (int n = 0; n < int(this->p) + 1 ; n++)
        {
          for (int m = 0; m < n + 1 ; m++)
            {
              P_n_m = this->assLegFunction->GetAssLegFunSph(n,abs(m),cos_alpha_);
              double realFact = P_n_m * pow(rho,double(n)) * strength;
              std::complex <double> a = exp(std::complex<double>(0.,-m*beta))*realFact;
              this->AddToCoeff(n,m,a);
            }
        }
    }

}

void MultipoleExpansion::AddNormDer(const double strength, const dealii::Point<3> &point, const dealii::Tensor<1, 3> &normal)

{
  if (strength == 0)
    {
    }
  else
    {
      this->is_zero = false;

      dealii::Point<3> pointRelPos = point + (-1.0*this->center);
      dealii::Tensor<1,3> normVersor = normal/sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);//normal.square());
      double rho = sqrt(pointRelPos.square());
      double dRhodN = (pointRelPos/rho)*normVersor;
      double beta = atan2(pointRelPos(1),pointRelPos(0));
      double dBetadN = (dealii::Point<3>( -pointRelPos(1), pointRelPos(0), 0.)/(pow(pointRelPos(0),2.)+pow(pointRelPos(1),2.)))*normVersor;

      double cos_alpha_ = pointRelPos(2)/rho;
      double dAlphadN = (dealii::Point<3>( cos_alpha_*cos(beta),cos_alpha_*sin(beta),-sqrt(1.-pow(cos_alpha_,2.)))/rho) * normVersor;


      double P_n_m;
      double dP_n_m_sin;
      for (int n = 0; n < int(this->p) + 1 ; n++)
        {
          for (int m = 0; m < n + 1 ; m++)
            {
              P_n_m =  this->assLegFunction->GetAssLegFunSph(n,abs(m),cos_alpha_);
              dP_n_m_sin = this->assLegFunction->GetAssLegFunSphDeriv(n,abs(m),cos_alpha_)*sqrt(1-pow(cos_alpha_,2.));
              std::complex <double> z = exp(std::complex <double>(0.,-double(m)*beta));
              z *= std::complex <double> (double(n)*pow(rho,double(n)-1.)*P_n_m*dRhodN -
                                          pow(rho,double(n))*dP_n_m_sin*dAlphadN, -double(m)*pow(rho,double(n))*P_n_m*dBetadN);
              z*=strength;

              this->AddToCoeff(n,m,z);
            }
        }
    }

}


void MultipoleExpansion::Add(const MultipoleExpansion &other) //translation of a multipole to its parent center

{
  if (other.is_zero)
    {
      this->is_zero = this->is_zero & other.is_zero;
    }
  else
    {
      this->is_zero = false;
      FullMatrix<double> &A_n_m = this->GetA_n_m();
      if (other.center.distance(this->center) > 1e-7)
        {

          dealii::Point<3> blockRelPos = other.center + (-1.0*this->center);
          double rho = sqrt(blockRelPos.square());
          double cos_alpha_ = blockRelPos(2)/rho;
          double beta = atan2(blockRelPos(1),blockRelPos(0));

          std::complex <double> imUnit(0.,1.);

          double P_nn_mm;

          for (int n = 0; n < int(this->p)+1 ; n++)
            {
              for (int m = 0; m < n+1 ; m++)
                {
                  std::complex <double> z(0.,0.);
                  for (int nn = 0; nn < n+1 ; nn++)
                    {
                      for (int mm = -1*nn; mm < nn+1 ; mm++)
                        {
                          if (abs(m-mm) >  n-nn)
                            {
                            }
                          else
                            {
                              std::complex <double> a = std::complex<double>((other.GetCoeff(abs(n-nn),abs(m-mm))).real(),GSL_SIGN(m-mm)*
                                                                             (other.GetCoeff(abs(n-nn),abs(m-mm))).imag());
                              P_nn_mm =  this->assLegFunction->GetAssLegFunSph(nn,abs(mm),cos_alpha_);
                              double realFact = P_nn_mm * pow(rho,double(nn)) * A_n_m(abs(nn),abs(mm)) *
                                                A_n_m(abs(n-nn),abs(m-mm)) / A_n_m(abs(n),abs(m));
                              realFact *= (pow(imUnit, double(abs(m)-abs(mm)-abs(m-mm)))).real();
                              z += realFact*(a*exp(std::complex <double>(0.,-mm*beta)));
                            }
                        }
                    }
                  this->AddToCoeff(n,m,z);
                }
            }
        }
      else
        {
          for (int n = 0; n < int(this->p)+1 ; n++)
            {
              for (int m = 0; m < n+1 ; m++)
                {
                  this->AddToCoeff(n,m,other.GetCoeff(n,m));
                }
            }
        }
    }
}


double MultipoleExpansion::Evaluate(const dealii::Point<3> &evalPoint)

{
  std::complex <double> fieldValue(0.,0.);
  if (this->is_zero)
    {
    }
  else
    {
      dealii::Point<3> blockRelPos = evalPoint + (-1.0*this->center);
      double rho = sqrt(blockRelPos.square());
      double cos_alpha_ = blockRelPos(2)/rho;
      double beta = atan2(blockRelPos(1),blockRelPos(0));

      double P_n_m;

      for (int n = 0; n < int(this->p)+1 ; n++)
        {
          P_n_m =  this->assLegFunction->GetAssLegFunSph(n,0,cos_alpha_);
          double realFact = P_n_m * pow(rho,double(-n-1));
          fieldValue += this->GetCoeff(n, 0)*realFact;
          for (int m = 1; m < n + 1 ; m++)
            {
              P_n_m =  this->assLegFunction->GetAssLegFunSph(n,abs(m),cos_alpha_);
              double realFact = P_n_m * pow(rho,double(-n-1));
              fieldValue += this->GetCoeff(n, m)*exp(std::complex <double>(0., m*beta))*2.*realFact;
            }
        }
    }
  return fieldValue.real();

}
