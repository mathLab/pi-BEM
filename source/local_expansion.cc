#include <iostream>
#include <math.h>


#include "local_expansion.h"

#define GSL_SIGN(x) (x<0 ? -1: (x>0 ? 1: 0))

FullMatrix<double> LocalExpansion::A_n_m = LocalExpansion::A_n_m_Matrix(20);

LocalExpansionCoeff
LocalExpansion::mExp_to_lExp_Coeff = LocalExpansion::mExp_to_lExp_Coeff_Build(LocalExpansion::A_n_m, 10);

std::vector <std::vector <std::map <int,std::map <int,double> > > >
LocalExpansion::lExp_to_lExp_Coeff = LocalExpansion::lExp_to_lExp_Coeff_Build(LocalExpansion::A_n_m, 10);

LocalExpansion::LocalExpansion()

{

  this->p = 0;
  this->center = Point<3>(0,0,0);
  this->assLegFunction = NULL;
  this->_L_n_m = NULL;
  this->is_zero = true;
}


LocalExpansion::LocalExpansion(const unsigned int order, const dealii::Point<3> &center, const AssLegFunction *assLegFunction)

{

  this->p = order;
  this->center = center;
  this->assLegFunction = assLegFunction;

  this->_L_n_m= new std::complex <double>[(this->p+1)*(this->p+2)/2];
  for (unsigned int i=0; i<(this->p+1)*(this->p+2)/2; ++i)
    this->_L_n_m[i] = std::complex <double>(0.0,0.0);
  this->is_zero = true;
}

LocalExpansion::LocalExpansion(const LocalExpansion &other)
{
  this->p=other.p;
  this->assLegFunction = other.assLegFunction;
  this->lExp_to_lExp_Coeff = other.lExp_to_lExp_Coeff;
  this->mExp_to_lExp_Coeff = other.mExp_to_lExp_Coeff;
  this->center = other.center;
  this->_L_n_m = new std::complex <double>[(this->p+1)*(this->p+2)/2];
  memcpy(this->_L_n_m, other.GetCoeffs(), sizeof(std::complex <double>)*(this->p+1)*(this->p+2)/2);
  this->is_zero = other.is_zero;
}

LocalExpansion::~LocalExpansion()
{
  if (_L_n_m != NULL)
    delete [] _L_n_m;

}

LocalExpansion &LocalExpansion::operator=( const LocalExpansion &other )
{
  this->p=other.p;
  this->assLegFunction = other.assLegFunction;
  this->lExp_to_lExp_Coeff = other.lExp_to_lExp_Coeff;
  this->mExp_to_lExp_Coeff = other.mExp_to_lExp_Coeff;
  this->center = other.center;
  this->_L_n_m = new std::complex <double>[(this->p+1)*(this->p+2)/2];
  memcpy(this->_L_n_m, other.GetCoeffs(), sizeof(std::complex <double>)*(this->p+1)*(this->p+2)/2);
  this->is_zero = other.is_zero;
  return *this;
}



void LocalExpansion::Add(const std::vector <double> &real, const std::vector <double> &imag)

{

  unsigned int count = 0;
  double sum = 0.0;
  for (unsigned int m = 0; m < this->p+1 ; m++)
    {
      for (unsigned int n = m; n < this->p + 1 ; n++)
        {
          std::complex <double> a(real.at(count),imag.at(count));
          sum += norm(a);
          this->AddToCoeff(n,n+m,a);
          a=std::conj(a);
          this->AddToCoeff(n,n-m,a);
          count++;
        }
    }
  if (sum > 1e-20)
    this->is_zero = false;


}



void LocalExpansion::Add(const LocalExpansion &other) // translation of local expansion

{
  if (other.is_zero)
    {
    }
  else
    {
      unsigned int p = this->p;
      if (other.center.distance(this->center) > 1e-7)
        {

          dealii::Point<3> blockRelPos = other.GetCenter() + (-1.0*this->center);
          double rho = sqrt(blockRelPos.square());
          double cos_alpha_ = blockRelPos(2)/rho;
          double beta = atan2(blockRelPos(1),blockRelPos(0));

          double P_nn_mm;

          for (int n = 0; n < int(p)+1 ; n++)
            {
              for (int m = 0; m < n+1 ; m++)
                {
                  std::complex <double> z = std::complex <double>(0.,0.);
                  for (int nn = n; nn < int(p)+1 ; nn++)
                    {
                      double rhoFact = pow(rho,double(nn-n));
                      for (int mm = -1*nn; mm < nn+1 ; mm++)
                        {
                          if (abs(mm-m) >  nn-n)
                            {
                            }
                          else
                            {
                              std::complex <double> a = std::complex <double>(other.GetCoeff(abs(nn),abs(mm)).real(),
                                                                              GSL_SIGN(mm)*other.GetCoeff(abs(nn),abs(mm)).imag());
                              P_nn_mm =  this->assLegFunction->GetAssLegFunSph(nn-n,abs(mm-m),cos_alpha_);
                              double realFact = P_nn_mm * rhoFact * lExp_to_lExp_Coeff[n][m][nn][mm];
                              z += a*std::complex<double>(cos((mm-m)*beta),sin((mm-m)*beta))*realFact;
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
      this->is_zero = false;
    }
}


void LocalExpansion::Add(const MultipoleExpansion &multipole) // multipole conversion into local expansion, and addition to the rest

{
//static unsigned int call_count=0;

  if (multipole.is_zero)
    {
    }
  else
    {
      //cout<<call_count<<" "<<std::setprecision(25)<<multipole.GetCoeff(0,0)<<endl;
      dealii::Point<3> blockRelPos = multipole.GetCenter() + (-1.0*this->center);
      double rho = sqrt(blockRelPos.square());
      double cos_alpha_ = blockRelPos(2)/rho;
      double beta = atan2(blockRelPos(1),blockRelPos(0));

      double P_nn_mm;
      std::complex <double> a;
      double realFact;
      for (int n = 0; n < int(this->p)+1 ; n++)
        {
          for (int m = 0; m < n+1 ; m++)
            {
              std::complex <double> z = std::complex <double>(0.,0.);
              for (int nn = 0; nn < int(this->p)+1 ; nn++)
                {
                  double rhoFact = pow(rho,double(-n-nn-1));
                  for (int mm = -1*nn; mm < 0 ; mm++)
                    {
                      a = multipole.GetCoeff(nn,abs(mm));
                      a = std::complex <double>(a.real(),-a.imag());
                      P_nn_mm =  this->assLegFunction->GetAssLegFunSph(nn+n,abs(mm-m),cos_alpha_);
                      realFact = P_nn_mm * rhoFact * mExp_to_lExp_Coeff.get(n,m,nn,mm);
                      z += a*std::complex <double>(cos((mm-m)*beta),sin((mm-m)*beta))*realFact;
                    }
                  for (int mm = 0; mm < nn+1 ; mm++)
                    {
                      a = multipole.GetCoeff(nn,abs(mm));
                      P_nn_mm =  this->assLegFunction->GetAssLegFunSph(nn+n,abs(mm-m),cos_alpha_);
                      double realFact = P_nn_mm * rhoFact * mExp_to_lExp_Coeff.get(n,m,nn,mm);
                      z += a*std::complex <double>(cos((mm-m)*beta),sin((mm-m)*beta))*realFact;
                    }
                }
              //cout<<call_count<<":   "<<z<<" ("<<a<<")"<<endl;
              this->AddToCoeff(n,m,z);
            }
        }
      //call_count++;
      this->is_zero = false;
    }
}


double LocalExpansion::Evaluate(const dealii::Point<3> &evalPoint)
{
  std::complex <double> fieldValue = std::complex <double>(0.,0.);
  if (this->is_zero)
    {
    }
  else
    {

      unsigned int p = this->p;
      dealii::Point<3> blockRelPos = evalPoint + (-1.0*this->center);
      double rho = sqrt(blockRelPos.square());
      double cos_alpha_ = blockRelPos(2)/rho;
      double beta = atan2(blockRelPos(1),blockRelPos(0));

      double P_n_m;
      for (int n = 0; n < int(p) + 1 ; n++)
        {
          P_n_m =  this->assLegFunction->GetAssLegFunSph(n,0,cos_alpha_);
          double realFact = P_n_m * pow(rho,double(n));
          fieldValue += this->GetCoeff(n, 0)*realFact;
          for (int m = 1; m < n + 1 ; m++)
            {
              P_n_m =  this->assLegFunction->GetAssLegFunSph(n,m,cos_alpha_);
              double realFact = P_n_m * pow(rho,double(n));
              //std::complex <double> complexFact = exp(std::complex <double>(0., 1.*m*beta))*2.*realFact;
              std::complex <double> complexFact = std::complex<double>(cos(m*beta),sin(m*beta))*2.*realFact;
              fieldValue += this->GetCoeff(n, m)*complexFact;
            }
        }
    }

  return fieldValue.real();

}


