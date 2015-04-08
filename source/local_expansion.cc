#include <iostream>
#include <gsl/gsl_sf_legendre.h>
#include <math.h>
#include <gsl/gsl_math.h>


#include "local_expansion.h"


gsl_matrix* LocalExpansion::A_n_m = LocalExpansion::A_n_m_Matrix(20);

std::vector <std::vector <std::vector <std::map <int,double> > > >
            LocalExpansion::mExp_to_lExp_Coeff = LocalExpansion::mExp_to_lExp_Coeff_Build(LocalExpansion::A_n_m, 10);
	    
std::vector <std::vector <std::map <int,std::map <int,double> > > >
            LocalExpansion::lExp_to_lExp_Coeff = LocalExpansion::lExp_to_lExp_Coeff_Build(LocalExpansion::A_n_m, 10);
	    
LocalExpansion::LocalExpansion()

{
	
this->p = 0;
this->center = Point<3>(0,0,0);
this->assLegFunction = NULL;

L_n_m.resize(this->p+1);

for (unsigned int i = 0; i < this->p+1; i++) {

(L_n_m.at(i)).resize(i+1);	
	
}



for (int n = 0; n < int(this->p)+1 ; n++) {

	for (int m = 0; m < n + 1 ; m++) {
		
		(L_n_m.at(n)).at(m) = std::complex <double>(0.,0.);
		
		}		
		
	}


}


LocalExpansion::LocalExpansion(const unsigned int order, const dealii::Point<3> center, const AssLegFunction *assLegFunction)

{
	
this->p = order;
this->center = center;
this->assLegFunction = assLegFunction;

L_n_m.resize(this->p+1);

for (unsigned int i = 0; i < this->p+1; i++) {

(L_n_m.at(i)).resize(i+1);	
	
}



for (int n = 0; n < int(this->p)+1 ; n++) {

	for (int m = 0; m < n + 1 ; m++) {
		
		(L_n_m.at(n)).at(m) = std::complex <double>(0.,0.);
		
		}		
		
	}


}


LocalExpansion::~LocalExpansion()
{


}

void LocalExpansion::Add(const std::vector <double> real, const std::vector <double> imag)

{

unsigned int count = 0;

for (unsigned int m = 0; m < this->p+1 ; m++) {

	for (unsigned int n = m; n < this->p + 1 ; n++) {

		(L_n_m.at(n)).at(n+m) += std::complex <double>(real.at(count),imag.at(count));
		(L_n_m.at(n)).at(n-m) = (L_n_m.at(n)).at(n+m)+std::complex <double>(real.at(count),-1*imag.at(count));

		count++;
		}		
	}

}



void LocalExpansion::Add(const LocalExpansion *parent) // translation of local expansion

{
		gsl_matrix *A_n_m = this->GetA_n_m();
		unsigned int p = this->p;
		dealii::Point<3> blockRelPos = parent->GetCenter() - this->center;
		double rho = sqrt(blockRelPos.square());
		double cos_alpha_ = blockRelPos(2)/rho;
		double beta = atan2(blockRelPos(1),blockRelPos(0));
		
		std::complex <double> imUnit = std::complex <double>(0.,1.);
		
		double P_nn_mm;
		
	 	for (int n = 0; n < int(p)+1 ; n++) {			

			for (int m = 0; m < n+1 ; m++) {

				std::complex <double> z = std::complex <double>(0.,0.);

				for (int nn = n; nn < int(p)+1 ; nn++) {
					double rhoFact = pow(rho,double(nn-n));
					for (int mm = -1*nn; mm < nn+1 ; mm++) {
						
						if (abs(mm-m) >  nn-n)
							{
							}

						else
							{
							std::complex <double> a = std::complex <double>(parent->GetL_n_m(abs(nn),abs(mm)).real(),GSL_SIGN(mm)*parent->GetL_n_m(abs(nn),abs(mm)).imag());
							P_nn_mm =  this->assLegFunction->GetAssLegFunSph(nn-n,abs(mm-m),cos_alpha_);
							double realFact = P_nn_mm * rhoFact * lExp_to_lExp_Coeff[n][m][nn][mm];
							z += a*exp(std::complex <double>(0.,(mm-m)*beta))*realFact;
							}
					
					}
					
				}

			(this->L_n_m.at(n)).at(m) += z;

			}
			

		}

}


void LocalExpansion::Add(const MultipoleExpansion *multipole) // multipole conversion into local expansion, and addition to the rest

{
		gsl_matrix *A_n_m = this->GetA_n_m();
		dealii::Point<3> blockRelPos = multipole->GetCenter() - this->center;
		double rho = sqrt(blockRelPos.square());
		double cos_alpha_ = blockRelPos(2)/rho;
		double beta = atan2(blockRelPos(1),blockRelPos(0));
		
		std::complex <double> imUnit = std::complex <double> (0.,1.);

		double P_nn_mm;    	
	 	for (int n = 0; n < int(this->p)+1 ; n++) {			

			for (int m = 0; m < n+1 ; m++) {

				std::complex <double> z = std::complex <double>(0.,0.);

				for (int nn = 0; nn < int(this->p)+1 ; nn++) {
					double rhoFact = pow(rho,double(-n-nn-1));
					for (int mm = -1*nn; mm < nn+1 ; mm++) {

						std::complex <double> a = std::complex <double>((multipole->GetM_n_m(abs(nn),abs(mm))).real(),GSL_SIGN(mm)*(multipole->GetM_n_m(abs(nn),abs(mm))).imag());
						P_nn_mm =  this->assLegFunction->GetAssLegFunSph(nn+n,abs(mm-m),cos_alpha_);
						double realFact = P_nn_mm * rhoFact * mExp_to_lExp_Coeff[n][m][nn][mm];
						z += a*exp(std::complex <double>(0.,(mm-m)*beta))*realFact;
						
							
					}
					
				}

			(this->L_n_m.at(n)).at(m) += z;
			}
			

		}


}


double LocalExpansion::Evaluate(const dealii::Point<3> evalPoint)

{


		unsigned int p = this->p;	
		std::complex <double> fieldValue = std::complex <double>(0.,0.);
		dealii::Point<3> blockRelPos = evalPoint - this->center;
		double rho = sqrt(blockRelPos.square());
		double cos_alpha_ = blockRelPos(2)/rho;
		double beta = atan2(blockRelPos(1),blockRelPos(0));

		double P_n_m;
		for (int n = 0; n < int(p) + 1 ; n++)
			{
			P_n_m =  this->assLegFunction->GetAssLegFunSph(n,0,cos_alpha_);
			double realFact = P_n_m * pow(rho,double(n));
			fieldValue += this->GetL_n_m(n, 0)*realFact;		
	  
			for (int m = 1; m < n + 1 ; m++)
				{
				P_n_m =  this->assLegFunction->GetAssLegFunSph(n,m,cos_alpha_);
				double realFact = P_n_m * pow(rho,double(n));
				std::complex <double> complexFact = exp(std::complex <double>(0., 1.*m*beta))*2.*realFact;
				fieldValue += this->GetL_n_m(n, m)*complexFact;
				}		
			}



		return fieldValue.real();



}


