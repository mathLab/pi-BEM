#include <iostream>

#include "multipole_expansion.h"
#include<deal.II/base/point.h>

#define GSL_SIGN(x) (x<0 ? -1: (x>0 ? 1: 0))

FullMatrix<double> MultipoleExpansion::A_n_m = MultipoleExpansion::A_n_m_Matrix(20);

MultipoleExpansion::MultipoleExpansion()

{
	
this->p = 0;
this->center = Point<3>(0,0,0);
this->assLegFunction = NULL;


M_n_m.resize(p+1);

for (unsigned int i = 0; i < p+1; i++) {

(M_n_m.at(i)).resize(i+1);	
	
}



for (int n = 0; n < int(p)+1 ; n++) {

	for (int m = 0; m < n + 1 ; m++) {
		
		(M_n_m.at(n)).at(m) = std::complex<double>(0.0,0.0);
		
		}		
		
	}

}


MultipoleExpansion::MultipoleExpansion(const unsigned int order, const dealii::Point<3> center, const AssLegFunction *assLegFunction)

{
	
this->p = order;
this->center = center;
this->assLegFunction = assLegFunction;


M_n_m.resize(p+1);

for (unsigned int i = 0; i < p+1; i++) {

(M_n_m.at(i)).resize(i+1);	
	
}



for (int n = 0; n < int(p)+1 ; n++) {

	for (int m = 0; m < n + 1 ; m++) {
		
		(M_n_m.at(n)).at(m) = std::complex<double>(0.0,0.0);
		
		}		
		
	}

}


MultipoleExpansion::~MultipoleExpansion()
{


}

void MultipoleExpansion::Add(const MultipoleExpansion *multipole,const double sol)

{
for (int n = 0; n < int(this->p)+1 ; n++)
	
	{
	for (int m = 0; m < n + 1 ; m++)
		
		{
		(this->M_n_m.at(n)).at(m) += multipole->GetM_n_m(n,m)*sol;
		
		}		
	}
			
}

void MultipoleExpansion::Add(const double strength, const dealii::Point<3> point)

{
  dealii::Tensor<1,3> pointRelPos = point - this->center;
		double rho = pointRelPos.norm();
		double cos_alpha_ = pointRelPos[2]/rho;
		double beta = atan2(pointRelPos[1],pointRelPos[0]);
		
 		double P_n_m;	
 		for (int n = 0; n < int(this->p) + 1 ; n++)
 			{
            		for (int m = 0; m < n + 1 ; m++)
				{
				//double P_n_m = gsl_matrix_get (P, n, abs(m));
            	P_n_m = this->assLegFunction->GetAssLegFunSph(n,abs(m),cos_alpha_);
            	double realFact = P_n_m * pow(rho,double(n)) * strength; 
				(this->M_n_m.at(n)).at(m) += exp(std::complex<double>(0.,-m*beta))*realFact;
				
				}
			}

}

void MultipoleExpansion::AddNormDer(const double strength, const dealii::Point<3> point, const dealii::Point<3> normal)

{
  dealii::Tensor<1,3> pointRelPos = point - this->center;
  dealii::Point<3> normVersor = normal/normal.norm();
		double rho = pointRelPos.norm();
		double dRhodN = (pointRelPos/rho)*normVersor;
		double beta = atan2(pointRelPos[1],pointRelPos[0]);
 		double dBetadN = (dealii::Point<3>( -pointRelPos[1], pointRelPos[0], 0.)/(pow(pointRelPos[0],2.)+pow(pointRelPos[1],2.)))*normVersor;

		double cos_alpha_ = pointRelPos[2]/rho;
		double dAlphadN = (dealii::Point<3>( cos_alpha_*cos(beta),cos_alpha_*sin(beta),-sqrt(1.-pow(cos_alpha_,2.)))/rho) * normVersor;								        
 		
 		
		double P_n_m;
		double dP_n_m_sin; 					
 		for (int n = 0; n < int(this->p) + 1 ; n++)
 			{
            		for (int m = 0; m < n + 1 ; m++)
				{
            			P_n_m =  this->assLegFunction->GetAssLegFunSph(n,abs(m),cos_alpha_);
            			dP_n_m_sin = this->assLegFunction->GetAssLegFunSphDeriv(n,abs(m),cos_alpha_)*sqrt(1-pow(cos_alpha_,2.));
            			std::complex <double> eim = exp(std::complex <double>(0.,-double(m)*beta));	
				std::complex <double> z = std::complex <double> (double(n)*pow(rho,double(n)-1.)*P_n_m*dRhodN - pow(rho,double(n))*dP_n_m_sin*dAlphadN, -double(m)*pow(rho,double(n))*P_n_m*dBetadN);
				(this->M_n_m.at(n)).at(m) += strength*z*eim;
				
				}
			}


}


void MultipoleExpansion::Add(const MultipoleExpansion *child) //translation of a multipole to its parent center

{
		FullMatrix<double> &A_n_m = this->GetA_n_m();
		dealii::Tensor<1,3> blockRelPos = child->center - this->center;
		double rho = blockRelPos.norm();
		double cos_alpha_ = blockRelPos[2]/rho;
		double beta = atan2(blockRelPos[1],blockRelPos[0]);

		std::complex <double> imUnit(0.,1.);
		
		double P_nn_mm;
		
	 	for (int n = 0; n < int(this->p)+1 ; n++) {			

			for (int m = 0; m < n+1 ; m++) {

				std::complex <double> z(0.,0.);

				for (int nn = 0; nn < n+1 ; nn++) {

					for (int mm = -1*nn; mm < nn+1 ; mm++) {

						if (abs(m-mm) >  n-nn)
							{
							}

						else
							{
							std::complex <double> a = std::complex<double>((child->GetM_n_m(abs(n-nn),abs(m-mm))).real(),GSL_SIGN(m-mm)*(child->GetM_n_m(abs(n-nn),abs(m-mm))).imag());
							P_nn_mm =  this->assLegFunction->GetAssLegFunSph(nn,abs(mm),cos_alpha_);
							double realFact = P_nn_mm * pow(rho,double(nn)) * A_n_m(abs(nn),abs(mm)) * A_n_m(abs(n-nn),abs(m-mm)) / A_n_m(abs(n),abs(m));
							realFact *= (pow(imUnit, double(abs(m)-abs(mm)-abs(m-mm)))).real();
							z += realFact*(a*exp(std::complex <double>(0.,-mm*beta)));

							}
					}
					
				}

			(this->M_n_m.at(n)).at(m) += z;

			}
			

		}

}


double MultipoleExpansion::Evaluate(const dealii::Point<3> evalPoint)

{

		std::complex <double> fieldValue(0.,0.);
		dealii::Tensor<1,3> blockRelPos = evalPoint - this->center;
		double rho = blockRelPos.norm();
		double cos_alpha_ = blockRelPos[2]/rho;
		double beta = atan2(blockRelPos[1],blockRelPos[0]);
		
		double P_n_m;

		for (int n = 0; n < int(this->p)+1 ; n++)	
			{
	
	 		P_n_m =  this->assLegFunction->GetAssLegFunSph(n,0,cos_alpha_);
	 		double realFact = P_n_m * pow(rho,double(-n-1));
	 		fieldValue += this->GetM_n_m(n, 0)*realFact;		
	
     			for (int m = 1; m < n + 1 ; m++)
				{
	
				P_n_m =  this->assLegFunction->GetAssLegFunSph(n,abs(m),cos_alpha_);
				double realFact = P_n_m * pow(rho,double(-n-1));		
				fieldValue += this->GetM_n_m(n, m)*exp(std::complex <double>(0., m*beta))*2.*realFact;
				}		
			}

		return fieldValue.real();

}


