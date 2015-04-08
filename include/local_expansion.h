#ifndef LOCALEXPANSION_H_
#define LOCALEXPANSION_H_

#include <math.h>
#include <string>
#include <vector>
#include <deal.II/lac/full_matrix.h>


#include "multipole_expansion.h"
#include "ass_leg_function.h"



class LocalExpansion
{
	
public:
		static FullMatrix<double> A_n_m;
		
		static std::vector <std::vector <std::vector <std::map <int,double > > > > mExp_to_lExp_Coeff;
		
		static std::vector <std::vector <std::map <int,std::map <int,double > > > > lExp_to_lExp_Coeff;	
private:

		mutable unsigned int p;
		
		mutable std::vector <std::vector <std::complex <double> > > L_n_m;
		
		mutable dealii::Point<3> center;
		
		mutable const AssLegFunction *assLegFunction;


public:
		LocalExpansion();
		
		LocalExpansion(const unsigned int order, const  dealii::Point<3> center, const AssLegFunction *assLegFunction);
		
		~LocalExpansion();

		void Add(const std::vector <double> real, const std::vector <double> imag);
		
		
		void Add(const LocalExpansion *parent);
		
		void Add(const MultipoleExpansion *multipole);
		
		double Evaluate(const dealii::Point<3> evalPoint);
		
		inline std::complex <double> GetL_n_m(const unsigned int n, const int m) const
			{return this->L_n_m.at(n).at(m);}
		
		inline dealii::Point<3> GetCenter() const
			{return this->center;}
		
		inline FullMatrix<double> &GetA_n_m() const 		
		    {return this->A_n_m;}	
		
	    static FullMatrix<double> & A_n_m_Matrix(unsigned int dimension)
	    	{
	    	A_n_m.reinit(dimension+1,dimension+1);
	    	for (unsigned int n = 0; n < dimension+1 ; n++)
	    			
	    			{	
	    				for (unsigned int m = 0; m < n+1 ; m++)		

	    					{
	    					double f1 = 1.;
	    					double f2 = 1.;

	    					for(int ii = n-m; ii > 0; ii-- )
	    						f1 *= ii;

	    					for(int ii = n+m; ii > 0; ii-- )
	    						f2 *= (ii);

	    					A_n_m(n, m) = pow(-1.,double(n))/sqrt(f1*f2); 
	    					}   
	    			}
	    	return A_n_m;
	    	}

	    static  std::vector <std::vector <std::vector <std::map <int,double > > > >
	    mExp_to_lExp_Coeff_Build(A_n_m(n,m), unsigned int p)
	    	{
	    	
		std::complex <double> imUnit = std::complex <double> (0.,1.);
    		std::vector <std::vector <std::vector <std::map <int,double > > > > realCoeff;
		realCoeff.resize(p+1);
	 	for (int n = 0; n < int(p)+1 ; n++) {			
			realCoeff[n].resize(n+1);
			for (int m = 0; m < n+1 ; m++) {
				realCoeff[n][m].resize(p+1);
				for (int nn = 0; nn < int(p)+1 ; nn++) {
					for (int mm = -1*nn; mm < nn+1 ; mm++) {
                                                
						double realFact = (*gsl_matrix_ptr(A_n_m,nn,abs(mm))) / (*gsl_matrix_ptr(A_n_m,n+nn,abs(m-mm))) * (*gsl_matrix_ptr(A_n_m,n,abs(m)));
						realFact *= (pow(imUnit, double(abs(m-mm)-abs(m)-abs(mm)))).real()/pow(-1.,nn);
						realCoeff[n][m][nn][mm] = realFact;
                                                	
					}
					
				}

			}
			

		}		
	    	
		return realCoeff;
	    	}

	    static  std::vector <std::vector <std::map <int,std::map <int,double > > > >
	    lExp_to_lExp_Coeff_Build(gsl_matrix * A_n_m, unsigned int p)
	    	{
	    	
		std::complex <double> imUnit = std::complex <double> (0.,1.);
    		std::vector <std::vector <std::map <int,std::map <int,double > > > > realCoeff;
		realCoeff.resize(p+1);
	 	for (int n = 0; n < int(p)+1 ; n++) {			
			realCoeff[n].resize(n+1);
			for (int m = 0; m < n+1 ; m++) {
				for (int nn = n; nn < int(p)+1 ; nn++) {
					for (int mm = -1*nn; mm < nn+1 ; mm++) {
                                                
						double realFact = (*gsl_matrix_ptr(A_n_m,nn-n,abs(mm-m))) / (*gsl_matrix_ptr(A_n_m,nn,abs(mm))) * (*gsl_matrix_ptr(A_n_m,n,abs(m)));
						realFact *= (pow(imUnit, double(abs(mm)-abs(mm-m)-abs(m)))).real()*pow(-1.,nn+n);
						realCoeff[n][m][nn][mm] = realFact;
                                                	
					}
					
				}

			}
			

		}
		
		
		
	    	return realCoeff;
	    	}
		
};





#endif /*LOCALEXPANSION_H_*/
