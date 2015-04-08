#ifndef MULTIPOLE_EXPANSION_H_
#define MULTIPOLE_EXPANSION_H_

#include <math.h>
#include <string>
#include <vector>
#include <complex>

#include<deal.II/base/point.h>

#include "ass_leg_function.h"

using namespace dealii;

typedef gsl_matrix double;

class MultipoleExpansion
{
public: 
	static gsl_matrix* A_n_m;
	
private:

		mutable unsigned int p;
		
		mutable std::vector <std::vector <std::complex<double> > > M_n_m;
		
		mutable dealii::Point<3> center;
		
		mutable const AssLegFunction *assLegFunction;
		

public:

	        MultipoleExpansion();
		
		MultipoleExpansion(const unsigned int order,const  dealii::Point<3> center,const AssLegFunction *assLegFunction);

		~MultipoleExpansion();

		void Add(const MultipoleExpansion *multipole,const double sol);
		
		void Add(const double strength, const dealii::Point<3> point);
		
		void Add(const MultipoleExpansion *child);
		
		void AddNormDer(const double strength, const dealii::Point<3> point, const dealii::Point<3> normal);
		
		double Evaluate(const dealii::Point<3> evalPoint);
		
		inline std::complex <double> GetM_n_m(const unsigned int n, const int m) const
			{return this->M_n_m.at(n).at(m);}
			
		inline dealii::Point<3> GetCenter() const
			{return this->center;}
		
		inline gsl_matrix * GetA_n_m() const		
		    {return this->A_n_m;}
				
	    	static gsl_matrix * A_n_m_Matrix(unsigned int dim)
	    		{
	    		gsl_matrix* A_n_m = gsl_matrix_calloc(dim+1,dim+1);
	    		for (unsigned int n = 0; n < dim+1 ; n++)
	    			
	    			{	
	    				for (unsigned int m = 0; m < n+1 ; m++)		

	    					{
	    					double f1 = 1.;
	    					double f2 = 1.;

	    					for(int ii = n-m; ii > 0; ii-- )
	    						f1 *= ii;

	    					for(int ii = n+m; ii > 0; ii-- )
	    						f2 *= (ii);

	    					*gsl_matrix_ptr (A_n_m, n, m) = pow(-1.,double(n))/sqrt(f1*f2); 
	    					}   
	    			}
	    		return A_n_m;
	    		}


};			



#endif /*MULTIPOLE_EXPANSION_H_*/
