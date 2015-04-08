//----------------------------  step-34.cc  ---------------------------
//    $Id: step-34.cc 18734 2009-04-25 13:36:48Z heltai $
//    Version: $Name$ 
//
//    Copyright (C) 2009 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//    Authors: Luca Heltai, Cataldo Manigrasso
//
//----------------------------  step-34.cc  ---------------------------


				 // @sect3{Include files}

				 // The program starts with including a bunch
				 // of include files that we will use in the
				 // various parts of the program. Most of them
				 // have been discussed in previous tutorials
				 // already:

#include<deal.II/base/smartpointer.h>
#include<deal.II/base/point.h>
#include<deal.II/base/utilities.h>

				 // And here are a few C++ standard header
				 // files that we will need:
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <set>


namespace LaplaceKernel
{
  template <int dim>
  double single_layer(const Point<dim> &R)
  {
    switch(dim)
      {
      case 2:
	return (-std::log(R.norm()) / (2*numbers::PI) );
	
      case 3:
	return (1./( R.norm()*4*numbers::PI ) );
	
      default:
	Assert(false, ExcInternalError());
	return 0.;
      }
  }
        


  template <int dim>
  Point<dim> double_layer(const Point<dim> &R)
  {
    switch(dim)
      {
      case 2:
	return R / ( -2*numbers::PI * R.square());
      case 3:
	return R / ( -4*numbers::PI * R.square() * R.norm() );
	
      default:
	Assert(false, ExcInternalError());
	return Point<dim>();
      }
  }
  
    template <int dim>
    void kernels(const Tensor<1, dim> &R, Tensor<1,dim> &D, double &d)
  {
    double r = R.norm();
    double r2 = r*r;
    switch(dim)
      {
      case 2:
	d = -std::log(r) / (2*numbers::PI);
	D = R / ( -2*numbers::PI * r2);
    break;
      case 3:
	d = (1./( r*4*numbers::PI ) );
	D = R / ( -4*numbers::PI * r2*r );
    break;
      default:
	Assert(false, ExcInternalError());
      }
  }
  
}
