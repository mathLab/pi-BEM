//-----------------------------------------------------------
//
//    Copyright (C) 2014 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------


#include "tests.h"
#include "multipole_expansion.h"

int main ()
{
  initlog();

  // we start creating a few charges located at the points P_i,
  // and having intensities a_i
  // the charges are roughly located around point P(3,3,3)
  // since we want to test the expansion of the neumann kernel
  // \partial G / \n = a * r * n / |r|^3
  // we also create the normals n_i corresponding to each charge
  std::vector<Point<3> > charges_locations;
  std::vector<double> charges_intensities;
  std::vector<Point<3> > normals_at_charges_locations;

  charges_locations.push_back(Point<3>(3.0-0.2,3.0,3.0));
  charges_intensities.push_back(0.1);
  normals_at_charges_locations.push_back(Point<3>(1.0,0.0,0.0));
  charges_locations.push_back(Point<3>(3.0,3.0-0.2,3.0));
  charges_intensities.push_back(0.4);
  normals_at_charges_locations.push_back(Point<3>(1.0,0.0,0.0));
  charges_locations.push_back(Point<3>(3.0+0.2,3.0,3.0-0.2));
  charges_intensities.push_back(1.0);
  normals_at_charges_locations.push_back(Point<3>(1.0,0.0,0.0));
  charges_locations.push_back(Point<3>(3.0+0.2,3.0,3.0));
  charges_intensities.push_back(0.2);
  normals_at_charges_locations.push_back(Point<3>(1.0,0.0,0.0));
  charges_locations.push_back(Point<3>(3.0-0.3,3.0,3.0+0.2));
  charges_intensities.push_back(0.7);
  normals_at_charges_locations.push_back(Point<3>(1.0,0.0,0.0));

  // we then avaluate the exact potential phi = sum_i (a_i*(P_i-Q)*n_i/|P_i-Q|^3) evaluated
  // at point Q(0,0,0)
  Point<3> zero(0.0,0.0,0.0);
  double exact_potential_normal_gradient = 0.0;
  for (unsigned int i=0; i<charges_locations.size(); ++i)
    {
      exact_potential_normal_gradient += -charges_intensities[i]*(charges_locations[i]*normals_at_charges_locations[i])/pow(charges_locations[i].distance(zero),3.0);
      deallog<<"Charge "<<i+1<<" of intensity "<<charges_intensities[i]<<" is located at point ("<<charges_locations[i]<<
             ")  and has normal ("<<normals_at_charges_locations[i]<<")"<<std::endl;
    }

  deallog<<"Exact potential normal gradient at point Q(0,0,0): "<<exact_potential_normal_gradient<<std::endl;


  // we now try to create a multipole expansion centered in P and use it
  // to use it to approximated the potential in Q due to the charges in P_i

  // we define the center, order of truncation and associated laplace functions
  // needed to construct the MultipoleExpansion class object
  Point<3> center(3.0,3.0,3.0);
  AssLegFunction *alf_ptr = new AssLegFunction();
  unsigned int truncation_order = 6;
  MultipoleExpansion multipole(truncation_order, center, alf_ptr);

  // the MultipoleExpansion created so far is ready to work, but it's empty
  // we use the AddNormDer(const double strength, const dealii::Point<3> &point, const dealii::Point<3> &normal)
  // member so that the expansion can account of the effect of each charge
  // considered
  for (unsigned int i=0; i<charges_locations.size(); ++i)
    multipole.AddNormDer(charges_intensities[i], charges_locations[i], normals_at_charges_locations[i]);

  // the only thing left to do is use the expansion to evaluate the potential
  // at point Q. This is done using the Evaluate(const dealii::Point<3> &evalPoint)
  // member function.
  double approx_potential = multipole.Evaluate(zero);

  deallog<<"Approx potential at point Q(0,0,0): "<<approx_potential<<std::endl;

  deallog<<"The absolute multipole approximation error is: "<<fabs(exact_potential_normal_gradient-approx_potential)<<std::endl;

  deallog<<"The relative multipole approximation error is: "<<fabs(exact_potential_normal_gradient-approx_potential)/fabs(exact_potential_normal_gradient)<<std::endl;

}
