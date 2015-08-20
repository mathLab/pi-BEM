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
#include "local_expansion.h"

int main ()
{
  initlog();

  // we start creating a few charges located at the points P_i,
  // and having intensities a_i
  // the charges are roughly located around point P(3,3,3)
  std::vector<Point<3> > charges_locations;
  std::vector<double> charges_intensities;

  charges_locations.push_back(Point<3>(3.0-0.2,3.0,3.0));
  charges_intensities.push_back(0.1);
  charges_locations.push_back(Point<3>(3.0,3.0-0.2,3.0));
  charges_intensities.push_back(0.4);
  charges_locations.push_back(Point<3>(3.0,3.0,3.0-0.2));
  charges_intensities.push_back(1.0);
  charges_locations.push_back(Point<3>(3.0+0.2,3.0,3.0));
  charges_intensities.push_back(0.2);
  charges_locations.push_back(Point<3>(3.0,3.0,3.0+0.2));
  charges_intensities.push_back(0.7);


  // we then avaluate the exact potential phi = sum_i (a_i/|P_i-Q|) evaluated
  // at point Q(0,0,0)
  Point<3> zero(0.0,0.0,0.0);
  double exact_potential = 0.0;
  for (unsigned int i=0; i<charges_locations.size(); ++i)
    {
      exact_potential += charges_intensities[i]/charges_locations[i].distance(zero);
      deallog<<"Charge "<<i+1<<" of intensity "<<charges_intensities[i]<<" is located at point ("<<charges_locations[i]<<")"<<std::endl;
    }


  // we now try to create a multipole expansion centered in P and use it
  // to use it compute a local expansion centered in Q'(2,2,2). The latter will be
  // used to approximate the potential in Q due to the charges in P_i

  // we define the center, order of truncation and associated laplace functions
  // needed to construct the MultipoleExpansion class object
  Point<3> center(3.0,3.0,3.0);
  AssLegFunction *alf_ptr = new AssLegFunction();
  unsigned int truncation_order = 6;
  MultipoleExpansion multipole(truncation_order, center, alf_ptr);

  // the MultipoleExpansion created so far is ready to work, but it's empty
  // we use the Add(const double strength, const dealii::Point<3> &point)
  // member so that the expansion can account of the effect of each charge
  // considered
  for (unsigned int i=0; i<charges_locations.size(); ++i)
    multipole.Add(charges_intensities[i],charges_locations[i]);


  // we now create a local expansion centered at Q'(-2,-2,-2)
  Point<3> new_center(-2.0,-2.0,-2.0);
  LocalExpansion local(truncation_order, new_center, alf_ptr);

  // we fill the empty local exapnsion with the contribution of the multipole
  // expansion previously created, using the Add(const MultipoleExpansion &multipole)
  // member function
  local.Add(multipole);

  // to test a multipole expansion copy, we create a new multipole and then copy
  // the original multipole content into the new expension. we test
  // both the assignement through = to an uninitialized expansion, and that
  // on an initialized one
  LocalExpansion new_local(truncation_order, center, alf_ptr);
  new_local = local;
  LocalExpansion new_new_local = new_local;

  // we get the values of the expensions
  std::complex <double> *old_local_values = local.GetCoeffs();
  std::complex <double> *new_local_values = new_new_local.GetCoeffs();

  deallog<<"Original and copied local expansions are created. "<<std::endl;
  std::complex <double> total_difference(0.0,0.0);
  for (unsigned int i=0; i<(truncation_order+1)*(truncation_order+2)/2; ++i)
    {
      total_difference += abs(old_local_values[i]-new_local_values[i]);
    }

  deallog<<std::endl;
  deallog<<"Total difference between original and copied local expansion entries: "<<total_difference<<std::endl;


}
