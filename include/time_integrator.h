#ifndef __thiwi__time_integrator_h
#define __thiwi__time_integrator_h

#include<deal.II/base/logstream.h>
#include<deal.II/base/exceptions.h>
#include<deal.II/base/parameter_handler.h>

#include "ode_argument.h"

using namespace std;
using namespace dealii;

class TimeIntegrator
{
public:
  /** Constructor for the TimeIntegrator class. The Solver class is
   * required to have a Solver.solve(Vector<double> &dst, const
   * Vector<double> &src) method that will be called by the time
   * integrator to find out about the solution to a given src. */
  TimeIntegrator(OdeArgument &solver);

  /** Declare parameters for this class to function properly. */
  static void declare_parameters(ParameterHandler &prm);

  /** Parse a parameter handler. */
  void parse_parameters(ParameterHandler &prm);

  /** Returns a list of implemented stepping schemes, suitable for
   * use with parameter handler class. */
  static string get_stepping_names();

  /** Reinit the solver with different types or values (for example,
   * after reading a paramter file). */
  void initialize();

  /** Evolve. */
  void start_ode(const double initial_data[]);

  /** To get the output frequency from other classes using the time integrator. */
  //inline unsigned int Output_frequency()
  //           {return output_frequency;}

  unsigned int output_frequency;

private:

  /** The bubble membrane poperties. */
  OdeArgument &solver;

  /** Type of time integration. RK, EI, EE. */
  const gsl_odeiv_step_type *type;

  /** Time integrator step control*/
  gsl_odeiv_step *step;

  /** Time integrator control. */
  gsl_odeiv_control *control;

  /** Evolution object.*/
  gsl_odeiv_evolve *evolution;

  /** Ode system. */
  gsl_odeiv_system system;

  /** Size of the ode system. */
  unsigned int n_dofs;

  /** Initial step size. */
  double initial_step_size;

  /** Absolute error tolerance for adaptive time stepping. */
  double abs_tol;

  /** Relative error tolerance for adaptive time stepping. */
  double rel_tol;

  /** Initial time. */
  double initial_time;

  /** Final time. */
  double final_time;

  /** Maximum number of time steps. */
  unsigned int max_n_steps;

  // /** Number of output times per second written. */
  //unsigned int output_frequency;

  /** Initialization flag.*/
  bool is_initialized;
};


#endif
