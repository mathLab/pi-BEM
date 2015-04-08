#ifndef __thiwi__ode_argument_h
#define __thiwi__ode_argument_h

/** Base class that needs to be inherited by any function that wants
 * to use the time integrator class. */

class OdeArgument {
public :
    /** Returns the number of degrees of freedom. Pure virtual function. */
    virtual unsigned int n_dofs() const = 0;
    
    /** Solve for the current time t, using src as source and dst as
     * destination.  This means we will be solving the following ode:
     *
     * \dot src = dst(t, src)
     *
     * Notice that any conversion from or to double pointers need to
     * be taken care of inside the calling function. 
     */
    virtual void solve(double t, double dst[],  const double src[]) = 0;
    
    /** This function is called at the end of each iteration step for
     * the ode solver. Once again, the conversion between pointers and
     * other forms of vectors need to be done inside the inheriting
     * class. */
    virtual void output_step(const double t, 
			     const unsigned int step_number,
			     const double h,
			     const double solution[]) const = 0;

    virtual ~OdeArgument();
};

#endif
