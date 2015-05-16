/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Luca Heltai, Cataldo Manigrasso, 2009
 */


#include "step_fma.h"





// @sect3{The main() function}

// This is the main function of this program. It is exactly like all previous
// tutorial programs:
int main (int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  try
    {
      using namespace dealii;
      using namespace MinFmm;

      const unsigned int degree = 1;

      // deallog.depth_console (3);
      // StepFMA<2> laplace_problem_2d(degree, mapping_degree);
      // laplace_problem_2d.run();

      StepFMA<3> laplace_problem_3d_fmm(degree, true);
      ParameterAcceptor::initialize("parameters_fma.prm", "used_parameters_fma.prm");
      laplace_problem_3d_fmm.run();

      StepFMA<3> laplace_problem_3d_direct(degree, false);
      ParameterAcceptor::initialize("parameters_fma.prm", "used_parameters_fma.prm");
      laplace_problem_3d_direct.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
