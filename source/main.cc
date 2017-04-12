

#include "driver.h"

int main (int argc, char *argv[])
{
  try
    {
      unsigned int threads;
      if (argc == 1)
        threads = numbers::invalid_unsigned_int;
      else
        threads = atoi(argv[1]);
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, threads);

      std::string pname = "parameters_bem_" + std::to_string(DIMENSION) + ".prm";
      std::string pname2 = "used_parameters_bem_" + std::to_string(DIMENSION) + ".prm";

      Driver<DIMENSION> driver;
      ParameterAcceptor::initialize(pname, pname2);
      driver.run();


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
