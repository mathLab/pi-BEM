// An owned matrix row needs boundary data for columns on other MPI ranks.
// Exercise real and complex preconditioners against an explicit small matrix.
#include "bem_problem.h"
#include "tests.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  const auto rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (rank == 0)
    initlog();
  AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 2,
              ExcMessage("This regression requires two MPI ranks"));
  ComputationalDomain<3> domain;
  GridGenerator::hyper_cube(domain.tria);
  BEMProblem<3> problem(domain);
  problem.fe = std::make_unique<FE_Q<2, 3>>(1);
  problem.dh.distribute_dofs(*problem.fe);
  AssertThrow(problem.dh.n_dofs() == 4, ExcInternalError());
  problem.this_cpu_set = IndexSet(4);
  problem.this_cpu_set.add_range(2 * rank, 2 * rank + 2);
  problem.this_cpu_set.compress();
  problem.preconditioner_type                   = "ILU";
  problem.preconditioner_band                   = 8;
  problem.is_preconditioner_initialized         = false;
  problem.is_preconditioner_complex_initialized = false;
  problem.constraints.close();
  problem.constraints_imag.close();
  problem.full_sparsity_pattern.reinit(problem.this_cpu_set, MPI_COMM_WORLD);
  problem.preconditioner_sparsity_pattern.reinit(problem.this_cpu_set,
                                                 MPI_COMM_WORLD,
                                                 4);
  IndexSet complex_owned(8);
  complex_owned.add_indices(problem.this_cpu_set);
  complex_owned.add_indices(problem.this_cpu_set, 4);
  complex_owned.compress();
  problem.preconditioner_complex_sparsity_pattern.reinit(complex_owned,
                                                         MPI_COMM_WORLD,
                                                         8);
  for (auto i : problem.this_cpu_set)
    for (unsigned int j = 0; j < 4; ++j)
      problem.full_sparsity_pattern.add(i, j);
  problem.full_sparsity_pattern.compress();
  problem.neumann_matrix.reinit(problem.full_sparsity_pattern);
  problem.dirichlet_matrix.reinit(problem.full_sparsity_pattern);
  for (auto *v : {&problem.dirichlet_nodes,
                  &problem.robin_nodes,
                  &problem.robin_matrix_diagonal,
                  &problem.robin_matrix_diagonal_imag,
                  &problem.alpha})
    v->reinit(problem.this_cpu_set, MPI_COMM_WORLD);
  for (auto i : problem.this_cpu_set)
    {
      problem.dirichlet_nodes(i)            = (i == 0 || i == 3) ? 1 : 0;
      problem.robin_nodes(i)                = i == 1 ? 1 : 0;
      problem.robin_matrix_diagonal(i)      = i == 1 ? 3 : 0;
      problem.robin_matrix_diagonal_imag(i) = 0;
      problem.alpha(i)                      = 4;
      for (unsigned int j = 0; j < 4; ++j)
        {
          problem.neumann_matrix.set(i, j, 1);
          problem.dirichlet_matrix.set(i, j, 2);
        }
    }
  for (auto *v : {&problem.dirichlet_nodes,
                  &problem.robin_nodes,
                  &problem.robin_matrix_diagonal,
                  &problem.robin_matrix_diagonal_imag,
                  &problem.alpha})
    v->compress(VectorOperation::insert);
  problem.neumann_matrix.compress(VectorOperation::insert);
  problem.dirichlet_matrix.compress(VectorOperation::insert);
  problem.assemble_preconditioner();
  problem.assemble_preconditioner_complex();
  const double expected[4][4] = {{-2, 7, 1, -2},
                                 {-2, 11, 1, -2},
                                 {-2, 7, 5, -2},
                                 {-2, 7, 1, -2}};
  for (auto i : problem.this_cpu_set)
    for (unsigned int j = 0; j < 4; ++j)
      {
        AssertThrow(std::abs(problem.band_system.el(i, j) - expected[i][j]) <
                      1e-14,
                    ExcInternalError());
        AssertThrow(std::abs(problem.band_system_complex.el(i, j) -
                             expected[i][j]) < 1e-14,
                    ExcInternalError());
        AssertThrow(std::abs(problem.band_system_complex.el(i + 4, j + 4) -
                             expected[i][j]) < 1e-14,
                    ExcInternalError());
      }
  if (rank == 0)
    deallog << "Real and complex preconditioners import remote column data"
            << std::endl;
}
