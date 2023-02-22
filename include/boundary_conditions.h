// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the pi-BEM authors.
//
// This file is part of the pi-BEM library.
//
// The pi-BEM is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License version 2.1 as published by the Free Software Foundation.
// The full text of the license can be found in the file LICENSE at
// the top level of the pi-BEM distribution.
//
// Authors: Nicola Giuliani, Andrea Mola, Luca Heltai

#ifndef boundary_conditions_h
#define boundary_conditions_h
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <string>

#include "../include/bem_problem.h"
#include "../include/computational_domain.h"

using namespace dealii;

enum BoundaryConditionType
{
  dirichlet,
  neumann,
  robin,
  invalid
};

/**
 * - BoundaryCondition. The class handles the boundary conditions. In particular
 *   - it reads the boundary conditions for the potential and its normal
 * derivative;
 *   - given the peculiarities of the BEM, the boundary conditions represent the
 * actual unknowns, thus it creates the vectors containing the variables and
 * fills them with the proper data;
 *   - it performs the error analysis on both unknowns.
 */
template <int dim>
class BoundaryConditions : public ParameterAcceptor
{
public:
  static constexpr unsigned int MAX_COMPONENTS      = 14;
  static constexpr unsigned int MAX_CONDITION_SLOTS = 2;

  BoundaryConditions(ComputationalDomain<dim> &comp_dom,
                     BEMProblem<dim> &         bem,
                     const MPI_Comm            comm         = MPI_COMM_WORLD,
                     unsigned int              n_components = 1)
    : n_components(n_components)
    , current_component(0)
    , winds(n_components * MAX_CONDITION_SLOTS)
    , potentials(n_components * MAX_CONDITION_SLOTS)
    , robin_coeffs(n_components * MAX_CONDITION_SLOTS)
    , comp_dom(comp_dom)
    , bem(bem)
    , phis(n_components)
    , dphi_dns(n_components)
    , mpi_communicator(comm)
    , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
    , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
    , pcout(std::cout, (this_mpi_process == 0))
  {
    dofs_number = 0, output_frequency = 1;
  }

  typedef typename DoFHandler<dim - 1, dim>::active_cell_iterator cell_it;

  virtual void
  declare_parameters(ParameterHandler &prm) override;

  virtual void
  parse_parameters(ParameterHandler &prm) override;

  void
  prepare_bem_vectors(TrilinosWrappers::MPI::Vector &rhs);

  void
  prepare_robin_datastructs(
    TrilinosWrappers::MPI::Vector &robin_matrix_diagonal,
    TrilinosWrappers::MPI::Vector &robin_rhs);

  void
  prepare_robin_datastructs(
    TrilinosWrappers::MPI::Vector &robin_matrix_diagonal,
    TrilinosWrappers::MPI::Vector &robin_matrix_diagonal_imag,
    TrilinosWrappers::MPI::Vector &robin_rhs,
    TrilinosWrappers::MPI::Vector &robin_rhs_imag);

  void
  solve_problem(bool reset_matrix = true);

  // solves the complex problem where the real parts of b.conditions and Robin
  // coefficients are given by the current component data structs, the imaginary
  // parts from the next component
  void
  solve_complex_problem(bool reset_matrix = true);

  void
  output_results(const std::string);

  void
  compute_errors(bool complex, bool current_is_real);

  // the components are activated one at a time
  unsigned int
  n_phi_components() const
  {
    return n_components;
  }

  void
  set_n_phi_components(unsigned int n_components)
  {
    AssertIndexRange(this->n_components, n_components + 1);
    if (this->n_components != n_components)
      {
        winds.resize(n_components * MAX_CONDITION_SLOTS);
        potentials.resize(n_components * MAX_CONDITION_SLOTS);
        robin_coeffs.resize(n_components * MAX_CONDITION_SLOTS);
        phis.resize(n_components);
        dphi_dns.resize(n_components);

        this->n_components = n_components;
      }
  }

  unsigned int
  get_current_phi_component() const
  {
    return current_component;
  }

  void
  set_current_phi_component(unsigned int component)
  {
    AssertIndexRange(component, n_components);
    current_component = component;
    bem.set_current_phi_component(component);
  }

  // wrap extraction of current component; public interface retrieves const&
  // objects
  const Functions::ParsedFunction<dim> &
  get_wind(unsigned int component, unsigned int slot) const
  {
    AssertIndexRange(component * MAX_CONDITION_SLOTS + slot,
                     n_components * MAX_CONDITION_SLOTS);
    return *winds[component * MAX_CONDITION_SLOTS + slot];
  }

  const Functions::ParsedFunction<dim> &
  get_potential(unsigned int component, unsigned int slot) const
  {
    AssertIndexRange(component * MAX_CONDITION_SLOTS + slot,
                     n_components * MAX_CONDITION_SLOTS);
    return *potentials[component * MAX_CONDITION_SLOTS + slot];
  }

  const Functions::ParsedFunction<dim> &
  get_robin_coeffs(unsigned int component, unsigned int slot) const
  {
    AssertIndexRange(component * MAX_CONDITION_SLOTS + slot,
                     n_components * MAX_CONDITION_SLOTS);
    return *robin_coeffs[component * MAX_CONDITION_SLOTS + slot];
  }

  // same as above, for the Vectors; however, retrieve non-const&
  TrilinosWrappers::MPI::Vector &
  get_phi()
  {
    return phis[current_component];
  }

  TrilinosWrappers::MPI::Vector &
  get_dphi_dn()
  {
    return dphi_dns[current_component];
  }

  TrilinosWrappers::MPI::Vector &
  get_phi(unsigned int component)
  {
    AssertIndexRange(component, n_components);
    return phis[component];
  }

  TrilinosWrappers::MPI::Vector &
  get_dphi_dn(unsigned int component)
  {
    AssertIndexRange(component, n_components);
    return dphi_dns[component];
  }

  /// get all components as blocks
  TrilinosWrappers::MPI::BlockVector
  get_phi_components(const std::vector<unsigned int> &component_idxs)
  {
    TrilinosWrappers::MPI::BlockVector ret;
    std::vector<IndexSet>              indexSets(n_components);
    for (auto comp : component_idxs)
      {
        for (auto i : this_cpu_set)
          {
            indexSets[comp].add_index(comp * this_cpu_set.size() + i);
          }
        indexSets[comp].compress();
      }

    ret.reinit(indexSets, mpi_communicator);

    for (auto comp : component_idxs)
      {
        for (auto i : this_cpu_set)
          {
            ret[comp * this_cpu_set.size() + i] = get_phi(comp)[i];
          }
      }

    ret.compress(VectorOperation::add);

    return ret;
  }

  TrilinosWrappers::MPI::BlockVector
  get_phi_components()
  {
    std::vector<unsigned int> component_idxs(n_components);
    for (unsigned int i = 0; i < n_components; ++i)
      {
        component_idxs[i] = i;
      }

    return get_phi_components(component_idxs);
  }

  std::string output_file_name;

protected:
  // non const extractors: needed for set_time
  Functions::ParsedFunction<dim> &
  get_wind(unsigned int component, unsigned int slot)
  {
    AssertIndexRange(component * MAX_CONDITION_SLOTS + slot,
                     n_components * MAX_CONDITION_SLOTS);
    return *winds[component * MAX_CONDITION_SLOTS + slot];
  }

  Functions::ParsedFunction<dim> &
  get_potential(unsigned int component, unsigned int slot)
  {
    AssertIndexRange(component * MAX_CONDITION_SLOTS + slot,
                     n_components * MAX_CONDITION_SLOTS);
    return *potentials[component * MAX_CONDITION_SLOTS + slot];
  }

  Functions::ParsedFunction<dim> &
  get_robin_coeffs(unsigned int component, unsigned int slot)
  {
    AssertIndexRange(component * MAX_CONDITION_SLOTS + slot,
                     n_components * MAX_CONDITION_SLOTS);
    return *robin_coeffs[component * MAX_CONDITION_SLOTS + slot];
  }

  unsigned int n_components;
  unsigned int current_component;
  // each vector holds MAX_COMPONENTS * MAX_CONDITION_SLOTS functions
  std::vector<std::unique_ptr<Functions::ParsedFunction<dim>>> winds;
  std::vector<std::unique_ptr<Functions::ParsedFunction<dim>>> potentials;
  std::vector<std::unique_ptr<Functions::ParsedFunction<dim>>> robin_coeffs;

  std::string node_displacement_type;

  SolverControl solver_control;

  ComputationalDomain<dim> &comp_dom;

  BEMProblem<dim> &bem;

  types::global_dof_index dofs_number;

  unsigned int output_frequency;

  TrilinosWrappers::MPI::Vector              tmp_rhs;
  std::vector<TrilinosWrappers::MPI::Vector> phis;
  std::vector<TrilinosWrappers::MPI::Vector> dphi_dns;

  MPI_Comm mpi_communicator;

  unsigned int n_mpi_processes;

  unsigned int this_mpi_process;

  bool can_determine_phi;

  IndexSet this_cpu_set;

  ConditionalOStream pcout;
};

#endif
