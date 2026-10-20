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

/** Boundary type codes used by ComputationalDomain and boundary-data
 * preparation. */
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
  /// Maximum number of real storage components declared in parameter files.
  static constexpr unsigned int MAX_COMPONENTS = 14;
  /// Number of alternative boundary functions per component (slots 0 and 1).
  static constexpr unsigned int MAX_CONDITION_SLOTS = 2;

  /// Bind boundary data to comp_dom and bem on comm. n_components counts
  /// real storage components; a complex potential uses two consecutive slots
  /// in the component array. Parameter functions are initialized during
  /// parsing.
  BoundaryConditions(ComputationalDomain<dim> &comp_dom,
                     BEMProblem<dim>          &bem,
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

  /// Populate known boundary values and rhs for the active component using
  /// the domain boundary types and function-slot map. Requires initialized
  /// DoFs, normals, parsed functions and distributed solution vectors.
  void
  prepare_bem_vectors(TrilinosWrappers::MPI::Vector &rhs);

  /// Evaluate a*phi+b*dphi_dn=g on Robin nodes and store a/b and g/b.
  /// The complex overload reads real coefficients from the active component
  /// and imaginary coefficients from the next one. b must be nonzero. Output
  /// vectors must be initialized with the problem partitioning.
  void
  prepare_robin_datastructs(
    TrilinosWrappers::MPI::Vector &robin_matrix_diagonal,
    TrilinosWrappers::MPI::Vector &robin_rhs);

  /// Evaluate a*phi+b*dphi_dn=g on Robin nodes and store a/b and g/b.
  /// The complex overload reads real coefficients from the active component
  /// and imaginary coefficients from the next one. b must be nonzero. Output
  /// vectors must be initialized with the problem partitioning.
  void
  prepare_robin_datastructs(
    TrilinosWrappers::MPI::Vector &robin_matrix_diagonal,
    TrilinosWrappers::MPI::Vector &robin_matrix_diagonal_imag,
    TrilinosWrappers::MPI::Vector &robin_rhs,
    TrilinosWrappers::MPI::Vector &robin_rhs_imag);

  /// Prepare and solve the active real component collectively. Reuse integral
  /// operators when reset_matrix=false and geometry/discretization are
  /// unchanged.
  void
  solve_problem(bool reset_matrix = true);

  // solves the complex problem where the real parts of b.conditions and Robin
  // coefficients are given by the current component data structs, the imaginary
  // parts from the next component
  /// Prepare and solve the active real component and its following imaginary
  /// component collectively; current_component+1 must exist. reset_matrix has
  /// the same operator-reuse contract as solve_problem().
  void
  solve_complex_problem(bool reset_matrix = true);

  /// Write the active component and its gradient fields using the supplied
  /// filename stem. Invoke on all ranks after the solution is available.
  void
  output_results(const std::string);

  /// Compute errors and gradients for the active component. For a real problem
  /// pass (false,true). For a complex pair pass (true,true) on the real part
  /// and (true,false) on the imaginary part; the neighboring component must
  /// exist.
  void
  compute_errors(bool complex, bool current_is_real);

  // the components are activated one at a time
  /// Return the number of real storage components, counting imaginary parts.
  unsigned int
  n_phi_components() const
  {
    return n_components;
  }

  /// Grow solution and function storage; shrinking is unsupported. Parse
  /// parameters and initialize distributed vectors before using added
  /// components.
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

  /// Return the zero-based active real storage component.
  unsigned int
  get_current_phi_component() const
  {
    return current_component;
  }

  /// Select a component smaller than n_phi_components() and synchronize the
  /// active component of the bound BEMProblem.
  void
  set_current_phi_component(unsigned int component)
  {
    AssertIndexRange(component, n_components);
    current_component = component;
    bem.set_current_phi_component(component);
  }

  // wrap extraction of current component; public interface retrieves const&
  // objects
  /// Return the parsed vector field for component and function slot (0 or 1).
  /// Its normal projection supplies Neumann data. Requires parsed functions.
  const Functions::ParsedFunction<dim> &
  get_wind(unsigned int component, unsigned int slot) const
  {
    AssertIndexRange(component * MAX_CONDITION_SLOTS + slot,
                     n_components * MAX_CONDITION_SLOTS);
    return *winds[component * MAX_CONDITION_SLOTS + slot];
  }

  /// Return the parsed scalar potential for component and function slot
  /// (0 or 1). Requires parsed functions.
  const Functions::ParsedFunction<dim> &
  get_potential(unsigned int component, unsigned int slot) const
  {
    AssertIndexRange(component * MAX_CONDITION_SLOTS + slot,
                     n_components * MAX_CONDITION_SLOTS);
    return *potentials[component * MAX_CONDITION_SLOTS + slot];
  }

  /// Return the parsed Robin coefficient function for component and slot
  /// (0 or 1). In 3D its entries are a, b, g in a*phi+b*dphi_dn=g.
  const Functions::ParsedFunction<dim> &
  get_robin_coeffs(unsigned int component, unsigned int slot) const
  {
    AssertIndexRange(component * MAX_CONDITION_SLOTS + slot,
                     n_components * MAX_CONDITION_SLOTS);
    return *robin_coeffs[component * MAX_CONDITION_SLOTS + slot];
  }

  // same as above, for the Vectors; however, retrieve non-const&
  /// Return mutable potential storage for the active or explicit component.
  /// Requires a solved/initialized component; resizing invalidates the
  /// reference.
  TrilinosWrappers::MPI::Vector &
  get_phi()
  {
    return phis[current_component];
  }

  /// Return mutable normal-derivative storage for the active or explicit
  /// component. Resizing component storage invalidates the reference.
  TrilinosWrappers::MPI::Vector &
  get_dphi_dn()
  {
    return dphi_dns[current_component];
  }

  /// Return mutable potential storage for the active or explicit component.
  /// Requires a solved/initialized component; resizing invalidates the
  /// reference.
  TrilinosWrappers::MPI::Vector &
  get_phi(unsigned int component)
  {
    AssertIndexRange(component, n_components);
    return phis[component];
  }

  /// Return mutable normal-derivative storage for the active or explicit
  /// component. Resizing component storage invalidates the reference.
  TrilinosWrappers::MPI::Vector &
  get_dphi_dn(unsigned int component)
  {
    AssertIndexRange(component, n_components);
    return dphi_dns[component];
  }

  /// get all components as blocks
  /// Build and return a distributed block-vector copy using component indices
  /// as block indices; the no-argument overload selects every component.
  /// Invoke collectively after solution vectors and ownership are initialized.
  /// The result is owned by the caller and is not a view into solution storage.
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

  /// Build and return a distributed block-vector copy using component indices
  /// as block indices; the no-argument overload selects every component.
  /// Invoke collectively after solution vectors and ownership are initialized.
  /// The result is owned by the caller and is not a view into solution storage.
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
  /// Return the parsed vector field for component and function slot (0 or 1).
  /// Its normal projection supplies Neumann data. Requires parsed functions.
  Functions::ParsedFunction<dim> &
  get_wind(unsigned int component, unsigned int slot)
  {
    AssertIndexRange(component * MAX_CONDITION_SLOTS + slot,
                     n_components * MAX_CONDITION_SLOTS);
    return *winds[component * MAX_CONDITION_SLOTS + slot];
  }

  /// Return the parsed scalar potential for component and function slot
  /// (0 or 1). Requires parsed functions.
  Functions::ParsedFunction<dim> &
  get_potential(unsigned int component, unsigned int slot)
  {
    AssertIndexRange(component * MAX_CONDITION_SLOTS + slot,
                     n_components * MAX_CONDITION_SLOTS);
    return *potentials[component * MAX_CONDITION_SLOTS + slot];
  }

  /// Return the parsed Robin coefficient function for component and slot
  /// (0 or 1). In 3D its entries are a, b, g in a*phi+b*dphi_dn=g.
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
