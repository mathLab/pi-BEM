// The program starts with including a bunch
// of include files that we will use in the
// various parts of the program. Most of them
// have been discussed in previous tutorials
// already:




#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/filtered_iterator.h>

#include "../include/boundary_conditions.h"

#include <deal.II/grid/filtered_iterator.h>

#include "../include/boundary_conditions.h"
#include "../include/vector_tools_integrate_difference.h"

template <int dim, class DH = DoFHandler<dim, dim + 1>>
class FilteredDataOut : public DataOut<dim, DH>
{
public:
  FilteredDataOut(const unsigned int subdomain_id)
    : subdomain_id(subdomain_id)
  {}

  virtual typename DataOut<dim, DH>::cell_iterator
  first_cell()
  {
    typename DataOut<dim, DH>::active_cell_iterator cell =
      this->dofs->begin_active();
    while ((cell != this->dofs->end()) &&
           (cell->subdomain_id() != subdomain_id))
      {
        ++cell;
      }
    return cell;
  }

  virtual typename DataOut<dim, DH>::cell_iterator
  next_cell(const typename DataOut<dim, DH>::cell_iterator &old_cell)
  {
    if (old_cell != this->dofs->end())
      {
        const IteratorFilters::SubdomainEqualTo predicate(subdomain_id);
        return ++(
          FilteredIterator<typename DataOut<dim, DH>::active_cell_iterator>(
            predicate, old_cell));
      }
    else
      {
        return old_cell;
      }
  }

private:
  const unsigned int subdomain_id;
};

#include "Teuchos_TimeMonitor.hpp"

using Teuchos::RCP;
using Teuchos::Time;
using Teuchos::TimeMonitor;

RCP<Time> PrepareTime = Teuchos::TimeMonitor::getNewTimer("PrepareBEMVectors");
RCP<Time> ErrorsTime  = Teuchos::TimeMonitor::getNewTimer("Errors");
RCP<Time> OutputTimer = Teuchos::TimeMonitor::getNewTimer("Output");
namespace
{
  template <typename VEC>
  void
  vector_shift(VEC &in_vec, double a_scalar)
  {
    for (auto i : in_vec.locally_owned_elements())
      in_vec[i] += a_scalar;
  }
} // namespace
template <int dim>
void
BoundaryConditions<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.declare_entry("Output file name", "result", Patterns::Anything());

  std::vector<std::string> bcond_names = {"Potential",
                                          "Wind function",
                                          "Robin coefficients"};
  std::vector<std::string> defaults_2d = {"x+y", "1; 1", "1; -1; 0"};
  std::vector<std::string> defaults_3d = {"x+y+z", "1; 1; 1", "0; 0; 0"};

  for (unsigned int d = 0; d < 2; ++d)
    {
      for (unsigned int comp = 0; comp < MAX_COMPONENTS; ++comp)
        {
          for (unsigned int cond = 0; cond < bcond_names.size(); ++cond)
            {
              for (unsigned int slot = 0; slot < MAX_CONDITION_SLOTS; ++slot)
                {
                  std::string label = bcond_names[cond] + " " +
                                      Utilities::int_to_string(comp) + " " +
                                      Utilities::int_to_string(slot) + " " +
                                      Utilities::int_to_string(d + 2) + "d";

                  prm.enter_subsection(label);
                  if (!cond)
                    {
                      // potentials are scalars
                      if (!d)
                        {
                          Functions::ParsedFunction<2>::declare_parameters(prm);
                          prm.set("Function expression", defaults_2d[cond]);
                        }
                      else
                        {
                          Functions::ParsedFunction<3>::declare_parameters(prm);
                          prm.set("Function expression", defaults_3d[cond]);
                        }
                    }
                  else
                    {
                      // other are vectors
                      if (!d)
                        {
                          Functions::ParsedFunction<2>::declare_parameters(prm,
                                                                           d);
                          prm.set("Function expression", defaults_2d[cond]);
                        }
                      else
                        {
                          Functions::ParsedFunction<3>::declare_parameters(prm,
                                                                           d);
                          prm.set("Function expression", defaults_3d[cond]);
                        }
                    }
                  prm.leave_subsection();
                }
            }
        }
    }
}

template <int dim>
void
BoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
{
  output_file_name = prm.get("Output file name");

  std::vector<std::string> bcond_names = {"Potential",
                                          "Wind function",
                                          "Robin coefficients"};

  for (unsigned int comp = 0; comp < n_components; ++comp)
    {
      for (unsigned int cond = 0; cond < bcond_names.size(); ++cond)
        {
          for (unsigned int slot = 0; slot < MAX_CONDITION_SLOTS; ++slot)
            {
              std::string label = bcond_names[cond] + " " +
                                  Utilities::int_to_string(comp) + " " +
                                  Utilities::int_to_string(slot) + " " +
                                  Utilities::int_to_string(dim) + "d";

              auto pos = comp * MAX_CONDITION_SLOTS + slot;

              prm.enter_subsection(label);
              switch (cond)
                {
                  case BoundaryConditionType::dirichlet:
                    potentials[pos].reset(
                      new Functions::ParsedFunction<dim>(1));
                    potentials[pos]->parse_parameters(prm);
                    break;
                  case BoundaryConditionType::neumann:
                    winds[pos].reset(new Functions::ParsedFunction<dim>(dim));
                    winds[pos]->parse_parameters(prm);
                    break;
                  case BoundaryConditionType::robin:
                    robin_coeffs[pos].reset(
                      new Functions::ParsedFunction<dim>(dim));
                    robin_coeffs[pos]->parse_parameters(prm);
                    break;
                  case BoundaryConditionType::invalid:
                  default:
                    break;
                }
              prm.leave_subsection();
            }
        }
    }
}

template <int dim>
void
BoundaryConditions<dim>::solve_problem(bool reset_matrix)
{
  for (unsigned int i = 0; i < MAX_CONDITION_SLOTS; ++i)
    {
      get_potential(current_component, i).set_time(0);
      get_wind(current_component, i).set_time(0);
      get_robin_coeffs(current_component, i).set_time(0);
    }

  const types::global_dof_index    n_dofs = bem.dh.n_dofs();
  std::vector<types::subdomain_id> dofs_domain_association(n_dofs);
  DoFTools::get_subdomain_association(bem.dh, dofs_domain_association);
  this_cpu_set.clear();
  this_cpu_set = bem.this_cpu_set;
  this_cpu_set.compress();

  get_phi().reinit(this_cpu_set, mpi_communicator);
  get_dphi_dn().reinit(this_cpu_set, mpi_communicator);
  tmp_rhs.reinit(this_cpu_set, mpi_communicator);
  bem.robin_matrix_diagonal.reinit(this_cpu_set, mpi_communicator);
  bem.robin_rhs.reinit(this_cpu_set, mpi_communicator);

  pcout << "Computing normal vector" << std::endl;
  if (reset_matrix)
    {
      bem.compute_normals();
    }
  prepare_bem_vectors(tmp_rhs);
  prepare_robin_datastructs(bem.robin_matrix_diagonal, bem.robin_rhs);

  bem.solve(get_phi(), get_dphi_dn(), tmp_rhs, reset_matrix);

  if (!bem.can_determine_phi)
    {
      pcout << "Computing phi shift" << std::endl;
      // TODO: it seems a bit wasteful to retrieve all n_dofs support pts
      std::vector<Point<dim>> support_points(n_dofs);
      DoFTools::map_dofs_to_support_points<dim - 1, dim>(*bem.mapping,
                                                         bem.dh,
                                                         support_points);
      double shift = 0.0;
      if (this_mpi_process == 0)
        {
          shift = get_potential(current_component, 0)
                    .value(support_points[*bem.this_cpu_set.begin()]) -
                  get_phi()(*bem.this_cpu_set.begin());
        }
      MPI_Bcast(&shift, 1, MPI_DOUBLE, 0, mpi_communicator);
      vector_shift(get_phi(), shift);

      pcout << "Phi shift of " << shift << std::endl;
    }
}

template <int dim>
void
BoundaryConditions<dim>::solve_complex_problem(bool reset_matrix)
{
  for (unsigned int i = 0; i < MAX_CONDITION_SLOTS; ++i)
    {
      // real parts - current component
      get_potential(current_component, i).set_time(0);
      get_wind(current_component, i).set_time(0);
      get_robin_coeffs(current_component, i).set_time(0);
      // imaginary parts - next component
      get_potential(current_component + 1, i).set_time(0);
      get_wind(current_component + 1, i).set_time(0);
      get_robin_coeffs(current_component + 1, i).set_time(0);
    }

  const types::global_dof_index    n_dofs = bem.dh.n_dofs();
  std::vector<types::subdomain_id> dofs_domain_association(n_dofs);
  DoFTools::get_subdomain_association(bem.dh, dofs_domain_association);
  this_cpu_set.clear();
  this_cpu_set = bem.this_cpu_set;
  this_cpu_set.compress();

  // real parts - current component
  get_phi().reinit(this_cpu_set, mpi_communicator);
  get_dphi_dn().reinit(this_cpu_set, mpi_communicator);
  tmp_rhs.reinit(this_cpu_set, mpi_communicator);
  bem.robin_matrix_diagonal.reinit(this_cpu_set, mpi_communicator);
  bem.robin_rhs.reinit(this_cpu_set, mpi_communicator);
  // imaginary parts - next component
  get_phi(current_component + 1).reinit(this_cpu_set, mpi_communicator);
  get_dphi_dn(current_component + 1).reinit(this_cpu_set, mpi_communicator);
  TrilinosWrappers::MPI::Vector tmp_rhs_imag;
  tmp_rhs_imag.reinit(this_cpu_set, mpi_communicator);
  bem.robin_matrix_diagonal_imag.reinit(this_cpu_set, mpi_communicator);
  bem.robin_rhs_imag.reinit(this_cpu_set, mpi_communicator);

  if (reset_matrix)
    {
      pcout << "Computing normal vector" << std::endl;
      bem.compute_normals();
    }

  // TODO: these calls waste the retrieval of the support points
  // real parts - current component
  prepare_bem_vectors(tmp_rhs);
  // imaginary parts - next component
  set_current_phi_component(current_component + 1);
  prepare_bem_vectors(tmp_rhs_imag);
  set_current_phi_component(current_component - 1);

  prepare_robin_datastructs(bem.robin_matrix_diagonal,
                            bem.robin_matrix_diagonal_imag,
                            bem.robin_rhs,
                            bem.robin_rhs_imag);

  bem.solve(get_phi(),
            get_phi(current_component + 1),
            get_dphi_dn(),
            get_dphi_dn(current_component + 1),
            tmp_rhs,
            tmp_rhs_imag,
            reset_matrix);

  if (!bem.can_determine_phi)
    {
      pcout << "Computing phi shift of real part" << std::endl;
      // TODO: it seems a bit wasteful to retrieve all n_dofs support pts
      std::vector<Point<dim>> support_points(n_dofs);
      DoFTools::map_dofs_to_support_points<dim - 1, dim>(*bem.mapping,
                                                         bem.dh,
                                                         support_points);
      double shift = 0.0;
      if (this_mpi_process == 0)
        {
          shift = get_potential(current_component, 0)
                    .value(support_points[*bem.this_cpu_set.begin()]) -
                  get_phi()(*bem.this_cpu_set.begin());
        }
      MPI_Bcast(&shift, 1, MPI_DOUBLE, 0, mpi_communicator);
      vector_shift(get_phi(), shift);

      pcout << "Phi shift of real part : " << shift << std::endl;
      pcout << "Computing phi shift of imaginary part" << std::endl;
      shift = 0.0;
      if (this_mpi_process == 0)
        {
          shift = get_potential(current_component + 1, 0)
                    .value(support_points[*bem.this_cpu_set.begin()]) -
                  get_phi(current_component + 1)(*bem.this_cpu_set.begin());
        }
      MPI_Bcast(&shift, 1, MPI_DOUBLE, 0, mpi_communicator);
      vector_shift(get_phi(current_component + 1), shift);

      pcout << "Phi shift of imaginary part : " << shift << std::endl;
    }
}

template <int dim>
void
BoundaryConditions<dim>::prepare_bem_vectors(TrilinosWrappers::MPI::Vector &rhs)
{
  Teuchos::TimeMonitor LocalTimer(*PrepareTime);

  // get_phi().reinit(this_cpu_set, mpi_communicator);
  // get_dphi_dn().reinit(this_cpu_set, mpi_communicator);
  rhs           = 0;
  get_phi()     = 0;
  get_dphi_dn() = 0;

  const types::global_dof_index n_dofs = bem.dh.n_dofs();
  std::vector<Point<dim>>       support_points(n_dofs);
  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*bem.mapping,
                                                     bem.dh,
                                                     support_points);

  std::vector<Point<dim>> vec_support_points(bem.gradient_dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*bem.mapping,
                                                     bem.gradient_dh,
                                                     vec_support_points);


  const unsigned int                   dofs_per_cell = bem.fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  FEValues<dim - 1, dim>               fe_v(*bem.mapping,
                              *bem.fe,
                              *bem.quadrature,
                              update_values | update_normal_vectors |
                                update_quadrature_points | update_JxW_values);

  // pcout << "local_dof_indices elements: " << local_dof_indices.size()
  //       << std::endl;
  // pcout << "phi elements: " << get_phi().size() << std::endl;
  // pcout << "dphi_dn elements: " << get_phi().size() << std::endl;
  // pcout << "support_points elements: " << support_points.size() << std::endl;
  // pcout << "vec_support_points elements: " << vec_support_points.size()
  //       << std::endl;

  // Vector<double> localized_dirichlet_nodes(bem.dirichlet_nodes);
  // Vector<double> localized_neumann_nodes(bem.neumann_nodes);

  Vector<double> coeffs(3);
  for (const auto &cell : bem.dh.active_cell_iterators())
    {
      fe_v.reinit(cell);
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int j = 0; j < bem.fe->dofs_per_cell; ++j)
        {
          if (this_cpu_set.is_element(local_dof_indices[j]))
            {
              auto slot =
                comp_dom.manifold2bcondition_slot_map[cell->manifold_id()];

              bool dirichlet = bem.dirichlet_nodes(local_dof_indices[j]) == 1;
              if (dirichlet)
                {
                  rhs(local_dof_indices[j]) =
                    get_potential(current_component, slot)
                      .value(support_points[local_dof_indices[j]]);
                  get_phi()(local_dof_indices[j]) = rhs(local_dof_indices[j]);
                }
              else
                {
                  bool neumann = bem.neumann_nodes(local_dof_indices[j]) == 1;
                  if (neumann)
                    {
                      Vector<double> imposed_pot_grad(dim);
                      get_wind(current_component, slot)
                        .vector_value(support_points[local_dof_indices[j]],
                                      imposed_pot_grad);

                      double tmp_dphi_dn = 0;
                      double normy       = 0;

                      for (unsigned int d = 0; d < dim; ++d)
                        {
                          types::global_dof_index dummy =
                            bem.sub_wise_to_original[local_dof_indices[j]];
                          types::global_dof_index vec_index =
                            bem.vec_original_to_sub_wise
                              [bem.gradient_dh.n_dofs() / dim * d + dummy];

                          Assert(
                            bem.vector_this_cpu_set.is_element(vec_index),
                            ExcMessage(
                              "vector cpu set and cpu set are inconsistent"));

                          tmp_dphi_dn += imposed_pot_grad[d] *
                                         bem.vector_normals_solution[vec_index];
                          normy += bem.vector_normals_solution[vec_index] *
                                   bem.vector_normals_solution[vec_index];
                        }

                      rhs(local_dof_indices[j])           = tmp_dphi_dn;
                      get_dphi_dn()(local_dof_indices[j]) = tmp_dphi_dn;
                    }
                  else
                    {
                      rhs(local_dof_indices[j])           = 0;
                      get_phi()(local_dof_indices[j])     = 0;
                      get_dphi_dn()(local_dof_indices[j]) = 0;
                    }
                }
            }
        }
    }
}

template <int dim>
void
BoundaryConditions<dim>::prepare_robin_datastructs(
  TrilinosWrappers::MPI::Vector &robin_matrix_diagonal,
  TrilinosWrappers::MPI::Vector &robin_rhs)
{
  Teuchos::TimeMonitor LocalTimer(*PrepareTime);

  const types::global_dof_index n_dofs = bem.dh.n_dofs();
  std::vector<Point<dim>>       support_points(n_dofs);
  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*bem.mapping,
                                                     bem.dh,
                                                     support_points);

  const unsigned int                   dofs_per_cell = bem.fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  Vector<double>                    coeffs(3);
  std::set<types::global_dof_index> processed;
  for (const auto &cell : bem.dh.active_cell_iterators())
    {
      if (comp_dom.manifold2bcondition_map[cell->manifold_id()] ==
          BoundaryConditionType::robin)
        {
          auto slot =
            comp_dom.manifold2bcondition_slot_map[cell->manifold_id()];
          cell->get_dof_indices(local_dof_indices);
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              if (this_cpu_set.is_element(local_dof_indices[j]) &&
                  !processed.count(local_dof_indices[j]))
                {
                  // evaluate robin coefficients
                  // coeffs(0) * phi + coeffs(1) * dphi_dn = coeffs(2)
                  get_robin_coeffs(current_component, slot)
                    .vector_value(support_points[local_dof_indices[j]], coeffs);
                  robin_matrix_diagonal(local_dof_indices[j]) =
                    coeffs(0) / coeffs(1);
                  robin_rhs(local_dof_indices[j]) = coeffs(2) / coeffs(1);

                  processed.insert(local_dof_indices[j]);
                }
            }
        }
    }
}

template <int dim>
void
BoundaryConditions<dim>::prepare_robin_datastructs(
  TrilinosWrappers::MPI::Vector &robin_matrix_diagonal,
  TrilinosWrappers::MPI::Vector &robin_matrix_diagonal_imag,
  TrilinosWrappers::MPI::Vector &robin_rhs,
  TrilinosWrappers::MPI::Vector &robin_rhs_imag)
{
  Teuchos::TimeMonitor LocalTimer(*PrepareTime);

  const types::global_dof_index n_dofs = bem.dh.n_dofs();
  std::vector<Point<dim>>       support_points(n_dofs);
  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*bem.mapping,
                                                     bem.dh,
                                                     support_points);

  const unsigned int                   dofs_per_cell = bem.fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  Vector<double>                    coeffs(3);
  Vector<double>                    coeffs_imag(3);
  std::set<types::global_dof_index> processed;
  for (const auto &cell : bem.dh.active_cell_iterators())
    {
      if (comp_dom.manifold2bcondition_map[cell->manifold_id()] ==
          BoundaryConditionType::robin)
        {
          auto slot =
            comp_dom.manifold2bcondition_slot_map[cell->manifold_id()];
          cell->get_dof_indices(local_dof_indices);
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              if (this_cpu_set.is_element(local_dof_indices[j]) &&
                  !processed.count(local_dof_indices[j]))
                {
                  // evaluate robin coefficients
                  // coeffs(0) * phi + coeffs(1) * dphi_dn = coeffs(2)
                  get_robin_coeffs(current_component, slot)
                    .vector_value(support_points[local_dof_indices[j]], coeffs);
                  get_robin_coeffs(current_component + 1, slot)
                    .vector_value(support_points[local_dof_indices[j]],
                                  coeffs_imag);
                  std::complex<double> c0(coeffs(0), coeffs_imag(0));
                  std::complex<double> c1(coeffs(1), coeffs_imag(1));
                  std::complex<double> c2(coeffs(2), coeffs_imag(2));

                  std::complex<double> diag = c0 / c1;
                  std::complex<double> rhs  = c2 / c1;

                  robin_matrix_diagonal(local_dof_indices[j]) = std::real(diag);
                  robin_matrix_diagonal_imag(local_dof_indices[j]) =
                    std::imag(diag);
                  robin_rhs(local_dof_indices[j])      = std::real(rhs);
                  robin_rhs_imag(local_dof_indices[j]) = std::imag(rhs);

                  processed.insert(local_dof_indices[j]);
                }
            }
        }
    }
}

template <int dim>
void
BoundaryConditions<dim>::compute_errors(bool complex, bool current_is_real)
{
  Teuchos::TimeMonitor LocalTimer(*ErrorsTime);

  AssertThrow(
    complex || current_is_real,
    ExcMessage(
      "Called compute errors for a real problem, but specified that the current component is imaginary"));

  // We still need to communicate our results to compute the errors.
  bem.compute_gradients(get_phi(), get_dphi_dn());
  Vector<double> localized_phi(get_phi());
  Vector<double> localized_dphi_dn(get_dphi_dn());
  Vector<double> localized_phi_other;
  Vector<double> localized_dphi_dn_other;
  if (complex)
    {
      if (current_is_real)
        {
          localized_phi_other     = get_phi(current_component + 1);
          localized_dphi_dn_other = get_dphi_dn(current_component + 1);
        }
      else
        {
          localized_phi_other     = get_phi(current_component - 1);
          localized_dphi_dn_other = get_dphi_dn(current_component - 1);
        }
    }
  Vector<double> localized_gradient_solution(
    bem.get_vector_gradients_solution());
  Vector<double> localized_surface_gradient_solution(
    bem.get_vector_surface_gradients_solution());
  Vector<double> localised_normals(bem.vector_normals_solution);

  Vector<double> localized_dirichlet_flags(bem.dirichlet_flags);
  Vector<double> localized_neumann_flags(bem.neumann_flags);

  // We let only the first processor do the error computations
  if (this_mpi_process == 0)
    {
      pcout << "computing errors on P0" << std::endl;

      Vector<double> phi_diff_cell(comp_dom.tria.n_active_cells());
      Vector<double> phi_diff_node(bem.dh.n_dofs());

      Vector<double> dphi_dn_diff_cell(comp_dom.tria.n_active_cells());
      Vector<double> dphi_dn_diff_node(bem.dh.n_dofs());

      Vector<double> gradphi_diff_cell(comp_dom.tria.n_active_cells());
      Vector<double> gradphi_diff_node(bem.gradient_dh.n_dofs());

      Vector<double> robin_diff_cell(comp_dom.tria.n_active_cells());
      Vector<double> robin_diff_node(bem.dh.n_dofs());
      Vector<double> coeffs(3);
      Vector<double> coeffs_other(3);

      std::vector<double>         phi_refval_node(bem.dh.n_dofs());
      std::vector<Vector<double>> gradphi_refval_node(bem.dh.n_dofs(),
                                                      Vector<double>(dim));

      std::vector<Point<dim>> support_points(bem.dh.n_dofs());
      DoFTools::map_dofs_to_support_points<dim - 1, dim>(*bem.mapping,
                                                         bem.dh,
                                                         support_points);

      std::vector<types::global_dof_index> dofs(bem.fe->dofs_per_cell);
      std::vector<types::global_dof_index> dofs_gradient(
        bem.gradient_fe->dofs_per_cell);

      // map material (an indirection step to boundary) to the actual expression
      std::map<types::material_id, const Function<dim, double> *>
        bcond_functions;

      // first, build the map for the potential itself
      for (const auto &pair : comp_dom.manifold2bcondition_map)
        {
          if (pair.second == BoundaryConditionType::dirichlet)
            {
              auto slot = comp_dom.manifold2bcondition_slot_map[pair.first];
              bcond_functions[pair.first] =
                &get_potential(current_component, slot);
            }
        }

      integrate_difference_based_on_material_id(*bem.mapping,
                                                bem.dh,
                                                localized_phi,
                                                bcond_functions,
                                                phi_diff_cell,
                                                QGauss<(dim - 1)>(
                                                  2 * (2 * bem.fe->degree + 1)),
                                                VectorTools::L2_norm);

      bcond_functions.clear();
      // now, build the map for the gradient
      for (const auto &pair : comp_dom.manifold2bcondition_map)
        {
          if (pair.second == BoundaryConditionType::neumann)
            {
              auto slot = comp_dom.manifold2bcondition_slot_map[pair.first];
              bcond_functions[pair.first] = &get_wind(current_component, slot);
            }
        }

      integrate_difference_based_on_material_id(*bem.mapping,
                                                bem.gradient_dh,
                                                localized_gradient_solution,
                                                bcond_functions,
                                                gradphi_diff_cell,
                                                QGauss<(dim - 1)>(
                                                  2 * (2 * bem.fe->degree + 1)),
                                                VectorTools::L2_norm);

      // collect the error on each dof evaluating the difference between the
      // actual solution and that imposed by the b.condition
      // consider that the same dof might be confronted with multiple functions;
      // sum their errors

      // process dirichlet nodes, store |phi - imposed_phi|
      std::vector<bool> covered(bem.dh.n_dofs());
      for (unsigned int slot = 0; slot < MAX_CONDITION_SLOTS; ++slot)
        {
          get_potential(current_component, slot)
            .value_list(support_points, phi_refval_node);

          std::fill(covered.begin(), covered.end(), false);
          for (const auto &cell : bem.dh.active_cell_iterators())
            {
              if (cell->is_locally_owned() &&
                  comp_dom.manifold2bcondition_map[cell->manifold_id()] ==
                    BoundaryConditionType::dirichlet &&
                  comp_dom.manifold2bcondition_slot_map[cell->manifold_id()] ==
                    slot)
                {
                  cell->get_dof_indices(dofs);
                  for (auto dof : dofs)
                    {
                      if (!covered[dof])
                        {
                          covered[dof] = true;

                          phi_diff_node[dof] +=
                            std::abs(phi_refval_node[dof] - localized_phi[dof]);
                        }
                    }
                }
            }
        }

      // process neumann nodes, store |dphi_dn - imposed_dphi_dn|
      // also, do |grad_phi - imposed_grad_phi|?
      for (unsigned int slot = 0; slot < MAX_CONDITION_SLOTS; ++slot)
        {
          get_wind(current_component, slot)
            .vector_value_list(support_points, gradphi_refval_node);

          std::fill(covered.begin(), covered.end(), false);
          for (const auto &cell : bem.dh.active_cell_iterators())
            {
              if (cell->is_locally_owned() &&
                  comp_dom.manifold2bcondition_map[cell->manifold_id()] ==
                    BoundaryConditionType::neumann &&
                  comp_dom.manifold2bcondition_slot_map[cell->manifold_id()] ==
                    slot)
                {
                  cell->get_dof_indices(dofs);
                  for (auto dof : dofs)
                    {
                      if (!covered[dof])
                        {
                          covered[dof] = true;

                          // compute dphi_dn from the vector_gradient_solution?
                          Tensor<1, dim> tmp_surface_gradphi, tmp_norm_gradphi,
                            tmp_gradphi, tmp_boundary_gradphi, tmp_norm;
                          double tmp_dphi_dn = 0;
                          for (unsigned int d = 0; d < dim; ++d)
                            {
                              tmp_boundary_gradphi[d] =
                                gradphi_refval_node[dof][d];

                              auto vec_dof = d * bem.dh.n_dofs() + dof;

                              tmp_norm[d] = localised_normals
                                [bem.vec_original_to_sub_wise[vec_dof]];

                              tmp_surface_gradphi[d] =
                                localized_surface_gradient_solution
                                  [bem.vec_original_to_sub_wise[vec_dof]];
                              tmp_norm_gradphi[d] =
                                localized_dphi_dn[dof] *
                                localised_normals
                                  [bem.vec_original_to_sub_wise[vec_dof]];
                              tmp_gradphi[d] =
                                tmp_surface_gradphi[d] + tmp_norm_gradphi[d];

                              // uses gradient_dh
                              gradphi_diff_node
                                [bem.vec_original_to_sub_wise[vec_dof]] +=
                                std::abs(
                                  gradphi_refval_node
                                    [bem.original_to_sub_wise[dof]][d] -
                                  localized_gradient_solution
                                    [bem.vec_original_to_sub_wise[vec_dof]]);

                              tmp_dphi_dn +=
                                localised_normals
                                  [bem.vec_original_to_sub_wise[vec_dof]] *
                                gradphi_refval_node
                                  [bem.original_to_sub_wise[dof]][d];
                            }

                          // if ((tmp_boundary_gradphi - tmp_gradphi).norm() >
                          //     1e-3)
                          //   {
                          //     pcout << "dof " << dof << std::endl;
                          //     pcout << "\tcoords " << support_points[dof]
                          //           << std::endl;
                          //     pcout << "\tnormal " << tmp_norm << std::endl;
                          //     pcout << "\tdphi_dn (boundary) " << tmp_dphi_dn
                          //           << std::endl;
                          //     pcout << "\tdphi_dn (solution) "
                          //           << localized_dphi_dn[dof] << std::endl;
                          //     pcout << "\tgradphi (boundary) "
                          //           << gradphi_refval_node[dof] << std::endl;
                          //     pcout << "\tgradphi (solution) " << tmp_gradphi
                          //           << std::endl;
                          //     pcout << "\tgradphi_norm " << tmp_norm_gradphi
                          //           << std::endl;
                          //     pcout << "\tgradphi_surface "
                          //           << tmp_surface_gradphi << std::endl;
                          //   }

                          dphi_dn_diff_node[bem.original_to_sub_wise[dof]] +=
                            std::abs(
                              tmp_dphi_dn -
                              localized_dphi_dn[bem.original_to_sub_wise[dof]]);
                        }
                    }
                }
            }
        }

      // process robin nodes, store |coeffs(0) * phi + coeffs(1) * dphi_dn -
      // coeffs(2)|
      for (unsigned int slot = 0; slot < MAX_CONDITION_SLOTS; ++slot)
        {
          std::fill(covered.begin(), covered.end(), false);
          for (const auto &cell : bem.dh.active_cell_iterators())
            {
              if (cell->is_locally_owned() &&
                  comp_dom.manifold2bcondition_map[cell->manifold_id()] ==
                    BoundaryConditionType::robin &&
                  comp_dom.manifold2bcondition_slot_map[cell->manifold_id()] ==
                    slot)
                {
                  cell->get_dof_indices(dofs);
                  for (auto dof : dofs)
                    {
                      if (!covered[dof])
                        {
                          get_robin_coeffs(current_component, slot)
                            .vector_value(support_points[dof], coeffs);
                          if (complex)
                            {
                              if (current_is_real)
                                {
                                  get_robin_coeffs(current_component + 1, slot)
                                    .vector_value(support_points[dof],
                                                  coeffs_other);
                                }
                              else
                                {
                                  get_robin_coeffs(current_component - 1, slot)
                                    .vector_value(support_points[dof],
                                                  coeffs_other);
                                }
                            }
                          else
                            {
                              coeffs_other = 0;
                            }
                          // assemble the complex values considering the current
                          // component's role
                          std::complex<double> c0(coeffs(0), coeffs_other(0));
                          std::complex<double> c1(coeffs(1), coeffs_other(1));
                          std::complex<double> c2(coeffs(2), coeffs_other(2));

                          double phi_other = 0, dphi_dn_other = 0;
                          if (complex)
                            {
                              phi_other     = localized_phi_other(dof);
                              dphi_dn_other = localized_dphi_dn_other(dof);
                            }
                          std::complex<double> ph(localized_phi(dof),
                                                  phi_other);
                          std::complex<double> dph_dn(localized_dphi_dn(dof),
                                                      dphi_dn_other);
                          if (!current_is_real)
                            {
                              double tmp = c0.real();
                              c0.real(c0.imag());
                              c0.imag(tmp);

                              tmp = c1.real();
                              c1.real(c1.imag());
                              c1.imag(tmp);

                              tmp = c2.real();
                              c2.real(c2.imag());
                              c2.imag(tmp);

                              tmp = ph.real();
                              ph.real(ph.imag());
                              ph.imag(tmp);

                              tmp = dph_dn.real();
                              dph_dn.real(dph_dn.imag());
                              dph_dn.imag(tmp);
                            }

                          auto residual = c0 * ph + c1 * dph_dn;
                          if (current_is_real)
                            {
                              robin_diff_node[dof] +=
                                std::abs(std::real(residual - c2));
                            }
                          else
                            {
                              robin_diff_node[dof] +=
                                std::abs(std::imag(residual - c2));
                            }
                        }
                    }
                }
            }
        }

      VectorTools::integrate_difference(*bem.mapping,
                                        bem.dh,
                                        dphi_dn_diff_node,
                                        ZeroFunction<dim, double>(1),
                                        dphi_dn_diff_cell,
                                        QGauss<(dim - 1)>(
                                          2 * (2 * bem.fe->degree + 1)),
                                        VectorTools::L2_norm);

      VectorTools::integrate_difference(*bem.mapping,
                                        bem.dh,
                                        robin_diff_node,
                                        ZeroFunction<dim, double>(1),
                                        robin_diff_cell,
                                        QGauss<(dim - 1)>(
                                          2 * (2 * bem.fe->degree + 1)),
                                        VectorTools::L2_norm);

      pcout << "   Number of active cells:       "
            << comp_dom.tria.n_active_cells() << std::endl;
      pcout << "   Number of degrees of freedom: " << bem.dh.n_dofs()
            << std::endl;

      // TODO: phi_max_error should probably use the phi_diff_node vector
      pcout << "Phi Nodes error L_inf norm: " << phi_diff_node.linfty_norm()
            << std::endl;
      pcout << "Phi Cells error L_2 norm: " << phi_diff_cell.l2_norm()
            << std::endl;

      pcout << "dPhidN Nodes error L_inf norm: "
            << dphi_dn_diff_node.linfty_norm() << std::endl;
      pcout << "dPhidN Cells error L_2 norm: " << dphi_dn_diff_cell.l2_norm()
            << std::endl;

      pcout << "Phi Nodes Gradient error L_inf norm: "
            << gradphi_diff_node.linfty_norm() << std::endl;
      pcout << "Phi Cells Gradient error L_2 norm: "
            << gradphi_diff_cell.l2_norm() << std::endl;

      pcout << "Robin Nodes Residual error L_inf norm: "
            << robin_diff_node.linfty_norm() << std::endl;
      pcout << "Robin Nodes Residual error L_2 norm: "
            << robin_diff_cell.l2_norm() << std::endl;

      std::string filename_vector =
        "error" +
        (current_component ? "_" + Utilities::int_to_string(current_component) :
                             "") +
        "_vector.vtu";
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation(
          dim, DataComponentInterpretation::component_is_part_of_vector);
      DataOut<dim - 1, dim> dataout_vector;
      dataout_vector.attach_dof_handler(bem.gradient_dh);
      dataout_vector.add_data_vector(
        gradphi_diff_node,
        std::vector<std::string>(dim, "phi_gradient_error"),
        DataOut<dim - 1, dim>::type_dof_data,
        data_component_interpretation);

      dataout_vector.build_patches(*bem.mapping,
                                   bem.mapping_degree,
                                   DataOut<dim - 1, dim>::curved_inner_cells);

      std::ofstream file_vector(filename_vector.c_str());

      dataout_vector.write_vtu(file_vector);

      std::string filename_scalar =
        "error" +
        (current_component ? "_" + Utilities::int_to_string(current_component) :
                             "") +
        "_scalar.vtu";
      DataOut<dim - 1, DoFHandler<dim - 1, dim>> dataout_scalar;
      dataout_scalar.attach_dof_handler(bem.dh);
      dataout_scalar.add_data_vector(
        phi_diff_node,
        std::vector<std::string>(1, "phi_error"),
        DataOut<dim - 1, DoFHandler<dim - 1, dim>>::type_dof_data);
      dataout_scalar.add_data_vector(
        dphi_dn_diff_node,
        std::vector<std::string>(1, "dphi_dn_error"),
        DataOut<dim - 1, DoFHandler<dim - 1, dim>>::type_dof_data);
      dataout_scalar.add_data_vector(
        robin_diff_node,
        std::vector<std::string>(1, "robin_c_error"),
        DataOut<dim - 1, DoFHandler<dim - 1, dim>>::type_dof_data);
      dataout_scalar.build_patches(
        *bem.mapping,
        bem.mapping_degree,
        DataOut<dim - 1, DoFHandler<dim - 1, dim>>::curved_inner_cells);

      std::ofstream file_scalar(filename_scalar.c_str());
      dataout_scalar.write_vtu(file_scalar);
    }
}

template <int dim>
void
BoundaryConditions<dim>::output_results(const std::string filename)
{
  Teuchos::TimeMonitor LocalTimer(*OutputTimer);

  // Even for the output we need to serialize the code and then perform the
  // output only on the first processor.
  const Vector<double> localized_phi(get_phi());
  const Vector<double> localized_dphi_dn(get_dphi_dn());
  const Vector<double> localized_alpha(bem.alpha);
  const Vector<double> localized_gradients(bem.get_vector_gradients_solution());
  const Vector<double> localized_surf_gradients(
    bem.get_vector_surface_gradients_solution());
  const Vector<double> localized_normals(bem.vector_normals_solution);

  if (this_mpi_process == 0)
    {
      std::string filename_scalar, filename_vector;
      filename_scalar = filename + "_scalar_results" + ".vtu";
      filename_vector = filename + "_vector_results" + ".vtu";

      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation(
          dim, DataComponentInterpretation::component_is_part_of_vector);

      DataOut<dim - 1, dim> dataout_scalar;
      DataOut<dim - 1, dim> dataout_vector;

      dataout_scalar.attach_dof_handler(bem.dh);
      dataout_vector.attach_dof_handler(bem.gradient_dh);



      dataout_scalar.add_data_vector(localized_phi,
                                     "phi",
                                     DataOut<dim - 1, dim>::type_dof_data);
      dataout_scalar.add_data_vector(localized_dphi_dn,
                                     "dphi_dn",
                                     DataOut<dim - 1, dim>::type_dof_data);
      dataout_scalar.add_data_vector(localized_alpha,
                                     "alpha",
                                     DataOut<dim - 1, dim>::type_dof_data);

      dataout_vector.add_data_vector(localized_gradients,
                                     std::vector<std::string>(dim,
                                                              "phi_gradient"),
                                     DataOut<dim - 1, dim>::type_dof_data,
                                     data_component_interpretation);
      dataout_vector.add_data_vector(
        localized_surf_gradients,
        std::vector<std::string>(dim, "phi_surf_gradient"),
        DataOut<dim - 1, dim>::type_dof_data,
        data_component_interpretation);
      dataout_vector.add_data_vector(
        localized_normals,
        std::vector<std::string>(dim, "normals_at_nodes"),
        DataOut<dim - 1, dim>::type_dof_data,
        data_component_interpretation);


      dataout_scalar.build_patches(*bem.mapping,
                                   bem.mapping_degree,
                                   DataOut<dim - 1, dim>::curved_inner_cells);

      std::ofstream file_scalar(filename_scalar.c_str());

      dataout_scalar.write_vtu(file_scalar);

      dataout_vector.build_patches(*bem.mapping,
                                   bem.mapping_degree,
                                   DataOut<dim - 1, dim>::curved_inner_cells);

      std::ofstream file_vector(filename_vector.c_str());

      dataout_vector.write_vtu(file_vector);
    }
}

template class BoundaryConditions<2>;
template class BoundaryConditions<3>;
