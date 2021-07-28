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

  prm.enter_subsection("Wind function 2d");
  {
    Functions::ParsedFunction<2>::declare_parameters(prm, 2);
    prm.set("Function expression", "1; 1");
  }
  prm.leave_subsection();

  prm.enter_subsection("Wind function 3d");
  {
    Functions::ParsedFunction<3>::declare_parameters(prm, 3);
    prm.set("Function expression", "1; 1; 1");
  }
  prm.leave_subsection();

  prm.enter_subsection("Potential 2d");
  {
    Functions::ParsedFunction<2>::declare_parameters(prm);
    prm.set("Function expression", "x+y");
  }
  prm.leave_subsection();

  prm.enter_subsection("Potential 3d");
  {
    Functions::ParsedFunction<3>::declare_parameters(prm);
    prm.set("Function expression", "x+y+z");
  }
  prm.leave_subsection();

  // Floor potentials (that is, the scalar phi at the bottom of the domain)
  prm.enter_subsection("Floor wind 2d");
  {
    Functions::ParsedFunction<2>::declare_parameters(prm, 2);
    prm.set("Function expression", "1; 1");
  }
  prm.leave_subsection();

  prm.enter_subsection("Floor wind 3d");
  {
    Functions::ParsedFunction<3>::declare_parameters(prm, 3);
    prm.set("Function expression", "1; 1; 1");
  }
  prm.leave_subsection();

  // Wall potentials (that is, the scalar phi at the walls)
  prm.enter_subsection("Wall wind 2d");
  {
    Functions::ParsedFunction<2>::declare_parameters(prm, 2);
    prm.set("Function expression", "1; 1");
  }
  prm.leave_subsection();

  prm.enter_subsection("Wall wind 3d");
  {
    Functions::ParsedFunction<3>::declare_parameters(prm, 3);
    prm.set("Function expression", "1; 1; 1");
  }
  prm.leave_subsection();

  prm.enter_subsection("Robin coefficients");
  {
    Functions::ParsedFunction<3>::declare_parameters(prm, 3);
    prm.set("Function expression", "1; 1; 0");
  }
  prm.leave_subsection();

  // hardcoded multiple components part - is there a way to dynamically query
  // the number of comps?
  for (unsigned int comp = 1; comp < MAX_COMPS; ++comp)
    {
      // Winds (that is, the phi gradient on the wetted)
      std::string section = std::string("Wind function ") +
                            Utilities::int_to_string(comp + 1) + " " +
                            std::string("2d");
      prm.enter_subsection(section);
      {
        Functions::ParsedFunction<2>::declare_parameters(prm, 2);
        prm.set("Function expression", "1; 1");
      }
      prm.leave_subsection();

      section = std::string("Wind function ") +
                Utilities::int_to_string(comp + 1) + " " + std::string("3d");
      prm.enter_subsection(section);
      {
        Functions::ParsedFunction<3>::declare_parameters(prm, 3);
        prm.set("Function expression", "1; 1; 1");
      }
      prm.leave_subsection();

      // Potentials (that is, the scalar phi on the domain surface)
      section = std::string("Potential ") + Utilities::int_to_string(comp + 1) +
                " " + std::string("2d");
      prm.enter_subsection(section);
      {
        Functions::ParsedFunction<2>::declare_parameters(prm);
        prm.set("Function expression", "x+y");
      }
      prm.leave_subsection();

      section = std::string("Potential ") + Utilities::int_to_string(comp + 1) +
                " " + std::string("3d");
      prm.enter_subsection(section);
      {
        Functions::ParsedFunction<3>::declare_parameters(prm);
        prm.set("Function expression", "x+y+z");
      }
      prm.leave_subsection();

      // Floor potentials (that is, the scalar phi at the bottom of the domain)
      section = std::string("Floor wind ") +
                Utilities::int_to_string(comp + 1) + " " + std::string("2d");
      prm.enter_subsection(section);
      {
        Functions::ParsedFunction<2>::declare_parameters(prm, 2);
        prm.set("Function expression", "1; 1");
      }
      prm.leave_subsection();

      section = std::string("Floor wind ") +
                Utilities::int_to_string(comp + 1) + " " + std::string("3d");
      prm.enter_subsection(section);
      {
        Functions::ParsedFunction<3>::declare_parameters(prm, 3);
        prm.set("Function expression", "1; 1; 1");
      }
      prm.leave_subsection();

      // Wall potentials (that is, the scalar phi at the walls)
      section = std::string("Wall wind ") + Utilities::int_to_string(comp + 1) +
                " " + std::string("2d");
      prm.enter_subsection(section);
      {
        Functions::ParsedFunction<2>::declare_parameters(prm, 2);
        prm.set("Function expression", "1; 1");
      }
      prm.leave_subsection();

      section = std::string("Wall wind ") + Utilities::int_to_string(comp + 1) +
                " " + std::string("3d");
      prm.enter_subsection(section);
      {
        Functions::ParsedFunction<1>::declare_parameters(prm, 3);
        prm.set("Function expression", "1; 1; 1");
      }
      prm.leave_subsection();

      prm.enter_subsection("Robin coefficients " +
                           Utilities::int_to_string(comp + 1));
      {
        Functions::ParsedFunction<dim>::declare_parameters(prm, 3);
        prm.set("Function expression", "1; 1; 0");
      }
      prm.leave_subsection();
    }
}

template <int dim>
void
BoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
{
  output_file_name = prm.get("Output file name");

  prm.enter_subsection(std::string("Wind function ") +
                       Utilities::int_to_string(dim) + std::string("d"));
  {
    winds[0].reset(new Functions::ParsedFunction<dim>(dim));
    winds[0]->parse_parameters(prm);
  }
  prm.leave_subsection();

  prm.enter_subsection(std::string("Potential ") +
                       Utilities::int_to_string(dim) + std::string("d"));
  {
    potentials[0].reset(new Functions::ParsedFunction<dim>(1));
    potentials[0]->parse_parameters(prm);
  }
  prm.leave_subsection();

  prm.enter_subsection(std::string("Wall wind ") +
                       Utilities::int_to_string(dim) + std::string("d"));
  {
    wallwinds[0].reset(new Functions::ParsedFunction<dim>(dim));
    wallwinds[0]->parse_parameters(prm);
  }
  prm.leave_subsection();

  prm.enter_subsection(std::string("Floor wind ") +
                       Utilities::int_to_string(dim) + std::string("d"));
  {
    floorwinds[0].reset(new Functions::ParsedFunction<dim>(dim));
    floorwinds[0]->parse_parameters(prm);
  }
  prm.leave_subsection();

  prm.enter_subsection(std::string("Robin coefficients"));
  {
    robin_coeffs[0].reset(new Functions::ParsedFunction<dim>(dim));
    robin_coeffs[0]->parse_parameters(prm);
  }
  prm.leave_subsection();

  for (unsigned int comp = 1; comp < n_components; ++comp)
    {
      prm.enter_subsection(std::string("Wind function ") +
                           Utilities::int_to_string(comp + 1) + " " +
                           Utilities::int_to_string(dim) + std::string("d"));
      {
        winds[comp].reset(new Functions::ParsedFunction<dim>(dim));
        winds[comp]->parse_parameters(prm);
      }
      prm.leave_subsection();

      prm.enter_subsection(std::string("Potential ") +
                           Utilities::int_to_string(comp + 1) + " " +
                           Utilities::int_to_string(dim) + std::string("d"));
      {
        potentials[comp].reset(new Functions::ParsedFunction<dim>(1));
        potentials[comp]->parse_parameters(prm);
      }
      prm.leave_subsection();

      prm.enter_subsection(std::string("Wall wind ") +
                           Utilities::int_to_string(comp + 1) + " " +
                           Utilities::int_to_string(dim) + std::string("d"));
      {
        wallwinds[comp].reset(new Functions::ParsedFunction<dim>(dim));
        wallwinds[comp]->parse_parameters(prm);
      }
      prm.leave_subsection();

      prm.enter_subsection(std::string("Floor wind ") +
                           Utilities::int_to_string(comp + 1) + " " +
                           Utilities::int_to_string(dim) + std::string("d"));
      {
        floorwinds[comp].reset(new Functions::ParsedFunction<dim>(dim));
        floorwinds[comp]->parse_parameters(prm);
      }
      prm.leave_subsection();

      prm.enter_subsection(std::string("Robin coefficients ") +
                           Utilities::int_to_string(comp + 1));
      {
        robin_coeffs[comp].reset(new Functions::ParsedFunction<dim>(3));
        robin_coeffs[comp]->parse_parameters(prm);
      }
      prm.leave_subsection();
    }
}

template <int dim>
void
BoundaryConditions<dim>::solve_problem(bool reset_matrix)
{
  get_potential().set_time(0);
  get_wind().set_time(0);
  get_wallwind().set_time(0);
  get_floorwind().set_time(0);
  get_robin_coeffs().set_time(0);

  const types::global_dof_index    n_dofs = bem.dh.n_dofs();
  std::vector<types::subdomain_id> dofs_domain_association(n_dofs);
  DoFTools::get_subdomain_association(bem.dh, dofs_domain_association);
  this_cpu_set.clear();
  this_cpu_set = bem.this_cpu_set;
  this_cpu_set.compress();

  get_phi().reinit(this_cpu_set, mpi_communicator);
  get_dphi_dn().reinit(this_cpu_set, mpi_communicator);
  tmp_rhs.reinit(this_cpu_set, mpi_communicator);

  pcout << "Computing normal vector" << std::endl;
  if (reset_matrix)
    {
      bem.compute_normals();
    }
  prepare_bem_vectors(tmp_rhs);
  prepare_robin_datastructs(bem.robin_matrix_diagonal, bem.robin_rhs);

  bem.solve(get_phi(), get_dphi_dn(), tmp_rhs, reset_matrix);
  have_dirichlet_bc = bem.have_dirichlet_bc;
  if (!have_dirichlet_bc)
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
          shift =
            get_potential().value(support_points[*bem.this_cpu_set.begin()]) -
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
  // real parts - current component
  get_potential().set_time(0);
  get_wind().set_time(0);
  get_wallwind().set_time(0);
  get_floorwind().set_time(0);
  get_robin_coeffs().set_time(0);
  // imaginary parts - next component
  get_potential(current_component + 1).set_time(0);
  get_wind(current_component + 1).set_time(0);
  get_wallwind(current_component + 1).set_time(0);
  get_floorwind(current_component + 1).set_time(0);
  get_robin_coeffs(current_component + 1).set_time(0);

  const types::global_dof_index    n_dofs = bem.dh.n_dofs();
  std::vector<types::subdomain_id> dofs_domain_association(n_dofs);
  DoFTools::get_subdomain_association(bem.dh, dofs_domain_association);
  this_cpu_set.clear();
  this_cpu_set = bem.this_cpu_set;
  this_cpu_set.compress();

  // real parts - current component
  // get_phi().reinit(this_cpu_set, mpi_communicator);
  // get_dphi_dn().reinit(this_cpu_set, mpi_communicator);
  tmp_rhs.reinit(this_cpu_set, mpi_communicator);
  // imaginary parts - next component
  // get_phi(current_component + 1).reinit(this_cpu_set, mpi_communicator);
  // get_dphi_dn(current_component + 1).reinit(this_cpu_set, mpi_communicator);
  TrilinosWrappers::MPI::Vector tmp_rhs_imag;
  tmp_rhs_imag.reinit(this_cpu_set, mpi_communicator);

  if (reset_matrix)
    {
      pcout << "Computing normal vector" << std::endl;
      bem.compute_normals();
    }
  pcout << "Preparing BEM vectors - real" << std::endl;
  // TODO: these calls waste the retrieval of the support points
  // real parts - current component
  prepare_bem_vectors(tmp_rhs);
  // imaginary parts - next component
  set_current_phi_component(current_component + 1);
  pcout << "Preparing BEM vectors - imaginary" << std::endl;
  prepare_bem_vectors(tmp_rhs_imag);
  set_current_phi_component(current_component - 1);

  pcout << "Preparing Robin data structures" << std::endl;
  prepare_robin_datastructs(bem.robin_matrix_diagonal,
                            bem.robin_matrix_diagonal_imag,
                            bem.robin_rhs,
                            bem.robin_rhs_imag);

  pcout << "Solve complex problem" << std::endl;
  bem.solve(get_phi(),
            get_phi(current_component + 1),
            get_dphi_dn(),
            get_dphi_dn(current_component + 1),
            tmp_rhs,
            tmp_rhs_imag,
            reset_matrix);
  have_dirichlet_bc = bem.have_dirichlet_bc;
  if (!have_dirichlet_bc)
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
          shift =
            get_potential().value(support_points[*bem.this_cpu_set.begin()]) -
            get_phi()(*bem.this_cpu_set.begin());
        }
      MPI_Bcast(&shift, 1, MPI_DOUBLE, 0, mpi_communicator);
      vector_shift(get_phi(), shift);

      pcout << "Phi shift of real part : " << shift << std::endl;
      pcout << "Computing phi shift of imaginary part" << std::endl;
      shift = 0.0;
      if (this_mpi_process == 0)
        {
          shift = get_potential(current_component + 1)
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

  get_phi().reinit(this_cpu_set, mpi_communicator);
  get_dphi_dn().reinit(this_cpu_set, mpi_communicator);

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

  Vector<double> coeffs(3);
  for (const auto &cell : bem.dh.active_cell_iterators())
    {
      fe_v.reinit(cell);
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int j = 0; j < bem.fe->dofs_per_cell; ++j)
        {
          if (this_cpu_set.is_element(local_dof_indices[j]))
            {
              bool dirichlet =
                std::find(comp_dom.dirichlet_boundary_ids.begin(),
                          comp_dom.dirichlet_boundary_ids.end(),
                          cell->boundary_id()) !=
                comp_dom.dirichlet_boundary_ids.end();

              if (dirichlet)
                {
                  // TODO: update boundary conditions check, evaluation
                  Assert(cell->boundary_id() == static_cast<types::boundary_id>(
                                                  BoundaryType::freesurface),
                         ExcInternalError());

                  // pcout << "Dirichlet node" << std::endl;
                  rhs(local_dof_indices[j]) =
                    get_potential().value(support_points[local_dof_indices[j]]);
                  get_phi()(local_dof_indices[j]) = rhs(local_dof_indices[j]);
                }
              else
                {
                  bool neumann =
                    std::find(comp_dom.neumann_boundary_ids.begin(),
                              comp_dom.neumann_boundary_ids.end(),
                              cell->boundary_id()) !=
                    comp_dom.neumann_boundary_ids.end();

                  if (neumann)
                    {
                      // TODO: update boundary conditions check and evaluation
                      Assert(
                        cell->boundary_id() == static_cast<types::boundary_id>(
                                                 BoundaryType::floor) ||
                          cell->boundary_id() ==
                            static_cast<types::boundary_id>(
                              BoundaryType::wall) ||
                          cell->boundary_id() ==
                            static_cast<types::boundary_id>(BoundaryType::hull),
                        ExcInternalError());

                      Vector<double> imposed_pot_grad(dim);
                      switch (static_cast<BoundaryType>(cell->boundary_id()))
                        {
                          case BoundaryType::floor:
                            get_floorwind().vector_value(
                              support_points[local_dof_indices[j]],
                              imposed_pot_grad);
                            break;
                          case BoundaryType::wall:
                            get_wallwind().vector_value(
                              support_points[local_dof_indices[j]],
                              imposed_pot_grad);
                            break;
                          case BoundaryType::hull:
                            get_wind().vector_value(
                              support_points[local_dof_indices[j]],
                              imposed_pot_grad);
                            break;
                          case BoundaryType::freesurface:
                          case BoundaryType::freesurface_robin:
                          case BoundaryType::invalid:
                          default:
                            break;
                        }

                      // pcout << "Neumann node" << std::endl;
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
                      // pcout << "Robin node" << std::endl;
                      // TODO: is there a good initial value for the Robin
                      // nodes? possibly, setting tmp_rhs=coeffs(2) and phi =
                      // coeffs(2)/coeffs(0)
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

  robin_matrix_diagonal.reinit(this_cpu_set, mpi_communicator);
  robin_rhs.reinit(this_cpu_set, mpi_communicator);

  Vector<double>                    coeffs(3);
  std::set<types::global_dof_index> processed;
  for (const auto &cell : bem.dh.active_cell_iterators())
    {
      bool dirichlet =
        std::find(comp_dom.dirichlet_boundary_ids.begin(),
                  comp_dom.dirichlet_boundary_ids.end(),
                  cell->boundary_id()) != comp_dom.dirichlet_boundary_ids.end();
      bool neumann =
        std::find(comp_dom.neumann_boundary_ids.begin(),
                  comp_dom.neumann_boundary_ids.end(),
                  cell->boundary_id()) != comp_dom.neumann_boundary_ids.end();

      if (!dirichlet && !neumann)
        {
          cell->get_dof_indices(local_dof_indices);
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              if (this_cpu_set.is_element(local_dof_indices[j]) &&
                  !processed.count(local_dof_indices[j]))
                {
                  // evaluate robin coefficients
                  // coeffs(0) * phi + coeffs(1) * dphi_dn = coeffs(2)
                  get_robin_coeffs().vector_value(
                    support_points[local_dof_indices[j]], coeffs);
                  robin_matrix_diagonal(local_dof_indices[j]) =
                    coeffs(0) / coeffs(1);
                  robin_rhs(local_dof_indices[j]) = coeffs(2) / coeffs(1);

                  // pcout << "set robin datastructs @ dof "
                  //       << local_dof_indices[j] << std::endl;
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

  robin_matrix_diagonal.reinit(this_cpu_set, mpi_communicator);
  robin_rhs.reinit(this_cpu_set, mpi_communicator);
  robin_matrix_diagonal_imag.reinit(this_cpu_set, mpi_communicator);
  robin_rhs_imag.reinit(this_cpu_set, mpi_communicator);

  Vector<double>                    coeffs(3);
  Vector<double>                    coeffs_imag(3);
  std::set<types::global_dof_index> processed;
  for (const auto &cell : bem.dh.active_cell_iterators())
    {
      bool dirichlet =
        std::find(comp_dom.dirichlet_boundary_ids.begin(),
                  comp_dom.dirichlet_boundary_ids.end(),
                  cell->boundary_id()) != comp_dom.dirichlet_boundary_ids.end();
      bool neumann =
        std::find(comp_dom.neumann_boundary_ids.begin(),
                  comp_dom.neumann_boundary_ids.end(),
                  cell->boundary_id()) != comp_dom.neumann_boundary_ids.end();

      if (!dirichlet && !neumann)
        {
          cell->get_dof_indices(local_dof_indices);
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              if (this_cpu_set.is_element(local_dof_indices[j]) &&
                  !processed.count(local_dof_indices[j]))
                {
                  // evaluate robin coefficients
                  // coeffs(0) * phi + coeffs(1) * dphi_dn = coeffs(2)
                  get_robin_coeffs().vector_value(
                    support_points[local_dof_indices[j]], coeffs);
                  get_robin_coeffs(current_component + 1)
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

                  // pcout << "set robin datastructs @ dof "
                  //       << local_dof_indices[j] << " diag: " << diag
                  //       << " rhs: " << rhs << std::endl;
                  processed.insert(local_dof_indices[j]);
                }
            }
        }
    }
}

template <int dim>
void
BoundaryConditions<dim>::compute_errors()
{
  Teuchos::TimeMonitor LocalTimer(*ErrorsTime);

  // We still need to communicate our results to compute the errors.
  bem.compute_gradients(get_phi(), get_dphi_dn());
  Vector<double> localized_phi(get_phi());
  Vector<double> localized_dphi_dn(get_dphi_dn());
  Vector<double> localized_gradient_solution(
    bem.get_vector_gradients_solution()); // vector_gradients_solution
  Vector<double> localised_normals(bem.vector_normals_solution);

  // We let only the first processor do the error computations
  if (this_mpi_process == 0)
    {
      pcout << "computing errors on P0" << std::endl;

      // pcout << "comp_dom.tria.n_active_cells() "
      //       << comp_dom.tria.n_active_cells() << std::endl;
      // pcout << "bem.dh.n_dofs() " << bem.dh.n_dofs() << std::endl;
      // pcout << "bem.gradient_dh.n_dofs() " << bem.gradient_dh.n_dofs()
      //       << std::endl;

      // pcout << "localized_phi.size() " << localized_phi.size() << std::endl;
      // pcout << "localized_dphi_dn.size() " << localized_dphi_dn.size()
      //       << std::endl;
      // pcout << "localized_gradient_solution.size() "
      //       << localized_gradient_solution.size() << std::endl;
      // pcout << "localised_normals.size() " << localised_normals.size()
      //       << std::endl;

      Vector<double>      phi_diff_cell(comp_dom.tria.n_active_cells());
      Vector<double>      phi_diff_node(bem.dh.n_dofs());
      std::vector<double> phi_refval_node(bem.dh.n_dofs());

      Vector<double> dphi_dn_diff_cell(comp_dom.tria.n_active_cells());
      Vector<double> dphi_dn_diff_node(bem.dh.n_dofs());
      std::vector<Vector<double>> dphi_dn_refval_node(bem.dh.n_dofs(),
                                                      Vector<double>(dim));

      Vector<double> gradphi_diff_cell(comp_dom.tria.n_active_cells());
      Vector<double> gradphi_diff_node(bem.gradient_dh.n_dofs());
      std::vector<Vector<double>> gradphi_refval_node(bem.dh.n_dofs(),
                                                      Vector<double>(dim));

      std::vector<Point<dim>> support_points(bem.dh.n_dofs());
      DoFTools::map_dofs_to_support_points<dim - 1, dim>(*bem.mapping,
                                                         bem.dh,
                                                         support_points);

      // map material (an indirection step to boundary) to the actual expression
      std::map<types::material_id, const Function<dim, double> *>
        bcond_functions;

      // first, build the map for the potential itself
      for (const auto &pair : comp_dom.manifold2boundary_map)
        {
          switch (static_cast<BoundaryType>(pair.second))
            {
              case BoundaryType::freesurface:
                bcond_functions[pair.first] = &get_potential();
                break;
              case BoundaryType::floor:
              case BoundaryType::wall:
              case BoundaryType::hull:
              case BoundaryType::freesurface_robin:
              case BoundaryType::invalid:
              default:
                // undefined -> no difference
                break;
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
      for (const auto &pair : comp_dom.manifold2boundary_map)
        {
          switch (static_cast<BoundaryType>(pair.second))
            {
              case BoundaryType::floor:
                bcond_functions[pair.first] = &get_floorwind();
                break;
              case BoundaryType::wall:
                bcond_functions[pair.first] = &get_wallwind();
                break;
              case BoundaryType::hull:
                bcond_functions[pair.first] = &get_wind();
                break;
              case BoundaryType::freesurface:
              case BoundaryType::freesurface_robin:
              case BoundaryType::invalid:
              default:
                // undefined -> no difference
                break;
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

      IndexSet                  covered;
      std::vector<BoundaryType> boundaries = {BoundaryType::floor,
                                              BoundaryType::wall,
                                              BoundaryType::hull};
      for (auto btype : boundaries)
        {
          std::set<types::boundary_id> selected = {
            static_cast<types::boundary_id>(btype)};
          IndexSet bdofs;
          // DoFTools::extract_boundary_dofs(bem.gradient_dh,
          //                                 ComponentMask(),
          //                                 bdofs,
          //                                 selected);
          std::vector<types::global_dof_index> dofs(
            bem.gradient_fe->dofs_per_cell);
          for (const auto &cell : bem.gradient_dh.active_cell_iterators())
            {
              if (cell->is_locally_owned() &&
                  selected.count(cell->boundary_id()))
                {
                  cell->get_dof_indices(dofs);

                  bdofs.add_indices(dofs.begin(), dofs.end());
                }
            }
          bdofs.compress();

          bdofs.subtract_set(covered);
          if (!bdofs.is_empty())
            {
              switch (btype)
                {
                  case BoundaryType::floor:
                    get_floorwind().vector_value_list(support_points,
                                                      gradphi_refval_node);
                    break;
                  case BoundaryType::wall:
                    get_wallwind().vector_value_list(support_points,
                                                     gradphi_refval_node);
                    break;
                  case BoundaryType::hull:
                    get_wind().vector_value_list(support_points,
                                                 gradphi_refval_node);
                    break;
                  case BoundaryType::freesurface:
                  case BoundaryType::freesurface_robin:
                  case BoundaryType::invalid:
                  default:
                    break;
                }

              for (auto i : bdofs)
                {
                  gradphi_diff_node[bem.vec_original_to_sub_wise[i]] =
                    gradphi_refval_node
                      [bem.original_to_sub_wise[i % bem.dh.n_dofs()]]
                      [i / bem.dh.n_dofs()];
                }
              covered.add_indices(bdofs);
            }
        }

      gradphi_diff_node *= -1.0;
      gradphi_diff_node.add(1., localized_gradient_solution);

      covered.clear();
      {
        std::set<types::boundary_id> selected = {
          static_cast<types::boundary_id>(BoundaryType::freesurface)};
        IndexSet bdofs;
        // DoFTools::extract_boundary_dofs(bem.dh,
        //                                 ComponentMask(),
        //                                 bdofs,
        //                                 selected);
        std::vector<types::global_dof_index> dofs(bem.fe->dofs_per_cell);
        for (const auto &cell : bem.dh.active_cell_iterators())
          {
            if (cell->is_locally_owned() && selected.count(cell->boundary_id()))
              {
                cell->get_dof_indices(dofs);

                bdofs.add_indices(dofs.begin(), dofs.end());
              }
          }
        bdofs.compress();

        bdofs.subtract_set(covered);
        if (!bdofs.is_empty())
          {
            get_potential().value_list(support_points, phi_refval_node);

            for (auto i : bdofs)
              {
                phi_diff_node[i] = phi_refval_node.at(i);
              }
            covered.add_indices(bdofs);
          }
      }

      phi_diff_node *= -1.0;
      phi_diff_node.add(1., localized_phi);

      covered.clear();
      // TODO: calc this when doing gradphi...
      boundaries = {BoundaryType::floor,
                    BoundaryType::wall,
                    BoundaryType::hull};
      for (auto btype : boundaries)
        {
          std::set<types::boundary_id> selected = {
            static_cast<types::boundary_id>(btype)};
          IndexSet bdofs;
          // DoFTools::extract_boundary_dofs(bem.gradient_dh,
          //                                 ComponentMask(),
          //                                 bdofs,
          //                                 selected);
          std::vector<types::global_dof_index> dofs(
            bem.gradient_fe->dofs_per_cell);
          for (const auto &cell : bem.gradient_dh.active_cell_iterators())
            {
              if (cell->is_locally_owned() &&
                  selected.count(cell->boundary_id()))
                {
                  cell->get_dof_indices(dofs);

                  bdofs.add_indices(dofs.begin(), dofs.end());
                }
            }
          bdofs.compress();

          bdofs.subtract_set(covered);
          if (!bdofs.is_empty())
            {
              switch (btype)
                {
                  case BoundaryType::floor:
                    get_floorwind().vector_value_list(support_points,
                                                      dphi_dn_refval_node);
                    break;
                  case BoundaryType::wall:
                    get_wallwind().vector_value_list(support_points,
                                                     dphi_dn_refval_node);
                    break;
                  case BoundaryType::hull:
                    get_wind().vector_value_list(support_points,
                                                 dphi_dn_refval_node);
                    break;
                  case BoundaryType::freesurface:
                  case BoundaryType::freesurface_robin:
                  case BoundaryType::invalid:
                  default:
                    break;
                }

              for (auto i : bdofs)
                {
                  dphi_dn_diff_node
                    [bem.original_to_sub_wise[i % bem.dh.n_dofs()]] +=
                    localised_normals[bem.vec_original_to_sub_wise[i]] *
                    dphi_dn_refval_node
                      [bem.original_to_sub_wise[i % bem.dh.n_dofs()]]
                      [i / bem.dh.n_dofs()];
                }

              covered.add_indices(bdofs);
            }
        }

      dphi_dn_diff_node *= -1.0;
      dphi_dn_diff_node.add(1., localized_dphi_dn);

      VectorTools::integrate_difference(*bem.mapping,
                                        bem.dh,
                                        dphi_dn_diff_node,
                                        ZeroFunction<dim, double>(1),
                                        dphi_dn_diff_cell,
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
      pcout << "Phi Cells Gradient  error L_2 norm: "
            << gradphi_diff_cell.l2_norm() << std::endl;

      std::string filename_vector = "vector_error.vtu";
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

      std::string           filename_scalar = "scalar_error.vtu";
      DataOut<dim - 1, dim> dataout_scalar;
      dataout_scalar.attach_dof_handler(bem.dh);
      dataout_scalar.add_data_vector(phi_node_error,
                                     std::vector<std::string>(1, "phi_error"),
                                     DataOut<dim - 1, dim>::type_dof_data);
      dataout_scalar.add_data_vector(dphi_dn_node_error,
                                     std::vector<std::string>(1,
                                                              "dphi_dn_error"),
                                     DataOut<dim - 1, dim>::type_dof_data);
      dataout_scalar.build_patches(*bem.mapping,
                                   bem.mapping_degree,
                                   DataOut<dim - 1, dim>::curved_inner_cells);

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
