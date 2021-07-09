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

  prm.declare_entry("Potential components", "1", Patterns::Integer());

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
  bem.compute_normals();
  prepare_bem_vectors();

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
BoundaryConditions<dim>::prepare_bem_vectors()
{
  Teuchos::TimeMonitor          LocalTimer(*PrepareTime);
  const types::global_dof_index n_dofs = bem.dh.n_dofs();

  get_phi().reinit(this_cpu_set, mpi_communicator);
  get_dphi_dn().reinit(this_cpu_set, mpi_communicator);

  std::vector<Point<dim>> support_points(n_dofs);
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

  for (const auto &cell : bem.dh.active_cell_iterators())
    {
      fe_v.reinit(cell);
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int j = 0; j < bem.fe->dofs_per_cell; ++j)
        {
          if (this_cpu_set.is_element(local_dof_indices[j]))
            {
              bool dirichlet = false;
              bool neumann   = false;
              for (auto dbound : comp_dom.dirichlet_boundary_ids)
                {
                  if (cell->boundary_id() == dbound)
                    {
                      Assert((cell->boundary_id() == dbound) ==
                               (cell->boundary_id() ==
                                static_cast<types::boundary_id>(
                                  BoundaryType::freesurface)),
                             ExcInternalError());
                      dirichlet = true;
                      break;
                    }
                }

              if (dirichlet)
                {
                  get_phi()(local_dof_indices[j]) =
                    get_potential().value(support_points[local_dof_indices[j]]);
                  tmp_rhs(local_dof_indices[j]) =
                    get_phi()(local_dof_indices[j]);
                }
              else
                {
                  for (auto nbound : comp_dom.neumann_boundary_ids)
                    {
                      if (cell->boundary_id() == nbound)
                        {
                          Assert((cell->boundary_id() == nbound) ==
                                   (cell->boundary_id() ==
                                      static_cast<types::boundary_id>(
                                        BoundaryType::floor) ||
                                    cell->boundary_id() ==
                                      static_cast<types::boundary_id>(
                                        BoundaryType::wall) ||
                                    cell->boundary_id() ==
                                      static_cast<types::boundary_id>(
                                        BoundaryType::hull)),
                                 ExcInternalError());
                          neumann = true;
                          break;
                        }
                    }

                  if (neumann)
                    {
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
                          case BoundaryType::invalid:
                          default:
                            break;
                        }
                      // get_wind().vector_value(
                      //   support_points[local_dof_indices[j]],
                      //   imposed_pot_grad);
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

                          tmp_dphi_dn +=
                            imposed_pot_grad[d] *
                            bem.get_vector_normals_solution()[vec_index];
                          normy +=
                            bem.get_vector_normals_solution()[vec_index] *
                            bem.get_vector_normals_solution()[vec_index];
                        }

                      tmp_rhs(local_dof_indices[j])       = tmp_dphi_dn;
                      get_dphi_dn()(local_dof_indices[j]) = tmp_dphi_dn;
                    }
                  else
                    {
                      tmp_rhs(local_dof_indices[j])       = 0;
                      get_dphi_dn()(local_dof_indices[j]) = 0;
                    }
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
  Vector<double> localized_gradient_solution(
    bem.get_vector_gradients_solution()); // vector_gradients_solution
  Vector<double> localized_phi(get_phi());
  Vector<double> localized_dphi_dn(get_dphi_dn());
  Vector<double> localised_normals(bem.get_vector_normals_solution());

  // We let only the first processor do the error computations
  if (this_mpi_process == 0)
    {
      pcout << "computing errors on P0" << std::endl;

      Vector<double>          phi_diff_cell(comp_dom.tria.n_active_cells());
      Vector<double>          gradphi_diff_cell(comp_dom.tria.n_active_cells());
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
              case BoundaryType::invalid:
              default:
                // undefined -> no difference
                break;
            }
        }

      VectorTools::integrate_difference(*bem.mapping,
                                        bem.dh,
                                        localized_phi,
                                        get_potential(),
                                        phi_diff_cell,
                                        QGauss<(dim - 1)>(
                                          2 * (2 * bem.fe->degree + 1)),
                                        VectorTools::L2_norm);
      // integrate_difference_based_on_material_id(*bem.mapping,
      //                                           bem.dh,
      //                                           localized_phi,
      //                                           bcond_functions,
      //                                           phi_diff_cell,
      //                                           QGauss<(dim - 1)>(
      //                                             2 * (2 * bem.fe->degree +
      //                                             1)),
      //                                           VectorTools::L2_norm);

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
              case BoundaryType::invalid:
              default:
                // undefined -> no difference
                break;
            }
        }

      VectorTools::integrate_difference(*bem.mapping,
                                        bem.gradient_dh,
                                        localized_gradient_solution,
                                        get_wind(),
                                        gradphi_diff_cell,
                                        QGauss<(dim - 1)>(
                                          2 * (2 * bem.fe->degree + 1)),
                                        VectorTools::L2_norm);
      // integrate_difference_based_on_material_id(*bem.mapping,
      //                                           bem.gradient_dh,
      //                                           localized_gradient_solution,
      //                                           bcond_functions,
      //                                           gradphi_diff_cell,
      //                                           QGauss<(dim - 1)>(
      //                                             2 * (2 * bem.fe->degree +
      //                                             1)),
      //                                           VectorTools::L2_norm);

      Vector<double>              gradphi_diff_node(bem.gradient_dh.n_dofs());
      std::vector<Vector<double>> gradphi_refval_node(bem.dh.n_dofs(),
                                                      Vector<double>(dim));
      // TODO: change to use mappping
      get_wind().vector_value_list(support_points, gradphi_refval_node);
      for (types::global_dof_index d = 0; d < dim; ++d)
        {
          for (types::global_dof_index i = 0; i < bem.dh.n_dofs(); ++i)
            {
              gradphi_diff_node
                [bem.vec_original_to_sub_wise[d * bem.dh.n_dofs() + i]] =
                  gradphi_refval_node[bem.original_to_sub_wise[i]][d];
            }
        }
      gradphi_diff_node *= -1.0;
      gradphi_diff_node.add(1., localized_gradient_solution);

      Vector<double>      phi_diff_node(bem.dh.n_dofs());
      std::vector<double> phi_refval_node(bem.dh.n_dofs());
      get_potential().value_list(support_points, phi_refval_node);
      for (types::global_dof_index i = 0; i < bem.dh.n_dofs(); ++i)
        {
          phi_diff_node[i] = phi_refval_node[i];
        }

      phi_diff_node *= -1.0;
      phi_diff_node.add(1., localized_phi);

      Vector<double>              dphi_dn_diff_node(bem.dh.n_dofs());
      std::vector<Vector<double>> dphi_dn_refval_node(bem.dh.n_dofs(),
                                                      Vector<double>(dim));
      // TODO: change to use mappping
      get_wind().vector_value_list(support_points, dphi_dn_refval_node);
      dphi_dn_diff_node = 0.;
      for (types::global_dof_index i = 0; i < bem.dh.n_dofs(); ++i)
        {
          for (unsigned int d = 0; d < dim; ++d)
            {
              dphi_dn_diff_node[bem.original_to_sub_wise[i]] +=
                localised_normals
                  [bem.vec_original_to_sub_wise[i + d * bem.dh.n_dofs()]] *
                dphi_dn_refval_node[bem.original_to_sub_wise[i]][d];
            }
        }

      dphi_dn_diff_node *= -1.0;
      dphi_dn_diff_node.add(1., localized_dphi_dn);


      Vector<double> dphi_dn_diff_cell(comp_dom.tria.n_active_cells());
      VectorTools::integrate_difference(*bem.mapping,
                                        bem.dh,
                                        dphi_dn_diff_node,
                                        ZeroFunction<dim, double>(1),
                                        dphi_dn_diff_cell,
                                        QGauss<(dim - 1)>(
                                          2 * (2 * bem.fe->degree + 1)),
                                        VectorTools::L2_norm);

      double       phi_max_error      = phi_diff_cell.linfty_norm();
      const double L2_error           = phi_diff_cell.l2_norm();
      const double dphi_dn_L2_error   = dphi_dn_diff_cell.l2_norm();
      const double grad_phi_max_error = gradphi_diff_node.linfty_norm();
      const double grad_L2_error      = gradphi_diff_cell.l2_norm();

      pcout << "   Number of active cells:       "
            << comp_dom.tria.n_active_cells() << std::endl;
      pcout << "   Number of degrees of freedom: " << bem.dh.n_dofs()
            << std::endl;

      // TODO: phi_max_error should probably use the phi_diff_node vector
      pcout << "Phi Nodes error L_inf norm: " << phi_max_error << std::endl;
      pcout << "Phi Cells error L_2 norm: " << L2_error << std::endl;

      pcout << "dPhidN Nodes error L_inf norm: "
            << dphi_dn_diff_node.linfty_norm() << std::endl;
      // TODO: this, too, is on cells
      pcout << "dPhidN Nodes error L_2 norm: " << dphi_dn_L2_error << std::endl;

      pcout << "Phi Nodes Gradient error L_inf norm: " << grad_phi_max_error
            << std::endl;
      pcout << "Phi Cells Gradient  error L_2 norm: " << grad_L2_error
            << std::endl;

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
  const Vector<double> localized_normals(bem.get_vector_normals_solution());

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
