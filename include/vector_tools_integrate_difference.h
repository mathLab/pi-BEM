#ifndef vector_tools_integrate_difference_h
#define vector_tools_integrate_difference_h

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/utilities.h>

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
// #include<deal.II/grid/tria_boundary_lib.h>

#include <deal.II/base/types.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_integrate_difference.templates.h>

// And here are a few C++ standard header
// files that we will need:
#include <deal2lkit/parsed_data_out.h>
#include <deal2lkit/utilities.h>

using namespace dealii;

// adapted from
// https://github.com/dealii/dealii/blob/master/include/deal.II/numerics/vector_tools_integrate_difference.templates.h
template <int dim, class InVector, class OutVector, int spacedim>
void
do_integrate_difference_based_on_material_id(
  const hp::MappingCollection<dim, spacedim> &mapping,
  const DoFHandler<dim, spacedim> &           dof,
  const InVector &                            fe_function,
  // const Function<spacedim, typename InVector::value_type> &exact_solution,
  const std::map<types::material_id,
                 const Function<spacedim, typename InVector::value_type> *>
    &                                 exact_solutions,
  OutVector &                         difference,
  const dealii::hp::QCollection<dim> &q,
  const VectorTools::NormType &       norm,
  const Function<spacedim> *          weight,
  const double                        exponent_1)
{
  using Number = typename InVector::value_type;
  // we mark the "exponent" parameter to this function "const" since it is
  // strictly incoming, but we need to set it to something different later
  // on, if necessary, so have a read-write version of it:
  double exponent = exponent_1;

  const unsigned int n_components = dof.get_fe(0).n_components();

#ifdef DEBUG
  for (const auto &pair : exact_solutions)
    {
      Assert(pair.second->n_components == n_components,
             ExcDimensionMismatch(pair.second->n_components, n_components));
    }
#endif

  if (weight != nullptr)
    {
      Assert((weight->n_components == 1) ||
               (weight->n_components == n_components),
             ExcDimensionMismatch(weight->n_components, n_components));
    }

  difference.reinit(dof.get_triangulation().n_active_cells());

  switch (norm)
    {
      case VectorTools::L2_norm:
      case VectorTools::H1_seminorm:
      case VectorTools::H1_norm:
      case VectorTools::Hdiv_seminorm:
        exponent = 2.;
        break;

      case VectorTools::L1_norm:
        exponent = 1.;
        break;

      default:
        break;
    }

  UpdateFlags update_flags =
    UpdateFlags(update_quadrature_points | update_JxW_values);
  switch (norm)
    {
      case VectorTools::H1_seminorm:
      case VectorTools::Hdiv_seminorm:
      case VectorTools::W1p_seminorm:
      case VectorTools::W1infty_seminorm:
        update_flags |= UpdateFlags(update_gradients);
        if (spacedim == dim + 1)
          update_flags |= UpdateFlags(update_normal_vectors);

        break;

      case VectorTools::H1_norm:
      case VectorTools::W1p_norm:
      case VectorTools::W1infty_norm:
        update_flags |= UpdateFlags(update_gradients);
        if (spacedim == dim + 1)
          update_flags |= UpdateFlags(update_normal_vectors);
        DEAL_II_FALLTHROUGH;

      default:
        update_flags |= UpdateFlags(update_values);
        break;
    }

  const dealii::hp::FECollection<dim, spacedim> &fe_collection =
    dof.get_fe_collection();
  VectorTools::internal::IDScratchData<dim, spacedim, Number> data(
    mapping, fe_collection, q, update_flags);

  // loop over all cells
  for (const auto &cell : dof.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        // initialize for this cell
        data.x_fe_values.reinit(cell);

        const dealii::FEValues<dim, spacedim> &fe_values =
          data.x_fe_values.get_present_fe_values();
        const unsigned int n_q_points = fe_values.n_quadrature_points;
        data.resize_vectors(n_q_points, n_components);

        if (update_flags & update_values)
          fe_values.get_function_values(fe_function, data.function_values);
        if (update_flags & update_gradients)
          fe_values.get_function_gradients(fe_function, data.function_grads);

        const auto found = exact_solutions.find(cell->material_id());
        if (found != exact_solutions.end())
          {
            difference(cell->active_cell_index()) = VectorTools::internal::
              integrate_difference_inner<dim, spacedim, Number>(*found->second,
                                                                norm,
                                                                weight,
                                                                update_flags,
                                                                exponent,
                                                                n_components,
                                                                data);
            // std::cout << "using boundary-selected integration" << std::endl;
          }
        else
          {
            difference(cell->active_cell_index()) = 0;
            // std::cout
            //   << "skipping cell difference due to missing function
            //   assignment"
            //   << std::endl;
          }
      }
    else
      // the cell is a ghost cell or is artificial. write a zero into the
      // corresponding value of the returned vector
      difference(cell->active_cell_index()) = 0;
}

template <int dim, class InVector, class OutVector, int spacedim>
void
integrate_difference_based_on_material_id(
  const Mapping<dim, spacedim> &   mapping,
  const DoFHandler<dim, spacedim> &dof,
  const InVector &                 fe_function,
  // const Function<spacedim, typename InVector::value_type> &exact_solution,
  const std::map<types::material_id,
                 const Function<spacedim, typename InVector::value_type> *>
    &                          exact_solutions,
  OutVector &                  difference,
  const Quadrature<dim> &      q,
  const VectorTools::NormType &norm,
  const Function<spacedim> *   weight   = nullptr,
  const double                 exponent = 2.0)
{
  do_integrate_difference_based_on_material_id(
    hp::MappingCollection<dim, spacedim>(mapping),
    dof,
    fe_function,
    exact_solutions,
    difference,
    hp::QCollection<dim>(q),
    norm,
    weight,
    exponent);
}

#endif