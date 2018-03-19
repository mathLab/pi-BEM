
// The program starts with including a bunch
// of include files that we will use in the
// various parts of the program. Most of them
// have been discussed in previous tutorials
// already:

#include <deal.II/grid/filtered_iterator.h>
#include "../include/boundary_conditions.h"


template<int dim, class DH=DoFHandler<dim,dim+1> >
class FilteredDataOut : public DataOut<dim, DH >
{
public:
  FilteredDataOut (const unsigned int subdomain_id)
    :
    subdomain_id (subdomain_id)
  {}
  virtual typename DataOut<dim, DH >::cell_iterator
  first_cell ()
  {
    typename DataOut<dim, DH >::active_cell_iterator
    cell = this->dofs->begin_active();
    while ((cell != this->dofs->end()) &&
           (cell->subdomain_id() != subdomain_id))
      ++cell;
    return cell;
  }
  virtual typename DataOut<dim, DH >::cell_iterator
  next_cell (const typename DataOut<dim, DH >::cell_iterator &old_cell)
  {
    if (old_cell != this->dofs->end())
      {
        const IteratorFilters::SubdomainEqualTo
        predicate(subdomain_id);
        return
          ++(FilteredIterator
             <typename DataOut<dim, DH >::active_cell_iterator>
             (predicate,old_cell));
      }
    else
      return old_cell;
  }
private:
  const unsigned int subdomain_id;
};

#include "Teuchos_TimeMonitor.hpp"

using Teuchos::Time;
using Teuchos::TimeMonitor;
using Teuchos::RCP;

RCP<Time> PrepareTime = Teuchos::TimeMonitor::getNewTimer("PrepareBEMVectors");
RCP<Time> ErrorsTime = Teuchos::TimeMonitor::getNewTimer("Errors");
RCP<Time> OutputTimer = Teuchos::TimeMonitor::getNewTimer("Output");

template <int dim>
void BoundaryConditions<dim>::declare_parameters(ParameterHandler &prm)
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

}

template <int dim>
void BoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
{

  output_file_name = prm.get("Output file name");


  prm.enter_subsection(std::string("Wind function ")+
                       Utilities::int_to_string(dim)+std::string("d"));
  {
    wind.parse_parameters(prm);
  }
  prm.leave_subsection();

  prm.enter_subsection(std::string("Potential ")+
                       Utilities::int_to_string(dim)+std::string("d"));
  {
    potential.parse_parameters(prm);
  }
  prm.leave_subsection();

}



/*template <int dim>
double* BoundaryConditions<dim>::initial_conditions() {

  initial_wave_shape.set_time(initial_time);
  initial_wave_potential.set_time(initial_time);
  wind.set_time(initial_time);

  Vector<double> instantWindValue(dim);
  Point<dim> zero(0,0,0);
  wind.vector_value(zero,instantWindValue);
  bem.pcout<<std::endl<<"Simulation time= "<<initial_time<<"   Vinf= ";
  instantWindValue.print(cout,4,false,true);
  bem.pcout<<std::endl;

  dofs_number = bem.dh.n_dofs()+bem.gradient_dh.n_dofs();

  phi.reinit(bem.dh.n_dofs());
  dphi_dn.reinit(bem.dh.n_dofs());
  tmp_rhs.reinit(bem.dh.n_dofs());

  DXDt_and_DphiDt_vector.resize(n_dofs());

  std::vector<Point<dim> > support_points(bem.dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( bem.mapping, bem.dh, support_points);
  std::vector<Point<dim> > gradient_support_points(bem.gradient_dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( bem.mapping, bem.gradient_dh, gradient_support_points);


  unsigned int j = dim-1;
  for(unsigned int i=j; i<bem.gradient_dh.n_dofs(); i=i+dim)
      {
      DXDt_and_DphiDt_vector[i] = initial_wave_shape.value(gradient_support_points[i]);
      //bem.pcout<<DXDt_and_DphiDt_vector[i]<<std::endl;
      }

   for (unsigned int i=0; i<bem.dh.n_dofs(); i++)
      {
      DXDt_and_DphiDt_vector[i+bem.gradient_dh.n_dofs()] = initial_wave_potential.value(gradient_support_points[i]);
      //bem.pcout<<DXDt_and_DphiDt_vector[i+bem.gradient_dh.n_dofs()]<<std::endl;
      }

   max_y_coor_value = 0;
   for (unsigned int i=0; i < bem.dh.n_dofs(); i++)
       {
       //for printout
       //bem.pcout<<"Node "<<i<< "["<<support_points[i]<<"] "<<std::endl;
       max_y_coor_value = std::max(max_y_coor_value,std::abs(support_points[i](1)));
       }

  std::string filename = ( output_file_name + "_" +
         Utilities::int_to_string(0) +
         ".vtk" );

  output_results(filename);

  return &DXDt_and_DphiDt_vector[0];
} //*/




template <int dim>
void BoundaryConditions<dim>::solve_problem()
{


  potential.set_time(0);
  wind.set_time(0);

  const types::global_dof_index n_dofs =  bem.dh.n_dofs();
  std::vector<types::subdomain_id> dofs_domain_association(n_dofs);
  DoFTools::get_subdomain_association   (bem.dh,dofs_domain_association);
  this_cpu_set.clear();
  this_cpu_set = bem.this_cpu_set;
  this_cpu_set.compress();



  phi.reinit(this_cpu_set,mpi_communicator);
  dphi_dn.reinit(this_cpu_set,mpi_communicator);
  tmp_rhs.reinit(this_cpu_set,mpi_communicator);
  pcout<<"Computing normal vector"<<std::endl;
  bem.compute_normals();
  prepare_bem_vectors();


  bem.solve(phi, dphi_dn, tmp_rhs);
  have_dirichlet_bc=bem.have_dirichlet_bc;
  if (!have_dirichlet_bc)
    {
      std::vector<Point<dim> > support_points(n_dofs);
      DoFTools::map_dofs_to_support_points<dim-1, dim>( *bem.mapping, bem.dh, support_points);
      double shift = 0.0;
      if (this_mpi_process == 0)
        shift = potential.value(support_points[*bem.this_cpu_set.begin()]) - phi(*bem.this_cpu_set.begin());
      MPI_Bcast(&shift,1,MPI_DOUBLE,0,mpi_communicator);
      vector_shift(phi,shift);
    }

  // bem.compute_gradients(phi, dphi_dn);



}

template <int dim>
const TrilinosWrappers::MPI::Vector &BoundaryConditions<dim>::get_phi()
{
  return phi;
}

template <int dim>
const TrilinosWrappers::MPI::Vector &BoundaryConditions<dim>::get_dphi_dn()
{
  return dphi_dn;
}



template <int dim>
void BoundaryConditions<dim>::prepare_bem_vectors()
{
  Teuchos::TimeMonitor LocalTimer(*PrepareTime);
  // bem.compute_normals();
  const types::global_dof_index n_dofs =  bem.dh.n_dofs();

  phi.reinit(this_cpu_set,mpi_communicator);
  dphi_dn.reinit(this_cpu_set,mpi_communicator);
  // tmp_rhs.reinit(this_cpu_set,mpi_communicator);


  std::vector<Point<dim> > support_points(n_dofs);
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *bem.mapping, bem.dh, support_points);

  std::vector<Point<dim> > vec_support_points(bem.gradient_dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *bem.mapping, bem.gradient_dh, vec_support_points);

  cell_it
  cell = bem.dh.begin_active(),
  endc = bem.dh.end();

  const unsigned int   dofs_per_cell   = bem.fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  FEValues<dim-1,dim> fe_v(*bem.mapping, *bem.fe, *bem.quadrature,
                           update_values |
                           update_cell_normal_vectors |
                           update_quadrature_points |
                           update_JxW_values);


  for (cell = bem.dh.begin_active(); cell != endc; ++cell)
    {
      fe_v.reinit(cell);
      cell->get_dof_indices(local_dof_indices);
      //const std::vector<Point<dim> > &node_normals = fe_v.get_cell_normal_vectors();////provv
      for (unsigned int j=0; j<bem.fe->dofs_per_cell; ++j)
        if (this_cpu_set.is_element(local_dof_indices[j]))
          {
            //bem.pcout<<cell<<" "<<cell->material_id()<<" ("<<node_normals[0]<<") "<<-normals_sys_solution(local_dof_indices[j])<<std::endl;
            bool dirichlet = false;
            bool neumann = false;
            for (auto dbound : comp_dom.dirichlet_boundary_ids)
              if (cell->material_id() == dbound)
                {
                  dirichlet = true;
                  break;
                }
            if (dirichlet)
              {
                //tmp_rhs(local_dof_indices[j]) = node_coors[j](0);
                phi(local_dof_indices[j]) = potential.value(support_points[local_dof_indices[j]]);
                tmp_rhs(local_dof_indices[j]) = potential.value(support_points[local_dof_indices[j]]);
                //bem.pcout<<"internalElse "<<local_dof_indices[j]<<" norm ("<<node_normals[j]<<")  "<<" pos ("<<node_coors[j]<<")    "<<node_coors[j](0)<<std::endl;
              }
            else
              {
                for (auto nbound : comp_dom.neumann_boundary_ids)
                  if (cell->material_id() == nbound)
                    {
                      neumann = true;
                      break;
                    }

                if (neumann)
                  {
                    //tmp_rhs(local_dof_indices[j]) = normals_sys_solution(local_dof_indices[j]);
                    //dphi_dn(local_dof_indices[j]) = normals_sys_solution(local_dof_indices[j]);
                    Vector<double> imposed_pot_grad(dim);
                    wind.vector_value(support_points[local_dof_indices[j]],imposed_pot_grad);
                    // Point<dim> imposed_potential_gradient;
                    double tmp_dphi_dn = 0;
                    double normy = 0;
                    // double tol = 1e-1;
                    for (unsigned int d=0; d<dim; ++d)
                      {
                        types::global_dof_index dummy = bem.sub_wise_to_original[local_dof_indices[j]];
                        types::global_dof_index vec_index = bem.vec_original_to_sub_wise[bem.gradient_dh.n_dofs()/dim*d+dummy];//bem.vector_start_per_process[this_mpi_process] + d*bem.this_cpu_set.n_elements() + local_dof_indices[j]-bem.start_per_process[this_mpi_process];//bem.gradient_dh.n_dofs()/dim*d+local_dof_indices[j];//bem.vector_start_per_process[this_mpi_process]+((local_dof_indices[j]-bem.start_per_process[this_mpi_process])*dim+d); //local_dof_indices[j]*dim+d;
                        // std::cout<<this_mpi_process<<" "<<support_points[local_dof_indices[j]]<<" "<<vec_support_points[vec_index]<<std::endl;
                        Assert(bem.vector_this_cpu_set.is_element(vec_index), ExcMessage("vector cpu set and cpu set are inconsistent"));
                        // Assert(support_points[local_dof_indices[j]]==vec_support_points[vec_index], ExcMessage("the support points of dh and gradient_dh are different"));
                        tmp_dphi_dn += imposed_pot_grad[d]*bem.vector_normals_solution[vec_index];
                        normy += bem.vector_normals_solution[vec_index] * bem.vector_normals_solution[vec_index];
                      }
                    // Assert(std::fabs(normy-1.)<tol, ExcMessage("you are using wrongly the normal vector"));
                    tmp_rhs(local_dof_indices[j]) = tmp_dphi_dn;
                    dphi_dn(local_dof_indices[j]) = tmp_dphi_dn;
                  }
                else
                  {
                    tmp_rhs(local_dof_indices[j]) = 0;
                    dphi_dn(local_dof_indices[j]) = 0;
                  }
              }

          }
    }



}

template <int dim>
void BoundaryConditions<dim>::compute_errors()
{
  Teuchos::TimeMonitor LocalTimer(*ErrorsTime);

  // We still need to communicate our results to compute the errors.
  bem.compute_gradients(phi,dphi_dn);
  Vector<double> localized_gradient_solution(bem.vector_gradients_solution);//vector_gradients_solution
  Vector<double> localized_phi(phi);
  Vector<double> localized_dphi_dn(dphi_dn);
  Vector<double> localised_normals(bem.vector_normals_solution);
  // Vector<double> localised_alpha(bem.alpha);
  // We let only the first processor do the error computations
  if (this_mpi_process == 0)
    {
      pcout<<"computing errors on P0"<<std::endl;
      // for(auto i : localised_alpha.locally_owned_elements())
      //   localised_alpha[i] -= 0.5;

      Vector<double> grad_difference_per_cell (comp_dom.tria.n_active_cells());
      std::vector<Point<dim> > support_points(bem.dh.n_dofs());
      double phi_max_error;// = localized_phi.linfty_norm();
      Vector<double> difference_per_cell (comp_dom.tria.n_active_cells());
      DoFTools::map_dofs_to_support_points<dim-1, dim>( *bem.mapping, bem.dh, support_points);

      // removing this because now solution is translated after BEM solution
      // is given back by bem_problem in solve_problem
      /*
            if (false)//(!have_dirichlet_bc)
              {
                std::vector<double> exact_sol(bem.dh.n_dofs());
                Vector<double> exact_sol_deal(bem.dh.n_dofs());
                potential.value_list(support_points,exact_sol);
                for (auto i : exact_sol_deal.locally_owned_elements())
                  exact_sol_deal[bem.original_to_sub_wise[i]] = exact_sol[bem.original_to_sub_wise[i]];
                auto exact_mean = VectorTools::compute_mean_value(*bem.mapping,bem.dh,QGauss<(dim-1)>(2*(2*bem.fe->degree+1)), exact_sol_deal, 0);
                exact_sol_deal.add(-exact_mean);
                auto my_mean = VectorTools::compute_mean_value(*bem.mapping,bem.dh,QGauss<(dim-1)>(2*(2*bem.fe->degree+1)), localized_phi, 0);
                localized_phi.add(-my_mean);
                localized_phi.sadd(1.,-1.,exact_sol_deal);
                std::cout<<"Extracting the mean value from phi solution"<<std::endl;
                VectorTools::integrate_difference (*bem.mapping, bem.dh, localized_phi,
                                                   ZeroFunction<dim, double> (1),
                                                   difference_per_cell,
                                                   QGauss<(dim-1)>(2*(2*bem.fe->degree+1)),
                                                   VectorTools::L2_norm);

                phi_max_error = localized_phi.linfty_norm();


              }
            else
              {
      */
      VectorTools::integrate_difference (*bem.mapping, bem.dh, localized_phi,
                                         potential,
                                         difference_per_cell,
                                         QGauss<(dim-1)>(2*(2*bem.fe->degree+1)),
                                         VectorTools::L2_norm);

      phi_max_error = difference_per_cell.linfty_norm();

//        }
      VectorTools::integrate_difference (*bem.mapping, bem.gradient_dh, localized_gradient_solution,
                                         wind,
                                         grad_difference_per_cell,
                                         QGauss<(dim-1)>(2*(2*bem.fe->degree+1)),
                                         VectorTools::L2_norm);
      const double grad_L2_error = grad_difference_per_cell.l2_norm();

      const double L2_error = difference_per_cell.l2_norm();


      Vector<double> vector_gradients_node_error(bem.gradient_dh.n_dofs());
      std::vector<Vector<double> > grads_nodes_errs(bem.dh.n_dofs(),Vector<double>(dim));
      wind.vector_value_list(support_points,grads_nodes_errs);
      for (types::global_dof_index d=0; d<dim; ++d)
        for (types::global_dof_index i=0; i<bem.dh.n_dofs(); ++i)
          vector_gradients_node_error(bem.vec_original_to_sub_wise[d*bem.dh.n_dofs()+i]) = grads_nodes_errs[bem.original_to_sub_wise[i]](d);
      vector_gradients_node_error*=-1.0;
      vector_gradients_node_error.add(1.,localized_gradient_solution);

      Vector<double> phi_node_error(bem.dh.n_dofs());
      std::vector<double> phi_nodes_errs(bem.dh.n_dofs());
      potential.value_list(support_points,phi_nodes_errs);
      for (types::global_dof_index i=0; i<bem.dh.n_dofs(); ++i)
        phi_node_error(i) = phi_nodes_errs[i];

      phi_node_error*=-1.0;
      phi_node_error.add(1.,localized_phi);

      Vector<double> dphi_dn_node_error(bem.dh.n_dofs());
      std::vector<Vector<double> > dphi_dn_nodes_errs(bem.dh.n_dofs(), Vector<double> (dim));
      wind.vector_value_list(support_points,dphi_dn_nodes_errs);
      dphi_dn_node_error = 0.;
      for (types::global_dof_index i=0; i<bem.dh.n_dofs(); ++i)
        {
          // dphi_dn_node_error[i] = 0.;
          for (unsigned int d=0; d<dim; ++d)
            {
              dphi_dn_node_error[bem.original_to_sub_wise[i]] += localised_normals[bem.vec_original_to_sub_wise[i+d*bem.dh.n_dofs()]] * dphi_dn_nodes_errs[bem.original_to_sub_wise[i]][d];


            }
        }
      dphi_dn_node_error*=-1.0;
      dphi_dn_node_error.add(1.,localized_dphi_dn);

      // dphi_dn_node_error.print(std::cout);
      Vector<double> difference_per_cell_2(comp_dom.tria.n_active_cells());
      VectorTools::integrate_difference (*bem.mapping, bem.dh, dphi_dn_node_error,
                                         ZeroFunction<dim, double> (1),
                                         difference_per_cell_2,
                                         QGauss<(dim-1)>(2*(2*bem.fe->degree+1)),
                                         VectorTools::L2_norm);
      const double dphi_dn_L2_error = difference_per_cell_2.l2_norm();

      const double grad_phi_max_error = vector_gradients_node_error.linfty_norm();
      const types::global_dof_index n_active_cells=comp_dom.tria.n_active_cells();
      const types::global_dof_index n_dofs=bem.dh.n_dofs();

      pcout << "   Number of active cells:       "
            << n_active_cells
            << std::endl
            << "   Number of degrees of freedom: "
            << n_dofs
            << std::endl
            ;

      pcout<<"Phi Nodes error L_inf norm: "<<phi_max_error<<std::endl;
      pcout<<"Phi Cells error L_2 norm: "<<L2_error<<std::endl;
      pcout<<"dPhidN Nodes error L_inf norm: "<<dphi_dn_node_error.linfty_norm()<<std::endl;
      pcout<<"dPhidN Nodes error L_2 norm: "<<dphi_dn_L2_error<<std::endl;
      pcout<<"Phi Nodes Gradient error L_inf norm: "<<grad_phi_max_error<<std::endl;
      pcout<<"Phi Cells Gradient  error L_2 norm: "<<grad_L2_error<<std::endl;
      // pcout<<"alpha Nodes error L_inf norm: "<<localised_alpha.linfty_norm()<<std::endl;
      // pcout<<"alpha Nodes error L_2 norm: "<<localised_alpha.l2_norm()<<std::endl;

      std::string filename_vector = "vector_error.vtu";
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation
      (dim, DataComponentInterpretation::component_is_part_of_vector);
      DataOut<dim-1, DoFHandler<dim-1, dim> > dataout_vector;
      dataout_vector.attach_dof_handler(bem.gradient_dh);
      dataout_vector.add_data_vector(vector_gradients_node_error, std::vector<std::string > (dim,"phi_gradient_error"), DataOut<dim-1, DoFHandler<dim-1, dim> >::type_dof_data, data_component_interpretation);

      dataout_vector.build_patches(*bem.mapping,
                                   bem.mapping_degree,
                                   DataOut<dim-1, DoFHandler<dim-1, dim> >::curved_inner_cells);

      std::ofstream file_vector(filename_vector.c_str());

      dataout_vector.write_vtu(file_vector);

      std::string filename_scalar = "scalar_error.vtu";
      DataOut<dim-1, DoFHandler<dim-1, dim> > dataout_scalar;
      dataout_scalar.attach_dof_handler(bem.dh);
      dataout_scalar.add_data_vector(phi_node_error, std::vector<std::string > (1,"phi_error"), DataOut<dim-1, DoFHandler<dim-1, dim> >::type_dof_data);
      dataout_scalar.add_data_vector(dphi_dn_node_error, std::vector<std::string > (1,"dphi_dn_error"), DataOut<dim-1, DoFHandler<dim-1, dim> >::type_dof_data);
      dataout_scalar.build_patches(*bem.mapping,
                                   bem.mapping_degree,
                                   DataOut<dim-1, DoFHandler<dim-1, dim> >::curved_inner_cells);

      std::ofstream file_scalar(filename_scalar.c_str());
      dataout_scalar.write_vtu(file_scalar);

    }
}

template <int dim>
void BoundaryConditions<dim>::output_results(const std::string filename)
{
  Teuchos::TimeMonitor LocalTimer(*OutputTimer);
  // At the time being the output is not running in parallel with saks
  // const Vector<double> localized_phi (phi);
  // const Vector<double> localized_dphi_dn (dphi_dn);
  // const Vector<double> localized_alpha (bem.alpha);
  // const Vector<double> localized_gradients (bem.vector_gradients_solution);
  // const Vector<double> localized_normals (bem.vector_normals_solution);
  //
  // if(this_mpi_process == 0)
  // {
  //   data_out_scalar.prepare_data_output(bem.dh);
  //   data_out_scalar.add_data_vector (localized_phi, "phi");
  //   data_out_scalar.add_data_vector (localized_dphi_dn, "dphidn");
  //   data_out_scalar.add_data_vector (localized_alpha, "alpha");
  //   data_out_scalar.write_data_and_clear("", bem.mapping);
  //
  //   data_out_vector.prepare_data_output(bem.gradient_dh);
  //   data_out_vector.add_data_vector (localized_gradients, "gradient");
  //   data_out_vector.add_data_vector (localized_normals, "normal");
  //   data_out_vector.write_data_and_clear("", bem.mapping);
  // }

  // Even for the output we need to serialize the code and then perform the
  // output only on the first processor.
  const Vector<double> localized_phi (phi);
  const Vector<double> localized_dphi_dn (dphi_dn);
  const Vector<double> localized_alpha (bem.alpha);
  const Vector<double> localized_gradients (bem.vector_gradients_solution);
  const Vector<double> localized_surf_gradients (bem.vector_surface_gradients_solution);
  const Vector<double> localized_normals (bem.vector_normals_solution);
  // localized_dphi_dn.print(std::cout);
  // localized_phi.print(std::cout);
  if (this_mpi_process == 0)
    {
      std::string filename_scalar, filename_vector;
      filename_scalar = filename + "_scalar_results" + ".vtu";
      filename_vector = filename + "_vector_results" + ".vtu";

      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation
      (dim, DataComponentInterpretation::component_is_part_of_vector);

      DataOut<dim-1, DoFHandler<dim-1, dim> > dataout_scalar;
      DataOut<dim-1, DoFHandler<dim-1, dim> > dataout_vector;

      dataout_scalar.attach_dof_handler(bem.dh);
      dataout_vector.attach_dof_handler(bem.gradient_dh);



      dataout_scalar.add_data_vector(localized_phi, "phi", DataOut<dim-1, DoFHandler<dim-1, dim> >::type_dof_data);
      dataout_scalar.add_data_vector(localized_dphi_dn, "dphi_dn", DataOut<dim-1, DoFHandler<dim-1, dim> >::type_dof_data);
      dataout_scalar.add_data_vector(localized_alpha, "alpha", DataOut<dim-1, DoFHandler<dim-1, dim> >::type_dof_data);

      dataout_vector.add_data_vector(localized_gradients, std::vector<std::string > (dim,"phi_gradient"), DataOut<dim-1, DoFHandler<dim-1, dim> >::type_dof_data, data_component_interpretation);
      dataout_vector.add_data_vector(localized_surf_gradients, std::vector<std::string > (dim,"phi_surf_gradient"), DataOut<dim-1, DoFHandler<dim-1, dim> >::type_dof_data, data_component_interpretation);
      dataout_vector.add_data_vector(localized_normals, std::vector<std::string > (dim,"normals_at_nodes"), DataOut<dim-1, DoFHandler<dim-1, dim> >::type_dof_data, data_component_interpretation);


      dataout_scalar.build_patches(*bem.mapping,
                                   bem.mapping_degree,
                                   DataOut<dim-1, DoFHandler<dim-1, dim> >::curved_inner_cells);

      std::ofstream file_scalar(filename_scalar.c_str());

      dataout_scalar.write_vtu(file_scalar);

      dataout_vector.build_patches(*bem.mapping,
                                   bem.mapping_degree,
                                   DataOut<dim-1, DoFHandler<dim-1, dim> >::curved_inner_cells);

      std::ofstream file_vector(filename_vector.c_str());

      dataout_vector.write_vtu(file_vector);

    }

}

template class BoundaryConditions<2>;
template class BoundaryConditions<3>;
