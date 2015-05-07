//----------------------------  step-34.cc  ---------------------------
//    $Id: step-34.cc 18734 2009-04-25 13:36:48Z heltai $
//    Version: $Name$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//    Authors: Luca Heltai, Cataldo Manigrasso
//
//----------------------------  step-34.cc  ---------------------------


// @sect3{Include files}

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

  dofs_number = comp_dom.dh.n_dofs()+comp_dom.gradient_dh.n_dofs();

  phi.reinit(comp_dom.dh.n_dofs());
  dphi_dn.reinit(comp_dom.dh.n_dofs());
  tmp_rhs.reinit(comp_dom.dh.n_dofs());

  DXDt_and_DphiDt_vector.resize(n_dofs());

  std::vector<Point<dim> > support_points(comp_dom.dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.dh, support_points);
  std::vector<Point<dim> > gradient_support_points(comp_dom.gradient_dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.gradient_dh, gradient_support_points);


  unsigned int j = dim-1;
  for(unsigned int i=j; i<comp_dom.gradient_dh.n_dofs(); i=i+dim)
      {
      DXDt_and_DphiDt_vector[i] = initial_wave_shape.value(gradient_support_points[i]);
      //bem.pcout<<DXDt_and_DphiDt_vector[i]<<std::endl;
      }

   for (unsigned int i=0; i<comp_dom.dh.n_dofs(); i++)
      {
      DXDt_and_DphiDt_vector[i+comp_dom.gradient_dh.n_dofs()] = initial_wave_potential.value(gradient_support_points[i]);
      //bem.pcout<<DXDt_and_DphiDt_vector[i+comp_dom.gradient_dh.n_dofs()]<<std::endl;
      }

   max_y_coor_value = 0;
   for (unsigned int i=0; i < comp_dom.dh.n_dofs(); i++)
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
void BoundaryConditions<dim>:: solve_problem()
{


  potential.set_time(0);
  wind.set_time(0);

  const unsigned int n_dofs =  comp_dom.dh.n_dofs();
  std::vector<types::subdomain_id> dofs_domain_association(n_dofs);
  DoFTools::get_subdomain_association   (comp_dom.dh,dofs_domain_association);
  this_cpu_set.set_size(n_dofs);
  for (unsigned int i=0; i<n_dofs; ++i)
    if (dofs_domain_association[i] == this_mpi_process)
      {
        this_cpu_set.add_index(i);
      }
  this_cpu_set.compress();

  phi.reinit(this_cpu_set,mpi_communicator);
  dphi_dn.reinit(this_cpu_set,mpi_communicator);
  tmp_rhs.reinit(this_cpu_set,mpi_communicator);

  prepare_bem_vectors();
  bem.solve(phi, dphi_dn, tmp_rhs);




}




template <int dim>
void BoundaryConditions<dim>::prepare_bem_vectors()
{

  comp_dom.compute_normals();

  const unsigned int n_dofs =  comp_dom.dh.n_dofs();

  phi.reinit(this_cpu_set,mpi_communicator);
  dphi_dn.reinit(this_cpu_set,mpi_communicator);
  tmp_rhs.reinit(this_cpu_set,mpi_communicator);


  std::vector<Point<dim> > support_points(n_dofs);
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.dh, support_points);

  cell_it
  cell = comp_dom.dh.begin_active(),
  endc = comp_dom.dh.end();

  const unsigned int   dofs_per_cell   = comp_dom.fe.dofs_per_cell;
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  FEValues<dim-1,dim> fe_v(*comp_dom.mapping, comp_dom.fe, *comp_dom.quadrature,
                           update_values |
                           update_cell_normal_vectors |
                           update_quadrature_points |
                           update_JxW_values);

  for (cell = comp_dom.dh.begin_active(); cell != endc; ++cell)
    {
      fe_v.reinit(cell);
      cell->get_dof_indices(local_dof_indices);
      //const std::vector<Point<dim> > &node_normals = fe_v.get_cell_normal_vectors();////provv
      for (unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j)
        if (this_cpu_set.is_element(local_dof_indices[j]))
          {
            //bem.pcout<<cell<<" "<<cell->material_id()<<" ("<<node_normals[0]<<") "<<-normals_sys_solution(local_dof_indices[j])<<std::endl;

            if (cell->material_id() != comp_dom.free_sur_ID1 &&
                cell->material_id() != comp_dom.free_sur_ID2 &&
                cell->material_id() != comp_dom.free_sur_ID3 )
              {
                if (cell->material_id() == comp_dom.wall_sur_ID1 ||
                    cell->material_id() == comp_dom.wall_sur_ID2 ||
                    cell->material_id() == comp_dom.wall_sur_ID3)
                  {
                    //tmp_rhs(local_dof_indices[j]) = normals_sys_solution(local_dof_indices[j]);
                    //dphi_dn(local_dof_indices[j]) = normals_sys_solution(local_dof_indices[j]);
                    Vector<double> imposed_pot_grad(dim);
                    wind.vector_value(support_points[local_dof_indices[j]],imposed_pot_grad);
                    Point<dim> imposed_potential_gradient;
                    for (unsigned int d=0; d<dim; ++d)
                      imposed_potential_gradient(d) = imposed_pot_grad(d);
                    tmp_rhs(local_dof_indices[j]) = imposed_potential_gradient*comp_dom.node_normals[local_dof_indices[j]];
                    dphi_dn(local_dof_indices[j]) = imposed_potential_gradient*comp_dom.node_normals[local_dof_indices[j]];
                  }
                else
                  {
                    tmp_rhs(local_dof_indices[j]) = 0;
                    dphi_dn(local_dof_indices[j]) = 0;
                  }
                //bem.pcout<<"internalIf   "<<local_dof_indices[j]<<" norm ("<<node_normals[j]<<")  "<<" pos ("<<node_coors[j]<<")    "<<node_normals[j]*datum_grad<<std::endl;
              }
            else
              {
                //tmp_rhs(local_dof_indices[j]) = node_coors[j](0);
                phi(local_dof_indices[j]) = potential.value(support_points[local_dof_indices[j]]);
                tmp_rhs(local_dof_indices[j]) = phi(local_dof_indices[j]);
                //bem.pcout<<"internalElse "<<local_dof_indices[j]<<" norm ("<<node_normals[j]<<")  "<<" pos ("<<node_coors[j]<<")    "<<node_coors[j](0)<<std::endl;
              }
            //bem.pcout<<tmp_rhs(local_dof_indices[j])<<"   phi "<<phi(local_dof_indices[j])<<std::endl;
          }

    }


}

template <int dim>
void BoundaryConditions<dim>::compute_errors()
{

  bem.compute_gradients(phi,dphi_dn);
  Vector<double> vector_gradients_solution(comp_dom.gradient_dh.n_dofs());
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
    for (unsigned int d=0; d<dim; ++d)
      vector_gradients_solution(3*i+d) = bem.node_gradients[i](d);

  Vector<float> grad_difference_per_cell (comp_dom.tria.n_active_cells());
  VectorTools::integrate_difference (*comp_dom.mapping, comp_dom.gradient_dh, vector_gradients_solution,
                                     wind,
                                     grad_difference_per_cell,
                                     QGauss<(dim-1)>(2*comp_dom.fe.degree+1),
                                     VectorTools::L2_norm);
  const double grad_L2_error = grad_difference_per_cell.l2_norm();

  Vector<float> difference_per_cell (comp_dom.tria.n_active_cells());
  VectorTools::integrate_difference (*comp_dom.mapping, comp_dom.dh, phi,
                                     potential,
                                     difference_per_cell,
                                     QGauss<(dim-1)>(2*comp_dom.fe.degree+1),
                                     VectorTools::L2_norm);
  const double L2_error = difference_per_cell.l2_norm();

  std::vector<Point<dim> > support_points(comp_dom.dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.dh, support_points);

  Vector<double> vector_gradients_node_error(comp_dom.gradient_dh.n_dofs());
  std::vector<Vector<double> > grads_nodes_errs(comp_dom.dh.n_dofs(),Vector<double>(dim));
  wind.vector_value_list(support_points,grads_nodes_errs);
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
    for (unsigned int d=0; d<dim; ++d)
      vector_gradients_node_error(3*i+d) = grads_nodes_errs[i](d);
  vector_gradients_node_error*=-1.0;
  vector_gradients_node_error.add(vector_gradients_solution);

  Vector<double> phi_node_error(comp_dom.dh.n_dofs());
  std::vector<double> phi_nodes_errs(comp_dom.dh.n_dofs());
  potential.value_list(support_points,phi_nodes_errs);
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
    phi_node_error(i) = phi_nodes_errs[i];

  phi_node_error*=-1.0;
  phi_node_error.add(phi);


  const double phi_max_error = phi_node_error.linfty_norm();
  const double grad_phi_max_error = vector_gradients_node_error.linfty_norm();
  const unsigned int n_active_cells=comp_dom.tria.n_active_cells();
  const unsigned int n_dofs=comp_dom.dh.n_dofs();

  std::cout << "   Number of active cells:       "
            << n_active_cells
            << std::endl
            << "   Number of degrees of freedom: "
            << n_dofs
            << std::endl
            ;

  std::cout<<"Phi Nodes error L_inf norm: "<<phi_max_error<<std::endl;
  std::cout<<"Phi Cells error L_2 norm: "<<L2_error<<std::endl;
  std::cout<<"Phi Nodes Gradient error L_inf norm: "<<grad_phi_max_error<<std::endl;
  std::cout<<"Phi Cells Gradient  error L_2 norm: "<<grad_L2_error<<std::endl;

}

template <int dim>
void BoundaryConditions<dim>::output_results(const std::string filename) const
{

  DataOut<dim-1, DoFHandler<dim-1, dim> > dataout;

  dataout.attach_dof_handler(comp_dom.gradient_dh);

  Vector<double> phi_dphidn_alpha(comp_dom.gradient_dh.n_dofs());

  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
    {
      phi_dphidn_alpha(3*i) = phi(i);
      phi_dphidn_alpha(3*i+1) = dphi_dn(i);
      phi_dphidn_alpha(3*i+2) = bem.alpha(i);
    }

  bem.compute_gradients(phi,dphi_dn);
  Vector<double> vector_gradients_solution(comp_dom.gradient_dh.n_dofs());
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
    for (unsigned int d=0; d<dim; ++d)
      vector_gradients_solution(3*i+d) = bem.node_gradients[i](d);

  comp_dom.compute_normals();
  Vector<double> vector_normals_solution(comp_dom.gradient_dh.n_dofs());
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
    for (unsigned int d=0; d<dim; ++d)
      vector_normals_solution(3*i+d) = comp_dom.node_normals[i](d);

  dataout.add_data_vector(phi_dphidn_alpha, "phi_dphidn_alpha");
  dataout.add_data_vector(vector_gradients_solution, "phi_gradient");
  dataout.add_data_vector(vector_normals_solution, "normals_at_nodes");
  dataout.build_patches(*comp_dom.mapping,
                        comp_dom.mapping->get_degree(),
                        DataOut<dim-1, DoFHandler<dim-1, dim> >::curved_inner_cells);

  std::ofstream file(filename.c_str());

  dataout.write_vtk(file);


}

//template class BoundaryConditions<3>;
