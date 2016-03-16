#include <preconditioner_bem.h>
#include <laplace_kernel.h>
// #include <bem_problem_access.h>
// template <int dim>
// void PreconditionerBEM<dim>::initialize(const BEMProblem<dim> * in_bem)
// {
//   bem_problem = in_bem;
// }

template <int dim>
void PreconditionerBEM<dim>::build_coarse_inverse(unsigned int level_mg)//, const TrilinosWrappers::MPI::Vector &dirichlet_nodes, const ConstraintMatrix &constraints, const Mapping<dim-1,dim> mapping)
{

  auto &mapping = this->get_bem_mapping();
  auto &dh = this->get_bem_dh();
  auto &quadrature = this->get_bem_quadrature();
  auto singular_quadrature_order = this->get_bem_singular_quadrature_order();
  auto &constraints = this->get_bem_constraint_matrix();
  auto &fe = this->get_bem_fe();
  auto &dirichlet_nodes = this->get_bem_dirichlet_nodes();
  auto &double_nodes_set = this->get_bem_double_nodes_set();

  constraints_coarse.clear();
  constraints_coarse.merge(constraints);
  constraints_coarse.close();
  tria_coarse.copy_triangulation(dh.get_triangulation());
  dh_coarse.initialize(tria_coarse, fe);
  // dh.distribute_mg_dofs(fe);
  // First We need to set up the sparsity pattern for the matrix.
  sp_coarse_full.reinit(dh_coarse.n_dofs(),dh_coarse.n_dofs(),dh_coarse.n_dofs());
  for(size_type i=0; i<dh_coarse.n_dofs(); ++i)
  {
    if(constraints.is_constrained(i))
    {
      const std::vector< std::pair < types::global_dof_index, double > >
      *entries = constraints.get_constraint_entries (i);
      sp_coarse_full.add(i,i);
      for (types::global_dof_index j=0; j< entries->size(); ++j)
        sp_coarse_full.add(i,(*entries)[j].first);

    }
    else
      for(size_type j=0; j<dh_coarse.n_dofs(); ++j)
        sp_coarse_full.add(i,j);

  }
  sp_coarse_full.compress();
  prec_matrix.reinit(sp_coarse_full);



  std::vector<QTelles<dim-1> > sing_quadratures;
  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    sing_quadratures.push_back
    (QTelles<dim-1>(singular_quadrature_order,
                    fe.get_unit_support_points()[i]));

  FEValues<dim-1,dim> fe_v(mapping,fe, quadrature,
                           update_values |
                           update_cell_normal_vectors |
                           update_quadrature_points |
                           update_JxW_values);

  const unsigned int n_q_points = fe_v.n_quadrature_points;

  std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);
  Vector<double>      local_matrix_row_i(fe.dofs_per_cell);
  std::vector<Point<dim> > support_points(dh_coarse.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( mapping, dh_coarse, support_points);

  typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;
  cell_it
  cell = dh_coarse.begin_active(),
  endc = dh_coarse.end();

  Point<dim> D;
  double s;
  std::cout<<"looping for the inverse"<<std::endl;

  for (cell = dh_coarse.begin_active(); cell != endc; ++cell)
    {
      fe_v.reinit(cell);
      cell->get_dof_indices(local_dof_indices);

      const std::vector<Point<dim> > &q_points = fe_v.get_quadrature_points();
      const std::vector<Tensor<1, dim> > &normals = fe_v.get_all_normal_vectors();

      for (types::global_dof_index i=0; i< dh_coarse.n_dofs() ; ++i) //these must now be the locally owned dofs. the rest should stay the same
      {
            local_matrix_row_i = 0;

            bool is_singular = false;
            unsigned int singular_index = numbers::invalid_unsigned_int;

            for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
              //if(local_dof_indices[j] == i)
              if (double_nodes_set[i].count(local_dof_indices[j]) > 0)
                {
                  singular_index = j;
                  is_singular = true;
                  break;
                }
            if(dirichlet_nodes[i] == 0)
            {
            if (is_singular == false)
              {
                for (unsigned int q=0; q<n_q_points; ++q)
                  {
                    const Tensor<1,dim> R = q_points[q] - support_points[i];
                    LaplaceKernel::kernels(R, D, s);
                    // if(support_points[i][0]==0.25&&support_points[i][1]==0.25)
                    //   pcout<<"D "<<D<<" s "<<s<<" , ";
                    for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
                      {
                        local_matrix_row_i(j) += ( ( D *
                                                             normals[q] ) *
                                                           fe_v.shape_value(j,q) *
                                                           fe_v.JxW(q)       );
                      }
                  }
              }
            else
              {
                Assert(singular_index != numbers::invalid_unsigned_int,
                       ExcInternalError());

                const Quadrature<dim-1> *
                singular_quadrature
                  = dynamic_cast<Quadrature<dim-1>*>(
                      &sing_quadratures[singular_index]);
                Assert(singular_quadrature, ExcInternalError());

                FEValues<dim-1,dim> fe_v_singular (mapping, fe, *singular_quadrature,
                                                   update_jacobians |
                                                   update_values |
                                                   update_cell_normal_vectors |
                                                   update_quadrature_points );

                fe_v_singular.reinit(cell);

                const std::vector<Tensor<1, dim> > &singular_normals = fe_v_singular.get_all_normal_vectors();
                const std::vector<Point<dim> > &singular_q_points = fe_v_singular.get_quadrature_points();

                for (unsigned int q=0; q<singular_quadrature->size(); ++q)
                  {
                    const Tensor<1,dim> R = singular_q_points[q] - support_points[i];
                    LaplaceKernel::kernels(R, D, s);

                    for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
                      {
                        local_matrix_row_i(j) += (( D *
                                                            singular_normals[q])                *
                                                          fe_v_singular.shape_value(j,q)        *
                                                          fe_v_singular.JxW(q)       );

                        }
                  }
              }

            // Finally, we need to add
            // the contributions of the
            // current cell to the
            // global matrix.
            }
            else
            {
              if (is_singular == false)
                {
                  for (unsigned int q=0; q<n_q_points; ++q)
                    {
                      const Tensor<1,dim> R = q_points[q] - support_points[i];
                      LaplaceKernel::kernels(R, D, s);
                      // if(support_points[i][0]==0.25&&support_points[i][1]==0.25)
                      //   pcout<<"D "<<D<<" s "<<s<<" , ";
                      for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
                        {
                          local_matrix_row_i(j) -= ( s *
                                                               fe_v.shape_value(j,q) *
                                                               fe_v.JxW(q) );

                        }
                    }
                }
              else
                {
                  Assert(singular_index != numbers::invalid_unsigned_int,
                         ExcInternalError());

                  const Quadrature<dim-1> *
                  singular_quadrature
                    = dynamic_cast<Quadrature<dim-1>*>(
                        &sing_quadratures[singular_index]);
                  Assert(singular_quadrature, ExcInternalError());

                  FEValues<dim-1,dim> fe_v_singular (mapping, fe, *singular_quadrature,
                                                     update_jacobians |
                                                     update_values |
                                                     update_cell_normal_vectors |
                                                     update_quadrature_points );

                  fe_v_singular.reinit(cell);

                  const std::vector<Tensor<1, dim> > &singular_normals = fe_v_singular.get_all_normal_vectors();
                  const std::vector<Point<dim> > &singular_q_points = fe_v_singular.get_quadrature_points();

                  for (unsigned int q=0; q<singular_quadrature->size(); ++q)
                    {
                      const Tensor<1,dim> R = singular_q_points[q] - support_points[i];
                      LaplaceKernel::kernels(R, D, s);

                      for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
                        {
                          local_matrix_row_i(j) -= ( s   *
                                                               fe_v_singular.shape_value(j,q)     *
                                                               fe_v_singular.JxW(q) );
                        }
                    }
                }

              // Finally, we need to add
              // the contributions of the
              // current cell to the
              // global matrix.

            }
            if(constraints.is_constrained(i))
            {
              if (constraints.is_constrained(i))
                {
                  prec_matrix.set(i,i,1);
                  const std::vector< std::pair < types::global_dof_index, double > >
                  *entries = constraints.get_constraint_entries (i);
                  for (types::global_dof_index j=0; j< entries->size(); ++j)
                    {
                      prec_matrix.set(i,(*entries)[j].first,(*entries)[j].second);
                      //pcout<<i<<" "<<(*entries)[j].first<<"  * "<<(*entries)[j].second<<std::endl;
                    }
                }

            }
            else
              for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
                {
                  prec_matrix.add(i,local_dof_indices[j],local_matrix_row_i(j));
                }
          }
    }

    // Computing alpha
    n_coarse = prec_matrix.n();
    std::cout<<"Computing alpha"<<std::endl;

    Vector<double> ones(n_coarse), alpha(n_coarse);

    vector_shift(ones, -1.);

    prec_matrix.vmult(alpha, ones);
    // prec_matrix.print_formatted(std::cout);
    std::cout<<"Adding alpha"<<std::endl;

    for(size_type i =0 ; i<n_coarse; ++i)
      if(dirichlet_nodes[i]==0 && !(constraints.is_constrained(i)))
        prec_matrix.add(i,i,alpha[i]);
    prec_matrix.compress(VectorOperation::add);
    // prec_matrix.print_formatted(std::cout);
    std::cout<<"initializing solver "<<n_coarse<<std::endl;
    // prec_solver.clear();
    prec_solver.initialize(prec_matrix);


}


template <int dim>
void PreconditionerBEM<dim>::build_projector()
{
  // auto mapping = this->get_bem_mapping();
  // auto dh = this->get_bem_dh();
  // auto quadrature = this->get_bem_quadrature();
  // auto fe = this->get_bem_fe();
  // FEValues<dim-1,dim> fe_v(mapping,fe, quadrature,
  //                          update_values );
  //
  // typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;
  // std::vector<Point<dim> > support_points(dh.n_dofs());
  // DoFTools::map_dofs_to_support_points<dim-1, dim>( mapping, dh, support_points);
  //
  // cell_it
  // cell = dh.begin_mg(level_mg),
  // endc = dh.end_mg(level_mg);
  // // which level??? (check sui support point per determinarlo)
  // for (cell = dh.begin_mg(level_mg); cell != endc; ++cell)
  // {
  //   fe_v.reinit(cell);
  //   cell->get_dof_indices(local_dof_indices);
  //   for(size_type i = 0; i<support_points.size(); ++i)
  //   {
  //
  //   }
  // }
  //
}
//
//
template <int dim>
void PreconditionerBEM<dim>::vmult(TrilinosWrappers::MPI::Vector &out, const TrilinosWrappers::MPI::Vector &in) const
{
  unsigned int level_mg = 1;
  Vector<double> loc(in);

  Vector<double> in_coarse(n_coarse);

  Vector<double> out_coarse(n_coarse);

  auto &dh = this->get_bem_dh();

  auto &constraints_fine = this->get_bem_constraint_matrix();

  std::cout<<loc.size()<<" "<<dh.n_dofs()<<" "<<in_coarse.size()<<" "<<dh_coarse.n_dofs()<<std::endl;
  VectorTools::interpolate_to_different_mesh(dh,loc,dh_coarse,constraints_coarse,in_coarse);

  prec_solver.vmult(out_coarse, in_coarse);

  loc = 0.;

  VectorTools::interpolate_to_different_mesh(dh_coarse,out_coarse,dh,constraints_fine,loc);

  for(auto i : out.locally_owned_elements())
    out[i] = loc[i];

  out.compress(VectorOperation::insert);

}

template class PreconditionerBEM<3>;
