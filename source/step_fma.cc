#include "step_fma.h"


// @sect4{StepFMA::StepFMA and StepFMA::read_parameters}

// The constructor initializes the various object in much the same way as
// done in the finite element programs such as step-4 or step-6. The only
// new ingredient here is the ParsedFunction object, which needs, at
// construction time, the specification of the number of components.
//
// For the exact solution the number of vector components is one, and no
// action is required since one is the default value for a ParsedFunction
// object. The wind, however, requires dim components to be
// specified. Notice that when declaring entries in a parameter file for the
// expression of the Functions::ParsedFunction, we need to specify the
// number of components explicitly, since the function
// Functions::ParsedFunction::declare_parameters is static, and has no
// knowledge of the number of components.

template <int gdim>
void test(BEMFMA<gdim> &fma)
{
  fma.trunc_order += 1;
}

template <int dim>
MinFmm::StepFMA<dim>::StepFMA(const unsigned int fe_degree, bool fmm_method, const MPI_Comm comm)
  :
  mpi_communicator (comm),
  fe(fe_degree),
  dh(tria),
  mapping(fe_degree),
  fmm_sol(fmm_method),
  eh("","u","L2, H1, Linfty"),
  n_mpi_processes (Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
  this_mpi_process (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
  pcout(std::cout)
{
  pcout.set_condition(this_mpi_process == 0);
}


template <int dim>
void MinFmm::StepFMA<dim>::declare_parameters (ParameterHandler &prm)
{


  prm.declare_entry("Initial refinement", "0",
                    Patterns::Integer());

  prm.declare_entry("Number of cycles", "2",
                    Patterns::Integer());

  prm.declare_entry("Number of octree cycles", "10",
                    Patterns::Integer());

  prm.declare_entry("External refinement", "5",
                    Patterns::Integer());
  prm.declare_entry("Extend solution on the -2,2 box", "false",
                    Patterns::Bool());
  prm.declare_entry("Run 2d simulation", "false",
                    Patterns::Bool());

  prm.declare_entry("Error from theory", "true",
                    Patterns::Bool());

  prm.declare_entry("Run 3d simulation", "true",
                    Patterns::Bool());

  prm.enter_subsection("Quadrature rules");
  {
    prm.declare_entry("Quadrature type", "gauss",
                      Patterns::Selection(QuadratureSelector<(dim-1)>::get_quadrature_names()));
    prm.declare_entry("Quadrature order", "4", Patterns::Integer());
    prm.declare_entry("Singular quadrature order", "8", Patterns::Integer());
  }
  prm.leave_subsection();

  prm.enter_subsection("Exact dPhi dn solution 2d");
  {
    Functions::ParsedFunction<2>::declare_parameters(prm);
    prm.set("Function expression", "x + y");
  }
  prm.leave_subsection();

  prm.enter_subsection("Exact dPhi dn solution 3d");
  {
    Functions::ParsedFunction<3>::declare_parameters(prm);
    prm.set("Function expression", "(x + y + z)");
  }
  prm.leave_subsection();

  prm.enter_subsection("Exact Phi solution 2d");
  {
    Functions::ParsedFunction<2>::declare_parameters(prm);
    prm.set("Function expression", "x+y");
  }
  prm.leave_subsection();

  prm.enter_subsection("Exact Phi solution 3d");
  {
    Functions::ParsedFunction<3>::declare_parameters(prm);
    prm.set("Function expression", "x+y+z");
  }
  prm.leave_subsection();
  prm.enter_subsection("Solver");
  SolverControl::declare_parameters(prm);
  prm.leave_subsection();
}

template <int dim>
void MinFmm::StepFMA<dim>::parse_parameters (ParameterHandler &prm)
{
  initial_ref = prm.get_integer("Initial refinement");
  n_cycles = prm.get_integer("Number of cycles");
  n_cycles_octree = prm.get_integer("Number of octree cycles");
  from_theory = prm.get_bool("Error from theory");
  external_refinement = prm.get_integer("External refinement");
  extend_solution = prm.get_bool("Extend solution on the -2,2 box");

  prm.enter_subsection("Quadrature rules");
  {
    quadrature =
      std_cxx11::shared_ptr<Quadrature<dim-1> >
      (new QuadratureSelector<dim-1> (prm.get("Quadrature type"),
                                      prm.get_integer("Quadrature order")));
    singular_quadrature_order = prm.get_integer("Singular quadrature order");
  }
  prm.leave_subsection();

  prm.enter_subsection(std::string("Exact dPhi dn solution ")+
                       Utilities::int_to_string(dim)+std::string("d"));
  {
    exact_dphi_dn_solution.parse_parameters(prm);
  }
  prm.leave_subsection();

  prm.enter_subsection(std::string("Exact Phi solution ")+
                       Utilities::int_to_string(dim)+std::string("d"));
  {
    exact_phi_solution.parse_parameters(prm);
  }
  prm.leave_subsection();

  prm.enter_subsection("Solver");
  solver_control.parse_parameters(prm);
  prm.leave_subsection();


  run_in_this_dimension = prm.get_bool("Run " +
                                       Utilities::int_to_string(dim) +
                                       "d simulation");
}



template <int dim>
void MinFmm::StepFMA<dim>::read_domain()
{
  static const Point<dim> center = Point<dim>();
  static const SphericalManifold<dim-1, dim> manifold(center);

  std::ifstream in;
  switch (dim)
    {
    case 2:
      in.open ("../../../coarse_circle.inp");
      break;

    case 3:
      in.open ("../../../coarse_sphere.inp");
      break;

    default:
      Assert (false, ExcNotImplemented());
    }

  GridIn<dim-1, dim> gi;
  gi.attach_triangulation (tria);
  gi.read_ucd (in);

  tria.set_all_manifold_ids(1);
  tria.set_manifold(1, manifold);
  if (initial_ref)
    tria.refine_global(initial_ref);
}


// @sect4{StepFMA::refine_and_resize}

// This function globally refines the mesh, distributes degrees of freedom,
// and resizes matrices and vectors.

template <int dim>
void MinFmm::StepFMA<dim>::refine_and_resize(bool ref)
{
  // BEMFMA(dh, double_nodes_set, dirichlet_nodes, mapping);
  if (ref)
    tria.refine_global(1);

  dh.distribute_dofs(fe);

  const types::global_dof_index n_dofs =  dh.n_dofs();

  std::vector<types::subdomain_id> dofs_domain_association(n_dofs);
  DoFTools::get_subdomain_association   (dh,dofs_domain_association);
  this_cpu_set.clear();
  this_cpu_set.set_size(n_dofs);

  for (types::global_dof_index i=0; i<n_dofs; ++i)
    if (dofs_domain_association[i] == this_mpi_process)
      {
        this_cpu_set.add_index(i);
      }
  this_cpu_set.compress();


  system_rhs.reinit(this_cpu_set,mpi_communicator);
  phi.reinit(this_cpu_set,mpi_communicator);
  dphi_dn.reinit(this_cpu_set,mpi_communicator);
  dirichlet_values.reinit(this_cpu_set,mpi_communicator);
  dirichlet_nodes.reinit(this_cpu_set,mpi_communicator);
  neumann_values.reinit(this_cpu_set,mpi_communicator);

  tril_sp.reinit(this_cpu_set,this_cpu_set,mpi_communicator);
  // full_sparsity_pattern.reinit(sol.vector_partitioner(), n_dofs);
  for (types::global_dof_index i=0; i<n_dofs; ++i)
    if (this_cpu_set.is_element(i))
      {
        for (types::global_dof_index j=0; j<n_dofs; ++j)
          tril_sp.add(i,j);
      }
  // full_sparsity_pattern.compress();
  tril_sp.compress();
  system_matrix.reinit(tril_sp);
  system_rhs.reinit(this_cpu_set,mpi_communicator);
  system_alpha.reinit(this_cpu_set,mpi_communicator);

  // if(fma == NULL)
  //   fma = new BEMFMA(dh, double_nodes_set, dirichlet_nodes, mapping);

}


template <int dim>
void MinFmm::StepFMA<dim>::assemble_direct_system()
{

  FEValues<dim-1,dim> fe_v(mapping, fe, *quadrature,
                           update_values |
                           update_cell_normal_vectors |
                           update_quadrature_points |
                           update_JxW_values);

  const unsigned int n_q_points = fe_v.n_quadrature_points;

  std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);

  std::vector<double> cell_dphi_dn(n_q_points);
  std::vector<double> cell_phi(n_q_points);

  Vector<double>      local_matrix_row_i(fe.dofs_per_cell);

  std::vector<Point<dim> > support_points(dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( mapping, dh, support_points);

  typename DoFHandler<dim-1,dim>::active_cell_iterator
  cell = dh.begin_active(),
  endc = dh.end();

  for (cell = dh.begin_active(); cell != endc; ++cell)
    {
      fe_v.reinit(cell);
      cell->get_dof_indices(local_dof_indices);

      const std::vector<Point<dim> > &q_points = fe_v.get_quadrature_points();
      const std::vector<Tensor<1, dim> > &normals = fe_v.get_all_normal_vectors();
      exact_phi_solution.value_list(q_points, cell_phi);
      exact_dphi_dn_solution.value_list(q_points, cell_dphi_dn);
      for (types::global_dof_index i=0; i<dh.n_dofs() ; ++i)
        {
          local_matrix_row_i = 0;

          bool is_singular = false;
          unsigned int singular_index = numbers::invalid_unsigned_int;

          for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
            if (local_dof_indices[j] == i)
              {
                singular_index = j;
                is_singular = true;
                break;
              }

          if (dirichlet_nodes[i]==0)
            {

              if (is_singular == false)
                {
                  for (unsigned int q=0; q<n_q_points; ++q)
                    {

                      const Tensor<1,dim> R = q_points[q] - support_points[i];

                      system_rhs(i) += ( LaplaceKernel::single_layer(R)   *
                                         cell_dphi_dn[q]                      *
                                         fe_v.JxW(q) );

                      for (unsigned int j=0; j<fe.dofs_per_cell; ++j)

                        local_matrix_row_i(j) -= ( ( LaplaceKernel::double_layer(R)     *
                                                     normals[q] )            *
                                                   fe_v.shape_value(j,q)     *
                                                   fe_v.JxW(q)       );
                    }
                }
              else
                {
                  Assert(singular_index != numbers::invalid_unsigned_int,
                         ExcInternalError());

                  const Quadrature<dim-1> & singular_quadrature =
                    get_singular_quadrature(singular_index);

                  FEValues<dim-1,dim> fe_v_singular (mapping, fe, singular_quadrature,
                                                     update_jacobians |
                                                     update_values |
                                                     update_cell_normal_vectors |
                                                     update_quadrature_points );

                  fe_v_singular.reinit(cell);

                  std::vector<double> singular_cell_dphi_dn( singular_quadrature.size());

                  const std::vector<Tensor<1,dim> > &singular_normals = fe_v_singular.get_all_normal_vectors();
                  const std::vector<Point<dim> > &singular_q_points = fe_v_singular.get_quadrature_points();

                  exact_phi_solution.value_list(singular_q_points, singular_cell_dphi_dn);

                  for (unsigned int q=0; q<singular_quadrature.size(); ++q)
                    {
                      const Tensor<1,dim> R = singular_q_points[q] - support_points[i];

                      system_rhs(i) += ( LaplaceKernel::single_layer(R)   *
                                         singular_cell_dphi_dn[q]                      *
                                         fe_v_singular.JxW(q) );

                      for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
                        {
                          local_matrix_row_i(j) -= (( LaplaceKernel::double_layer(R) *
                                                      singular_normals[q])                *
                                                    fe_v_singular.shape_value(j,q)        *
                                                    fe_v_singular.JxW(q)       );
                        }
                    }
                }

              for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
                system_matrix.add(i,local_dof_indices[j],local_matrix_row_i(j));
            }
          else
            {

              if (is_singular == false)
                {
                  for (unsigned int q=0; q<n_q_points; ++q)
                    {

                      const Tensor<1,dim> R = q_points[q] - support_points[i];

                      system_rhs(i) += ( ( LaplaceKernel::double_layer(R)     *
                                           normals[q] )            *
                                         cell_phi[q]              *
                                         fe_v.JxW(q)       );

                      for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
                        local_matrix_row_i(j) -= ( LaplaceKernel::single_layer(R)   *
                                                   fe_v.shape_value(j,q)      *
                                                   fe_v.JxW(q) );
                    }
                }
              else
                {
                  Assert(singular_index != numbers::invalid_unsigned_int,
                         ExcInternalError());

                  const Quadrature<dim-1> & singular_quadrature =
                    get_singular_quadrature(singular_index);

                  FEValues<dim-1,dim> fe_v_singular (mapping, fe, singular_quadrature,
                                                     update_jacobians |
                                                     update_values |
                                                     update_cell_normal_vectors |
                                                     update_quadrature_points );

                  fe_v_singular.reinit(cell);

                  std::vector<double> singular_cell_phi( singular_quadrature.size());

                  const std::vector<Tensor<1, dim> > &singular_normals = fe_v_singular.get_all_normal_vectors();
                  const std::vector<Point<dim> > &singular_q_points = fe_v_singular.get_quadrature_points();

                  exact_phi_solution.value_list(singular_q_points, singular_cell_phi);

                  for (unsigned int q=0; q<singular_quadrature.size(); ++q)
                    {
                      const Tensor<1,dim> R = singular_q_points[q] - support_points[i];

                      system_rhs(i) += (( LaplaceKernel::double_layer(R) *
                                          singular_normals[q])                *
                                        singular_cell_phi[q]        *
                                        fe_v_singular.JxW(q)       );

                      for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
                        {
                          local_matrix_row_i(j) -= ( LaplaceKernel::single_layer(R)   *
                                                     fe_v_singular.shape_value(j,q)        *
                                                     fe_v_singular.JxW(q) );
                        }
                    }
                }

              for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
                system_matrix.add(i,local_dof_indices[j],local_matrix_row_i(j));


            }
        }
    }

  TrilinosWrappers::MPI::Vector ones(this_cpu_set, mpi_communicator);
  vector_shift(ones, -1.);

  system_matrix.vmult(system_alpha, ones);
  vector_shift(system_alpha, 1.);
  std::ofstream ofs;
  ofs.open ("direct_alpha.txt", std::ofstream::out | std::ofstream::app);
  ofs << system_alpha.linfty_norm() << " "<< system_alpha.l2_norm()<< std::endl;;
  ofs.close();

  for (types::global_dof_index i = 0; i<dh.n_dofs(); ++i)
    system_matrix.add(i,i,system_alpha[i]);
}

template <int dim>
void MinFmm::StepFMA<dim>::compute_boundary_condition()
{
  std::vector<Point<dim> > support_points(dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( mapping, dh, support_points);
  for (types::global_dof_index i=0; i<dh.n_dofs(); ++i)
    {
      if (support_points[i][0] < -4)
        dirichlet_nodes[i] = 1;
      else
        dirichlet_nodes[i] = 0;

      dirichlet_values[i]=exact_phi_solution.value(support_points[i]);
      neumann_values[i]=exact_dphi_dn_solution.value(support_points[i]);



    }

}

template <int dim>
void MinFmm::StepFMA<dim>::compute_double_nodes_set()
{
  double tol=1e-4;
  double_nodes_set.clear();
  double_nodes_set.resize(dh.n_dofs());
  std::vector<Point<dim> > support_points(dh.n_dofs());

  DoFTools::map_dofs_to_support_points<dim-1, dim>( mapping,
                                                    dh, support_points);

  for (types::global_dof_index i=0; i<dh.n_dofs(); ++i)
    {
      for (types::global_dof_index j=0; j<dh.n_dofs(); ++j)
        {
          if (support_points[i].distance(support_points[j]) < tol)
            {
              double_nodes_set[i].insert(j);
            }
        }

    }


}
// @sect4{StepFMA::solve_system}

// The next function simply solves the linear system.
template <int dim>
void MinFmm::StepFMA<dim>::solve_system()
{
  if (fmm_sol)
    {
      Operator::MinBEMOperator<dim> oppy(fma, mpi_communicator, this_cpu_set, this_mpi_process);
      fma.generate_octree_blocking();
      fma.direct_integrals();
      fma.multipole_integrals();
      std::cout<<"setting alpha"<<std::endl;
      oppy.set_alpha();
      system_alpha = oppy.get_alpha();
      std::cout<<"computing rhs"<<std::endl;
      oppy.compute_rhs(system_rhs, dirichlet_values, neumann_values);
      std::cout<<"solving"<<std::endl;
      SolverGMRES<TrilinosWrappers::MPI::Vector > solver (solver_control);
      solver.solve (oppy, phi, system_rhs, PreconditionIdentity());

      for (types::global_dof_index i=0; i<dh.n_dofs(); ++i)
        {
          if (dirichlet_nodes[i]==0)
            dphi_dn[i] = neumann_values[i];
          else
            {
              dphi_dn[i] = phi[i];
              phi[i] = dirichlet_values[i];
            }
        }
    }
  else
    {
      SolverGMRES<TrilinosWrappers::MPI::Vector > solver (solver_control);
      solver.solve (system_matrix, phi, system_rhs, PreconditionIdentity());
      for (types::global_dof_index i=0; i<dh.n_dofs(); ++i)
        {
          if (dirichlet_nodes[i]==0)
            {
              dphi_dn[i] = neumann_values[i];
              phi[i] *= 2;
              // phi[i] = dirichlet_values[i];
              // if(phi[i]<1e-4)
              //   std::cout<<"errore"<<std::endl;
            }
          else
            {
              dphi_dn[i] = phi[i];
              phi[i] = dirichlet_values[i];
              // dphi_dn[i] = neumann_values[i];
              // if(phi[i]<1e-4)
              //   std::cout<<"errore"<<std::endl;

            }
        }


    }
}


template <int dim>
void MinFmm::StepFMA<dim>::compute_errors(const unsigned int cycle)
{
  if (from_theory)
    {
      std::cout<<"using sak to compute the errors"<<std::endl;
      eh.error_from_exact(mapping, dh, phi, exact_phi_solution,0);
      eh.error_from_exact(mapping, dh, dphi_dn, exact_dphi_dn_solution,1);
    }
  else
    {
      if (!fmm_sol)
        {
          save_direct_solution(n_cycles);
          std::cout<<"using sak to compute the errors"<<std::endl;
          eh.error_from_exact(mapping, dh, phi, exact_phi_solution,0);
          eh.error_from_exact(mapping, dh, dphi_dn, exact_dphi_dn_solution,1);

        }
      else
        {
          std::cout<<"using sak to compute the errors"<<std::endl;
          Vector<double> phi_dummy(phi);
          Vector<double> alpha_dummy(system_alpha);
          Vector<double> phi_direct(phi_dummy.size());
          Vector<double> alpha_dir(alpha_dummy.size());
          read_direct_solution(n_cycles, phi_direct, alpha_dir);
          phi_direct -= phi_dummy;
          alpha_dir -= alpha_dummy;
          eh.error_from_exact(mapping, dh, phi, exact_phi_solution,0);
          eh.error_from_exact(mapping, dh, alpha_dir, ZeroFunction<dim>(),1);

        }
    }


}

template <int dim>
void MinFmm::StepFMA<dim>::save_direct_solution(const unsigned int cycle)
{
  std::string file_name1, file_name2;
  Vector<double> phi_dir(phi);
  Vector<double> alpha_dir(system_alpha);
  vector_shift(alpha_dir, -1.);
  alpha_dir*=-1.;
  file_name1 = "direct_solution_" + Utilities::int_to_string(cycle) + ".bin";
  std::ofstream dir_sol (file_name1.c_str());
  phi_dir.block_write(dir_sol);
  file_name2 = "alpha_direct_solution_" + Utilities::int_to_string(cycle) + ".bin";
  std::ofstream alpha_dir_sol (file_name2.c_str());
  alpha_dir.block_write(alpha_dir_sol);

}

template <int dim>
void MinFmm::StepFMA<dim>::read_direct_solution(const unsigned int cycle, Vector<double> &phi_direct, Vector<double> &alpha_dir)
{
  std::string file_name1, file_name2;
  file_name1 = "direct_solution_" + Utilities::int_to_string(cycle) + ".bin";
  std::ifstream dir_sol (file_name1.c_str());
  phi_direct.block_read(dir_sol);
  file_name2 = "alpha_direct_solution_" + Utilities::int_to_string(cycle) + ".bin";
  std::ifstream alpha_dir_sol (file_name2.c_str());
  alpha_dir.block_read(alpha_dir_sol);

}

// Singular integration requires a careful selection of the quadrature
// rules. In particular the deal.II library provides quadrature rules which
// are tailored for logarithmic singularities (QGaussLog, QGaussLogR), as
// well as for 1/R singularities (QGaussOneOverR).
//
// Singular integration is typically obtained by constructing weighted
// quadrature formulas with singular weights, so that it is possible to
// write
//
// \f[ \int_K f(x) s(x) dx = \sum_{i=1}^N w_i f(q_i) \f]
//
// where $s(x)$ is a given singularity, and the weights and quadrature
// points $w_i,q_i$ are carefully selected to make the formula above an
// equality for a certain class of functions $f(x)$.
//
// In all the finite element examples we have seen so far, the weight of the
// quadrature itself (namely, the function $s(x)$), was always constantly
// equal to 1.  For singular integration, we have two choices: we can use
// the definition above, factoring out the singularity from the integrand
// (i.e., integrating $f(x)$ with the special quadrature rule), or we can
// ask the quadrature rule to "normalize" the weights $w_i$ with $s(q_i)$:
//
// \f[ \int_K f(x) s(x) dx = \int_K g(x) dx = \sum_{i=1}^N
//   \frac{w_i}{s(q_i)} g(q_i) \f]
//
// We use this second option, through the @p factor_out_singularity
// parameter of both QGaussLogR and QGaussOneOverR.
//
// These integrals are somewhat delicate, especially in two dimensions, due
// to the transformation from the real to the reference cell, where the
// variable of integration is scaled with the determinant of the
// transformation.
//
// In two dimensions this process does not result only in a factor appearing
// as a constant factor on the entire integral, but also on an additional
// integral altogether that needs to be evaluated:
//
// \f[ \int_0^1 f(x)\ln(x/\alpha) dx = \int_0^1 f(x)\ln(x) dx - \int_0^1
//  f(x) \ln(\alpha) dx.  \f]
//
// This process is taken care of by the constructor of the QGaussLogR class,
// which adds additional quadrature points and weights to take into
// consideration also the second part of the integral.
//
// A similar reasoning should be done in the three dimensional case, since
// the singular quadrature is tailored on the inverse of the radius $r$ in
// the reference cell, while our singular function lives in real space,
// however in the three dimensional case everything is simpler because the
// singularity scales linearly with the determinant of the
// transformation. This allows us to build the singular two dimensional
// quadrature rules only once and, reuse them over all cells.
//
// In the one dimensional singular integration this is not possible, since
// we need to know the scaling parameter for the quadrature, which is not
// known a priori. Here, the quadrature rule itself depends also on the size
// of the current cell. For this reason, it is necessary to create a new
// quadrature for each singular integration.
//
// The different quadrature rules are built inside the
// get_singular_quadrature, which is specialized for dim=2 and dim=3, and
// they are retrieved inside the assemble_direct_system function. The index given
// as an argument is the index of the unit support point where the
// singularity is located.

template<int dim>
const Quadrature<dim-1> & MinFmm::StepFMA<dim>::get_singular_quadrature(const unsigned int index) const
{
  Assert(index < fe.dofs_per_cell,
         ExcIndexRange(0, fe.dofs_per_cell, index));

  //static std::vector<QGaussOneOverR<2> > quadratures;
  static std::vector<QTelles<dim-1> > quadratures;
  if (quadratures.size() == 0)
    for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      quadratures.push_back(QTelles<dim-1>(singular_quadrature_order,
                                           fe.get_unit_support_points()[i]));

  return quadratures[index];
}

//
// template<>
// const Quadrature<1> & MinFmm::StepFMA<2>::get_singular_quadrature(const unsigned int index) const
// {
//     Assert(index < fe.dofs_per_cell,
//            ExcIndexRange(0, fe.dofs_per_cell, index));
//
//     static std::vector<QTelles<1> > quadratures;
//     if (quadratures.size() == 0)
//         for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
//             quadratures.push_back(QTelles<1>(singular_quadrature_order,
//                                              fe.get_unit_support_points()[i]));
//
//     return quadratures[index];
// }

//
// template<>
// const Quadrature<1> &MinFmm::StepFMA<2>::get_singular_quadrature(
//   const DoFHandler<1,2>::active_cell_iterator &cell,
//   const unsigned int index) const
// {
//   Assert(index < fe.dofs_per_cell,
//          ExcIndexRange(0, fe.dofs_per_cell, index));
//
//   static Quadrature<1> *q_pointer = NULL;
//   if (q_pointer) delete q_pointer;
//
//   q_pointer = new QGaussLogR<1>(singular_quadrature_order,
//                                 fe.get_unit_support_points()[index],
//                                 1./cell->measure(), true);
//   return (*q_pointer);
// }



// @sect4{StepFMA::compute_exterior_solution}

// We'd like to also know something about the value of the potential $\phi$
// in the exterior domain: after all our motivation to consider the boundary
// integral problem was that we wanted to know the velocity in the exterior
// domain!
//
// To this end, let us assume here that the boundary element domain is
// contained in the box $[-2,2]^{\text{dim}}$, and we extrapolate the actual
// solution inside this box using the convolution with the fundamental
// solution. The formula for this is given in the introduction.
//
// The reconstruction of the solution in the entire space is done on a
// continuous finite element grid of dimension dim. These are the usual
// ones, and we don't comment any further on them. At the end of the
// function, we output this exterior solution in, again, much the usual way.
// template <int dim>
// void MinFmm::StepFMA<dim>::compute_exterior_solution()
// {
//   Triangulation<dim>  external_tria;
//   GridGenerator::hyper_cube(external_tria, -2, 2);
//
//   FE_Q<dim>           external_fe(1);
//   DoFHandler<dim>     external_dh (external_tria);
//   TrilinosWrappers::MPI::Vector      external_phi;
//
//   external_tria.refine_global(external_refinement);
//   external_dh.distribute_dofs(external_fe);
//   external_phi.reinit(external_dh.n_dofs());
//
//   typename DoFHandler<dim-1,dim>::active_cell_iterator
//   cell = dh.begin_active(),
//   endc = dh.end();
//
//
//   FEValues<dim-1,dim> fe_v(mapping, fe, *quadrature,
//                            update_values |
//                            update_cell_normal_vectors |
//                            update_quadrature_points |
//                            update_JxW_values);
//
//   const unsigned int n_q_points = fe_v.n_quadrature_points;
//
//   std::vector<types::global_dof_index> dofs(fe.dofs_per_cell);
//
//   std::TrilinosWrappers::MPI::Vector local_phi(n_q_points);
//   std::TrilinosWrappers::MPI::Vector normal_wind(n_q_points);
//   std::vector<TrilinosWrappers::MPI::Vector > local_wind(n_q_points, TrilinosWrappers::MPI::Vector(dim) );
//
//   std::vector<Point<dim> > external_support_points(external_dh.n_dofs());
//   DoFTools::map_dofs_to_support_points<dim>(StaticMappingQ1<dim>::mapping,
//                                             external_dh, external_support_points);
//
//   for (cell = dh.begin_active(); cell != endc; ++cell)
//     {
//       fe_v.reinit(cell);
//
//       const std::vector<Point<dim> > &q_points = fe_v.get_quadrature_points();
//       const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors();
//
//       cell->get_dof_indices(dofs);
//       fe_v.get_function_values(phi, local_phi);
//
//       wind.vector_value_list(q_points, local_wind);
//
//       for (unsigned int q=0; q<n_q_points; ++q)
//         {
//           normal_wind[q] = 0;
//           for (unsigned int d=0; d<dim; ++d)
//             normal_wind[q] += normals[q][d]*local_wind[q](d);
//         }
//
//       for (unsigned int i=0; i<external_dh.n_dofs(); ++i)
//         for (unsigned int q=0; q<n_q_points; ++q)
//           {
//
//             const Tensor<1,dim> R = q_points[q] - external_support_points[i];
//
//             external_phi(i) += ( ( LaplaceKernel::single_layer(R) *
//                                    normal_wind[q]
//                                    +
//                                    (LaplaceKernel::double_layer(R) *
//                                     normals[q] )            *
//                                    local_phi[q] )           *
//                                  fe_v.JxW(q) );
//           }
//     }
//
//   DataOut<dim> data_out;
//
//   data_out.attach_dof_handler(external_dh);
//   data_out.add_data_vector(external_phi, "external_phi");
//   data_out.build_patches();
//
//   const std::string
//   filename = Utilities::int_to_string(dim) + "d_external.vtk";
//   std::ofstream file(filename.c_str());
//
//   data_out.write_vtk(file);
// }


// @sect4{StepFMA::output_results}

// Outputting the results of our computations is a rather mechanical
// tasks. All the components of this function have been discussed before.
template <int dim>
void MinFmm::StepFMA<dim>::output_results(const unsigned int cycle)
{
  DataOut<dim-1, DoFHandler<dim-1, dim> > dataout;

  dataout.attach_dof_handler(dh);
  dataout.add_data_vector(phi, "phi",
                          DataOut<dim-1, DoFHandler<dim-1, dim> >::type_dof_data);
  dataout.add_data_vector(dphi_dn, "dphi_dn",
                          DataOut<dim-1, DoFHandler<dim-1, dim> >::type_dof_data);
  dataout.build_patches(mapping,
                        mapping.get_degree(),
                        DataOut<dim-1, DoFHandler<dim-1, dim> >::curved_inner_cells);
  std::string multipole_str;

  if (fmm_sol)
    multipole_str = "_fmm";
  else
    multipole_str = "_direct";
  std::string filename = ( Utilities::int_to_string(dim) +
                           "d_boundary_solution_" +
                           Utilities::int_to_string(cycle) +
                           multipole_str +
                           ".vtu" );
  std::ofstream file(filename.c_str());

  dataout.write_vtu(file);

}


// @sect4{StepFMA::run}

// This is the main function. It should be self explanatory in its
// briefness:
template <int dim>
void MinFmm::StepFMA<dim>::run()
{
  // ParameterHandler prm;
  //
  // declare_parameters(prm);
  // parse_parameters(prm);

  //read_parameters("parameters.prm");

  if (run_in_this_dimension == false)
    {
      deallog << "Run in dimension " << dim
              << " explicitly disabled in parameter file. "
              << std::endl;
      return;
    }

  read_domain();

  for (unsigned int cycle=0; cycle<n_cycles; ++cycle)
    {
      refine_and_resize();
      compute_double_nodes_set();
      compute_boundary_condition();
      if (cycle == 0)
        fma.init_fma(dh, double_nodes_set, dirichlet_nodes, mapping);
      if (!fmm_sol)
        assemble_direct_system();
      solve_system();
      compute_errors(cycle);
      output_results(cycle);
    }
  std::cout<<"using sak to print the output"<<std::endl;
  eh.output_table(std::cout,0);
  eh.output_table(std::cout,1);
  tria.set_manifold(1);
  // if (extend_solution == true)
  //   compute_exterior_solution();
}


template <int dim>
void MinFmm::StepFMA<dim>::run_for_octree()
{
  // ParameterHandler prm;
  //
  // declare_parameters(prm);
  // parse_parameters(prm);

  //read_parameters("parameters.prm");

  if (run_in_this_dimension == false)
    {
      deallog << "Run in dimension " << dim
              << " explicitly disabled in parameter file. "
              << std::endl;
      return;
    }

  read_domain();

  for (unsigned int cycle=0; cycle<n_cycles_octree; ++cycle)
    {
      if (cycle==0)
        refine_and_resize();
      else
        refine_and_resize(false);
      compute_double_nodes_set();
      compute_boundary_condition();
      if (cycle == 0)
        fma.init_fma(dh, double_nodes_set, dirichlet_nodes, mapping);
      if (!fmm_sol)
        assemble_direct_system();
      solve_system();
      compute_errors(cycle);
      output_results(cycle);
      test(fma);
    }
  std::cout<<"using sak to print the output"<<std::endl;
  eh.output_table(std::cout,0);
  eh.output_table(std::cout,1);
  tria.set_manifold(1);
  // if (extend_solution == true)
  //   compute_exterior_solution();
}

template class MinFmm::StepFMA<3>;
