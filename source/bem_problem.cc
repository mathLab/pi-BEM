

#include "../include/bem_problem.h"
#include "../include/laplace_kernel.h"
#include <deal.II/numerics/error_estimator.h>


#include <iostream>
#include <iomanip>

#include "Teuchos_TimeMonitor.hpp"

using Teuchos::Time;
using Teuchos::TimeMonitor;
using Teuchos::RCP;

#define ENTRY EntryRaiiObject obj ## LINE (__FUNCTION__);

struct EntryRaiiObject
{
  EntryRaiiObject(const char *f) : f_( f )
  {
    printf("Entered into %s\n", f_);
  }
  ~EntryRaiiObject()
  {
    printf("Exited from %s\n", f_);
  }
  const char *f_;
};

RCP<Time> ConstraintsTime = Teuchos::TimeMonitor::getNewTimer("Compute Constraints Time");
RCP<Time> AssembleTime = Teuchos::TimeMonitor::getNewTimer("Assemble Time");
RCP<Time> NormalsTime = Teuchos::TimeMonitor::getNewTimer("Normals Time");
RCP<Time> SurfaceGradientTime = Teuchos::TimeMonitor::getNewTimer("SurfaceGradientTime Time");
RCP<Time> GradientTime = Teuchos::TimeMonitor::getNewTimer("Gradient Time");
RCP<Time> LacSolveTime = Teuchos::TimeMonitor::getNewTimer("LAC Solve Time");
RCP<Time> ReinitTime = Teuchos::TimeMonitor::getNewTimer("BEM Reinitialisation Time");

// @sect4{BEMProblem::BEMProblem and
// BEMProblem::read_parameters}
// The constructor initializes the
// variuous object in much the same
// way as done in the finite element
// programs such as step-4 or
// step-6. The only new ingredient
// here is the ParsedFunction object,
// which needs, at construction time,
// the specification of the number of
// components.
//
// For the exact solution the number
// of vector components is one, and
// no action is required since one is
// the default value for a
// ParsedFunction object. The wind,
// however, requires dim components
// to be specified. Notice that when
// declaring entries in a parameter
// file for the expression of the
// Functions::ParsedFunction, we need
// to specify the number of
// components explicitly, since the
// function
// Functions::ParsedFunction::declare_parameters
// is static, and has no knowledge of
// the number of components.
template <>
BEMProblem<3>::BEMProblem(ComputationalDomain<3> &comp_dom,
                          // const unsigned int fe_degree,
                          MPI_Comm comm)
  :
  pcout(std::cout),
  comp_dom(comp_dom),
  parsed_fe("Scalar FE", "FE_Q(1)"),
  parsed_gradient_fe("Vector FE", "FESystem[FE_Q(1)^3]","u,u,u",3),
  dh(comp_dom.tria),
  gradient_dh(comp_dom.tria),
  mpi_communicator (comm),
  n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
  this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator))
{
  // Only output on first processor.
  pcout.set_condition(this_mpi_process == 0);
}
template <>
BEMProblem<2>::BEMProblem(ComputationalDomain<2> &comp_dom,
                          // const unsigned int fe_degree,
                          MPI_Comm comm)
  :
  pcout(std::cout),
  comp_dom(comp_dom),
  parsed_fe("Scalar FE", "FE_Q(1)"),
  parsed_gradient_fe("Vector FE", "FESystem[FE_Q(1)^2]","u,u",2),
  dh(comp_dom.tria),
  gradient_dh(comp_dom.tria),
  mpi_communicator (comm),
  n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
  this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator))
{
  // Only output on first processor.
  pcout.set_condition(this_mpi_process == 0);
}


template <int dim>
void BEMProblem<dim>::reinit()
{
  // ENTRY
  Teuchos::TimeMonitor LocalTimer(*ReinitTime);

  fe = parsed_fe();
  gradient_fe = parsed_gradient_fe();
  // fe = new FE_DGQArbitraryNodes<dim-1, dim>(QGauss<1> (2));
  // gradient_fe = new FESystem<dim-1,dim>(FE_DGQArbitraryNodes<dim-1,dim>(QGauss<1> (2)),dim);
  // // auto hhh = new FE_DGQArbitraryNodes<dim-1, dim>(QGauss<1> (2));
  std::string foo = fe->get_name();
  std::cout<<foo<<std::endl;
  // FiniteElement<dim-1,dim> * pippo = FETools::get_fe_by_name<dim-1, dim>(foo);
  // std::cout<<pippo->get_name()<<std::endl;

  dh.distribute_dofs(*fe);
  gradient_dh.distribute_dofs(*gradient_fe);

  // we should choose the appropriate renumbering strategy and then stick with it.
  // in step 32 they use component_wise which is very straight-forward but maybe the quickest
  // is subdomain_wise (step 17, 18)
  DoFRenumbering::component_wise (dh);
  DoFRenumbering::component_wise (gradient_dh);

  pcout<<"re-ordering vector"<<std::endl;

  compute_reordering_vectors();

  DoFRenumbering::subdomain_wise (dh);
  DoFRenumbering::subdomain_wise (gradient_dh);

  vector_constraints.reinit();
  DoFTools::make_hanging_node_constraints (gradient_dh,vector_constraints);
  vector_constraints.close();
  if (mapping_type == "FE")
    {
      map_vector.reinit(gradient_dh.n_dofs());
      // Fills the euler vector with information from the Triangulation
      VectorTools::get_position_vector(gradient_dh, map_vector);
      vector_constraints.distribute(map_vector);
    }
  // mapping_degree = fe->get_degree();
  if (!mapping)
    {
      if (mapping_type == "FE")
        mapping = SP(new MappingFEField<dim-1, dim> (gradient_dh, map_vector));
      else
        mapping = SP(new MappingQ<dim-1, dim> (mapping_degree));
    }



  const types::global_dof_index n_dofs =  dh.n_dofs();

  pcout<<dh.n_dofs()<<" "<<gradient_dh.n_dofs()<<std::endl;
  std::vector<types::subdomain_id> dofs_domain_association(n_dofs);

  DoFTools::get_subdomain_association   (dh,dofs_domain_association);
  std::vector<types::subdomain_id> vector_dofs_domain_association(gradient_dh.n_dofs());

  DoFTools::get_subdomain_association   (gradient_dh,vector_dofs_domain_association);


  this_cpu_set.clear();
  vector_this_cpu_set.clear();
  this_cpu_set.set_size(n_dofs);
  vector_this_cpu_set.set_size(gradient_dh.n_dofs());


  // We compute this two vector in order to use an eventual DoFRenumbering::subdomain_wise
  // At the time being we don't. We need to decide the better strategy.


  // We need to enforce consistency between the non-ghosted IndexSets.
  // To be changed accordingly with the DoFRenumbering strategy.
  pcout<<"you are using "<<sizeof(dh.n_dofs())<<" bytes indices"<<std::endl;
  pcout<<"setting cpu_sets"<<std::endl;

  for (types::global_dof_index i=0; i<n_dofs; ++i)
    if (dofs_domain_association[i] == this_mpi_process)
      {

        this_cpu_set.add_index(i);
        types::global_dof_index dummy=sub_wise_to_original[i];
        for (unsigned int idim=0; idim<dim; ++idim)
          {
            vector_this_cpu_set.add_index(vec_original_to_sub_wise[gradient_dh.n_dofs()/dim*idim+dummy]);
          }
      }


  // for (unsigned int i=0; i<gradient_dh.n_dofs(); ++i)
  //   if (vector_dofs_domain_association[i] == this_mpi_process)
  //     {
  //       vector_this_cpu_set.add_index(i);
  //       // for(unsigned int idim=0; idim<dim; ++idim)
  //       // {
  //       //   vector_this_cpu_set.add_index(i*dim+idim);
  //       // }
  //     }

  this_cpu_set.compress();
  vector_this_cpu_set.compress();
  // std::cout<<"set the cpu sets"<<std::endl;
  // std::vector<types::global_dof_index> localized_ndfos(n_mpi_processes);
  // std::vector<types::global_dof_index> localized_vector_ndfos(n_mpi_processes);
  // start_per_process.resize (n_mpi_processes);
  // vector_start_per_process.resize (n_mpi_processes);
  //
  // localized_ndfos[this_mpi_process] = this_cpu_set.n_elements();
  // localized_vector_ndfos[this_mpi_process] = vector_this_cpu_set.n_elements();
  //
  // Utilities::MPI::sum (localized_ndfos, mpi_communicator, start_per_process);
  // Utilities::MPI::sum (localized_vector_ndfos, mpi_communicator, vector_start_per_process);
  //
  // for(unsigned int i=start_per_process.size()-1; i>0; --i)
  // {
  //   start_per_process[i] = start_per_process[i-1];
  //   vector_start_per_process[i] = vector_start_per_process[i-1];
  // }
  // start_per_process[0] = 0;
  // vector_start_per_process[0] = 0;
  // for(unsigned int i=2; i<start_per_process.size(); ++i)
  // {
  //   start_per_process[i] += start_per_process[i-1];
  //   vector_start_per_process[i] += vector_start_per_process[i-1];
  // }
  // start_per_process[0] = 0;
  // vector_start_per_process[0] = 0;

  // std::cout<<this_mpi_process<<" "<<start_per_process[this_mpi_process]<<" "<<vector_start_per_process[this_mpi_process]<<std::endl;


  // At this point we just need to create a ghosted IndexSet for the scalar
  // DoFHandler. This can be through the builtin dealii function.
  // this_cpu_set.print(std::cout);
  MPI_Barrier(mpi_communicator);
  ghosted_set.clear();
  ghosted_set.set_size(dh.n_dofs());
  ghosted_set = DoFTools::dof_indices_with_subdomain_association(dh, this_mpi_process);
  ghosted_set.compress();
  // std::cout<<"set ghosted set"<<std::endl;

  // standard TrilinosWrappers::MPI::Vector reinitialization.
  system_rhs.reinit(this_cpu_set,mpi_communicator);
  sol.reinit(this_cpu_set,mpi_communicator);
  alpha.reinit(this_cpu_set,mpi_communicator);
  serv_phi.reinit(this_cpu_set,mpi_communicator);
  serv_dphi_dn.reinit(this_cpu_set,mpi_communicator);
  serv_tmp_rhs.reinit(this_cpu_set,mpi_communicator);


  // TrilinosWrappers::SparsityPattern for the BEM matricesreinitialization
  pcout<<"re-initializing sparsity patterns and matrices"<<std::endl;
  if (solution_method == "Direct")
    {
      full_sparsity_pattern.reinit(this_cpu_set, mpi_communicator);

      for (auto i : this_cpu_set)
        {
          for (types::global_dof_index j=0; j<dh.n_dofs(); ++j)
            full_sparsity_pattern.add(i,j);

        }

      full_sparsity_pattern.compress();
      neumann_matrix.reinit(full_sparsity_pattern);
      dirichlet_matrix.reinit(full_sparsity_pattern);
    }
  pcout<<"re-initialized sparsity patterns and matrices"<<std::endl;
  preconditioner_band = 100;
  preconditioner_sparsity_pattern.reinit(this_cpu_set, mpi_communicator, (types::global_dof_index) preconditioner_band);
  is_preconditioner_initialized = false;

  dirichlet_nodes.reinit(this_cpu_set,mpi_communicator);
  neumann_nodes.reinit(this_cpu_set,mpi_communicator);
  compute_dirichlet_and_neumann_dofs_vectors();
  compute_double_nodes_set();

  fma.init_fma(dh, double_nodes_set, dirichlet_nodes, *mapping, quadrature_order, singular_quadrature_order);



  // We need a TrilinosWrappers::MPI::Vector to reinit the SparsityPattern for
  // the parallel mass matrices.
  TrilinosWrappers::MPI::Vector helper(vector_this_cpu_set, mpi_communicator);
  //These are just for test
  // IndexSet vector_active_dofs;
  // IndexSet vector_relevant_dofs;
  IndexSet trial_index_set;
  // vector_active_dofs.clear();
  // vector_relevant_dofs.clear();
  trial_index_set.clear();
  // DoFTools::extract_locally_active_dofs(gradient_dh, vector_active_dofs);//, vector_active_dofs);
  trial_index_set = DoFTools::dof_indices_with_subdomain_association(gradient_dh, this_mpi_process);
  // Assert(trial_index_set == vector_this_cpu_set, ExcNotImplemented());
  // // The following functions returns the entire dof set.
  // DoFTools::extract_locally_relevant_dofs(gradient_dh, vector_relevant_dofs);
  // pcout<<vector_active_dofs.n_elements()<<"  "<<vector_relevant_dofs.n_elements()<<"   "<<vector_this_cpu_set.n_elements()<<std::endl;


  // This is the only way we could create the SparsityPattern, through the Epetramap of an
  // existing vector.
  // vector_sparsity_pattern.reinit(helper.vector_partitioner(), helper.vector_partitioner());
  vector_sparsity_pattern.reinit(vector_this_cpu_set, vector_this_cpu_set, mpi_communicator);
  DoFTools::make_sparsity_pattern (gradient_dh, vector_sparsity_pattern, vector_constraints, true, this_mpi_process);
  vector_sparsity_pattern.compress();



}


template<>
const Quadrature<2> &BEMProblem<3>::get_singular_quadrature(const unsigned int index) const
{
  Assert(index < fe->dofs_per_cell,
         ExcIndexRange(0, fe->dofs_per_cell, index));



  static std::vector<Quadrature<2> > quadratures;
  {
    if (quadratures.size() == 0)
      for (unsigned int i=0; i<fe->dofs_per_cell; ++i)
        {
          quadratures.push_back(QSplit<2> (QDuffy (singular_quadrature_order,1.),fe->get_unit_support_points()[i]));
        }
  }

  return quadratures[index];

}

template<>
const Quadrature<1> &BEMProblem<2>::get_singular_quadrature(const unsigned int index) const
{
  Assert(index < fe->dofs_per_cell,
         ExcIndexRange(0, fe->dofs_per_cell, index));

  static std::vector<Quadrature<1> > quadratures;
  if (quadratures.size() == 0)
    for (unsigned int i=0; i<fe->dofs_per_cell; ++i)
      {
        quadratures.push_back(QTelles<1>(singular_quadrature_order,
                                         fe->get_unit_support_points()[i]));
      }
  return quadratures[index];
}

template <int dim>
void BEMProblem<dim>::declare_parameters (ParameterHandler &prm)
{

  // In the solver section, we set
  // all SolverControl
  // parameters. The object will then
  // be fed to the GMRES solver in
  // the solve_system() function.

  prm.enter_subsection("Solver");
  SolverControl::declare_parameters(prm);
  prm.leave_subsection();

  prm.declare_entry("Preconditioner","ILU",
                    Patterns::Selection("ILU|AMG"));

  prm.declare_entry("Solution method", "Direct",
                    Patterns::Selection("Direct|FMA"));

  prm.enter_subsection("Quadrature rules");
  {
    prm.declare_entry("Quadrature type", "gauss",
                      Patterns::Selection(QuadratureSelector<(dim-1)>::get_quadrature_names()));
    prm.declare_entry("Quadrature order", "4", Patterns::Integer());
    prm.declare_entry("Singular quadrature order", "5", Patterns::Integer());
  }
  prm.leave_subsection();

  prm.declare_entry("Mapping Type","FE",
                    Patterns::Selection("FE|Q"));

  prm.declare_entry("Mapping Q Degree","1",Patterns::Integer());

  prm.declare_entry("Continuos gradient across edges","true",Patterns::Bool());

}

template <int dim>
void BEMProblem<dim>::parse_parameters (ParameterHandler &prm)
{

  prm.enter_subsection("Solver");
  solver_control.parse_parameters(prm);
  prm.leave_subsection();

  preconditioner_type = prm.get("Preconditioner");

  solution_method = prm.get("Solution method");


  prm.enter_subsection("Quadrature rules");
  {
    quadrature =
      std_cxx1x::shared_ptr<Quadrature<dim-1> >
      (new QuadratureSelector<dim-1> (prm.get("Quadrature type"),
                                      prm.get_integer("Quadrature order")));
    quadrature_order = prm.get_integer("Quadrature order");
    singular_quadrature_order = prm.get_integer("Singular quadrature order");
  }
  prm.leave_subsection();

  mapping_type = prm.get("Mapping Type");
  mapping_degree = prm.get_integer("Mapping Q Degree");
  continuos_gradient = prm.get_bool("Continuos gradient across edges");


}


template <int dim>
void BEMProblem<dim>::compute_dirichlet_and_neumann_dofs_vectors()
{
  have_dirichlet_bc = false;


  Vector<double> non_partitioned_dirichlet_nodes(dh.n_dofs());
  Vector<double> non_partitioned_neumann_nodes(dh.n_dofs());



  cell_it
  cell = dh.begin_active(),
  endc = dh.end();


  vector_shift(non_partitioned_neumann_nodes, 1.);
  std::vector<types::global_dof_index> dofs(fe->dofs_per_cell);
  std::vector<types::global_dof_index> gradient_dofs(gradient_fe->dofs_per_cell);
  unsigned int helper_dirichlet=0;
  for (; cell != endc; ++cell)
    {
      if (cell->subdomain_id() == this_mpi_process)
        {
          bool dirichlet = false;
          for (auto dummy : comp_dom.dirichlet_boundary_ids)
            {
              if (dummy == cell->material_id())
                {
                  cell->get_dof_indices(dofs);
                  for (unsigned int i=0; i<fe->dofs_per_cell; ++i)
                    {
                      non_partitioned_dirichlet_nodes(dofs[i]) = 1;
                      non_partitioned_neumann_nodes(dofs[i]) = 0;
                      //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"  surfNodes: "<<dirichlet_nodes(dofs[i])<<"  otherNodes: "<<neumann_nodes(dofs[i])<<std::endl;
                    }
                  dirichlet = true;
                  helper_dirichlet = 1.;
                  break;
                }
            }
          if (!dirichlet)
            {
              cell->get_dof_indices(dofs);
              // for(unsigned int i=0; i<fe->dofs_per_cell; ++i)
              // {
              //   non_partitioned_neumann_nodes(dofs[i]) = 1;
              //   non_partitioned_dirichlet_nodes(dofs[i]) = 0;
              // }

            }

          // if (cell->material_id() == comp_dom.dirichlet_sur_ID1 ||
          //     cell->material_id() == comp_dom.dirichlet_sur_ID2 ||
          //     cell->material_id() == comp_dom.dirichlet_sur_ID3)
          //   {
          //     // This is a free surface node.
          //     cell->get_dof_indices(dofs);
          //     for (unsigned int i=0; i<fe->dofs_per_cell; ++i)
          //       {
          //         non_partitioned_dirichlet_nodes(dofs[i]) = 1;
          //         non_partitioned_neumann_nodes(dofs[i]) = 0;
          //         //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"  surfNodes: "<<dirichlet_nodes(dofs[i])<<"  otherNodes: "<<neumann_nodes(dofs[i])<<std::endl;
          //       }
          //   }
          // else
          //   {
          //     for (unsigned int i=0; i<fe->dofs_per_cell; ++i)
          //       {
          //         cell->get_dof_indices(dofs);
          //         //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"  surfNodes: "<<dirichlet_nodes(dofs[i])<<"  otherNodes: "<<neumann_nodes(dofs[i])<<std::endl;
          //       }
          //
          //   }
        }


    }

  for (types::global_dof_index i=0; i<dh.n_dofs(); ++i)
    if (this_cpu_set.is_element(i))
      {
        dirichlet_nodes(i)=non_partitioned_dirichlet_nodes(i);
        neumann_nodes(i)=non_partitioned_neumann_nodes(i);
      }
  // dirichlet_nodes.add(non_partitioned_dirichlet_nodes, true);// = non_partitioned_dirichlet_nodes;
  // neumann_nodes.add(non_partitioned_neumann_nodes, true);// = non_partitioned_neumann_nodes;
  unsigned int helper_dirichlet_2;
  // std::cout<<this_mpi_process<<" , "<<helper_dirichlet<<std::endl;
  MPI_Allreduce(&helper_dirichlet,&helper_dirichlet_2,1,MPI_UNSIGNED,MPI_MAX,mpi_communicator);
  // std::cout<<this_mpi_process<<" , "<<helper_dirichlet<<" , "<<helper_dirichlet_2<<std::endl;
  if (helper_dirichlet_2>0)
    have_dirichlet_bc=true;
  // std::cout<<this_mpi_process<<" , "<<have_dirichlet_bc<<std::endl;
  //for (unsigned int i=0; i<dh.n_dofs(); ++i)
  //    if (this_mpi_process == 1)
  //       pcout<<i<<" "<<dirichlet_nodes(i)<<" "<<neumann_nodes(i)<<std::endl;

}

template <int dim>
void BEMProblem<dim>::compute_double_nodes_set()
{
  double tol=1e-10;
  double_nodes_set.clear();
  double_nodes_set.resize(dh.n_dofs());
  std::vector<Point<dim> > support_points(dh.n_dofs());

  DoFTools::map_dofs_to_support_points<dim-1, dim>( *mapping,
                                                    dh, support_points);

  typename DoFHandler<dim-1,dim>::active_cell_iterator
  cell = dh.begin_active(),
  endc = dh.end();
  std::vector<types::global_dof_index> face_dofs(fe->dofs_per_face);

  edge_set.clear();
  edge_set.set_size(dh.n_dofs());

  for (cell=dh.begin_active(); cell!=endc; ++cell)
    {
      for (unsigned int f=0; f<GeometryInfo<dim-1>::faces_per_cell; ++f)
        if ( cell->face(f)->at_boundary() )
          {
            cell->face(f)->get_dof_indices(face_dofs);
            for (unsigned int k=0; k<face_dofs.size(); ++k)
              edge_set.add_index(face_dofs[k]);
          }
    }
  edge_set.compress();

  for (types::global_dof_index i=0; i<dh.n_dofs(); ++i)
    double_nodes_set[i].insert(i);
  for (auto i : edge_set)//(types::global_dof_index i=0; i<dh.n_dofs(); ++i)
    {
      for (auto j : edge_set)
        {
          if (support_points[i].distance(support_points[j]) < tol)
            {
              double_nodes_set[i].insert(j);
            }
        }

    }



}

template <int dim>
void BEMProblem<dim>::compute_reordering_vectors()
{
  original_to_sub_wise.resize(dh.n_dofs());
  sub_wise_to_original.resize(dh.n_dofs());
  vec_original_to_sub_wise.resize(gradient_dh.n_dofs());
  vec_sub_wise_to_original.resize(gradient_dh.n_dofs());

  DoFRenumbering::compute_subdomain_wise(original_to_sub_wise, dh);
  DoFRenumbering::compute_subdomain_wise(vec_original_to_sub_wise, gradient_dh);

  for (types::global_dof_index i=0; i<gradient_dh.n_dofs(); ++i)
    {
      if (i<dh.n_dofs())
        {
          sub_wise_to_original[original_to_sub_wise[i]]=i;
        }
      vec_sub_wise_to_original[vec_original_to_sub_wise[i]]=i;
    }

}
template <int dim>
void BEMProblem<dim>::assemble_system()
{
  Teuchos::TimeMonitor LocalTimer(*AssembleTime);
  pcout<<"(Directly) Assembling system matrices"<<std::endl;

  neumann_matrix = 0;
  dirichlet_matrix = 0;



  // Next, we initialize an FEValues
  // object with the quadrature
  // formula for the integration of
  // the kernel in non singular
  // cells. This quadrature is
  // selected with the parameter
  // file, and needs to be quite
  // precise, since the functions we
  // are integrating are not
  // polynomial functions.
  FEValues<dim-1,dim> fe_v(*mapping,*fe, *quadrature,
                           update_values |
                           update_cell_normal_vectors |
                           update_quadrature_points |
                           update_JxW_values);

  const unsigned int n_q_points = fe_v.n_quadrature_points;

  std::vector<types::global_dof_index> local_dof_indices(fe->dofs_per_cell);
  pcout<<fe->dofs_per_cell<<" "<<std::endl;
  // Unlike in finite element
  // methods, if we use a collocation
  // boundary element method, then in
  // each assembly loop we only
  // assemble the information that
  // refers to the coupling between
  // one degree of freedom (the
  // degree associated with support
  // point $i$) and the current
  // cell. This is done using a
  // vector of fe->dofs_per_cell
  // elements, which will then be
  // distributed to the matrix in the
  // global row $i$. The following
  // object will hold this
  // information:
  Vector<double>      local_neumann_matrix_row_i(fe->dofs_per_cell);
  Vector<double>      local_dirichlet_matrix_row_i(fe->dofs_per_cell);

  // Now that we have checked that
  // the number of vertices is equal
  // to the number of degrees of
  // freedom, we construct a vector
  // of support points which will be
  // used in the local integrations:
  std::vector<Point<dim> > support_points(dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *mapping, dh, support_points);


  // After doing so, we can start the
  // integration loop over all cells,
  // where we first initialize the
  // FEValues object and get the
  // values of $\mathbf{\tilde v}$ at
  // the quadrature points (this
  // vector field should be constant,
  // but it doesn't hurt to be more
  // general):


  cell_it
  cell = dh.begin_active(),
  endc = dh.end();

  Point<dim> D;
  double s;

  for (cell = dh.begin_active(); cell != endc; ++cell)
    {
      fe_v.reinit(cell);
      cell->get_dof_indices(local_dof_indices);

      const std::vector<Point<dim> > &q_points = fe_v.get_quadrature_points();
      const std::vector<Tensor<1, dim> > &normals = fe_v.get_normal_vectors();

      // We then form the integral over
      // the current cell for all
      // degrees of freedom (note that
      // this includes degrees of
      // freedom not located on the
      // current cell, a deviation from
      // the usual finite element
      // integrals). The integral that
      // we need to perform is singular
      // if one of the local degrees of
      // freedom is the same as the
      // support point $i$. A the
      // beginning of the loop we
      // therefore check wether this is
      // the case, and we store which
      // one is the singular index:
      for (types::global_dof_index i=0; i<dh.n_dofs() ; ++i) //these must now be the locally owned dofs. the rest should stay the same
        {
          if (this_cpu_set.is_element(i))
            {
              local_neumann_matrix_row_i = 0;
              local_dirichlet_matrix_row_i = 0;

              bool is_singular = false;
              unsigned int singular_index = numbers::invalid_unsigned_int;

              for (unsigned int j=0; j<fe->dofs_per_cell; ++j)
                //if(local_dof_indices[j] == i)
                if (double_nodes_set[i].count(local_dof_indices[j]) > 0)
                  {
                    singular_index = j;
                    is_singular = true;
                    break;
                  }

              // We then perform the
              // integral. If the index $i$
              // is not one of the local
              // degrees of freedom, we
              // simply have to add the
              // single layer terms to the
              // right hand side, and the
              // double layer terms to the
              // matrix:
              if (is_singular == false)
                {
                  for (unsigned int q=0; q<n_q_points; ++q)
                    {
                      const Tensor<1,dim> R = q_points[q] - support_points[i];
                      LaplaceKernel::kernels(R, D, s);
                      // if(support_points[i][0]==0.25&&support_points[i][1]==0.25)
                      //   pcout<<"D "<<D<<" s "<<s<<" , ";
                      for (unsigned int j=0; j<fe->dofs_per_cell; ++j)
                        {
                          local_neumann_matrix_row_i(j) += ( ( D *
                                                               normals[q] ) *
                                                             fe_v.shape_value(j,q) *
                                                             fe_v.JxW(q)       );
                          local_dirichlet_matrix_row_i(j) += ( s *
                                                               fe_v.shape_value(j,q) *
                                                               fe_v.JxW(q) );

                        }
                    }
                }
              else
                {
                  // Now we treat the more
                  // delicate case. If we
                  // are here, this means
                  // that the cell that
                  // runs on the $j$ index
                  // contains
                  // support_point[i]. In
                  // this case both the
                  // single and the double
                  // layer potential are
                  // singular, and they
                  // require special
                  // treatment.
                  //
                  // Whenever the
                  // integration is
                  // performed with the
                  // singularity inside the
                  // given cell, then a
                  // special quadrature
                  // formula is used that
                  // allows one to
                  // integrate arbitrary
                  // functions against a
                  // singular weight on the
                  // reference cell.
                  // Notice that singular
                  // integration requires a
                  // careful selection of
                  // the quadrature
                  // rules. In particular
                  // the deal.II library
                  // provides quadrature
                  // rules which are
                  // taylored for
                  // logarithmic
                  // singularities
                  // (QGaussLog,
                  // QGaussLogR), as well
                  // as for 1/R
                  // singularities
                  // (QGaussOneOverR).
                  //
                  // Singular integration
                  // is typically obtained
                  // by constructing
                  // weighted quadrature
                  // formulas with singular
                  // weights, so that it is
                  // possible to write
                  //
                  // \f[
                  //   \int_K f(x) s(x) dx = \sum_{i=1}^N w_i f(q_i)
                  // \f]
                  //
                  // where $s(x)$ is a given
                  // singularity, and the weights
                  // and quadrature points
                  // $w_i,q_i$ are carefully
                  // selected to make the formula
                  // above an equality for a
                  // certain class of functions
                  // $f(x)$.
                  //
                  // In all the finite
                  // element examples we
                  // have seen so far, the
                  // weight of the
                  // quadrature itself
                  // (namely, the function
                  // $s(x)$), was always
                  // constantly equal to 1.
                  // For singular
                  // integration, we have
                  // two choices: we can
                  // use the definition
                  // above, factoring out
                  // the singularity from
                  // the integrand (i.e.,
                  // integrating $f(x)$
                  // with the special
                  // quadrature rule), or
                  // we can ask the
                  // quadrature rule to
                  // "normalize" the
                  // weights $w_i$ with
                  // $s(q_i)$:
                  //
                  // \f[
                  //   \int_K f(x) s(x) dx =
                  //   \int_K g(x) dx = \sum_{i=1}^N \frac{w_i}{s(q_i)} g(q_i)
                  // \f]
                  //
                  // We use this second
                  // option, through the @p
                  // factor_out_singularity
                  // parameter of both
                  // QGaussLogR and
                  // QGaussOneOverR.
                  //
                  // These integrals are
                  // somewhat delicate,
                  // especially in two
                  // dimensions, due to the
                  // transformation from
                  // the real to the
                  // reference cell, where
                  // the variable of
                  // integration is scaled
                  // with the determinant
                  // of the transformation.
                  //
                  // In two dimensions this
                  // process does not
                  // result only in a
                  // factor appearing as a
                  // constant factor on the
                  // entire integral, but
                  // also on an additional
                  // integral alltogether
                  // that needs to be
                  // evaluated:
                  //
                  // \f[
                  //  \int_0^1 f(x)\ln(x/\alpha) dx =
                  //  \int_0^1 f(x)\ln(x) dx - \int_0^1 f(x) \ln(\alpha) dx.
                  // \f]
                  //
                  // This process is taken care of by
                  // the constructor of the QGaussLogR
                  // class, which adds additional
                  // quadrature points and weights to
                  // take into consideration also the
                  // second part of the integral.
                  //
                  // A similar reasoning
                  // should be done in the
                  // three dimensional
                  // case, since the
                  // singular quadrature is
                  // taylored on the
                  // inverse of the radius
                  // $r$ in the reference
                  // cell, while our
                  // singular function
                  // lives in real space,
                  // however in the three
                  // dimensional case
                  // everything is simpler
                  // because the
                  // singularity scales
                  // linearly with the
                  // determinant of the
                  // transformation. This
                  // allows us to build the
                  // singular two
                  // dimensional quadrature
                  // rules once and for all
                  // outside the loop over
                  // all cells, using only
                  // a pointer where needed.
                  //
                  // Notice that in one
                  // dimensional
                  // integration this is
                  // not possible, since we
                  // need to know the
                  // scaling parameter for
                  // the quadrature, which
                  // is not known a
                  // priori. Here, the
                  // quadrature rule itself
                  // depends also on the
                  // size of the current
                  // cell. For this reason,
                  // it is necessary to
                  // create a new
                  // quadrature for each
                  // singular
                  // integration. Since we
                  // create it using the
                  // new operator of C++,
                  // we also need to
                  // destroy it using the
                  // dual of new:
                  // delete. This is done
                  // at the end, and only
                  // if dim == 2.
                  //
                  // Putting all this into a
                  // dimension independent
                  // framework requires a little
                  // trick. The problem is that,
                  // depending on dimension, we'd
                  // like to either assign a
                  // QGaussLogR<1> or a
                  // QGaussOneOverR<2> to a
                  // Quadrature<dim-1>. C++
                  // doesn't allow this right
                  // away, and neither is a
                  // static_cast
                  // possible. However, we can
                  // attempt a dynamic_cast: the
                  // implementation will then
                  // look up at run time whether
                  // the conversion is possible
                  // (which we <em>know</em> it
                  // is) and if that isn't the
                  // case simply return a null
                  // pointer. To be sure we can
                  // then add a safety check at
                  // the end:
                  Assert(singular_index != numbers::invalid_unsigned_int,
                         ExcInternalError());

                  const Quadrature<dim-1> *
                  singular_quadrature
                    = &(get_singular_quadrature(singular_index));
                  Assert(singular_quadrature, ExcInternalError());

                  FEValues<dim-1,dim> fe_v_singular (*mapping, *fe, *singular_quadrature,
                                                     update_jacobians |
                                                     update_values |
                                                     update_cell_normal_vectors |
                                                     update_quadrature_points );

                  fe_v_singular.reinit(cell);

                  const std::vector<Tensor<1, dim> > &singular_normals = fe_v_singular.get_normal_vectors();
                  const std::vector<Point<dim> > &singular_q_points = fe_v_singular.get_quadrature_points();

                  for (unsigned int q=0; q<singular_quadrature->size(); ++q)
                    {
                      const Tensor<1,dim> R = singular_q_points[q] - support_points[i];
                      LaplaceKernel::kernels(R, D, s);

                      for (unsigned int j=0; j<fe->dofs_per_cell; ++j)
                        {
                          local_neumann_matrix_row_i(j) += (( D *
                                                              singular_normals[q])                *
                                                            fe_v_singular.shape_value(j,q)        *
                                                            fe_v_singular.JxW(q)       );

                          local_dirichlet_matrix_row_i(j) += ( s   *
                                                               fe_v_singular.shape_value(j,q)     *
                                                               fe_v_singular.JxW(q) );
                        }
                    }
                }

              // Finally, we need to add
              // the contributions of the
              // current cell to the
              // global matrix.
              for (unsigned int j=0; j<fe->dofs_per_cell; ++j)
                {
                  neumann_matrix.add(i,local_dof_indices[j],local_neumann_matrix_row_i(j));
                  dirichlet_matrix.add(i,local_dof_indices[j],local_dirichlet_matrix_row_i(j));
                }
            }
        }
    }

  // The second part of the integral
  // operator is the term
  // $\alpha(\mathbf{x}_i)
  // \phi_j(\mathbf{x}_i)$. Since we
  // use a collocation scheme,
  // $\phi_j(\mathbf{x}_i)=\delta_{ij}$
  // and the corresponding matrix is
  // a diagonal one with entries
  // equal to $\alpha(\mathbf{x}_i)$.

  // One quick way to compute this
  // diagonal matrix of the solid
  // angles, is to use the Neumann
  // matrix itself. It is enough to
  // multiply the matrix with a
  // vector of elements all equal to
  // -1, to get the diagonal matrix
  // of the alpha angles, or solid
  // angles (see the formula in the
  // introduction for this). The
  // result is then added back onto
  // the system matrix object to
  // yield the final form of the
  // matrix:

  /*
    pcout<<"Neumann"<<std::endl;
    for (unsigned int i = 0; i < dh.n_dofs(); i++)
        {
        if (this_cpu_set.is_element(i))
           {
           pcout<<this_mpi_process<<" *** ";
           for (unsigned int j = 0; j < dh.n_dofs(); j++)
               {
               pcout<<neumann_matrix(i,j)<<" ";
               }
           pcout<<std::endl;
           }
        }



    pcout<<"Dirichlet"<<std::endl;
    for (unsigned int i = 0; i < dh.n_dofs(); i++)
        {
        if (this_cpu_set.is_element(i))
           {
           pcout<<this_mpi_process<<" *** ";
           for (unsigned int j = 0; j < dh.n_dofs(); j++)
               {
               pcout<<dirichlet_matrix(i,j)<<" ";
               }
           pcout<<std::endl;
           }
        }
        //*/
  pcout<<"done assembling system matrices"<<std::endl;
  // std::cout<<"printing Neumann Matrix"<<std::endl;
  // for(unsigned int i=0; i<dh.n_dofs(); ++i)
  // {
  //   for(unsigned int j=0; j<dh.n_dofs(); ++j)
  //     std::cout<<neumann_matrix(i,j)<<" ";
  //   std::cout<<std::endl;
  // }
  // std::cout<<"printing Dirichlet Matrix"<<std::endl;
  // for(unsigned int i=0; i<dh.n_dofs(); ++i)
  // {
  //   for(unsigned int j=0; j<dh.n_dofs(); ++j)
  //     std::cout<<dirichlet_matrix(i,j)<<" ";
  //   std::cout<<std::endl;
  // }
}



template <int dim>
void BEMProblem<dim>::compute_alpha()
{
  static TrilinosWrappers::MPI::Vector ones, zeros, dummy;
  if (ones.size() != dh.n_dofs())
    {
      ones.reinit(this_cpu_set,mpi_communicator);
      vector_shift(ones, -1.);
      zeros.reinit(this_cpu_set,mpi_communicator);
      dummy.reinit(this_cpu_set,mpi_communicator);
    }


  if (solution_method == "Direct")
    {
      neumann_matrix.vmult(alpha, ones);
    }
  else
    {
      AssertThrow(dim == 3,
                  ExcMessage("FMA only works in 3D"));

      fma.generate_multipole_expansions(ones,zeros);
      fma.multipole_matr_vect_products(ones,zeros,alpha,dummy);
    }

  //alpha.print(pcout);
  //for (unsigned int i=0; i<alpha.size(); ++i)
  //    {
  //    cout<<std::setprecision(20)<<alpha(i)<<endl;
  //    }


}

template <int dim>
void BEMProblem<dim>::vmult(TrilinosWrappers::MPI::Vector &dst, const TrilinosWrappers::MPI::Vector &src) const
{
  serv_phi = src;
  if (!have_dirichlet_bc)
    {
      vector_shift(serv_phi, -serv_phi.l2_norm());
    }
  serv_dphi_dn = src;



  TrilinosWrappers::MPI::Vector matrVectProdN;
  TrilinosWrappers::MPI::Vector matrVectProdD;

  matrVectProdN.reinit(this_cpu_set,mpi_communicator);
  matrVectProdD.reinit(this_cpu_set,mpi_communicator);

  dst = 0;


  serv_phi.scale(neumann_nodes);
  serv_dphi_dn.scale(dirichlet_nodes);

  if (solution_method == "Direct")
    {
      dirichlet_matrix.vmult(dst, serv_dphi_dn);
      dst *= -1;
      neumann_matrix.vmult_add(dst, serv_phi);
      serv_phi.scale(alpha);
      dst += serv_phi;

    }
  else
    {
      AssertThrow(dim == 3,
                  ExcMessage("FMA only works in 3D"));

      fma.generate_multipole_expansions(serv_phi,serv_dphi_dn);
      fma.multipole_matr_vect_products(serv_phi,serv_dphi_dn,matrVectProdN,matrVectProdD);
      serv_phi.scale(alpha);
      dst += matrVectProdD;
      dst *= -1;
      dst += matrVectProdN;
      dst += serv_phi;
    }

  //std::cout<<"*** "<<serv_phi(0)<<" or "<<serv_dphi_dn(0)<<"   src: "<<src(0)<<"  dst: "<<dst(0)<<std::endl;
  // in fully neumann bc case, we have to rescale the vector to have a zero mean
  // one
  if (!have_dirichlet_bc)
    vector_shift(dst, -dst.l2_norm());
}


template <int dim>
void BEMProblem<dim>::compute_rhs(TrilinosWrappers::MPI::Vector &dst, const TrilinosWrappers::MPI::Vector &src) const
{

  serv_phi = src;
  serv_dphi_dn = src;

  static TrilinosWrappers::MPI::Vector matrVectProdN;
  static TrilinosWrappers::MPI::Vector matrVectProdD;


  matrVectProdN.reinit(this_cpu_set,mpi_communicator);
  matrVectProdD.reinit(this_cpu_set,mpi_communicator);


  serv_phi.scale(dirichlet_nodes);
  serv_dphi_dn.scale(neumann_nodes);

  if (solution_method == "Direct")
    {
      neumann_matrix.vmult(dst, serv_phi);
      serv_phi.scale(alpha);
      dst += serv_phi;
      dst *= -1;
      dirichlet_matrix.vmult_add(dst, serv_dphi_dn);
    }
  else
    {
      AssertThrow(dim == 3,
                  ExcMessage("FMA only works in 3D"));

      fma.generate_multipole_expansions(serv_phi,serv_dphi_dn);
      fma.multipole_matr_vect_products(serv_phi,serv_dphi_dn,matrVectProdN,matrVectProdD);
      serv_phi.scale(alpha);
      dst += matrVectProdN;
      dst += serv_phi;
      dst *= -1;
      dst += matrVectProdD;
    }




}





// @sect4{BEMProblem::solve_system}

// The next function simply solves
// the linear system.
template <int dim>
void BEMProblem<dim>::solve_system(TrilinosWrappers::MPI::Vector &phi, TrilinosWrappers::MPI::Vector &dphi_dn,
                                   const TrilinosWrappers::MPI::Vector &tmp_rhs)
{
  Teuchos::TimeMonitor LocalTimer(*LacSolveTime);
  SolverGMRES<TrilinosWrappers::MPI::Vector > solver (solver_control,
                                                      SolverGMRES<TrilinosWrappers::MPI::Vector >::AdditionalData(100));

  system_rhs = 0;
  sol = 0;
  alpha = 0;


  compute_alpha();

  //for (unsigned int i = 0; i < alpha.size(); i++)
  //    if (this_cpu_set.is_element(i))
  //       pcout<<std::setprecision(20)<<alpha(i)<<std::endl;




  compute_rhs(system_rhs, tmp_rhs);


  compute_constraints(constr_cpu_set, constraints, tmp_rhs);
  ConstrainedOperator<TrilinosWrappers::MPI::Vector, BEMProblem<dim> >
  cc(*this, constraints,constr_cpu_set,mpi_communicator);


  cc.distribute_rhs(system_rhs);
  system_rhs.compress(VectorOperation::insert);
  // vmult(sol,system_rhs);
  // Assert(sol.vector_partitioner().SameAs(system_rhs.vector_partitioner()),ExcMessage("Schizofrenia???"));
  // cc.vmult(sol,system_rhs);
  // Assert(sol.locally_owned_elements()==system_rhs.locally_owned_elements(),ExcMessage("IndexSet a muzzo..."));
  // Assert(sol.vector_partitioner().SameAs(system_rhs.vector_partitioner()),ExcMessage("Ma boh..."));


  if (solution_method == "Direct")
    {
      //SparseDirectUMFPACK &inv = fma.FMA_preconditioner(alpha);
      //solver.solve (*this, sol, system_rhs, inv);
      assemble_preconditioner();
      //solver.solve (cc, sol, system_rhs, PreconditionIdentity());
      sol.sadd(1., 0., system_rhs);
      solver.solve (cc, sol, system_rhs, preconditioner);
    }
  else
    {
      AssertThrow(dim == 3,
                  ExcMessage("FMA only works in 3D"));

      TrilinosWrappers::PreconditionILU &fma_preconditioner = fma.FMA_preconditioner(alpha,constraints);
      solver.solve (cc, sol, system_rhs, fma_preconditioner);
      // solver.solve (cc, sol, system_rhs, PreconditionIdentity());
    }

  // cc.apply_constraint(sol);
  //pcout<<"sol = [";
  //for (unsigned int i = 0; i < dh.n_dofs(); i++)
  //    pcout<<sol(i)<<"; ";
  //pcout<<"];"<<std::endl;

  // for (unsigned int i = 0; i < sol.size(); i++)
  //   if (this_cpu_set.is_element(i))
  //      pcout<<std::setprecision(20)<<sol(i)<<std::endl;




///////////////////////////////////
  /*
    std::vector<Point<dim> > support_points(dh.n_dofs());
    DoFTools::map_dofs_to_support_points<dim-1, dim>( mapping, dh, support_points);
      pcout<<"**solution "<<std::endl;
     for (unsigned int i = 0; i < alpha.size(); i++)
         if (this_cpu_set.is_element(i))
            pcout<<i<<" ("<<this_mpi_process<<")  "<<support_points[i](0)+support_points[i](1)+support_points[i](2)<<"   "<<sol(i)<<std::endl;

     pcout<<"SOLUTION "<<std::endl;
     for (unsigned int i = 0; i < alpha.size(); i++)
         if (this_cpu_set.is_element(i))
            pcout<<i<<" ("<<this_mpi_process<<")  "<<sol(i)<<std::endl;
  */
//////////////////////////////////

  for (types::global_dof_index i=0; i <dirichlet_nodes.size(); i++)
    {
      if (this_cpu_set.is_element(i))
        {
          if (dirichlet_nodes(i) == 0)
            {
              phi(i) = sol(i);
            }
          else
            {
              dphi_dn(i) = sol(i);
            }
        }
    }
  phi(this_cpu_set.nth_index_in_set(0))=phi(this_cpu_set.nth_index_in_set(0));
  dphi_dn(this_cpu_set.nth_index_in_set(0))=dphi_dn(this_cpu_set.nth_index_in_set(0));
  phi.compress(VectorOperation::insert);
  dphi_dn.compress(VectorOperation::insert);

  //if (!have_dirichlet_bc)
  //   vector_shift(phi,-phi.l2_norm());

  //for (unsigned int i=0;i<dh.n_dofs();++i)
  // std::cout<<i<<" "<<tmp_rhs(i)<<" "<<dphi_dn(i)<<" "<<phi(i)<<" "<<dirichlet_nodes(i)<<std::endl;

  //pcout<<"sol "<<std::endl;
  //for (unsigned int i = 0; i < sol.size(); i++)
  //    {
  //pcout<<i<<" "<<sol(i)<<" ";
  //std::set<unsigned int> doubles = double_nodes_set[i];
  // for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
  //    pcout<<*it<<"("<<dirichlet_nodes(*it)<<") ";
  //pcout<<"phi "<<phi(i)<<"  dphi_dn "<<dphi_dn(i);
  //pcout<<std::endl;
  //    }
}



// This method performs a Bem resolution,
// either in a direct or multipole method
template <int dim>
void BEMProblem<dim>::solve(TrilinosWrappers::MPI::Vector &phi, TrilinosWrappers::MPI::Vector &dphi_dn,
                            const TrilinosWrappers::MPI::Vector &tmp_rhs)
{

  if (solution_method == "Direct")
    {
      assemble_system();
      // neumann_matrix.print(std::cout);
      // dirichlet_matrix.print(std::cout);
    }
  else
    {
      AssertThrow(dim == 3,
                  ExcMessage("FMA only works in 3D"));

      fma.generate_octree_blocking();
      // fma.compute_m2l_flags();
      fma.direct_integrals();
      fma.multipole_integrals();
    }

  solve_system(phi,dphi_dn,tmp_rhs);
}


template <int dim>
void BEMProblem<dim>::compute_constraints(IndexSet &c_cpu_set, ConstraintMatrix &c, const TrilinosWrappers::MPI::Vector &tmp_rhs)

{
  Teuchos::TimeMonitor LocalTimer(*ConstraintsTime);
  // We need both the normal vector and surface gradients to apply correctly dirichlet-dirichlet
  // double node constraints.
  // compute_normals();
  compute_surface_gradients(tmp_rhs);

  // communication is needed here: there is one matrix per process: thus the vector needed to set
  // inhomogeneities has to be copied locally
  Vector<double> localized_surface_gradients(vector_surface_gradients_solution);
  Vector<double> localized_normals(vector_normals_solution);
  Vector<double> localized_dirichlet_nodes(dirichlet_nodes);
  Vector<double> loc_tmp_rhs(tmp_rhs.size());
  loc_tmp_rhs = tmp_rhs;

  // we start clearing the constraint matrix
  c.clear();

  // here we prepare the constraint matrix so as to account for the presence hanging
  // nodes

  ConstraintMatrix c_hn;
  DoFTools::make_hanging_node_constraints (dh,c_hn);
  c_hn.close();

  std::vector<types::subdomain_id> dofs_domain_association(dh.n_dofs());

  DoFTools::get_subdomain_association   (dh,dofs_domain_association);
  // here we prepare the constraint matrix so as to account for the presence of double and
  // triple dofs

  // we start looping on the dofs
  for (types::global_dof_index i=0; i <tmp_rhs.size(); i++)
    {
      // if (this_cpu_set.is_element(i))
      // {
      // in the next line we compute the "first" among the set of double nodes: this node
      // is the first dirichlet node in the set, and if no dirichlet node is there, we get the
      // first neumann node

      std::set<types::global_dof_index> doubles = double_nodes_set[i];
      types::global_dof_index firstOfDoubles = *doubles.begin();
      for (std::set<types::global_dof_index>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
        {
          // if(this_cpu_set.is_element(*it))
          if (localized_dirichlet_nodes(*it) == 1)
            {
              firstOfDoubles = *it;
              break;
            }
        }

      // for each set of double nodes, we will perform the correction only once, and precisely
      // when the current node is the first of the set
      if (i == firstOfDoubles)
        {
          // the vector entry corresponding to the first node of the set does not need modification,
          // thus we erase ti form the set
          doubles.erase(i);

          // if the current (first) node is a dirichlet node, for all its neumann doubles we will
          // impose that the potential is equal to that of the first node: this means that in the
          // matrix vector product we will put the potential value of the double node
          if (localized_dirichlet_nodes(i) == 1)
            {
              for (std::set<types::global_dof_index>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
                {
                  // if(this_cpu_set.is_element(*it))
                  {
                    if (localized_dirichlet_nodes(*it) == 1)
                      {
                        // this is the dirichlet-dirichlet case on flat edges: here we impose that
                        // dphi_dn on the two (or more) sides is equal.
                        double normal_distance = 0;

                        // types::global_dof_index owner_el_1 = DoFTools::count_dofs_with_subdomain_association (dh, dofs_domain_association[i]);
                        // types::global_dof_index owner_el_2 = DoFTools::count_dofs_with_subdomain_association (dh, dofs_domain_association[*it]);

                        for (unsigned int idim=0; idim < dim; ++idim)
                          {
                            types::global_dof_index dummy_1 = sub_wise_to_original[i];
                            types::global_dof_index dummy_2 = sub_wise_to_original[*it];
                            types::global_dof_index index1 = vec_original_to_sub_wise[gradient_dh.n_dofs()/dim*idim+dummy_1];//vector_start_per_process[dofs_domain_association[i]] + idim * owner_el_1 + (i - start_per_process[dofs_domain_association[i]]); //gradient_dh.n_dofs()/dim*idim+i;//vector_start_per_process[this_mpi_process] + (i - start_per_process[this_mpi_process]) * dim + idim; //i*dim+idim
                            types::global_dof_index index2 = vec_original_to_sub_wise[gradient_dh.n_dofs()/dim*idim+dummy_2];//vector_start_per_process[dofs_domain_association[*it]] + idim * owner_el_2 + ((*it) - start_per_process[dofs_domain_association[*it]]);//gradient_dh.n_dofs()/dim*idim+(*it); //vector_start_per_process[this_mpi_process] + ((*it) - start_per_process[this_mpi_process]) * dim + idim;//(*it)*dim+idim
                            normal_distance += localized_normals[index1] * localized_normals[index2];
                          }
                        normal_distance /= normal_distance;
                        if ( normal_distance < 1e-4 )
                          {
                            c.add_line(*it);
                            c.add_entry(*it,i,1);
                          }
                        // this is the dirichlet-dirichlet case on sharp edges: both normal gradients
                        // can be computed from surface gradients of phi and assingned as BC
                        else if (continuos_gradient)
                          {
                            c.add_line(*it);
                            double norm_i_norm_it = 0;
                            double surf_it_norm_i = 0;
                            double surf_i_norm_it = 0;

                            // types::global_dof_index owner_el_1 = DoFTools::count_dofs_with_subdomain_association (dh, dofs_domain_association[i]);
                            // types::global_dof_index owner_el_2 = DoFTools::count_dofs_with_subdomain_association (dh, dofs_domain_association[*it]);

                            // We no longer have a std::vector of Point<dim> so we need to perform the scalar product
                            for (unsigned int idim=0; idim < dim; ++idim)
                              {
                                types::global_dof_index dummy_1 = sub_wise_to_original[i];
                                types::global_dof_index dummy_2 = sub_wise_to_original[*it];

                                types::global_dof_index index1 = vec_original_to_sub_wise[gradient_dh.n_dofs()/dim*idim+dummy_1];//vector_start_per_process[dofs_domain_association[i]] + idim * owner_el_1 + (i - start_per_process[dofs_domain_association[i]]);//gradient_dh.n_dofs()/dim*idim+i;//vector_start_per_process[this_mpi_process] + (i - start_per_process[this_mpi_process]) * dim + idim;
                                types::global_dof_index index2 = vec_original_to_sub_wise[gradient_dh.n_dofs()/dim*idim+dummy_2];//vector_start_per_process[dofs_domain_association[*it]] + idim * owner_el_2 + ((*it) - start_per_process[dofs_domain_association[*it]]);//gradient_dh.n_dofs()/dim*idim+(*it);//vector_start_per_process[this_mpi_process] + ((*it) - start_per_process[this_mpi_process]) * dim + idim;
                                norm_i_norm_it += localized_normals[index1]*localized_normals[index2];
                                surf_it_norm_i += localized_surface_gradients[index2]*localized_normals[index1];
                                surf_i_norm_it += localized_surface_gradients[index1]*localized_normals[index2];
                              }
                            double this_normal_gradient = (1.0/(1.0-pow(norm_i_norm_it,2))) *
                                                          (surf_it_norm_i+
                                                           (surf_i_norm_it)*(norm_i_norm_it));
                            double other_normal_gradient = (1.0/(1.0-pow(norm_i_norm_it,2))) *
                                                           (surf_i_norm_it+
                                                            (surf_it_norm_i)*(norm_i_norm_it));
                            //std::cout<<"i="<<i<<" j="<<*it<<std::endl;
                            //std::cout<<"ni=("<<node_normals[i]<<")  nj=("<<node_normals[*it]<<")"<<std::endl;
                            //std::cout<<"grad_s_phi_i=("<<node_surface_gradients[i]<<")  grad_s_phi_j=("<<node_surface_gradients[*it]<<")"<<std::endl;
                            //std::cout<<"dphi_dn_i="<<this_normal_gradient<<" dphi_dn_j="<<other_normal_gradient<<std::endl;
                            //Point<3> this_full_gradient = node_normals[i]*this_normal_gradient + node_surface_gradients[i];
                            //Point<3> other_full_gradient = node_normals[*it]*other_normal_gradient + node_surface_gradients[*it];
                            //std::cout<<"grad_phi_i=("<<this_full_gradient<<")  grad_phi_j=("<<other_full_gradient<<")"<<std::endl;
                            c.add_line(i);
                            c.set_inhomogeneity(i,this_normal_gradient);
                            c.add_line(*it);
                            c.set_inhomogeneity(*it,other_normal_gradient);
                          }
                      }
                    else
                      {
                        c.add_line(*it);
                        c.set_inhomogeneity(*it,loc_tmp_rhs(i));
                        //dst(*it) = phi(*it)/alpha(*it);
                      }
                  }
                }
            }

          // if the current (first) node is a neumann node, for all its doubles we will impose that
          // the potential is equal to that of the first node: this means that in the matrix vector
          // product we will put the difference between the potential at the fist node in the doubles
          // set, and the current double node
          if (localized_dirichlet_nodes(i) == 0)
            {
              for (std::set<types::global_dof_index>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
                {
                  c.add_line(*it);
                  c.add_entry(*it,i,1);
                  //dst(*it) = phi(*it)/alpha(*it)-phi(i)/alpha(i);
                }
            }
        }
      // else if(firstOfDoubles == *doubles.begin())
      // {
      //   for(std::set<types::global_dof_index>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
      //     if(*it!=firstOfDoubles)
      //       {
      //         c.add_line(*it);
      //         c.add_entry(*it,firstOfDoubles,1);
      //       }
      // }
      // }
    }

  c.merge(c_hn);
  c.close();

  c_cpu_set.clear();
  c_cpu_set.set_size(this_cpu_set.size());
  for (types::global_dof_index i=0; i<dh.n_dofs(); ++i)
    {
      if (this_cpu_set.is_element(i))
        {
          c_cpu_set.add_index(i);
          if (c.is_constrained(i))
            {
              const std::vector< std::pair < types::global_dof_index, double > >
              *entries = c.get_constraint_entries (i);
              for (types::global_dof_index j=0; j< entries->size(); ++j)
                c_cpu_set.add_index((*entries)[j].first);

            }

        }
    }
  c_cpu_set.compress();

  /*
  pcout<<"CONSTAINT MATRIX CHECK "<<std::endl;
    for (types::global_dof_index i=0; i<dh.n_dofs(); ++i)
        {
        std::set <types::global_dof_index> duplicates = double_nodes_set[i];
        if (duplicates.size()>1)
           {
           pcout<<"Proc: "<<this_mpi_process<<" i= "<<i<<" ("<<localized_dirichlet_nodes(i)<<") duplicates: ";
           for (std::set<types::global_dof_index>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
               pcout<<" "<<*pos;
           pcout<<std::endl;
           }
        }

    for(unsigned int i=0; i<dh.n_dofs(); ++i)
      if( (constraints.is_constrained(i)) )
        {pcout<<"Proc: "<<this_mpi_process<<" i= "<<i<<" (";
    const std::vector< std::pair < types::global_dof_index, double > >
      * entries = constraints.get_constraint_entries (i);
          pcout<<entries->size()<<")  Entries:";
    for(unsigned int j=0; j< entries->size(); ++j)
       pcout<<" "<<(*entries)[j].first<<" ("<<(*entries)[j].second<<") ";
         pcout<<" Inomogeneities: "<<constraints.get_inhomogeneity(i)<<std::endl;
        }
  */
}

template<int dim>
void BEMProblem<dim>::assemble_preconditioner()
{


  if (is_preconditioner_initialized == false)
    {
      // pcout<<"Initialising preconditioner"<<std::endl;
      for (types::global_dof_index i=0; i<dh.n_dofs(); ++i)
        if (this_cpu_set.is_element( i))
          {
            // types::global_dof_index start_helper, end_helper;
            // if(i>preconditioner_band/2)
            //   start_helper = i-preconditioner_band/2;
            // else
            //   start_helper = (types::global_dof_index) 0;
            // if(i+preconditioner_band/2 < dh.n_dofs())
            //   end_helper = i+preconditioner_band/2;
            // else
            //   end_helper = dh.n_dofs();
            //   for(types::global_dof_index j=start_helper; j<end_helper; ++j)
            // pcout<<start_helper<<" "<<std::min((types::global_dof_index)(i+preconditioner_band/2),(types::global_dof_index)dh.n_dofs())<<std::endl;
            types::global_dof_index start_helper= ((i) > preconditioner_band/2) ? (i-preconditioner_band/2): ((types::global_dof_index)0);
            for (types::global_dof_index j=start_helper; j<std::min((types::global_dof_index)(i+preconditioner_band/2),(types::global_dof_index)dh.n_dofs()); ++j)
              preconditioner_sparsity_pattern.add(i,j);
          }
      preconditioner_sparsity_pattern.compress();
      band_system.reinit(preconditioner_sparsity_pattern);
      is_preconditioner_initialized = true;
    }
  else
    band_system = 0;


  for (types::global_dof_index i=0; i<dh.n_dofs(); ++i)
    {
      if (this_cpu_set.is_element(i))
        {
          if (constraints.is_constrained(i))
            band_system.add(i, i, 1);
          // types::global_dof_index start_helper, end_helper;
          // if(i>preconditioner_band/2)
          //   start_helper = i-preconditioner_band/2;
          // else
          //   start_helper = (types::global_dof_index) 0;
          // if(i+preconditioner_band/2 < dh.n_dofs())
          //   end_helper = i+preconditioner_band/2;
          // else
          //   end_helper = dh.n_dofs();
          // for(types::global_dof_index j=start_helper; j<end_helper; ++j)
          types::global_dof_index start_helper= ((i) > preconditioner_band/2) ? (i-preconditioner_band/2): ((types::global_dof_index)0);

          for (types::global_dof_index j=start_helper; j<std::min((types::global_dof_index)i+preconditioner_band/2,(types::global_dof_index)dh.n_dofs()); ++j)
            {
              if (constraints.is_constrained(i) == false)
                {
                  if (dirichlet_nodes(i) == 0)
                    {
                      // Nodo di Dirichlet
                      band_system.add(i,j,neumann_matrix(i,j));

                      if (i == j)
                        band_system.add(i, j, alpha(i));

                    }
                  else
                    band_system.add(i,j,-dirichlet_matrix(i,j));
                }
            }
        }
    }




  preconditioner.initialize(band_system);

  /*
  band_system.vmult(sol,alpha);
  pcout<<"**solution "<<std::endl;
   for (unsigned int i = 0; i < alpha.size(); i++)
       if (this_cpu_set.is_element(i))
          pcout<<i<<" ("<<this_mpi_process<<")  "<<sol(i)<<"   "<<sol(i)<<std::endl;
  */
}



template <int dim>
void BEMProblem<dim>::compute_gradients(const TrilinosWrappers::MPI::Vector &glob_phi, const TrilinosWrappers::MPI::Vector &glob_dphi_dn)
{
  Teuchos::TimeMonitor LocalTimer(*GradientTime);

  // We need the solution to be stored on a parallel vector with ghost elements. We let
  // Trilinos take care of it.

  TrilinosWrappers::MPI::Vector phi(ghosted_set);
  phi.reinit(glob_phi,false,true);
  TrilinosWrappers::MPI::Vector dphi_dn(ghosted_set);
  dphi_dn.reinit(glob_dphi_dn,false,true);



  // We reinit the gradient solution
  vector_gradients_solution.reinit(vector_this_cpu_set,mpi_communicator);

  typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;


  // The matrix and rhs of our problem. We must decide if compute the mass matrix just once and for all or not.
  TrilinosWrappers::SparseMatrix vector_gradients_matrix;
  TrilinosWrappers::MPI::Vector vector_gradients_rhs(vector_this_cpu_set,mpi_communicator);
  vector_gradients_matrix.reinit (vector_sparsity_pattern);


  // The vector FEValues to used in the assemblage
  FEValues<dim-1,dim> vector_fe_v(*mapping, *gradient_fe, *quadrature,
                                  update_values | update_gradients |
                                  update_cell_normal_vectors |
                                  update_quadrature_points |
                                  update_JxW_values);

  // The scalar FEValues to interpolate the known value of phi
  FEValues<dim-1,dim> fe_v(*mapping, *fe, *quadrature,
                           update_values | update_gradients |
                           update_cell_normal_vectors |
                           update_quadrature_points |
                           update_JxW_values);

  const unsigned int vector_n_q_points = vector_fe_v.n_quadrature_points;
  const unsigned int   vector_dofs_per_cell   = gradient_fe->dofs_per_cell;
  std::vector<types::global_dof_index> vector_local_dof_indices (vector_dofs_per_cell);


  std::vector< Tensor<1,dim> > phi_surf_grads(vector_n_q_points);
  std::vector<double> phi_norm_grads(vector_n_q_points);
  std::vector<Vector<double> > q_vector_normals_solution(vector_n_q_points,
                                                         Vector<double>(dim));

  FullMatrix<double>   local_gradients_matrix (vector_dofs_per_cell, vector_dofs_per_cell);
  Vector<double>       local_gradients_rhs (vector_dofs_per_cell);



  std::vector<Point<dim> > support_points(dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *mapping, dh, support_points);
  std::vector<types::global_dof_index> face_dofs(fe->dofs_per_face);

  Quadrature <dim-1> dummy_quadrature(fe->get_unit_support_points());
  FEValues<dim-1,dim> dummy_fe_v(*mapping, *fe, dummy_quadrature,
                                 update_values | update_gradients |
                                 update_cell_normal_vectors |
                                 update_quadrature_points);

  const unsigned int   dofs_per_cell = fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  const unsigned int n_q_points = dummy_fe_v.n_quadrature_points;
  std::vector< Tensor<1,dim> > dummy_phi_surf_grads(n_q_points);

  cell_it
  vector_cell = gradient_dh.begin_active();

  cell_it
  cell = dh.begin_active(),
  endc = dh.end();


  for (; cell!=endc; ++cell,++vector_cell)
    {
      Assert(cell->index() == vector_cell->index(), ExcInternalError());
      Assert(cell->subdomain_id() == vector_cell->subdomain_id(), ExcInternalError());

      if (cell->subdomain_id() == this_mpi_process)
        {
          fe_v.reinit (cell);
          vector_fe_v.reinit (vector_cell);
          local_gradients_matrix = 0;
          local_gradients_rhs = 0;
          const std::vector<Tensor<1, dim> > &vector_node_normals = vector_fe_v.get_normal_vectors();
          fe_v.get_function_gradients(phi, phi_surf_grads);
          fe_v.get_function_values(dphi_dn, phi_norm_grads);
          unsigned int comp_i, comp_j;




          for (unsigned int q=0; q<vector_n_q_points; ++q)
            {
              Tensor<1, dim> node_normal_grad_dir;
              for (unsigned int i=0; i<dim; ++i)
                node_normal_grad_dir[i] = q_vector_normals_solution[q][i];
              Tensor<1, dim> gradient = vector_node_normals[q]*phi_norm_grads[q] + phi_surf_grads[q];
              for (unsigned int i=0; i<vector_dofs_per_cell; ++i)
                {
                  comp_i = gradient_fe->system_to_component_index(i).first;
                  for (unsigned int j=0; j<vector_dofs_per_cell; ++j)
                    {
                      comp_j = gradient_fe->system_to_component_index(j).first;
                      if (comp_i == comp_j)
                        {
                          local_gradients_matrix(i,j) += vector_fe_v.shape_value(i,q)*
                                                         vector_fe_v.shape_value(j,q)*
                                                         vector_fe_v.JxW(q);
                        }
                    }
                  local_gradients_rhs(i) += (vector_fe_v.shape_value(i, q)) *
                                            gradient[comp_i] * vector_fe_v.JxW(q);
                }
            }
          vector_cell->get_dof_indices (vector_local_dof_indices);

          vector_constraints.distribute_local_to_global
          (local_gradients_matrix,
           local_gradients_rhs,
           vector_local_dof_indices,
           vector_gradients_matrix,
           vector_gradients_rhs);
        }
    }

  // At this point we can compress anything and solve via GMRES.
  vector_gradients_matrix.compress(VectorOperation::add);
  vector_gradients_rhs.compress(VectorOperation::add);

  SolverGMRES<TrilinosWrappers::MPI::Vector > solver (solver_control,
                                                      SolverGMRES<TrilinosWrappers::MPI::Vector >::AdditionalData(1000));

  solver.solve (vector_gradients_matrix, vector_gradients_solution, vector_gradients_rhs, PreconditionIdentity());

  vector_constraints.distribute(vector_gradients_solution);
}

template <int dim>
void BEMProblem<dim>::compute_surface_gradients(const TrilinosWrappers::MPI::Vector &tmp_rhs)
{
  Teuchos::TimeMonitor LocalTimer(*SurfaceGradientTime);
  TrilinosWrappers::MPI::Vector phi(ghosted_set);
  phi.reinit(tmp_rhs,false,true);

  vector_surface_gradients_solution.reinit(vector_this_cpu_set,mpi_communicator);


  typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;


  TrilinosWrappers::SparseMatrix vector_surface_gradients_matrix;
  TrilinosWrappers::MPI::Vector vector_surface_gradients_rhs(vector_this_cpu_set,mpi_communicator);


  vector_surface_gradients_matrix.reinit (vector_sparsity_pattern);



  FEValues<dim-1,dim> vector_fe_v(*mapping, *gradient_fe, *quadrature,
                                  update_values | update_gradients |
                                  update_cell_normal_vectors |
                                  update_quadrature_points |
                                  update_JxW_values);

  FEValues<dim-1,dim> fe_v(*mapping, *fe, *quadrature,
                           update_values | update_gradients |
                           update_cell_normal_vectors |
                           update_quadrature_points |
                           update_JxW_values);

  const unsigned int vector_n_q_points = vector_fe_v.n_quadrature_points;
  const unsigned int   vector_dofs_per_cell   = gradient_fe->dofs_per_cell;
  std::vector<types::global_dof_index> vector_local_dof_indices (vector_dofs_per_cell);

  std::vector< Tensor<1,dim> > phi_surf_grads(vector_n_q_points);
  std::vector<double> phi_norm_grads(vector_n_q_points);
  std::vector<Vector<double> > q_vector_normals_solution(vector_n_q_points,
                                                         Vector<double>(dim));

  FullMatrix<double>   local_gradients_matrix (vector_dofs_per_cell, vector_dofs_per_cell);
  Vector<double>       local_gradients_rhs (vector_dofs_per_cell);



  std::vector<Point<dim> > support_points(dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *mapping, dh, support_points);
  std::vector<types::global_dof_index> face_dofs(fe->dofs_per_face);

  Quadrature <dim-1> dummy_quadrature(fe->get_unit_support_points());
  FEValues<dim-1,dim> dummy_fe_v(*mapping, *fe, dummy_quadrature,
                                 update_values | update_gradients |
                                 update_cell_normal_vectors |
                                 update_quadrature_points);

  const unsigned int   dofs_per_cell = fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  const unsigned int n_q_points = dummy_fe_v.n_quadrature_points;
  std::vector< Tensor<1,dim> > dummy_phi_surf_grads(n_q_points);

  cell_it
  vector_cell = gradient_dh.begin_active();

  cell_it
  cell = dh.begin_active(),
  endc = dh.end();


  for (; cell!=endc; ++cell,++vector_cell)
    {
      Assert(cell->index() == vector_cell->index(), ExcInternalError());
      Assert(cell->subdomain_id() == vector_cell->subdomain_id(), ExcInternalError());

      if (cell->subdomain_id() == this_mpi_process)
        {
          fe_v.reinit (cell);
          vector_fe_v.reinit (vector_cell);
          local_gradients_matrix = 0;
          local_gradients_rhs = 0;
          fe_v.get_function_gradients(phi, phi_surf_grads);
          unsigned int comp_i, comp_j;




          for (unsigned int q=0; q<vector_n_q_points; ++q)
            {
              Tensor<1,dim> gradient = phi_surf_grads[q];
              for (unsigned int i=0; i<vector_dofs_per_cell; ++i)
                {
                  comp_i = gradient_fe->system_to_component_index(i).first;
                  for (unsigned int j=0; j<vector_dofs_per_cell; ++j)
                    {
                      comp_j = gradient_fe->system_to_component_index(j).first;
                      if (comp_i == comp_j)
                        {
                          local_gradients_matrix(i,j) += vector_fe_v.shape_value(i,q)*
                                                         vector_fe_v.shape_value(j,q)*
                                                         vector_fe_v.JxW(q);
                        }
                    }
                  local_gradients_rhs(i) += (vector_fe_v.shape_value(i, q)) *
                                            gradient[comp_i] * vector_fe_v.JxW(q);
                }
            }
          vector_cell->get_dof_indices (vector_local_dof_indices);

          vector_constraints.distribute_local_to_global
          (local_gradients_matrix,
           local_gradients_rhs,
           vector_local_dof_indices,
           vector_surface_gradients_matrix,
           vector_surface_gradients_rhs);

        }
    }

  vector_surface_gradients_matrix.compress(VectorOperation::add);
  vector_surface_gradients_rhs.compress(VectorOperation::add);

  SolverGMRES<TrilinosWrappers::MPI::Vector > solver (solver_control,
                                                      SolverGMRES<TrilinosWrappers::MPI::Vector >::AdditionalData(1000));

  solver.solve (vector_surface_gradients_matrix, vector_surface_gradients_solution, vector_surface_gradients_rhs, PreconditionIdentity());

  vector_constraints.distribute(vector_surface_gradients_solution);

}


template <int dim>
void BEMProblem<dim>::compute_normals()
{
  Teuchos::TimeMonitor LocalTimer(*NormalsTime);
  vector_normals_solution.reinit(vector_this_cpu_set,mpi_communicator);

  typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;


  TrilinosWrappers::SparseMatrix vector_normals_matrix;
  TrilinosWrappers::MPI::Vector vector_normals_rhs(vector_this_cpu_set,mpi_communicator);


  vector_normals_matrix.reinit (vector_sparsity_pattern);



  FEValues<dim-1,dim> vector_fe_v(*mapping, *gradient_fe, *quadrature,
                                  update_values | update_gradients |
                                  update_cell_normal_vectors |
                                  update_quadrature_points |
                                  update_JxW_values);

  const unsigned int vector_n_q_points = vector_fe_v.n_quadrature_points;

  const unsigned int   vector_dofs_per_cell   = gradient_fe->dofs_per_cell;

  std::vector<types::global_dof_index> vector_local_dof_indices (vector_dofs_per_cell);

  std::vector<Vector<double> > q_vector_normals_solution(vector_n_q_points,
                                                         Vector<double>(dim));

  FullMatrix<double>   local_normals_matrix (vector_dofs_per_cell, vector_dofs_per_cell);
  Vector<double>       local_normals_rhs (vector_dofs_per_cell);


  cell_it
  vector_cell = gradient_dh.begin_active(),
  vector_endc = gradient_dh.end();


  for (; vector_cell!=vector_endc; ++vector_cell)
    {

      if (vector_cell->subdomain_id() == this_mpi_process)
        {
          vector_fe_v.reinit (vector_cell);
          local_normals_matrix = 0;
          local_normals_rhs = 0;
          const std::vector<Tensor<1, dim> > &vector_node_normals = vector_fe_v.get_normal_vectors();
          unsigned int comp_i, comp_j;

          for (unsigned int q=0; q<vector_n_q_points; ++q)
            for (unsigned int i=0; i<vector_dofs_per_cell; ++i)
              {
                comp_i = gradient_fe->system_to_component_index(i).first;
                for (unsigned int j=0; j<vector_dofs_per_cell; ++j)
                  {
                    comp_j = gradient_fe->system_to_component_index(j).first;
                    if (comp_i == comp_j)
                      {
                        local_normals_matrix(i,j) += vector_fe_v.shape_value(i,q)*
                                                     vector_fe_v.shape_value(j,q)*
                                                     vector_fe_v.JxW(q);
                      }
                  }

                local_normals_rhs(i) += (vector_fe_v.shape_value(i, q)) *
                                        vector_node_normals[q][comp_i] * vector_fe_v.JxW(q);
              }

          vector_cell->get_dof_indices (vector_local_dof_indices);

          vector_constraints.distribute_local_to_global
          (local_normals_matrix,
           local_normals_rhs,
           vector_local_dof_indices,
           vector_normals_matrix,
           vector_normals_rhs);
        }
    }

  vector_normals_matrix.compress(VectorOperation::add);
  vector_normals_rhs.compress(VectorOperation::add);

  SolverGMRES<TrilinosWrappers::MPI::Vector > solver (solver_control,
                                                      SolverGMRES<TrilinosWrappers::MPI::Vector >::AdditionalData(1000));

  solver.solve (vector_normals_matrix, vector_normals_solution, vector_normals_rhs, PreconditionIdentity());

  vector_constraints.distribute(vector_normals_solution);

}

template<int dim>
void BEMProblem<dim>::adaptive_refinement(const TrilinosWrappers::MPI::Vector &error_vector)
{

  Vector<float> estimated_error_per_cell (comp_dom.tria.n_active_cells());
  Vector<double> helper(error_vector);

  KellyErrorEstimator<dim-1, dim>::estimate (*mapping,
                                             dh,
                                             QGauss<dim-2> (3),
                                             typename FunctionMap<dim>::type(),
                                             helper,
                                             estimated_error_per_cell);

  pgr.mark_cells(estimated_error_per_cell, comp_dom.tria);
  //  GridRefinement::refine_and_coarsen_fixed_number (comp_dom.tria,
  //                                                  estimated_error_per_cell,
  //                                                  refinement_threshold, coarsening_threshold);
  comp_dom.tria.prepare_coarsening_and_refinement();
  comp_dom.tria.execute_coarsening_and_refinement ();


}





template class BEMProblem<2>;
template class BEMProblem<3>;
