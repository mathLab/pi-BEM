//----------------------------  step-34.cc  ---------------------------
//    $Id: step-34.cc 18734 2009-04-25 13:36:48Z heltai $
//    Version: $Name$
//
//    Copyright (C) 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//    Authors: Luca Heltai, Cataldo Manigrasso
//
//----------------------------  step-34.cc  ---------------------------


#include "../include/bem_problem.h"
#include "../include/laplace_kernel.h"

#include <iostream>
#include <iomanip>

#include "Teuchos_TimeMonitor.hpp"

using Teuchos::Time;
using Teuchos::TimeMonitor;
using Teuchos::RCP;

RCP<Time> AssembleTime = TimeMonitor::getNewTimer("Assemble Time");
RCP<Time> LacSolveTime = TimeMonitor::getNewTimer("LAC Solve Time");

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
template <int dim>
BEMProblem<dim>::BEMProblem(ComputationalDomain<dim> &comp_dom,
                            BEMFMA<dim> &fma,
                            MPI_Comm comm)
  :
  pcout(std::cout),
  comp_dom(comp_dom), fma(fma),
  mpi_communicator (MPI_COMM_WORLD),
  n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
  this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator))
{
  // Only output on first processor.
  pcout.set_condition(this_mpi_process == 0);
}

template <int dim>
void BEMProblem<dim>::reinit()
{
  const unsigned int n_dofs =  comp_dom.dh.n_dofs();
  //const types::global_dof_index n_local_dofs = DoFTools::count_dofs_with_subdomain_association(comp_dom.dh,this_mpi_process);
  //this_cpu_set = DoFTools::dof_indices_with_subdomain_association(comp_dom.dh,this_mpi_process);
  std::vector<types::subdomain_id> dofs_domain_association(n_dofs);
  DoFTools::get_subdomain_association   (comp_dom.dh,dofs_domain_association);
  this_cpu_set.set_size(n_dofs);
  for (unsigned int i=0; i<n_dofs; ++i)
    if (dofs_domain_association[i] == this_mpi_process)
      {
        this_cpu_set.add_index(i);
      }
  this_cpu_set.compress();
  //pcout<<"****  Proc: "<<this_mpi_process<<"  Dofs: "<<n_local_dofs<<"  of "<<n_dofs<<std::endl;

  system_rhs.reinit(this_cpu_set,mpi_communicator);
  sol.reinit(this_cpu_set,mpi_communicator);
  alpha.reinit(this_cpu_set,mpi_communicator);
  serv_phi.reinit(this_cpu_set,mpi_communicator);
  serv_dphi_dn.reinit(this_cpu_set,mpi_communicator);
  serv_tmp_rhs.reinit(this_cpu_set,mpi_communicator);


  full_sparsity_pattern.reinit(sol.vector_partitioner(), n_dofs);
  for (unsigned int i=0; i<n_dofs; ++i)
    if (this_cpu_set.is_element(i))
      {
        for (unsigned int j=0; j<n_dofs; ++j)
          full_sparsity_pattern.add(i,j);
      }
  full_sparsity_pattern.compress();
  neumann_matrix.reinit(full_sparsity_pattern);
  dirichlet_matrix.reinit(full_sparsity_pattern);
  preconditioner_band = 100;
  preconditioner_sparsity_pattern.reinit(sol.vector_partitioner(), (unsigned int) preconditioner_band);
  is_preconditioner_initialized = false;

  surface_nodes.reinit(this_cpu_set,mpi_communicator);
  other_nodes.reinit(this_cpu_set,mpi_communicator);
  comp_dom.compute_phi_nodes();
  surface_nodes = comp_dom.surface_nodes;
  other_nodes = comp_dom.other_nodes;


  /*
  for (unsigned int i=0; i<surface_nodes.size(); ++i)
      if (this_cpu_set.is_element(i))
         {
         if (true)//(this_mpi_process == 0)
            {
            pcout<<"Proc: "<<this_mpi_process<<"   i: "<<i<<" Sn: "<<surface_nodes(i)<<" "<<comp_dom.surface_nodes(i) <<std::endl;
            pcout<<"Proc: "<<this_mpi_process<<"   i: "<<i<<" On: "<<other_nodes(i)<<" "<<comp_dom.other_nodes(i) <<std::endl;
            }
         }
  //*/
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

  prm.declare_entry("Solution method", "Direct",
                    Patterns::Selection("Direct|FMA"));
}

template <int dim>
void BEMProblem<dim>::parse_parameters (ParameterHandler &prm)
{

  prm.enter_subsection("Solver");
  solver_control.parse_parameters(prm);
  prm.leave_subsection();

  solution_method = prm.get("Solution method");


}



template <int dim>
void BEMProblem<dim>::assemble_system()
{
  TimeMonitor LocalTimer(*AssembleTime);
  pcout<<"(Directly) Assembling system matrices"<<std::endl;

  neumann_matrix = 0;
  dirichlet_matrix = 0;



  std::vector<QGaussOneOverR<2> > sing_quadratures_3d;
  for (unsigned int i=0; i<comp_dom.fe.dofs_per_cell; ++i)
    sing_quadratures_3d.push_back
    (QGaussOneOverR<2>(comp_dom.singular_quadrature_order,
                       comp_dom.fe.get_unit_support_points()[i], true));


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
  FEValues<dim-1,dim> fe_v(*comp_dom.mapping,comp_dom.fe, *comp_dom.quadrature,
                           update_values |
                           update_cell_normal_vectors |
                           update_quadrature_points |
                           update_JxW_values);

  const unsigned int n_q_points = fe_v.n_quadrature_points;

  std::vector<unsigned int> local_dof_indices(comp_dom.fe.dofs_per_cell);

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
  // vector of fe.dofs_per_cell
  // elements, which will then be
  // distributed to the matrix in the
  // global row $i$. The following
  // object will hold this
  // information:
  Vector<double>      local_neumann_matrix_row_i(comp_dom.fe.dofs_per_cell);
  Vector<double>      local_dirichlet_matrix_row_i(comp_dom.fe.dofs_per_cell);

  // Now that we have checked that
  // the number of vertices is equal
  // to the number of degrees of
  // freedom, we construct a vector
  // of support points which will be
  // used in the local integrations:
  std::vector<Point<dim> > support_points(comp_dom.dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.dh, support_points);


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
  cell = comp_dom.dh.begin_active(),
  endc = comp_dom.dh.end();

  Point<dim> D;
  double s;

  for (cell = comp_dom.dh.begin_active(); cell != endc; ++cell)
    {
      fe_v.reinit(cell);
      cell->get_dof_indices(local_dof_indices);

      const std::vector<Point<dim> > &q_points = fe_v.get_quadrature_points();
      const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors();

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
      for (unsigned int i=0; i<comp_dom.dh.n_dofs() ; ++i) //these must now be the locally owned dofs. the rest should stay the same
        {
          if (this_cpu_set.is_element(i))
            {
              local_neumann_matrix_row_i = 0;
              local_dirichlet_matrix_row_i = 0;

              bool is_singular = false;
              unsigned int singular_index = numbers::invalid_unsigned_int;

              for (unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j)
                //if(local_dof_indices[j] == i)
                if (comp_dom.double_nodes_set[i].count(local_dof_indices[j]) > 0)
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

                      for (unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j)
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
                    = (dim == 2
                       ?
                       dynamic_cast<Quadrature<dim-1>*>(
                         new QGaussLogR<1>(comp_dom.singular_quadrature_order,
                                           Point<1>((double)singular_index),
                                           1./cell->measure(), true))
                       :
                       (dim == 3
                        ?
                        dynamic_cast<Quadrature<dim-1>*>(
                          &sing_quadratures_3d[singular_index])
                        :
                        0));
                  Assert(singular_quadrature, ExcInternalError());

                  FEValues<dim-1,dim> fe_v_singular (*comp_dom.mapping, comp_dom.fe, *singular_quadrature,
                                                     update_jacobians |
                                                     update_values |
                                                     update_cell_normal_vectors |
                                                     update_quadrature_points );

                  fe_v_singular.reinit(cell);

                  const std::vector<Point<dim> > &singular_normals = fe_v_singular.get_normal_vectors();
                  const std::vector<Point<dim> > &singular_q_points = fe_v_singular.get_quadrature_points();

                  for (unsigned int q=0; q<singular_quadrature->size(); ++q)
                    {
                      const Tensor<1,dim> R = singular_q_points[q] - support_points[i];
                      LaplaceKernel::kernels(R, D, s);

                      for (unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j)
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
                  if (dim==2)
                    delete singular_quadrature;
                }

              // Finally, we need to add
              // the contributions of the
              // current cell to the
              // global matrix.
              for (unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j)
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
    for (unsigned int i = 0; i < comp_dom.dh.n_dofs(); i++)
        {
        if (this_cpu_set.is_element(i))
           {
           pcout<<this_mpi_process<<" *** ";
           for (unsigned int j = 0; j < comp_dom.dh.n_dofs(); j++)
               {
               pcout<<neumann_matrix(i,j)<<" ";
               }
           pcout<<std::endl;
           }
        }



    pcout<<"Dirichlet"<<std::endl;
    for (unsigned int i = 0; i < comp_dom.dh.n_dofs(); i++)
        {
        if (this_cpu_set.is_element(i))
           {
           pcout<<this_mpi_process<<" *** ";
           for (unsigned int j = 0; j < comp_dom.dh.n_dofs(); j++)
               {
               pcout<<dirichlet_matrix(i,j)<<" ";
               }
           pcout<<std::endl;
           }
        }
        //*/
  pcout<<"done assembling system matrices"<<std::endl;
}



template <int dim>
void BEMProblem<dim>::compute_alpha()
{
  static TrilinosWrappers::MPI::Vector ones, zeros, dummy;
  if (ones.size() != comp_dom.dh.n_dofs())
    {
      ones.reinit(this_cpu_set,mpi_communicator);
      ones.add(-1.);
      zeros.reinit(this_cpu_set,mpi_communicator);
      dummy.reinit(this_cpu_set,mpi_communicator);
    }


  if (solution_method == "Direct")
    {
      neumann_matrix.vmult(alpha, ones);
    }
  else
    {
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
  serv_dphi_dn = src;



  TrilinosWrappers::MPI::Vector matrVectProdN;
  TrilinosWrappers::MPI::Vector matrVectProdD;

  matrVectProdN.reinit(this_cpu_set,mpi_communicator);
  matrVectProdD.reinit(this_cpu_set,mpi_communicator);

  dst = 0;


  serv_phi.scale(other_nodes);
  serv_dphi_dn.scale(surface_nodes);

  //phi.print(pcout);
  //dphi_dn.print(pcout);
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
      fma.generate_multipole_expansions(serv_phi,serv_dphi_dn);
      fma.multipole_matr_vect_products(serv_phi,serv_dphi_dn,matrVectProdN,matrVectProdD);
      serv_phi.scale(alpha);
      dst += matrVectProdD;
      dst *= -1;
      dst += matrVectProdN;
      dst += serv_phi;
    }
  // in fully neumann bc case, we have to rescale the vector to have a zero mean
  // one
// DA RIVEDERE!!! PER ORA QUESTO CASO IN MPI NON GIREREBBE
  //if (comp_dom.surface_nodes.linfty_norm() < 1e-10)
  //   dst.add(-dst.l2_norm());

}


template <int dim>
void BEMProblem<dim>::compute_rhs(TrilinosWrappers::MPI::Vector &dst, const TrilinosWrappers::MPI::Vector &src) const
{



  serv_phi = src;
  serv_dphi_dn = src;

  static TrilinosWrappers::MPI::Vector matrVectProdN;
  static TrilinosWrappers::MPI::Vector matrVectProdD;

  //const types::global_dof_index n_local_dofs = DoFTools::count_dofs_with_subdomain_association(comp_dom.dh,this_mpi_process);
  IndexSet this_cpu_set = comp_dom.dh.locally_owned_dofs();

  matrVectProdN.reinit(this_cpu_set,mpi_communicator);
  matrVectProdD.reinit(this_cpu_set,mpi_communicator);


  serv_phi.scale(surface_nodes);
  serv_dphi_dn.scale(other_nodes);

  //for (unsigned int ppp = 0; ppp < src.size(); ++ppp)
  //pcout<<ppp<<"  phi(ppp)="<<phi(ppp)<<"  |||   dphi_dn(ppp)="<<dphi_dn(ppp)<<std::endl;

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
  TimeMonitor LocalTimer(*LacSolveTime);
  SolverGMRES<TrilinosWrappers::MPI::Vector > solver (solver_control,
                                                      SolverGMRES<TrilinosWrappers::MPI::Vector >::AdditionalData(50));


  system_rhs = 0;
  sol = 0;
  alpha = 0;


  compute_alpha();

  //pcout<<"alpha "<<std::endl;
  //for (unsigned int i = 0; i < alpha.size(); i++)
  //    if (this_cpu_set.is_element(i))
  //       pcout<<i<<" ("<<this_mpi_process<<")  "<<alpha(i)<<std::endl;




  compute_rhs(system_rhs, tmp_rhs);

  compute_constraints(constraints, tmp_rhs);
  ConstrainedOperator<TrilinosWrappers::MPI::Vector, BEMProblem<dim> >
  cc(*this, constraints);

  cc.distribute_rhs(system_rhs);




  if (solution_method == "Direct")
    {
      //SparseDirectUMFPACK &inv = fma.FMA_preconditioner(alpha);
      //solver.solve (*this, sol, system_rhs, inv);
      assemble_preconditioner();
      //solver.solve (cc, sol, system_rhs, PreconditionIdentity());
      solver.solve (cc, sol, system_rhs, preconditioner);
    }
  else
    {
      TrilinosWrappers::PreconditionILU &fma_preconditioner = fma.FMA_preconditioner(alpha,constraints);
      solver.solve (cc, sol, system_rhs, fma_preconditioner);
      //solver.solve (cc, sol, system_rhs, PreconditionIdentity());
    }

  //pcout<<"sol = [";
  //for (unsigned int i = 0; i < comp_dom.dh.n_dofs(); i++)
  //    pcout<<sol(i)<<"; ";
  //pcout<<"];"<<std::endl;




///////////////////////////////////
  /*
    std::vector<Point<dim> > support_points(comp_dom.dh.n_dofs());
    DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.dh, support_points);
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


  for (unsigned int i=0; i <surface_nodes.size(); i++)
    {
      if (this_cpu_set.is_element(i))
        {
          if (surface_nodes(i) == 0)
            {
              phi(i) = sol(i);
            }
          else
            {
              dphi_dn(i) = sol(i);
            }
        }
    }


  //for (unsigned int i=0;i<comp_dom.dh.n_dofs();++i)
  // cout<<i<<" "<<tmp_rhs(i)<<" "<<dphi_dn(i)<<" "<<phi(i)<<" "<<surface_nodes(i)<<endl;

  //pcout<<"sol "<<std::endl;
  //for (unsigned int i = 0; i < sol.size(); i++)
  //    {
  //pcout<<i<<" "<<sol(i)<<" ";
  //std::set<unsigned int> doubles = comp_dom.double_nodes_set[i];
  // for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
  //    pcout<<*it<<"("<<surface_nodes(*it)<<") ";
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
    }
  else
    {
      comp_dom.generate_octree_blocking();
      fma.direct_integrals();
      fma.multipole_integrals();
    }

  solve_system(phi,dphi_dn,tmp_rhs);

}


template <int dim>
void BEMProblem<dim>::compute_constraints(ConstraintMatrix &c, const TrilinosWrappers::MPI::Vector &tmp_rhs)

{
  // communication is needed here: there is one matrix per process: thus the vector needed to set
  // inhomogeneities has to be copied locally
  comp_dom.compute_normals();
  compute_surface_gradients(tmp_rhs);

  Vector<double> loc_tmp_rhs(tmp_rhs.size());
  loc_tmp_rhs = tmp_rhs;

  // we start clearing the constraint matrix
  c.clear();

  // here we prepare the constraint matrix so as to account for the presence hanging
  // nodes

  DoFTools::make_hanging_node_constraints (comp_dom.dh,c);


  // here we prepare the constraint matrix so as to account for the presence of double and
  // triple dofs

  // we start looping on the dofs
  for (unsigned int i=0; i <tmp_rhs.size(); i++)
    {
      //if (this_cpu_set.is_element(i))
      {
        // in the next line we compute the "first" among the set of double nodes: this node
        // is the first dirichlet node in the set, and if no dirichlet node is there, we get the
        // first neumann node

        std::set<unsigned int> doubles = comp_dom.double_nodes_set[i];
        unsigned int firstOfDoubles = *doubles.begin();
        for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
          {
            if (comp_dom.surface_nodes(*it) == 1)
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
            if (comp_dom.surface_nodes(i) == 1)
              {
                for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
                  {
                    if (comp_dom.surface_nodes(*it) == 1)
                      {
                        // this is the dirichlet-dirichlet case on flat edges: here we impose that
                        // dphi_dn on the two (or more) sides is equal.
                        if ( comp_dom.node_normals[*it].distance(comp_dom.node_normals[i]) < 1e-4 )
                          {
                            c.add_line(*it);
                            c.add_entry(*it,i,1);
                          }
                        // this is the dirichlet-dirichlet case on sharp edges: both normal gradients
                        // can be computed from surface gradients of phi and assingned as BC
                        else
                          {
                            c.add_line(*it);
                            double this_normal_gradient = (1.0/(1.0-pow(comp_dom.node_normals[i]*comp_dom.node_normals[*it],2))) *
                                                          (node_surface_gradients[*it]*comp_dom.node_normals[i]+
                                                           (node_surface_gradients[i]*comp_dom.node_normals[*it])*(comp_dom.node_normals[i]*comp_dom.node_normals[*it]));
                            double other_normal_gradient = (1.0/(1.0-pow(comp_dom.node_normals[*it]*comp_dom.node_normals[i],2))) *
                                                           (node_surface_gradients[i]*comp_dom.node_normals[*it]+
                                                            (node_surface_gradients[*it]*comp_dom.node_normals[i])*(comp_dom.node_normals[*it]*comp_dom.node_normals[i]));
                            //std::cout<<"i="<<i<<" j="<<*it<<std::endl;
                            //std::cout<<"ni=("<<comp_dom.node_normals[i]<<")  nj=("<<comp_dom.node_normals[*it]<<")"<<std::endl;
                            //std::cout<<"grad_s_phi_i=("<<node_surface_gradients[i]<<")  grad_s_phi_j=("<<node_surface_gradients[*it]<<")"<<std::endl;
                            //std::cout<<"dphi_dn_i="<<this_normal_gradient<<" dphi_dn_j="<<other_normal_gradient<<std::endl;
                            //Point<3> this_full_gradient = comp_dom.node_normals[i]*this_normal_gradient + node_surface_gradients[i];
                            //Point<3> other_full_gradient = comp_dom.node_normals[*it]*other_normal_gradient + node_surface_gradients[*it];
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

            // if the current (first) node is a neumann node, for all its doubles we will impose that
            // the potential is equal to that of the first node: this means that in the matrix vector
            // product we will put the difference between the potential at the fist node in the doubles
            // set, and the current double node
            if (comp_dom.surface_nodes(i) == 0)
              {
                for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
                  {
                    c.add_line(*it);
                    c.add_entry(*it,i,1);
                    //dst(*it) = phi(*it)/alpha(*it)-phi(i)/alpha(i);
                  }
              }
          }
      }
    }

  c.close();

  /*
  pcout<<"CONSTAINT MATRIX CHECK "<<std::endl;
    for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
        {
        std::set <unsigned int> duplicates = comp_dom.double_nodes_set[i];
        if (duplicates.size()>1)
           {
           pcout<<"Proc: "<<this_mpi_process<<" i= "<<i<<" ("<<comp_dom.surface_nodes(i)<<") duplicates: ";
           for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
               pcout<<" "<<*pos;
           pcout<<std::endl;
           }
        }

    for(unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      if( (constraints.is_constrained(i)) )
        {pcout<<"Proc: "<<this_mpi_process<<" i= "<<i<<" (";
    const std::vector< std::pair < unsigned int, double > >
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
      for (int i=0; (unsigned int) i<comp_dom.dh.n_dofs(); ++i)
        if (this_cpu_set.is_element((unsigned int) i))
          for (int j=std::max(i-preconditioner_band/2,0); j<std::min(i+preconditioner_band/2,(int)comp_dom.dh.n_dofs()); ++j)
            preconditioner_sparsity_pattern.add((unsigned int) i,(unsigned int) j);
      preconditioner_sparsity_pattern.compress();
      band_system.reinit(preconditioner_sparsity_pattern);
      is_preconditioner_initialized = true;
    }
  else
    band_system = 0;


  for (int i=0; (unsigned int) i<comp_dom.dh.n_dofs(); ++i)
    {
      if (this_cpu_set.is_element((unsigned int) i))
        {
          if (constraints.is_constrained(i))
            band_system.add((unsigned int)i,(unsigned int) i, 1);
          for (int j=std::max(i-preconditioner_band/2,0); j<std::min(i+preconditioner_band/2,(int)comp_dom.dh.n_dofs()); ++j)
            {
              if (constraints.is_constrained(i) == false)
                {
                  if (surface_nodes(i) == 0)
                    {
                      // Nodo di Dirichlet
                      band_system.add(i,j,neumann_matrix(i,j));

                      if (i == j)
                        band_system.add((unsigned int) i,(unsigned int) j, alpha(i));

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
  Vector<double> phi(glob_phi.size());
  phi = glob_phi;
  Vector<double> dphi_dn(glob_dphi_dn.size());
  dphi_dn = glob_dphi_dn;


  typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;

  SparsityPattern      gradients_sparsity_pattern;
  gradients_sparsity_pattern.reinit(comp_dom.gradient_dh.n_dofs(),
                                    comp_dom.gradient_dh.n_dofs(),
                                    comp_dom.gradient_dh.max_couplings_between_dofs());
  ConstraintMatrix  vector_constraints;
  vector_constraints.clear();
  DoFTools::make_hanging_node_constraints (comp_dom.gradient_dh,vector_constraints);
  vector_constraints.close();
  DoFTools::make_sparsity_pattern (comp_dom.gradient_dh, gradients_sparsity_pattern, vector_constraints);
  gradients_sparsity_pattern.compress();


  SparseMatrix<double> vector_gradients_matrix;
  Vector<double> vector_gradients_rhs;

  vector_gradients_matrix.reinit (gradients_sparsity_pattern);
  vector_gradients_rhs.reinit(comp_dom.gradient_dh.n_dofs());
  Vector<double> vector_gradients_solution(comp_dom.gradient_dh.n_dofs());


  FEValues<dim-1,dim> vector_fe_v(*comp_dom.mapping, comp_dom.gradient_fe, *comp_dom.quadrature,
                                  update_values | update_gradients |
                                  update_cell_normal_vectors |
                                  update_quadrature_points |
                                  update_JxW_values);

  FEValues<dim-1,dim> fe_v(*comp_dom.mapping, comp_dom.fe, *comp_dom.quadrature,
                           update_values | update_gradients |
                           update_cell_normal_vectors |
                           update_quadrature_points |
                           update_JxW_values);

  const unsigned int vector_n_q_points = vector_fe_v.n_quadrature_points;
  const unsigned int   vector_dofs_per_cell   = comp_dom.gradient_fe.dofs_per_cell;
  std::vector<unsigned int> vector_local_dof_indices (vector_dofs_per_cell);

  std::vector< Tensor<1,dim> > phi_surf_grads(vector_n_q_points);
  std::vector<double> phi_norm_grads(vector_n_q_points);
  std::vector<Vector<double> > q_vector_normals_solution(vector_n_q_points,
                                                         Vector<double>(dim));

  FullMatrix<double>   local_gradients_matrix (vector_dofs_per_cell, vector_dofs_per_cell);
  Vector<double>       local_gradients_rhs (vector_dofs_per_cell);


  cell_it
  vector_cell = comp_dom.gradient_dh.begin_active(),
  vector_endc = comp_dom.gradient_dh.end();

  cell_it
  cell = comp_dom.dh.begin_active(),
  endc = comp_dom.dh.end();

  std::vector<Point<dim> > support_points(comp_dom.dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.dh, support_points);
  std::vector<unsigned int> face_dofs(comp_dom.fe.dofs_per_face);

  Quadrature <dim-1> dummy_quadrature(comp_dom.fe.get_unit_support_points());
  FEValues<dim-1,dim> dummy_fe_v(*comp_dom.mapping, comp_dom.fe, dummy_quadrature,
                                 update_values | update_gradients |
                                 update_cell_normal_vectors |
                                 update_quadrature_points);

  const unsigned int   dofs_per_cell = comp_dom.fe.dofs_per_cell;
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  const unsigned int n_q_points = dummy_fe_v.n_quadrature_points;
  std::vector< Tensor<1,dim> > dummy_phi_surf_grads(n_q_points);


  for (; cell!=endc; ++cell,++vector_cell)
    {
      Assert(cell->index() == vector_cell->index(), ExcInternalError());

      fe_v.reinit (cell);
      vector_fe_v.reinit (vector_cell);
      local_gradients_matrix = 0;
      local_gradients_rhs = 0;
      const std::vector<Point<dim> > &vector_node_normals = vector_fe_v.get_normal_vectors();
      fe_v.get_function_gradients(phi, phi_surf_grads);
      fe_v.get_function_values(dphi_dn, phi_norm_grads);
      unsigned int comp_i, comp_j;




      for (unsigned int q=0; q<vector_n_q_points; ++q)
        {
          Point<dim> node_normal_grad_dir(q_vector_normals_solution[q](0),
                                          q_vector_normals_solution[q](1),
                                          q_vector_normals_solution[q](2));
          Point<dim> gradient = vector_node_normals[q]*phi_norm_grads[q] + phi_surf_grads[q];
          for (unsigned int i=0; i<vector_dofs_per_cell; ++i)
            {
              comp_i = comp_dom.gradient_fe.system_to_component_index(i).first;
              for (unsigned int j=0; j<vector_dofs_per_cell; ++j)
                {
                  comp_j = comp_dom.gradient_fe.system_to_component_index(j).first;
                  if (comp_i == comp_j)
                    {
                      local_gradients_matrix(i,j) += vector_fe_v.shape_value(i,q)*
                                                     vector_fe_v.shape_value(j,q)*
                                                     vector_fe_v.JxW(q);
                    }
                }
              local_gradients_rhs(i) += (vector_fe_v.shape_value(i, q)) *
                                        gradient(comp_i) * vector_fe_v.JxW(q);
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

  SparseDirectUMFPACK gradients_inverse;
  gradients_inverse.initialize(vector_gradients_matrix);
  gradients_inverse.vmult(vector_gradients_solution, vector_gradients_rhs);
  vector_constraints.distribute(vector_gradients_solution);


  node_gradients.resize(comp_dom.dh.n_dofs());

  for (unsigned int i=0; i<comp_dom.gradient_dh.n_dofs()/dim; ++i)
    {
      for (unsigned int d=0; d<dim; d++)
        node_gradients[i](d) = vector_gradients_solution(3*i+d);
      //cout<<i<<" Gradient: "<<node_gradients[i]<<endl;
    }


}



template <int dim>
void BEMProblem<dim>::compute_surface_gradients(const TrilinosWrappers::MPI::Vector &tmp_rhs)
{
  Vector<double> phi(tmp_rhs.size());
  phi = tmp_rhs;
  phi.scale(comp_dom.surface_nodes);


  typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;

  SparsityPattern      gradients_sparsity_pattern;
  gradients_sparsity_pattern.reinit(comp_dom.gradient_dh.n_dofs(),
                                    comp_dom.gradient_dh.n_dofs(),
                                    comp_dom.gradient_dh.max_couplings_between_dofs());
  ConstraintMatrix  vector_constraints;
  vector_constraints.clear();
  DoFTools::make_hanging_node_constraints (comp_dom.gradient_dh,vector_constraints);
  vector_constraints.close();
  DoFTools::make_sparsity_pattern (comp_dom.gradient_dh, gradients_sparsity_pattern, vector_constraints);
  gradients_sparsity_pattern.compress();


  SparseMatrix<double> vector_gradients_matrix;
  Vector<double> vector_gradients_rhs;

  vector_gradients_matrix.reinit (gradients_sparsity_pattern);
  vector_gradients_rhs.reinit(comp_dom.gradient_dh.n_dofs());
  Vector<double> vector_gradients_solution(comp_dom.gradient_dh.n_dofs());


  FEValues<dim-1,dim> vector_fe_v(*comp_dom.mapping, comp_dom.gradient_fe, *comp_dom.quadrature,
                                  update_values | update_gradients |
                                  update_cell_normal_vectors |
                                  update_quadrature_points |
                                  update_JxW_values);

  FEValues<dim-1,dim> fe_v(*comp_dom.mapping, comp_dom.fe, *comp_dom.quadrature,
                           update_values | update_gradients |
                           update_cell_normal_vectors |
                           update_quadrature_points |
                           update_JxW_values);

  const unsigned int vector_n_q_points = vector_fe_v.n_quadrature_points;
  const unsigned int   vector_dofs_per_cell   = comp_dom.gradient_fe.dofs_per_cell;
  std::vector<unsigned int> vector_local_dof_indices (vector_dofs_per_cell);

  std::vector< Tensor<1,dim> > phi_surf_grads(vector_n_q_points);
  std::vector<double> phi_norm_grads(vector_n_q_points);
  std::vector<Vector<double> > q_vector_normals_solution(vector_n_q_points,
                                                         Vector<double>(dim));

  FullMatrix<double>   local_gradients_matrix (vector_dofs_per_cell, vector_dofs_per_cell);
  Vector<double>       local_gradients_rhs (vector_dofs_per_cell);


  cell_it
  vector_cell = comp_dom.gradient_dh.begin_active(),
  vector_endc = comp_dom.gradient_dh.end();

  cell_it
  cell = comp_dom.dh.begin_active(),
  endc = comp_dom.dh.end();

  std::vector<Point<dim> > support_points(comp_dom.dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.dh, support_points);
  std::vector<unsigned int> face_dofs(comp_dom.fe.dofs_per_face);

  Quadrature <dim-1> dummy_quadrature(comp_dom.fe.get_unit_support_points());
  FEValues<dim-1,dim> dummy_fe_v(*comp_dom.mapping, comp_dom.fe, dummy_quadrature,
                                 update_values | update_gradients |
                                 update_cell_normal_vectors |
                                 update_quadrature_points);

  const unsigned int   dofs_per_cell = comp_dom.fe.dofs_per_cell;
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  const unsigned int n_q_points = dummy_fe_v.n_quadrature_points;
  std::vector< Tensor<1,dim> > dummy_phi_surf_grads(n_q_points);


  for (; cell!=endc; ++cell,++vector_cell)
    {
      Assert(cell->index() == vector_cell->index(), ExcInternalError());

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
              comp_i = comp_dom.gradient_fe.system_to_component_index(i).first;
              for (unsigned int j=0; j<vector_dofs_per_cell; ++j)
                {
                  comp_j = comp_dom.gradient_fe.system_to_component_index(j).first;
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

  SparseDirectUMFPACK gradients_inverse;
  gradients_inverse.initialize(vector_gradients_matrix);
  gradients_inverse.vmult(vector_gradients_solution, vector_gradients_rhs);
  vector_constraints.distribute(vector_gradients_solution);


  node_surface_gradients.resize(comp_dom.dh.n_dofs());

  for (unsigned int i=0; i<comp_dom.gradient_dh.n_dofs()/dim; ++i)
    {
      for (unsigned int d=0; d<dim; d++)
        node_surface_gradients[i](d) = vector_gradients_solution(3*i+d);
      //cout<<i<<" Gradient: "<<node_gradients[i]<<endl;
    }


}


template class BEMProblem<3>;
