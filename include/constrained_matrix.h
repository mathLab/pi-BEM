// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the pi-BEM authors.
//
// This file is part of the pi-BEM library.
//
// The pi-BEM is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License version 2.1 as published by the Free Software Foundation.
// The full text of the license can be found in the file LICENSE at
// the top level of the pi-BEM distribution.
//
// Authors: Nicola Giuliani, Andrea Mola, Luca Heltai

#ifndef __deal2__constrained_matrix_h
#  define __deal2__constrained_matrix_h



#  include <deal.II/base/config.h>

#  include <deal.II/base/memory_consumption.h>
#  include <deal.II/base/smartpointer.h>
#  include <deal.II/base/thread_management.h>
#  include <deal.II/base/types.h>

#  include <deal.II/lac/affine_constraints.h>
#  include <deal.II/lac/vector_memory.h>

#  include <algorithm>
#  include <vector>

DEAL_II_NAMESPACE_OPEN

template <typename number>
class Vector;
template <class VECTOR>
class ConstrainedMatrixBlock;


/*! @addtogroup Matrix2
 *@{
 */


/**
 *
 * @author Luca Heltai 2011
 */

template <class VEC, class MATRIX>
class ConstrainedOperator
{
public:
  ConstrainedOperator(const MATRIX &                   m,
                      const AffineConstraints<double> &c,
                      const IndexSet &                 c_cpu_set,
                      MPI_Comm                         comm = MPI_COMM_WORLD)
    : constraints(c)
    , matrix(m)
    , constr_cpu_set(c_cpu_set)
    , mpi_communicator(comm)
    , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
    , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  {}



  void
  vmult(VEC &dst, const VEC &src) const;

  void
  distribute_rhs(VEC &rhs) const;

private:
  const AffineConstraints<double> &constraints;
  const MATRIX &                   matrix;
  const IndexSet                   constr_cpu_set;
  MPI_Comm                         mpi_communicator;
  unsigned int                     n_mpi_processes;
  unsigned int                     this_mpi_process;
};

/*@}*/
/*---------------------- Inline functions -----------------------------------*/


//--------------------------------Iterators--------------------------------------//


template <class VEC, class MATRIX>
void
ConstrainedOperator<VEC, MATRIX>::vmult(VEC &dst, const VEC &src) const
{
  // Vector<double> loc_src(src);
  VEC loc_src(constr_cpu_set);
  loc_src.reinit(src, false, true);
  // loc_src = src;

  // std::cout<<"in vector "<<std::endl;
  // for (unsigned int i = 0; i < src.size(); i++)
  //     if (src.locally_owned_elements().is_element(i))
  //        std::cout<<i<<" ("<<this_mpi_process<<")  "<<loc_src(i)<<std::endl;

  matrix.vmult(dst, src);
  // dst.compress(VectorOperation::insert);
  // IndexSet dummy(dst.locally_owned_elements());

  for (unsigned int i = 0; i < src.size(); ++i)
    if ((constraints.is_constrained(i)) &&
        (src.locally_owned_elements().is_element(i)))
      {
        // dst(i) -= dst(i);
        dst(i) += src(i) - dst(i);
        const std::vector<std::pair<types::global_dof_index, double>> *entries =
          constraints.get_constraint_entries(i);
        for (unsigned int j = 0; j < entries->size(); ++j)
          dst(i) -= (*entries)[j].second * loc_src((*entries)[j].first);
      }
    else
      dst(i) += 0.;

  dst.compress(VectorOperation::add);
  // constraints.condense(dst);
  // std::cout<<"out vector "<<std::endl;
  // for (unsigned int i = 0; i < dst.size(); i++)
  //     if (dst.locally_owned_elements().is_element(i))
  //        std::cout<<i<<" ("<<this_mpi_process<<")  "<<dst(i)<<std::endl;

  // std::cout<<"check vector "<<std::endl;
  // for (unsigned int i = 0; i < dst.size(); i++)
  //     if (dst.locally_owned_elements().is_element(i) &&
  //     !(dummy.is_element(i)))
  //        std::cout<<i<<" ("<<this_mpi_process<<")  "<<dst(i)<<std::endl;
}

template <class VEC, class MATRIX>
void
ConstrainedOperator<VEC, MATRIX>::distribute_rhs(VEC &rhs) const
{
  for (auto i : rhs.locally_owned_elements())
    if ((constraints.is_constrained(i)) &&
        (rhs.locally_owned_elements().is_element(i)))
      rhs(i) = constraints.get_inhomogeneity(i);
    else
      rhs(i) = rhs(i);
}


DEAL_II_NAMESPACE_CLOSE

#endif
/*----------------------------   filtered_matrix.h ---------------------------*/
