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

#ifndef __deal2__constrained_matrix_complex_h
#  define __deal2__constrained_matrix_complex_h

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
class ConstrainedComplexOperator
{
public:
  ConstrainedComplexOperator(const MATRIX &                   m,
                             const AffineConstraints<double> &c,
                             const AffineConstraints<double> &c_imag,
                             const IndexSet &                 c_cpu_set,
                             MPI_Comm comm = MPI_COMM_WORLD)
    : constraints(c)
    , constraints_imag(c_imag)
    , matrix(m)
    , constr_cpu_set(c_cpu_set)
    , mpi_communicator(comm)
    , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
    , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  {
    constr_cpu_set_complex = IndexSet(2 * constr_cpu_set.size());
    constr_cpu_set_complex.add_indices(constr_cpu_set);
    constr_cpu_set_complex.add_indices(constr_cpu_set, constr_cpu_set.size());
    constr_cpu_set_complex.compress();

    this_cpu_set_imag = IndexSet(2 * matrix.this_cpu_set.size());
    this_cpu_set_imag.add_indices(matrix.this_cpu_set,
                                  matrix.this_cpu_set.size());
    this_cpu_set_imag.compress();

    c_dst.reinit(matrix.this_cpu_set);
    c_dst_imag.reinit(matrix.this_cpu_set);
    c_src.reinit(matrix.this_cpu_set);
    c_src_imag.reinit(matrix.this_cpu_set);
  }

  void
  vmult(VEC &dst, const VEC &src) const;

  void
  distribute_rhs(VEC &rhs) const;

private:
  const AffineConstraints<double> &constraints, constraints_imag;
  const MATRIX &                   matrix;
  IndexSet constr_cpu_set, constr_cpu_set_complex, this_cpu_set_imag;
  // this class is meant to be const when used by GMRES, but these caches are
  // kept allocated
  mutable VEC  c_dst, c_dst_imag, c_src, c_src_imag;
  MPI_Comm     mpi_communicator;
  unsigned int n_mpi_processes;
  unsigned int this_mpi_process;
};

/*@}*/
/*---------------------- Inline functions -----------------------------------*/


//--------------------------------Iterators--------------------------------------//


template <class VEC, class MATRIX>
void
ConstrainedComplexOperator<VEC, MATRIX>::vmult(VEC &dst, const VEC &src) const
{
  auto imag_offset = matrix.this_cpu_set.size();
  // store localized constrained values
  VEC loc_src(constr_cpu_set_complex);
  loc_src.reinit(src, false, true);

  // std::cout<<"in vector "<<std::endl;
  // for (unsigned int i = 0; i < src.size(); i++)
  //     if (src.locally_owned_elements().is_element(i))
  //        std::cout<<i<<" ("<<this_mpi_process<<")  "<<loc_src(i)<<std::endl;

  // split into real, imag
  c_src      = 0;
  c_src_imag = 0;
  c_dst      = 0;
  c_dst_imag = 0;
  for (auto i : matrix.this_cpu_set)
    {
      c_src(i)      = src(i);
      c_src_imag(i) = src(i + imag_offset);
    }
  c_src.compress(VectorOperation::add);
  c_src_imag.compress(VectorOperation::add);

  matrix.vmult(c_dst, c_dst_imag, c_src, c_src_imag);
  c_dst.compress(VectorOperation::add);
  c_dst_imag.compress(VectorOperation::add);

  for (auto i : matrix.this_cpu_set)
    {
      dst(i)               = c_dst(i);
      dst(i + imag_offset) = c_dst_imag(i);
    }

  matrix.vmult(c_dst, c_src);
  c_dst.compress(VectorOperation::add);

  // for (unsigned int i = 0; i < src.size(); ++i)
  for (auto i : matrix.this_cpu_set)
    {
      // if ((constraints.is_constrained(i)) &&
      //     (src.locally_owned_elements().is_element(i)))
      if (constraints.is_constrained(i))
        {
          // dst(i) -= dst(i);
          dst(i) += src(i) - dst(i);
          dst(i + imag_offset) += src(i + imag_offset) - dst(i + imag_offset);

          for (const auto &entry : *constraints.get_constraint_entries(i))
            {
              dst(i) -= entry.second * loc_src(entry.first);
            }
          // can the coefficients be different? or, this is just a repeat of the
          // first loop with an offset?
          for (const auto &entry : *constraints_imag.get_constraint_entries(i))
            {
              dst(i + imag_offset) -=
                entry.second * loc_src(entry.first + imag_offset);
            }
          // const std::vector<std::pair<types::global_dof_index, double>>
          //   *entries = constraints.get_constraint_entries(i);
          // for (unsigned int j = 0; j < entries->size(); ++j)
          //   dst(i) -= (*entries)[j].second * loc_src((*entries)[j].first);
        }
      else
        {
          dst(i) += 0.;
          dst(i + imag_offset) += 0.;
        }
    }

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
ConstrainedComplexOperator<VEC, MATRIX>::distribute_rhs(VEC &rhs) const
{
  for (auto i : matrix.this_cpu_set)
    {
      if (constraints.is_constrained(i))
        {
          rhs(i) = constraints.get_inhomogeneity(i);
          rhs(i + matrix.this_cpu_set.size()) =
            constraints_imag.get_inhomogeneity(i);
        }
      else
        {
          rhs(i) = rhs(i);
          rhs(i + matrix.this_cpu_set.size()) =
            rhs(i + matrix.this_cpu_set.size());
        }
    }
}


DEAL_II_NAMESPACE_CLOSE

#endif
/*----------------------------   filtered_matrix.h ---------------------------*/
