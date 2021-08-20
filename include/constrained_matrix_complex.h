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
#define __deal2__constrained_matrix_complex_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/types.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector_memory.h>

#include <algorithm>
#include <vector>

DEAL_II_NAMESPACE_OPEN

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
  IndexSet                         constr_cpu_set, constr_cpu_set_complex;
  // this class is meant to be const when used by GMRES, but these caches are
  // kept allocated
  mutable VEC  c_dst, c_dst_imag, c_src, c_src_imag;
  MPI_Comm     mpi_communicator;
  unsigned int n_mpi_processes;
  unsigned int this_mpi_process;
};

template <class VEC, class MATRIX>
void
ConstrainedComplexOperator<VEC, MATRIX>::vmult(VEC &dst, const VEC &src) const
{
  auto offset = matrix.this_cpu_set.size();
  // store localized constrained values
  VEC loc_src(constr_cpu_set_complex);
  loc_src.reinit(src, false, true);

  // split into real, imag
  c_src      = 0;
  c_src_imag = 0;
  c_dst      = 0;
  c_dst_imag = 0;
  for (auto i : matrix.this_cpu_set)
    {
      c_src(i)      = src(i);
      c_src_imag(i) = src(i + offset);
    }
  c_src.compress(VectorOperation::add);
  c_src_imag.compress(VectorOperation::add);

  matrix.vmult(c_dst, c_dst_imag, c_src, c_src_imag);
  c_dst.compress(VectorOperation::add);
  c_dst_imag.compress(VectorOperation::add);

  for (auto i : matrix.this_cpu_set)
    {
      dst(i)          = c_dst(i);
      dst(i + offset) = c_dst_imag(i);
    }

  for (auto i : matrix.this_cpu_set)
    {
      if (constraints.is_constrained(i))
        {
          dst(i) += src(i) - dst(i);
          dst(i + offset) += src(i + offset) - dst(i + offset);

          for (const auto &entry : *constraints.get_constraint_entries(i))
            {
              dst(i) -= entry.second * loc_src(entry.first);
            }

          for (const auto &entry : *constraints_imag.get_constraint_entries(i))
            {
              dst(i + offset) -= entry.second * loc_src(entry.first + offset);
            }
        }
      else
        {
          dst(i) += 0.;
          dst(i + offset) += 0.;
        }
    }

  dst.compress(VectorOperation::add);
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

// using BlockVector
template <class PreconditionerType>
class BlockPreconditioner
{
public:
  BlockPreconditioner(PreconditionerType &preconditioner)
    : preconditioner(preconditioner)
  {}

  virtual void
  vmult(TrilinosWrappers::MPI::BlockVector &      dst,
        const TrilinosWrappers::MPI::BlockVector &src) const
  {
    for (unsigned int i = 0; i < src.n_blocks(); ++i)
      {
        preconditioner.vmult(dst.block(i), src.block(i));
      }
  }

  virtual void
  Tvmult(TrilinosWrappers::MPI::BlockVector &      dst,
         const TrilinosWrappers::MPI::BlockVector &src) const
  {
    for (unsigned int i = 0; i < src.n_blocks(); ++i)
      {
        preconditioner.Tvmult(dst.block(i), src.block(i));
      }
  }

protected:
  PreconditionerType &preconditioner;
};

template <class VEC, class MATRIX>
class BlockConstrainedComplexOperator
{
public:
  BlockConstrainedComplexOperator(const MATRIX &                   m,
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

    // c_dst.reinit(matrix.this_cpu_set);
    // c_dst_imag.reinit(matrix.this_cpu_set);
    // c_src.reinit(matrix.this_cpu_set);
    // c_src_imag.reinit(matrix.this_cpu_set);
  }

  void
  vmult(VEC &dst, const VEC &src) const;

  void
  distribute_rhs(VEC &rhs) const;

private:
  const AffineConstraints<double> &constraints, constraints_imag;
  const MATRIX &                   matrix;
  IndexSet                         constr_cpu_set, constr_cpu_set_complex;
  // this class is meant to be const when used by GMRES, but these caches are
  // kept allocated
  // mutable VEC  c_dst, c_dst_imag, c_src, c_src_imag;
  MPI_Comm     mpi_communicator;
  unsigned int n_mpi_processes;
  unsigned int this_mpi_process;
};

template <class VEC, class MATRIX>
void
BlockConstrainedComplexOperator<VEC, MATRIX>::vmult(VEC &      dst,
                                                    const VEC &src) const
{
  // store localized constrained values
  std::vector<IndexSet> parallel_partitioning{constr_cpu_set, constr_cpu_set};
  VEC                   loc_src(parallel_partitioning);
  // loc_src.block(0).reinit(src.block(0), false, true);
  // loc_src.block(1).reinit(src.block(1), false, true);
  for (unsigned int i = 0; i < src.n_blocks(); ++i)
    {
      loc_src.block(i).reinit(src.block(i), false, true);
    }

  // split into real, imag
  // c_src      = 0;
  // c_src_imag = 0;
  // c_dst      = 0;
  // c_dst_imag = 0;
  // for (auto i : matrix.this_cpu_set)
  //   {
  //     c_src(i)      = src(i);
  //     c_src_imag(i) = src(i + offset);
  //   }
  // c_src.compress(VectorOperation::add);
  // c_src_imag.compress(VectorOperation::add);

  matrix.vmult(dst.block(0), dst.block(1), src.block(0), src.block(1));
  // c_dst.compress(VectorOperation::add);
  // c_dst_imag.compress(VectorOperation::add);

  // for (auto i : matrix.this_cpu_set)
  //   {
  //     dst(i)          = c_dst(i);
  //     dst(i + offset) = c_dst_imag(i);
  //   }

  for (auto i : matrix.this_cpu_set)
    {
      if (constraints.is_constrained(i))
        {
          dst.block(0)(i) += src.block(0)(i) - dst.block(0)(i);
          dst.block(1)(i) += src.block(1)(i) - dst.block(1)(i);

          for (const auto &entry : *constraints.get_constraint_entries(i))
            {
              dst.block(0)(i) -= entry.second * loc_src.block(0)(entry.first);
            }

          for (const auto &entry : *constraints_imag.get_constraint_entries(i))
            {
              dst.block(1)(i) -= entry.second * loc_src.block(1)(entry.first);
            }
        }
      else
        {
          // noop
          dst.block(0)(i) = dst.block(0)(i);
          dst.block(1)(i) = dst.block(1)(i);
        }
    }

  dst.compress(VectorOperation::add);
}

template <class VEC, class MATRIX>
void
BlockConstrainedComplexOperator<VEC, MATRIX>::distribute_rhs(VEC &rhs) const
{
  for (auto i : matrix.this_cpu_set)
    {
      if (constraints.is_constrained(i))
        {
          rhs.block(0)(i) = constraints.get_inhomogeneity(i);
          rhs.block(1)(i) = constraints_imag.get_inhomogeneity(i);
        }
      else
        {
          // noop
          rhs.block(0)(i) = rhs.block(0)(i);
          rhs.block(1)(i) = rhs.block(1)(i);
        }
    }

  rhs.compress(VectorOperation::insert);
}

DEAL_II_NAMESPACE_CLOSE

#endif
