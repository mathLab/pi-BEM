//---------------------------------------------------------------------------
//    $Id: filtered_matrix.h 23248 2011-01-23 06:03:57Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__constrained_matrix_h
#define __deal2__constrained_matrix_h



#include<deal.II/base/config.h>
#include<deal.II/base/smartpointer.h>
#include<deal.II/base/thread_management.h>
#include<deal.II/base/memory_consumption.h>
#include<deal.II/lac/pointer_matrix.h>
#include<deal.II/lac/constraint_matrix.h>
#include<deal.II/lac/vector_memory.h>
#include <vector>
#include <algorithm>

DEAL_II_NAMESPACE_OPEN

template <typename number> class Vector;
template <class VECTOR> class ConstrainedMatrixBlock;


/*! @addtogroup Matrix2
 *@{
 */


/**
 *
 * @author Luca Heltai 2011
 */

template<class VEC, class MATRIX>
class ConstrainedOperator
{
public:
  ConstrainedOperator(const MATRIX &m,
                      const ConstraintMatrix &c,
                      const IndexSet &c_cpu_set,
                      MPI_Comm comm = MPI_COMM_WORLD) :
    cm(c),
    matrix(m),
    constr_cpu_set(c_cpu_set),
    mpi_communicator (comm),
    n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
    this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator))
  {}


  // void compute_constraint_transpose(const unsigned int problem_size);

  void constraint_vmult(VEC &dst, const VEC &src) const;

  void constraint_tvmult(VEC &dst, const VEC &src) const;

  void vmult(VEC &dst, const VEC &src) const;

  void distribute_rhs(VEC &rhs) const;

  void finalise(VEC &u) const;

private:
  const ConstraintMatrix &cm;
  const MATRIX &matrix;
  const IndexSet constr_cpu_set;
  ConstraintMatrix cmt;
  MPI_Comm mpi_communicator;
  unsigned int n_mpi_processes;
  unsigned int this_mpi_process;
};

/*@}*/
/*---------------------- Inline functions -----------------------------------*/


//--------------------------------Iterators--------------------------------------//



// template<class VEC, class MATRIX>
// void ConstrainedOperator<VEC,MATRIX>::compute_constraint_transpose(const unsigned int problem_size)
// {
//   for(unsigned int i=0; i<problem_size; ++i)
//     if(cm.is_constrained(i))
//     {
//
//     }
// }


template<class VEC, class MATRIX>
void ConstrainedOperator<VEC,MATRIX>::constraint_vmult(VEC &dst, const VEC &src) const
{
  VEC loc_src(src.locally_owned_elements(), constr_cpu_set, mpi_communicator);
  loc_src.reinit(src,false,true);

  for (unsigned int i=0; i<dst.size(); ++i)
    {
      if (dst.locally_owned_elements().is_element(i))
        {
          if (cm.is_constrained(i))
            {
              dst(i) = 0;
              const std::vector< std::pair < unsigned int, double > >
              *entries = cm.get_constraint_entries (i);
              for (unsigned int j=0; j < entries->size(); ++j)
                {
                  unsigned int pos = (*entries)[j].first;
                  dst(i) +=  loc_src(pos) * (*entries)[j].second;
                }
            }
          else
            dst(i)=loc_src(i);
        }
    }

}

template<class VEC, class MATRIX>
void ConstrainedOperator<VEC,MATRIX>::constraint_tvmult(VEC &dst, const VEC &src) const
{

  VEC loc_dst(constr_cpu_set);
  for (unsigned int i=0; i<src.size(); ++i)
    {
      if (src.locally_owned_elements().is_element(i))
        {
          if (cm.is_constrained(i))
            {
              loc_dst(i)=0;
              const std::vector< std::pair < unsigned int, double > >
              *entries = cm.get_constraint_entries (i);
              for (unsigned int j=0; j < entries->size(); ++j)
                {
                  unsigned int pos = (*entries)[j].first;
                  loc_dst(pos) += src(i) * (*entries)[j].second;
                }
            }
          else
            loc_dst(i)=src(i);
        }
    }

  std::cout<<"Tvmult comm"<<std::endl;
  dst.add(loc_dst, true);

}

template<class VEC, class MATRIX>
void ConstrainedOperator<VEC,MATRIX>::distribute_rhs(VEC &rhs) const
{
  VEC constraint_term(rhs.locally_owned_elements(), mpi_communicator);
  VEC rhs_tmp(rhs.locally_owned_elements(), mpi_communicator);

  for (unsigned int i=0; i<rhs.size();  ++i)
    constraint_term(i)=-cm.get_inhomogeneity(i);

  matrix.vmult(rhs_tmp, constraint_term);
  constraint_term = rhs;
  constraint_term.add(rhs_tmp);
  rhs = 0.;
  constraint_tvmult(rhs,constraint_term);

}

template<class VEC, class MATRIX>
void ConstrainedOperator<VEC,MATRIX>::vmult(VEC &dst, const VEC &src) const
{


  auto dst_tmp = dst;
  auto src_tmp = src;

  constraint_vmult(src_tmp,src);
  matrix.vmult(dst_tmp,src_tmp);
  constraint_tvmult(dst,dst_tmp);


}

template<class VEC, class MATRIX>
void ConstrainedOperator<VEC,MATRIX>::finalise(VEC &u) const
{
  cm.distribute(u);
}
// template<class VEC, class MATRIX>
// void ConstrainedOperator<VEC,MATRIX>::distribute_rhs(VEC &rhs) const
// {
//   for (unsigned int i=0; i<rhs.size(); ++i)
//     if ( (constraints.is_constrained(i)) &&
//          (rhs.locally_owned_elements().is_element(i)) )
//       rhs(i) = constraints.get_inhomogeneity(i);
// }


DEAL_II_NAMESPACE_CLOSE

#endif
/*----------------------------   filtered_matrix.h     ---------------------------*/
