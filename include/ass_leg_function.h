//----------------------------  step-34.cc  ---------------------------
//    $Id: step-34.cc 18734 2009-04-25 13:36:48Z heltai $
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors
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

#ifndef __polimi__ass_leg_function_h
#define __polimi__ass_leg_function_h

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/matrix_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1_eulerian.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <vector>

using namespace dealii;



class AssLegFunction
{
public:

  AssLegFunction();

  ~AssLegFunction();


  void AssLegFunSph(const unsigned int p, const unsigned int m, double x, double y[]);

  void AssLegFunSphDeriv(const unsigned int p, const unsigned int m, double x, double y[], double dy[]);

  double GetAssLegFunSph(const unsigned int n, const unsigned int m, double x) const;

  double GetAssLegFunSphDeriv(const unsigned int n, const unsigned int m, double x) const;


private:

  std::vector <std::vector <double (AssLegFunction::*)(const double ) const> > leg_pointers;
  std::vector <std::vector <double (AssLegFunction::*)(const double ) const> > leg_der_pointers;

  /** Declaration of each of the functions. */
  double P_0_0 (const double x) const;
  double P_0_0_Deriv (const double x) const;
  double P_1_0 (const double x) const;
  double P_1_0_Deriv (const double x) const;
  double P_1_1 (const double x) const;
  double P_1_1_Deriv (const double x) const;
  double P_2_0 (const double x) const;
  double P_2_0_Deriv (const double x) const;
  double P_2_1 (const double x) const;
  double P_2_1_Deriv (const double x) const;
  double P_2_2 (const double x) const;
  double P_2_2_Deriv (const double x) const;
  double P_3_0 (const double x) const;
  double P_3_0_Deriv (const double x) const;
  double P_3_1 (const double x) const;
  double P_3_1_Deriv (const double x) const;
  double P_3_2 (const double x) const;
  double P_3_2_Deriv (const double x) const;
  double P_3_3 (const double x) const;
  double P_3_3_Deriv (const double x) const;
  double P_4_0 (const double x) const;
  double P_4_0_Deriv (const double x) const;
  double P_4_1 (const double x) const;
  double P_4_1_Deriv (const double x) const;
  double P_4_2 (const double x) const;
  double P_4_2_Deriv (const double x) const;
  double P_4_3 (const double x) const;
  double P_4_3_Deriv (const double x) const;
  double P_4_4 (const double x) const;
  double P_4_4_Deriv (const double x) const;
  double P_5_0 (const double x) const;
  double P_5_0_Deriv (const double x) const;
  double P_5_1 (const double x) const;
  double P_5_1_Deriv (const double x) const;
  double P_5_2 (const double x) const;
  double P_5_2_Deriv (const double x) const;
  double P_5_3 (const double x) const;
  double P_5_3_Deriv (const double x) const;
  double P_5_4 (const double x) const;
  double P_5_4_Deriv (const double x) const;
  double P_5_5 (const double x) const;
  double P_5_5_Deriv (const double x) const;
  double P_6_0 (const double x) const;
  double P_6_0_Deriv (const double x) const;
  double P_6_1 (const double x) const;
  double P_6_1_Deriv (const double x) const;
  double P_6_2 (const double x) const;
  double P_6_2_Deriv (const double x) const;
  double P_6_3 (const double x) const;
  double P_6_3_Deriv (const double x) const;
  double P_6_4 (const double x) const;
  double P_6_4_Deriv (const double x) const;
  double P_6_5 (const double x) const;
  double P_6_5_Deriv (const double x) const;
  double P_6_6 (const double x) const;
  double P_6_6_Deriv (const double x) const;
  double P_7_0 (const double x) const;
  double P_7_0_Deriv (const double x) const;
  double P_7_1 (const double x) const;
  double P_7_1_Deriv (const double x) const;
  double P_7_2 (const double x) const;
  double P_7_2_Deriv (const double x) const;
  double P_7_3 (const double x) const;
  double P_7_3_Deriv (const double x) const;
  double P_7_4 (const double x) const;
  double P_7_4_Deriv (const double x) const;
  double P_7_5 (const double x) const;
  double P_7_5_Deriv (const double x) const;
  double P_7_6 (const double x) const;
  double P_7_6_Deriv (const double x) const;
  double P_7_7 (const double x) const;
  double P_7_7_Deriv (const double x) const;
  double P_8_0 (const double x) const;
  double P_8_0_Deriv (const double x) const;
  double P_8_1 (const double x) const;
  double P_8_1_Deriv (const double x) const;
  double P_8_2 (const double x) const;
  double P_8_2_Deriv (const double x) const;
  double P_8_3 (const double x) const;
  double P_8_3_Deriv (const double x) const;
  double P_8_4 (const double x) const;
  double P_8_4_Deriv (const double x) const;
  double P_8_5 (const double x) const;
  double P_8_5_Deriv (const double x) const;
  double P_8_6 (const double x) const;
  double P_8_6_Deriv (const double x) const;
  double P_8_7 (const double x) const;
  double P_8_7_Deriv (const double x) const;
  double P_8_8 (const double x) const;
  double P_8_8_Deriv (const double x) const;
  double P_9_0 (const double x) const;
  double P_9_0_Deriv (const double x) const;
  double P_9_1 (const double x) const;
  double P_9_1_Deriv (const double x) const;
  double P_9_2 (const double x) const;
  double P_9_2_Deriv (const double x) const;
  double P_9_3 (const double x) const;
  double P_9_3_Deriv (const double x) const;
  double P_9_4 (const double x) const;
  double P_9_4_Deriv (const double x) const;
  double P_9_5 (const double x) const;
  double P_9_5_Deriv (const double x) const;
  double P_9_6 (const double x) const;
  double P_9_6_Deriv (const double x) const;
  double P_9_7 (const double x) const;
  double P_9_7_Deriv (const double x) const;
  double P_9_8 (const double x) const;
  double P_9_8_Deriv (const double x) const;
  double P_9_9 (const double x) const;
  double P_9_9_Deriv (const double x) const;
  double P_10_0 (const double x) const;
  double P_10_0_Deriv (const double x) const;
  double P_10_1 (const double x) const;
  double P_10_1_Deriv (const double x) const;
  double P_10_2 (const double x) const;
  double P_10_2_Deriv (const double x) const;
  double P_10_3 (const double x) const;
  double P_10_3_Deriv (const double x) const;
  double P_10_4 (const double x) const;
  double P_10_4_Deriv (const double x) const;
  double P_10_5 (const double x) const;
  double P_10_5_Deriv (const double x) const;
  double P_10_6 (const double x) const;
  double P_10_6_Deriv (const double x) const;
  double P_10_7 (const double x) const;
  double P_10_7_Deriv (const double x) const;
  double P_10_8 (const double x) const;
  double P_10_8_Deriv (const double x) const;
  double P_10_9 (const double x) const;
  double P_10_9_Deriv (const double x) const;
  double P_10_10 (const double x) const;
  double P_10_10_Deriv (const double x) const;
  double P_11_0 (const double x) const;
  double P_11_0_Deriv (const double x) const;
  double P_11_1 (const double x) const;
  double P_11_1_Deriv (const double x) const;
  double P_11_2 (const double x) const;
  double P_11_2_Deriv (const double x) const;
  double P_11_3 (const double x) const;
  double P_11_3_Deriv (const double x) const;
  double P_11_4 (const double x) const;
  double P_11_4_Deriv (const double x) const;
  double P_11_5 (const double x) const;
  double P_11_5_Deriv (const double x) const;
  double P_11_6 (const double x) const;
  double P_11_6_Deriv (const double x) const;
  double P_11_7 (const double x) const;
  double P_11_7_Deriv (const double x) const;
  double P_11_8 (const double x) const;
  double P_11_8_Deriv (const double x) const;
  double P_11_9 (const double x) const;
  double P_11_9_Deriv (const double x) const;
  double P_11_10 (const double x) const;
  double P_11_10_Deriv (const double x) const;
  double P_11_11 (const double x) const;
  double P_11_11_Deriv (const double x) const;
  double P_12_0 (const double x) const;
  double P_12_0_Deriv (const double x) const;
  double P_12_1 (const double x) const;
  double P_12_1_Deriv (const double x) const;
  double P_12_2 (const double x) const;
  double P_12_2_Deriv (const double x) const;
  double P_12_3 (const double x) const;
  double P_12_3_Deriv (const double x) const;
  double P_12_4 (const double x) const;
  double P_12_4_Deriv (const double x) const;
  double P_12_5 (const double x) const;
  double P_12_5_Deriv (const double x) const;
  double P_12_6 (const double x) const;
  double P_12_6_Deriv (const double x) const;
  double P_12_7 (const double x) const;
  double P_12_7_Deriv (const double x) const;
  double P_12_8 (const double x) const;
  double P_12_8_Deriv (const double x) const;
  double P_12_9 (const double x) const;
  double P_12_9_Deriv (const double x) const;
  double P_12_10 (const double x) const;
  double P_12_10_Deriv (const double x) const;
  double P_12_11 (const double x) const;
  double P_12_11_Deriv (const double x) const;
  double P_12_12 (const double x) const;
  double P_12_12_Deriv (const double x) const;
  double P_13_0 (const double x) const;
  double P_13_0_Deriv (const double x) const;
  double P_13_1 (const double x) const;
  double P_13_1_Deriv (const double x) const;
  double P_13_2 (const double x) const;
  double P_13_2_Deriv (const double x) const;
  double P_13_3 (const double x) const;
  double P_13_3_Deriv (const double x) const;
  double P_13_4 (const double x) const;
  double P_13_4_Deriv (const double x) const;
  double P_13_5 (const double x) const;
  double P_13_5_Deriv (const double x) const;
  double P_13_6 (const double x) const;
  double P_13_6_Deriv (const double x) const;
  double P_13_7 (const double x) const;
  double P_13_7_Deriv (const double x) const;
  double P_13_8 (const double x) const;
  double P_13_8_Deriv (const double x) const;
  double P_13_9 (const double x) const;
  double P_13_9_Deriv (const double x) const;
  double P_13_10 (const double x) const;
  double P_13_10_Deriv (const double x) const;
  double P_13_11 (const double x) const;
  double P_13_11_Deriv (const double x) const;
  double P_13_12 (const double x) const;
  double P_13_12_Deriv (const double x) const;
  double P_13_13 (const double x) const;
  double P_13_13_Deriv (const double x) const;
  double P_14_0 (const double x) const;
  double P_14_0_Deriv (const double x) const;
  double P_14_1 (const double x) const;
  double P_14_1_Deriv (const double x) const;
  double P_14_2 (const double x) const;
  double P_14_2_Deriv (const double x) const;
  double P_14_3 (const double x) const;
  double P_14_3_Deriv (const double x) const;
  double P_14_4 (const double x) const;
  double P_14_4_Deriv (const double x) const;
  double P_14_5 (const double x) const;
  double P_14_5_Deriv (const double x) const;
  double P_14_6 (const double x) const;
  double P_14_6_Deriv (const double x) const;
  double P_14_7 (const double x) const;
  double P_14_7_Deriv (const double x) const;
  double P_14_8 (const double x) const;
  double P_14_8_Deriv (const double x) const;
  double P_14_9 (const double x) const;
  double P_14_9_Deriv (const double x) const;
  double P_14_10 (const double x) const;
  double P_14_10_Deriv (const double x) const;
  double P_14_11 (const double x) const;
  double P_14_11_Deriv (const double x) const;
  double P_14_12 (const double x) const;
  double P_14_12_Deriv (const double x) const;
  double P_14_13 (const double x) const;
  double P_14_13_Deriv (const double x) const;
  double P_14_14 (const double x) const;
  double P_14_14_Deriv (const double x) const;
  double P_15_0 (const double x) const;
  double P_15_0_Deriv (const double x) const;
  double P_15_1 (const double x) const;
  double P_15_1_Deriv (const double x) const;
  double P_15_2 (const double x) const;
  double P_15_2_Deriv (const double x) const;
  double P_15_3 (const double x) const;
  double P_15_3_Deriv (const double x) const;
  double P_15_4 (const double x) const;
  double P_15_4_Deriv (const double x) const;
  double P_15_5 (const double x) const;
  double P_15_5_Deriv (const double x) const;
  double P_15_6 (const double x) const;
  double P_15_6_Deriv (const double x) const;
  double P_15_7 (const double x) const;
  double P_15_7_Deriv (const double x) const;
  double P_15_8 (const double x) const;
  double P_15_8_Deriv (const double x) const;
  double P_15_9 (const double x) const;
  double P_15_9_Deriv (const double x) const;
  double P_15_10 (const double x) const;
  double P_15_10_Deriv (const double x) const;
  double P_15_11 (const double x) const;
  double P_15_11_Deriv (const double x) const;
  double P_15_12 (const double x) const;
  double P_15_12_Deriv (const double x) const;
  double P_15_13 (const double x) const;
  double P_15_13_Deriv (const double x) const;
  double P_15_14 (const double x) const;
  double P_15_14_Deriv (const double x) const;
  double P_15_15 (const double x) const;
  double P_15_15_Deriv (const double x) const;
  double P_16_0 (const double x) const;
  double P_16_0_Deriv (const double x) const;
  double P_16_1 (const double x) const;
  double P_16_1_Deriv (const double x) const;
  double P_16_2 (const double x) const;
  double P_16_2_Deriv (const double x) const;
  double P_16_3 (const double x) const;
  double P_16_3_Deriv (const double x) const;
  double P_16_4 (const double x) const;
  double P_16_4_Deriv (const double x) const;
  double P_16_5 (const double x) const;
  double P_16_5_Deriv (const double x) const;
  double P_16_6 (const double x) const;
  double P_16_6_Deriv (const double x) const;
  double P_16_7 (const double x) const;
  double P_16_7_Deriv (const double x) const;
  double P_16_8 (const double x) const;
  double P_16_8_Deriv (const double x) const;
  double P_16_9 (const double x) const;
  double P_16_9_Deriv (const double x) const;
  double P_16_10 (const double x) const;
  double P_16_10_Deriv (const double x) const;
  double P_16_11 (const double x) const;
  double P_16_11_Deriv (const double x) const;
  double P_16_12 (const double x) const;
  double P_16_12_Deriv (const double x) const;
  double P_16_13 (const double x) const;
  double P_16_13_Deriv (const double x) const;
  double P_16_14 (const double x) const;
  double P_16_14_Deriv (const double x) const;
  double P_16_15 (const double x) const;
  double P_16_15_Deriv (const double x) const;
  double P_16_16 (const double x) const;
  double P_16_16_Deriv (const double x) const;
  double P_17_0 (const double x) const;
  double P_17_0_Deriv (const double x) const;
  double P_17_1 (const double x) const;
  double P_17_1_Deriv (const double x) const;
  double P_17_2 (const double x) const;
  double P_17_2_Deriv (const double x) const;
  double P_17_3 (const double x) const;
  double P_17_3_Deriv (const double x) const;
  double P_17_4 (const double x) const;
  double P_17_4_Deriv (const double x) const;
  double P_17_5 (const double x) const;
  double P_17_5_Deriv (const double x) const;
  double P_17_6 (const double x) const;
  double P_17_6_Deriv (const double x) const;
  double P_17_7 (const double x) const;
  double P_17_7_Deriv (const double x) const;
  double P_17_8 (const double x) const;
  double P_17_8_Deriv (const double x) const;
  double P_17_9 (const double x) const;
  double P_17_9_Deriv (const double x) const;
  double P_17_10 (const double x) const;
  double P_17_10_Deriv (const double x) const;
  double P_17_11 (const double x) const;
  double P_17_11_Deriv (const double x) const;
  double P_17_12 (const double x) const;
  double P_17_12_Deriv (const double x) const;
  double P_17_13 (const double x) const;
  double P_17_13_Deriv (const double x) const;
  double P_17_14 (const double x) const;
  double P_17_14_Deriv (const double x) const;
  double P_17_15 (const double x) const;
  double P_17_15_Deriv (const double x) const;
  double P_17_16 (const double x) const;
  double P_17_16_Deriv (const double x) const;
  double P_17_17 (const double x) const;
  double P_17_17_Deriv (const double x) const;
  double P_18_0 (const double x) const;
  double P_18_0_Deriv (const double x) const;
  double P_18_1 (const double x) const;
  double P_18_1_Deriv (const double x) const;
  double P_18_2 (const double x) const;
  double P_18_2_Deriv (const double x) const;
  double P_18_3 (const double x) const;
  double P_18_3_Deriv (const double x) const;
  double P_18_4 (const double x) const;
  double P_18_4_Deriv (const double x) const;
  double P_18_5 (const double x) const;
  double P_18_5_Deriv (const double x) const;
  double P_18_6 (const double x) const;
  double P_18_6_Deriv (const double x) const;
  double P_18_7 (const double x) const;
  double P_18_7_Deriv (const double x) const;
  double P_18_8 (const double x) const;
  double P_18_8_Deriv (const double x) const;
  double P_18_9 (const double x) const;
  double P_18_9_Deriv (const double x) const;
  double P_18_10 (const double x) const;
  double P_18_10_Deriv (const double x) const;
  double P_18_11 (const double x) const;
  double P_18_11_Deriv (const double x) const;
  double P_18_12 (const double x) const;
  double P_18_12_Deriv (const double x) const;
  double P_18_13 (const double x) const;
  double P_18_13_Deriv (const double x) const;
  double P_18_14 (const double x) const;
  double P_18_14_Deriv (const double x) const;
  double P_18_15 (const double x) const;
  double P_18_15_Deriv (const double x) const;
  double P_18_16 (const double x) const;
  double P_18_16_Deriv (const double x) const;
  double P_18_17 (const double x) const;
  double P_18_17_Deriv (const double x) const;
  double P_18_18 (const double x) const;
  double P_18_18_Deriv (const double x) const;
  double P_19_0 (const double x) const;
  double P_19_0_Deriv (const double x) const;
  double P_19_1 (const double x) const;
  double P_19_1_Deriv (const double x) const;
  double P_19_2 (const double x) const;
  double P_19_2_Deriv (const double x) const;
  double P_19_3 (const double x) const;
  double P_19_3_Deriv (const double x) const;
  double P_19_4 (const double x) const;
  double P_19_4_Deriv (const double x) const;
  double P_19_5 (const double x) const;
  double P_19_5_Deriv (const double x) const;
  double P_19_6 (const double x) const;
  double P_19_6_Deriv (const double x) const;
  double P_19_7 (const double x) const;
  double P_19_7_Deriv (const double x) const;
  double P_19_8 (const double x) const;
  double P_19_8_Deriv (const double x) const;
  double P_19_9 (const double x) const;
  double P_19_9_Deriv (const double x) const;
  double P_19_10 (const double x) const;
  double P_19_10_Deriv (const double x) const;
  double P_19_11 (const double x) const;
  double P_19_11_Deriv (const double x) const;
  double P_19_12 (const double x) const;
  double P_19_12_Deriv (const double x) const;
  double P_19_13 (const double x) const;
  double P_19_13_Deriv (const double x) const;
  double P_19_14 (const double x) const;
  double P_19_14_Deriv (const double x) const;
  double P_19_15 (const double x) const;
  double P_19_15_Deriv (const double x) const;
  double P_19_16 (const double x) const;
  double P_19_16_Deriv (const double x) const;
  double P_19_17 (const double x) const;
  double P_19_17_Deriv (const double x) const;
  double P_19_18 (const double x) const;
  double P_19_18_Deriv (const double x) const;
  double P_19_19 (const double x) const;
  double P_19_19_Deriv (const double x) const;
  double P_20_0 (const double x) const;
  double P_20_0_Deriv (const double x) const;
  double P_20_1 (const double x) const;
  double P_20_1_Deriv (const double x) const;
  double P_20_2 (const double x) const;
  double P_20_2_Deriv (const double x) const;
  double P_20_3 (const double x) const;
  double P_20_3_Deriv (const double x) const;
  double P_20_4 (const double x) const;
  double P_20_4_Deriv (const double x) const;
  double P_20_5 (const double x) const;
  double P_20_5_Deriv (const double x) const;
  double P_20_6 (const double x) const;
  double P_20_6_Deriv (const double x) const;
  double P_20_7 (const double x) const;
  double P_20_7_Deriv (const double x) const;
  double P_20_8 (const double x) const;
  double P_20_8_Deriv (const double x) const;
  double P_20_9 (const double x) const;
  double P_20_9_Deriv (const double x) const;
  double P_20_10 (const double x) const;
  double P_20_10_Deriv (const double x) const;
  double P_20_11 (const double x) const;
  double P_20_11_Deriv (const double x) const;
  double P_20_12 (const double x) const;
  double P_20_12_Deriv (const double x) const;
  double P_20_13 (const double x) const;
  double P_20_13_Deriv (const double x) const;
  double P_20_14 (const double x) const;
  double P_20_14_Deriv (const double x) const;
  double P_20_15 (const double x) const;
  double P_20_15_Deriv (const double x) const;
  double P_20_16 (const double x) const;
  double P_20_16_Deriv (const double x) const;
  double P_20_17 (const double x) const;
  double P_20_17_Deriv (const double x) const;
  double P_20_18 (const double x) const;
  double P_20_18_Deriv (const double x) const;
  double P_20_19 (const double x) const;
  double P_20_19_Deriv (const double x) const;
  double P_20_20 (const double x) const;
  double P_20_20_Deriv (const double x) const;
};


#endif /*ASSLEGFUNCTION_H_*/
