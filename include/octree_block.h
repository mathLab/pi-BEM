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

#ifndef octree_block_h
#define octree_block_h

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>

using namespace dealii;

template <int dim>
class OctreeBlock
{
public:
  typedef typename DoFHandler<dim - 1, dim>::active_cell_iterator cell_it;

private:
  unsigned int level;

  types::global_dof_index parentId;

  unsigned int numChildren;

  types::global_dof_index childrenId[8];

  std::vector<std::set<types::global_dof_index>> nearNeigh;

  std::vector<std::set<types::global_dof_index>> intList;

  std::vector<std::set<types::global_dof_index>> nonIntList;

  Point<dim> pMin;

  double delta;

  std::vector<types::global_dof_index> nodesId;

  std::map<cell_it, std::vector<types::global_dof_index>> quadPointsId;

public:
  OctreeBlock();

  OctreeBlock(unsigned int            level,
              types::global_dof_index parent,
              Point<dim>              pMin,
              double                  delta);

  OctreeBlock(const OctreeBlock<dim> &other);

  ~OctreeBlock();

  void
  CopyContent(const OctreeBlock *other);

  void
  AddNode(types::global_dof_index nodeId);

  void
  AddQuadPoint(cell_it elemPointer, types::global_dof_index quadPointId);

  std::vector<types::global_dof_index>
  GetBlockNodeList() const;

  void
  DelNodeList();

  std::map<cell_it, std::vector<types::global_dof_index>>
  GetBlockQuadPointsList() const;

  void
  DelQuadPointsList();

  types::global_dof_index
  GetBlockNodesNum() const;

  unsigned int
  GetBlockChildrenNum() const;

  types::global_dof_index
  GetParentId() const;

  void
  AddChild(types::global_dof_index childId);

  types::global_dof_index
  GetChildId(unsigned int idInList) const;

  Point<dim>
  GetPMin() const;

  double
  GetDelta() const;

  void
  AddNearNeigh(unsigned int sublevel, const types::global_dof_index nnBlockId);

  unsigned int
  NumNearNeigh(unsigned int sublevel) const;

  unsigned int
  NumNearNeighLevels() const;

  std::set<types::global_dof_index>
  GetNearNeighs(unsigned int sublevel) const;

  void
  AddBlockToIntList(unsigned int                  sublevel,
                    const types::global_dof_index intListBlockId);

  types::global_dof_index
  NumIntList(unsigned int sublevel) const;

  std::set<types::global_dof_index>
  GetIntList(unsigned int sublevel) const;

  std::vector<std::set<types::global_dof_index>>
  GetIntList() const;

  void
  AddBlockToNonIntList(unsigned int                  sublevel,
                       const types::global_dof_index intListBlockId);

  types::global_dof_index
  NumNonIntList(unsigned int sublevel) const;

  std::set<types::global_dof_index>
  GetNonIntList(unsigned int sublevel) const;

  void
  SetNearNeighSize(unsigned int sublevels);

  void
  SetIntListSize(unsigned int sublevels);

  void
  SetNonIntListSize(unsigned int sublevels);

  unsigned int
  GetNearNeighSize() const;

  types::global_dof_index
  GetIntListSize() const;

  types::global_dof_index
  GetNonIntListSize() const;
};

#endif
