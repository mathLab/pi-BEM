#ifndef octree_block_h
#define octree_block_h

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/matrix_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>

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

// And here are a few C++ standard header
// files that we will need:
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>


using namespace dealii;

template <int dim>
class OctreeBlock
{

public :

  typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;


private :

  unsigned int level;

  unsigned int parentId;

  unsigned int numChildren;

  unsigned int childrenId[8];

  std::vector <std::set <unsigned int> > nearNeigh;

  std::vector <std::set <unsigned int> > intList;

  std::vector <std::set <unsigned int> > nonIntList;

  Point<dim> pMin;

  double delta;

  std::vector <unsigned int> nodesId;

  std::map <cell_it, std::vector<unsigned int> > quadPointsId;


public:

  OctreeBlock();

  OctreeBlock(unsigned int level, unsigned int parent, Point<dim> pMin, double delta);

  OctreeBlock(const OctreeBlock<dim> &other);

  ~OctreeBlock();

  void CopyContent(const  OctreeBlock *other);

  void AddNode(unsigned int nodeId);

  void AddQuadPoint(cell_it elemPointer, unsigned int quadPointId);

  std::vector <unsigned int> GetBlockNodeList() const;

  void DelNodeList();

  std::map <cell_it, std::vector<unsigned int> > GetBlockQuadPointsList() const;

  void DelQuadPointsList();

  unsigned int GetBlockNodesNum() const;

  unsigned int GetBlockChildrenNum() const;

  unsigned int GetParentId() const;

  void AddChild(unsigned int childId);

  unsigned int GetChildId(unsigned int idInList) const;

  Point<dim> GetPMin() const;

  double GetDelta() const;

  void AddNearNeigh(unsigned int sublevel, const unsigned int nnBlockId);

  unsigned int NumNearNeigh(unsigned int sublevel) const;

  unsigned int NumNearNeighLevels() const;

  std::set <unsigned int> GetNearNeighs(unsigned int sublevel) const;

  void AddBlockToIntList(unsigned int sublevel, const unsigned int intListBlockId);

  unsigned int NumIntList(unsigned int sublevel) const;

  unsigned int NumIntListLevels() const;

  std::set <unsigned int> GetIntList(unsigned int sublevel) const;

  std::vector<std::set <unsigned int> > GetIntList() const;

  void AddBlockToNonIntList(unsigned int sublevel, const unsigned int intListBlockId);

  unsigned int NumNonIntList(unsigned int sublevel) const;

  unsigned int NumNonIntListLevels() const;

  std::set <unsigned int> GetNonIntList(unsigned int sublevel) const;

  void SetNearNeighSize(unsigned int sublevels);

  void SetIntListSize(unsigned int sublevels);

  void SetNonIntListSize(unsigned int sublevels);

  unsigned int GetNearNeighSize() const;

  unsigned int GetIntListSize() const;

  unsigned int GetNonIntListSize() const;

};


#endif
