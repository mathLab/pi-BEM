#include "octree_block.h"

template <int dim>
OctreeBlock<dim>::OctreeBlock()
  : level(0)
  , parentId(0)
  , numChildren(0)
  , nearNeigh(1)
  , intList(1)
  , nonIntList(1)
  , pMin()
  , delta(0)
  , nodesId(0)
  , quadPointsId()
{
  for (int i = 0; i < dim; i++)
    {
      this->pMin(i) = 0.0;
    }
}

template <int dim>
OctreeBlock<dim>::OctreeBlock(unsigned int            level,
                              types::global_dof_index parent,
                              Point<dim>              pMin,
                              double                  delta)
  : level(level)
  , parentId(parent)
  , numChildren(0)
  , nearNeigh(1)
  , intList(1)
  , nonIntList(1)
  , pMin(pMin)
  , delta(delta)
  , nodesId(0)
  , quadPointsId()
{}

// TODO: not declaring will enable move semantics
/*
OctreeBlock<dim>::OctreeBlock(const OctreeBlock<dim> &other)
{
  this->parentId = other.parentId;
  this->level    = other.level;
  this->pMin     = other.pMin;
  this->delta    = other.delta;
  for (types::global_dof_index i = 0; i < other.nodesId.size(); i++)
    {
      this->nodesId.push_back(other.nodesId[i]);
    }
  this->quadPointsId = other.quadPointsId;
  this->numChildren  = other.numChildren;
  for (unsigned int i = 0; i < this->numChildren; i++)
    {
      this->childrenId[i] = other.childrenId[i];
    }
  this->nearNeigh  = other.nearNeigh;
  this->intList    = other.intList;
  this->nonIntList = other.nonIntList;
}

template <int dim>
OctreeBlock<dim>::~OctreeBlock()
{
  this->nearNeigh.clear();
  this->intList.clear();
  this->nonIntList.clear();
  this->nodesId.clear();
  this->quadPointsId.clear();
}

template <int dim>
void
OctreeBlock<dim>::CopyContent(const OctreeBlock<dim> *other)
{
  this->parentId = other->parentId;
  this->level    = other->level;
  this->pMin     = other->pMin;
  this->delta    = other->delta;
  for (types::global_dof_index i = 0; i < other->nodesId.size(); i++)
    {
      this->nodesId.push_back(other->nodesId[i]);
    }
  this->quadPointsId = other->quadPointsId;
  this->numChildren  = other->numChildren;
  for (unsigned int i = 0; i < this->numChildren; i++)
    {
      this->childrenId[i] = other->childrenId[i];
    }
  this->nearNeigh  = other->nearNeigh;
  this->intList    = other->intList;
  this->nonIntList = other->nonIntList;
}
*/

template <int dim>
void
OctreeBlock<dim>::AddNode(types::global_dof_index nodeId)
{
  this->nodesId.push_back(nodeId);
}

template <int dim>
void
OctreeBlock<dim>::AddQuadPoint(cell_it                 elemPointer,
                               types::global_dof_index quadPointId)
{
  this->quadPointsId[elemPointer].push_back(quadPointId);
}

template <int dim>
inline const std::vector<types::global_dof_index> &
OctreeBlock<dim>::GetBlockNodeList() const
{
  return this->nodesId;
}

template <int dim>
inline void
OctreeBlock<dim>::DelNodeList()
{
  this->nodesId.clear();
}

template <int dim>
inline const std::map<typename DoFHandler<dim - 1, dim>::active_cell_iterator,
                      std::vector<types::global_dof_index>> &
OctreeBlock<dim>::GetBlockQuadPointsList() const
{
  return this->quadPointsId;
}

template <int dim>
inline void
OctreeBlock<dim>::DelQuadPointsList()
{
  this->quadPointsId.clear();
}


template <int dim>
inline types::global_dof_index
OctreeBlock<dim>::GetBlockNodesNum() const
{
  return this->nodesId.size();
}

template <int dim>
inline unsigned int
OctreeBlock<dim>::GetBlockChildrenNum() const
{
  return this->numChildren;
}

template <int dim>
inline types::global_dof_index
OctreeBlock<dim>::GetParentId() const
{
  return this->parentId;
}

template <int dim>
inline void
OctreeBlock<dim>::AddChild(types::global_dof_index childId)
{
  this->childrenId[numChildren] = childId;
  this->numChildren += 1;
}

template <int dim>
inline types::global_dof_index
OctreeBlock<dim>::GetChildId(unsigned int idInList) const
{
  return this->childrenId[idInList];
}

template <int dim>
inline Point<dim>
OctreeBlock<dim>::GetPMin() const
{
  return this->pMin;
}

template <int dim>
inline double
OctreeBlock<dim>::GetDelta() const
{
  return this->delta;
}

template <int dim>
inline void
OctreeBlock<dim>::AddNearNeigh(unsigned int                  sublevel,
                               const types::global_dof_index nnBlockId)
{
  AssertIndexRange(sublevel, nearNeigh.size());
  this->nearNeigh[sublevel].insert(nnBlockId);
}

template <int dim>
inline unsigned int
OctreeBlock<dim>::NumNearNeigh(unsigned int sublevel) const
{
  AssertIndexRange(sublevel, nearNeigh.size());
  return this->nearNeigh[sublevel].size();
}

template <int dim>
inline unsigned int
OctreeBlock<dim>::NumNearNeighLevels() const
{
  return this->nearNeigh.size();
}

template <int dim>
inline const std::set<types::global_dof_index> &
OctreeBlock<dim>::GetNearNeighs(unsigned int sublevel) const
{
  AssertIndexRange(sublevel, nearNeigh.size());
  return this->nearNeigh[sublevel];
}

template <int dim>
inline void
OctreeBlock<dim>::AddBlockToIntList(
  unsigned int                  sublevel,
  const types::global_dof_index intListBlockId)
{
  AssertIndexRange(sublevel, intList.size());
  this->intList[sublevel].insert(intListBlockId);
}

template <int dim>
inline types::global_dof_index
OctreeBlock<dim>::NumIntList(unsigned int sublevel) const
{
  AssertIndexRange(sublevel, intList.size());
  return this->intList[sublevel].size();
}

template <int dim>
inline const std::set<types::global_dof_index> &
OctreeBlock<dim>::GetIntList(unsigned int sublevel) const
{
  AssertIndexRange(sublevel, intList.size());
  return this->intList[sublevel];
}

template <int dim>
inline const std::vector<std::set<types::global_dof_index>> &
OctreeBlock<dim>::GetIntList() const
{
  return this->intList;
}

template <int dim>
inline void
OctreeBlock<dim>::AddBlockToNonIntList(
  unsigned int                  sublevel,
  const types::global_dof_index intListBlockId)
{
  AssertIndexRange(sublevel, nonIntList.size());
  this->nonIntList[sublevel].insert(intListBlockId);
}

template <int dim>
inline types::global_dof_index
OctreeBlock<dim>::NumNonIntList(unsigned int sublevel) const
{
  AssertIndexRange(sublevel, nonIntList.size());
  return this->nonIntList[sublevel].size();
}

template <int dim>
inline const std::set<types::global_dof_index> &
OctreeBlock<dim>::GetNonIntList(unsigned int sublevel) const
{
  AssertIndexRange(sublevel, nonIntList.size());
  return this->nonIntList[sublevel];
}

template <int dim>
inline void
OctreeBlock<dim>::SetNearNeighSize(unsigned int sublevels)
{
  this->nearNeigh.resize(sublevels);
}

template <int dim>
inline void
OctreeBlock<dim>::SetIntListSize(unsigned int sublevels)
{
  this->intList.resize(sublevels);
}

template <int dim>
inline void
OctreeBlock<dim>::SetNonIntListSize(unsigned int sublevels)
{
  this->nonIntList.resize(sublevels);
}

template <int dim>
inline unsigned int
OctreeBlock<dim>::GetNearNeighSize() const
{
  return this->nearNeigh.size();
}

template <int dim>
inline types::global_dof_index
OctreeBlock<dim>::GetIntListSize() const
{
  return this->intList.size();
}

template <int dim>
inline types::global_dof_index
OctreeBlock<dim>::GetNonIntListSize() const
{
  return this->nonIntList.size();
}

template class OctreeBlock<2>;
template class OctreeBlock<3>;
