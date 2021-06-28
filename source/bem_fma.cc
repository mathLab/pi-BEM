#define TOLL 0.001

#include <omp.h>


#include <functional>

#include "../include/bem_fma.h"
#include "../include/laplace_kernel.h"
#include "Teuchos_TimeMonitor.hpp"
using Teuchos::RCP;
using Teuchos::Time;
using Teuchos::TimeMonitor;
using namespace tbb;

RCP<Time> MatrVec =
  Teuchos::TimeMonitor::getNewCounter("Multipole MatrVec Products Time");
RCP<Time> MultGen =
  Teuchos::TimeMonitor::getNewCounter("Multipole Generation Time");
RCP<Time> MultInt =
  Teuchos::TimeMonitor::getNewCounter("Multipole Integral Time");
RCP<Time> ListCreat =
  Teuchos::TimeMonitor::getNewCounter("Octree Generation Time");
RCP<Time> DirInt = Teuchos::TimeMonitor::getNewCounter("Direct Integral Time");
RCP<Time> PrecondTime =
  Teuchos::TimeMonitor::getNewCounter("FMA_preconditioner Time");
RCP<Time> LocEval =
  Teuchos::TimeMonitor::getNewCounter("Local Evaluation Time");

template <int dim>
BEMFMA<dim>::BEMFMA(MPI_Comm mpi_commy)
  : mpi_communicator(mpi_commy)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout(std::cout, (this_mpi_process == 0))
{}

template <int dim>
BEMFMA<dim>::~BEMFMA()
{
  if (blocks.size() > 0)
    {
      for (types::global_dof_index ii = 0; ii < num_blocks; ii++)
        {
          delete blocks[ii];
        }
    }
}

template <int dim>
void
BEMFMA<dim>::init_fma(
  const DoFHandler<dim - 1, dim>                       &input_dh,
  const std::vector<std::set<types::global_dof_index>> &db_in,
  const TrilinosWrappers::MPI::Vector                  &input_sn,
  const Mapping<dim - 1, dim>                          &input_mapping,
  unsigned int                                          quad_order,
  unsigned int                                          sing_quad_order)
{
  quadrature_order          = quad_order;
  singular_quadrature_order = sing_quad_order;
  fma_dh                    = &input_dh;
  dirichlet_nodes           = new const Vector<double>(
    input_sn); // per quadratura singolare e octree generator
  this_cpu_set.clear();
  this_cpu_set.set_size(fma_dh->n_dofs());
  this_cpu_set.add_indices(input_sn.locally_owned_elements());
  this_cpu_set.compress();
  double_nodes_set = &db_in; // da passare al metodo che fa il precondizionatore
  fma_mapping      = &input_mapping;
}

template <int dim>
void
BEMFMA<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Octree Params");
  {
    prm.declare_entry("Number of Octree Levels", "10", Patterns::Integer());

    prm.declare_entry(
      "Maximum Number of Collocation Points per Childless Block",
      "20",
      Patterns::Integer());
  }
  prm.leave_subsection();

  prm.enter_subsection("FMA Params");
  {
    prm.declare_entry("FMA Truncation Order", "6", Patterns::Integer());
  }
  prm.leave_subsection();

  prm.enter_subsection("FMA Params");
  {
    prm.declare_entry("Granularity for TBB simple for cycles",
                      "10",
                      Patterns::Integer());
  }
  prm.leave_subsection();
}

template <int dim>
void
BEMFMA<dim>::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Octree Params");
  {
    num_octree_levels       = prm.get_integer("Number of Octree Levels");
    max_num_nodes_per_block = prm.get_integer(
      "Maximum Number of Collocation Points per Childless Block");
  }
  prm.leave_subsection();

  prm.enter_subsection("FMA Params");
  {
    trunc_order = prm.get_integer("FMA Truncation Order");
  }
  prm.leave_subsection();

  prm.enter_subsection("FMA Params");
  {
    tbb_granularity = prm.get_integer("Granularity for TBB simple for cycles");
  }
  prm.leave_subsection();
}

template <int dim>
void
BEMFMA<dim>::direct_integrals()
{
#ifdef _OPENMP
  direct_integrals_omp();
#else
  direct_integrals_tbb();
#endif
}

template <int dim>
void
BEMFMA<dim>::direct_integrals_tbb()
{
  pcout << "Computing direct integrals (TBB)..." << std::endl;
  Teuchos::TimeMonitor LocalTimer(*DirInt);
  // The following function performs
  // the direct integrals
  // for the fast multipole algorithm
  // and saves the results into two
  // sparse matrices which will be
  // also used for precondictioning:
  // the actual preconditioner is a
  // third sparse matrix
  // declaration of the 3d singular
  // quadrature to be used

  // We compute a vector containing all the possible
  // singular quadratures. We need them to properly treat the direct
  // contributions among nodes on the same block.
  // TODO: is this the same as BEMProblem::get_singular_quadrature()?
  // potentially, there's a different dh
  std::vector<QTelles<dim - 1>> sing_quadratures;
  for (unsigned int i = 0; i < fma_dh->get_fe().dofs_per_cell; ++i)
    {
      sing_quadratures.push_back(
        QTelles<dim - 1>(singular_quadrature_order,
                         fma_dh->get_fe().get_unit_support_points()[i]));
    }
  const types::global_dof_index dofs_per_cell = fma_dh->get_fe().dofs_per_cell;

  // vector containing the ids of the dofs
  // of each cell: it will be used to transfer
  // the computed local rows of the matrices
  // into the global matrices
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // vector to store parts of rows of neumann
  // and dirichlet matrix obtained in local
  // operations
  Vector<double> local_neumann_matrix_row_i(fma_dh->get_fe().dofs_per_cell);
  Vector<double> local_dirichlet_matrix_row_i(fma_dh->get_fe().dofs_per_cell);

  // Now that we have checked that
  // the number of vertices is equal
  // to the number of degrees of
  // freedom, we construct a vector
  // of support points which will be
  // used in the local integrations:
  std::vector<Point<dim>> support_points(fma_dh->n_dofs());
  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*fma_mapping,
                                                     *fma_dh,
                                                     support_points);

  // After doing so, we can start the
  // integration loop over all cells,
  // where we first initialize the
  // FEValues object and get the
  // values of $\mathbf{\tilde v}$ at
  // the quadrature points (this
  // vector field should be constant,
  // but it doesn't hurt to be more
  // general):

  // first, we (re)initialize the
  // preconditioning matricies by
  // generating the corresponding
  // sparsity pattern, obtained by
  // means of the octree blocking
  // the idea here is that we take
  // each childless block containing
  // at least a dof (such dof index is i)
  // and its interaction list.
  // for each of the dofs i
  // we create a set with all the
  // dofs (this dof index is j) of the
  // elements having
  // at least one quad point in the
  // interaction list blocks: these
  // dofs will determine the non null
  // elements ij of the precondition
  // matrix

  /// TODO understand the bandwith of the preconditioner
  types::global_dof_index preconditioner_band =
    125 * fma_dh->get_fe().dofs_per_cell;
  TrilinosWrappers::MPI::Vector helper(this_cpu_set, mpi_communicator);
  init_prec_sparsity_pattern.reinit(this_cpu_set,
                                    mpi_communicator,
                                    preconditioner_band);

  // In the following we use WorkStream to parallelise, through TBB, the setting
  // up of the initial preconditioner that does not consider any constraint. We
  // define two structs that are needed: the first one is empty since we have
  // decided to use the capture of lambda functions to let the worker know what
  // it needs. The second one instead is filled by each worker and passed down
  // by reference to the copier that manage any racing conditions copying
  // properly the computed data where they belong.
  struct InitPrecScratch
  {};

  // Every copier needs to thing, the global indices of the row associated with
  // each block and the indices of the coloumns to be added to each row of the
  // sparsity pattern.
  struct InitPrecCopy
  {
    std::vector<types::global_dof_index>              block_indices;
    std::vector<std::vector<types::global_dof_index>> col_indices;
  };

  // The worker function uses the capture to know the actual state of the
  // BEMFMA<dim> class. In this way we can perform the computation of the column
  // to be added at each row quite straigtforwardly. Since all the workers must
  // be able to run in parallel we must be sure that no racing condition occurs.
  // We use the global IndexSet this_cpu_set to know if we the computation
  // belogs to the actual processor or not, thus using a MPI strategy. In this
  // function we compute the sparisity pattern due to the contributions of the
  // blocks that are in the childlessList.
  auto f_init_prec_childless_worker = [this](types::global_dof_index kk,
                                             InitPrecScratch &,
                                             InitPrecCopy &copy_data) {
    // We resize everything to be sure to compute, and then copy only the needed
    // data.
    copy_data.block_indices.clear();
    copy_data.col_indices.clear();
    // for each block in the childless
    // list we get the list of nodes and
    // we check if it contains nodes:
    // if no nodes are contained there is
    // nothing to do

    types::global_dof_index blockId = this->childlessList[kk];

    OctreeBlock<dim> *block1 = this->blocks[blockId];

    // std::vector<types::global_dof_index>
    const auto &block1Nodes = block1->GetBlockNodeList();

    std::vector<types::global_dof_index> local_dof_indices(
      this->fma_dh->get_fe().dofs_per_cell);
    if (block1Nodes.size() > 0)
      {
        // if block1 contains nodes,
        // we need to get all the quad points
        // in the intList blocks of block1
        // (such quad points will be used for
        // direct integrals)

        unsigned int intListSubLevs = block1->GetIntListSize();
        // std::set<types::global_dof_index>
        const auto &block1IntList = block1->GetIntList(intListSubLevs - 1);

        // in this set we will put all the
        // dofs of the cell to whom
        // the quad points belong
        std::set<types::global_dof_index> directNodes;

        // start looping on the intList
        // blocks (block2 here)
        for (auto pos = block1IntList.begin(); pos != block1IntList.end();
             pos++)
          {
            OctreeBlock<dim> *block2 = this->blocks[*pos];
            // std::map<cell_it, std::vector<types::global_dof_index>>
            const auto &blockQuadPointsList = block2->GetBlockQuadPointsList();

            // get the list of quad points
            // in block2 and loop on it
            for (auto it = blockQuadPointsList.begin();
                 it != blockQuadPointsList.end();
                 it++)
              {
                // the key of the map (*it.first pointer) is
                // the cell of the quad point: we will
                // get its dofs and put them in the set
                // of direct nodes
                cell_it cell = (*it).first;
                cell->get_dof_indices(local_dof_indices);
                directNodes.insert(local_dof_indices.cbegin(),
                                   local_dof_indices.cend());
              }
          }

        // the loop over blocks in intList
        // is over: for all the nodes in
        // block1, we know nodes in directNodes
        // list have direct integrals, so
        // we use them to create the
        // direct contributions matrices
        // sparsity pattern
        for (types::global_dof_index i = 0; i < block1Nodes.size(); i++)
          {
            if (this_cpu_set.is_element(block1Nodes[i]))
              {
                copy_data.block_indices.push_back(block1Nodes[i]);
                copy_data.col_indices.emplace_back(directNodes.cbegin(),
                                                   directNodes.cend());
              }
          }
      }
  };

  // The copier function uses the InitPrecCopy structure to know the global
  // indices to add to the global initial sparsity pattern. We use once again
  // the capture to access the global memory.
  auto f_init_prec_copier = [this](const InitPrecCopy &copy_data) {
    for (types::global_dof_index i = 0; i < copy_data.col_indices.size(); ++i)
      {
        this->init_prec_sparsity_pattern.add_entries(
          copy_data.block_indices[i],
          copy_data.col_indices[i].begin(),
          copy_data.col_indices[i].end(),
          false);
      }
  };

  // We need to create two empty structures that will be copied by WorkStream
  // and passed to each worker-copier to compute the sparsity pattern for blocks
  // in the childlessList.
  InitPrecScratch foo_scratch;
  InitPrecCopy    foo_copy;
  WorkStream::run(0,
                  childlessList.size(),
                  f_init_prec_childless_worker,
                  f_init_prec_copier,
                  foo_scratch,
                  foo_copy);

  // unfortunately, the direct integrals must not be computed only for the
  // quadPoints in the intList: if a bigger block is in the nonIntList of
  // another block, the bound for the multipole expansion application does
  // not hold, and so we must compute direct integrals. Here we scan the
  // nonIntlists of each block at each level to look for bigger blocks and
  // initialize the prec matrices sparsity pattern with the corresponding nodes
  for (unsigned int level = 1; level < num_octree_levels + 1;
       level++) // loop over levels
    {
      types::global_dof_index startBlockLevel = startLevel[level];

      // For each level we need again WorkStream to compute the entries in the
      // sparisity pattern that belong to block that are in the nonIntList of
      // the current block but that are of greater size. Once again we use
      // IndexSet to know if we the node in the block belong to the current
      // proccesor.
      auto f_init_prec_level_worker = [this,
                                       &startBlockLevel,
                                       &level](types::global_dof_index jj,
                                               InitPrecScratch &,
                                               InitPrecCopy &copy_data) {
        copy_data.block_indices.clear();
        copy_data.col_indices.clear();

        OctreeBlock<dim> *block1 =
          this->blocks[this->dofs_filled_blocks[level][jj]];
        // std::vector<types::global_dof_index>
        const auto &nodesBlk1Ids = block1->GetBlockNodeList();
        std::vector<types::global_dof_index> local_dof_indices(
          this->fma_dh->get_fe().dofs_per_cell);

        // again, no need to perform next operations if block has no nodes
        if (nodesBlk1Ids.size() > 0) // !!!CHECK, IT SEEMS TO BE USELESS
          {
            // for each block containing nodes, loop over all sublevels in his
            // NN list (this is because if a block remains childless BEFORE
            // the last level, at this point we need to compute all its
            // contributions up to the bottom level)

            for (unsigned int subLevel = 0;
                 subLevel < block1->NumNearNeighLevels();
                 subLevel++)
              {
                // in this vectors we are saving the nodes needing direct
                // integrals
                std::set<types::global_dof_index> directNodes;
                // std::set<types::global_dof_index>
                const auto &nonIntList = block1->GetNonIntList(subLevel);

                // loop over well separated blocks of higher size (level): in
                // this case we must use direct evaluation: for each block we
                // get the quad points list
                for (auto pos = nonIntList.begin();
                     pos != nonIntList.lower_bound(startBlockLevel);
                     pos++)
                  {
                    OctreeBlock<dim> *block2 = this->blocks[*pos];
                    // std::map<cell_it, std::vector<types::global_dof_index>>
                    const auto &blockQuadPointsList =
                      block2->GetBlockQuadPointsList();

                    // we loop on the cells of the quad blocks (*it.first
                    // pointer) and put their dofs in the direct list
                    for (auto it = blockQuadPointsList.begin();
                         it != blockQuadPointsList.end();
                         it++)
                      {
                        cell_it cell = (*it).first;
                        cell->get_dof_indices(local_dof_indices);
                        directNodes.insert(local_dof_indices.cbegin(),
                                           local_dof_indices.cend());
                      }
                  } // end loop over blocks of a sublevel of nonIntList

                // we use the nodes in directList, to create the sparsity
                // pattern
                for (types::global_dof_index i = 0; i < nodesBlk1Ids.size();
                     i++)
                  {
                    if (this_cpu_set.is_element(nodesBlk1Ids[i]))
                      {
                        copy_data.block_indices.push_back(nodesBlk1Ids[i]);
                        copy_data.col_indices.emplace_back(directNodes.begin(),
                                                           directNodes.end());
                      }
                  }
              } // end loop over sublevels
          }     // end if: is there any node in the block?
      };

      // we loop over blocks of each level. We call WorkStream with the same
      // copier as before.
      InitPrecScratch level_scratch;
      InitPrecCopy    level_copy;
      WorkStream::run(0,
                      dofs_filled_blocks[level].size(),
                      f_init_prec_level_worker,
                      f_init_prec_copier,
                      level_scratch,
                      level_copy);
    } // end loop over octree levels

  // sparsity pattern is ready and can be compressed; the direct matrices
  // and the preconditioner one are the initialized with the sparsity
  // pattern just computed
  init_prec_sparsity_pattern.compress();
  double filling_percentage =
    double(init_prec_sparsity_pattern.n_nonzero_elements()) /
    double(fma_dh->n_dofs() * fma_dh->n_dofs()) * 100.;

  pcout << "Preconditioner: " << init_prec_sparsity_pattern.n_nonzero_elements()
        << " Nonzeros out of " << fma_dh->n_dofs() * fma_dh->n_dofs() << ":  "
        << filling_percentage << "%" << std::endl;

  prec_neumann_matrix.reinit(init_prec_sparsity_pattern);
  prec_dirichlet_matrix.reinit(init_prec_sparsity_pattern);
  init_preconditioner.reinit(init_prec_sparsity_pattern);

  // We need to set up the parallel assembling of the direct contributions in
  // our FMA. Once again we don't use any scratch data. We use the capture of
  // lambda functions or additional parameter to let the worker know what it
  // needs to properly compute the near field interactions.
  struct DirectScratchData
  {};

  // Basically we are applying a local to global operation so we need only great
  // care in the copy. We need to memorise any contribution associated with any
  // node in the block. For this reason we use vectors of vectors that we
  // dynamically build up. We need also to memorise the dof_indices associated
  // with such nodes in order to copy these contributions properly in the global
  // memory.
  struct DirectCopyData
  {
    DirectCopyData(){
      // each thread will hold a local copy of Multipole expansions. here they
      // are initialized in a very dumb way, but they're always overwritten
      // so...
    };

    // The working copy constructor for the copy structure
    DirectCopyData(const DirectCopyData &in_vec)
    {
      vec_local_neumann_matrix_row_i = in_vec.vec_local_neumann_matrix_row_i;
      vec_local_dirichlet_matrix_row_i =
        in_vec.vec_local_dirichlet_matrix_row_i;
      vec_local_dof_indices = in_vec.vec_local_dof_indices;
      vec_node_index        = in_vec.vec_node_index;
    };

    // The Destructor needs to make foo_fma to point to NULL (for this reason it
    // is mutable const)
    ~DirectCopyData(){};

    // The pointer we use to copy everything back.
    std::vector<Vector<double>> vec_local_neumann_matrix_row_i;
    std::vector<Vector<double>> vec_local_dirichlet_matrix_row_i;
    std::vector<std::vector<types::global_dof_index>> vec_local_dof_indices;
    std::vector<types::global_dof_index>              vec_node_index;
    std::vector<types::global_dof_index>              vec_start_helper;
  };

  // The worker function, it computes the direct integral checking that the dofs
  // belong to the IndexSet of the processor. As we did for the preconditioner
  // we firstly compute all the direct contributions associated with the blocks
  // in the childlessList and we secondly we will take care of all the blocks in
  // the nonIntList that does not respect the bounds for the multipole
  // expansion.
  auto f_worker_direct_childless_non_int_list =
    [this](typename std::vector<types::global_dof_index>::iterator block_it,
           DirectScratchData &,
           DirectCopyData                &copy_data,
           const std::vector<Point<dim>> &support_points,
           std::vector<QTelles<dim - 1>> &sing_quadratures) {
      // this is the Id of the block
      copy_data.vec_local_dof_indices.clear();
      copy_data.vec_local_neumann_matrix_row_i.clear();
      copy_data.vec_local_dirichlet_matrix_row_i.clear();
      copy_data.vec_node_index.clear();
      copy_data.vec_start_helper.clear();
      types::global_dof_index blockId = *block_it;

      // and this is the block pointer
      OctreeBlock<dim> *block1 = this->blocks[blockId];
      // we get the block node list
      // std::vector<types::global_dof_index>
      const auto &block1Nodes = block1->GetBlockNodeList();

      // if a block has no nodes (if it only contains quad points), there is
      // nothing to do if instead there are nodes, we start integrating
      if (block1Nodes.size() > 0)
        {
          // we first get all the blocks in
          // the intList of the current block (block1) and loop over these
          // blocks, to create a list of ALL the quadrature points that lie in
          // the interaction list blocks: these quad points have to be
          // integrated directly. the list of direct quad points has to be a
          // std::map of std::set of integers, meaning that to each cell, we
          // associate a std::set containing all the direct quad point ids
          unsigned int intListNumLevs = block1->GetIntListSize();
          // std::set<types::global_dof_index>
          const auto &block1IntList = block1->GetIntList(intListNumLevs - 1);

          std::map<cell_it, std::set<types::global_dof_index>> directQuadPoints;
          for (auto pos = block1IntList.begin(); pos != block1IntList.end();
               pos++)
            {
              // now for each block block2 we get the list of quad points
              OctreeBlock<dim> *block2 = this->blocks[*pos];
              // std::map<cell_it, std::vector<types::global_dof_index>>
              const auto &blockQuadPointsList =
                block2->GetBlockQuadPointsList();

              // we now loop among the cells of the list and for each cell we
              // loop among its quad points, to copy them into the direct quad
              // points list
              for (const auto &pair : blockQuadPointsList)
                {
                  directQuadPoints[pair.first].insert(pair.second.begin(),
                                                      pair.second.end());
                }
            }

          // we are now ready to go: for each node, we know which quad points
          // are to be treated directly, and for each node, we will now perform
          // the integral. we then start looping on the nodes of the block
          types::global_dof_index helper_index = 0;
          for (types::global_dof_index i = 0; i < block1Nodes.size(); i++)
            {
              types::global_dof_index nodeIndex = block1Nodes[i];

              if (this->this_cpu_set.is_element(nodeIndex))
                {
                  copy_data.vec_node_index.push_back(nodeIndex);
                  copy_data.vec_start_helper.push_back(helper_index);

                  // we loop on the list of quad points to be treated directly
                  for (auto it = directQuadPoints.begin();
                       it != directQuadPoints.end();
                       it++)
                    {
                      // the vectors with the local integrals for the cell
                      // must first be zeroed
                      copy_data.vec_local_neumann_matrix_row_i.emplace_back(
                        this->fma_dh->get_fe().dofs_per_cell);
                      copy_data.vec_local_dirichlet_matrix_row_i.emplace_back(
                        this->fma_dh->get_fe().dofs_per_cell);

                      // we get the first entry of the map, i.e. the cell
                      // pointer and we check if the cell contains the current
                      // node, to decide if singular of regular quadrature is
                      // to be used
                      cell_it cell = (*it).first;
                      copy_data.vec_local_dof_indices.emplace_back(
                        this->fma_dh->get_fe().dofs_per_cell);
                      cell->get_dof_indices(
                        copy_data.vec_local_dof_indices.back());

                      // we copy the cell quad points in this set
                      // std::set<types::global_dof_index>
                      const auto & cellQuadPoints = (*it).second;
                      bool         is_singular    = false;
                      unsigned int singular_index =
                        numbers::invalid_unsigned_int;

                      for (unsigned int j = 0;
                           j < this->fma_dh->get_fe().dofs_per_cell;
                           ++j)
                        {
                          if ((*(this->double_nodes_set))[nodeIndex].count(
                                copy_data.vec_local_dof_indices.back()[j]) > 0)
                            {
                              singular_index = j;
                              is_singular    = true;
                              break;
                            }
                        }

                      // first case: the current node does not belong to the
                      // current cell: we use regular quadrature
                      if (!is_singular)
                        {
                          // we start looping on the quad points of the cell:
                          // *pos will be the index of the quad point
                          for (auto pos = cellQuadPoints.begin();
                               pos != cellQuadPoints.end();
                               pos++)
                            {
                              // here we compute the distance R between the node
                              // and the quad point

                              // MAGARI USARE FEVALUES CON IL DOFHANDLER CRETINO
                              // DISCONTINUO E IL MAPPING bem_fma
                              Point<dim> D;
                              double     s = 0.;

                              Assert(this->quadPoints.count(cell) > 0,
                                     StandardExceptions::ExcInvalidIterator());
                              const Tensor<1, dim> R =
                                this->quadPoints[cell][*pos] -
                                support_points[nodeIndex];
                              LaplaceKernel::kernels(R, D, s);

                              // and here are the integrals for each of the
                              // degrees of freedom of the cell: note how the
                              // quadrature values (position, normals,
                              // jacobianXweight, shape functions) are taken
                              // from the precomputed ones in
                              // ComputationalDomain class
                              for (unsigned int j = 0;
                                   j < this->fma_dh->get_fe().dofs_per_cell;
                                   ++j)
                                {
                                  Assert(
                                    this->quadNormals.count(cell) > 0,
                                    StandardExceptions::ExcInvalidIterator());
                                  Assert(
                                    this->quadShapeFunValues.count(cell) > 0,
                                    StandardExceptions::ExcInvalidIterator());
                                  Assert(
                                    this->quadJxW.count(cell) > 0,
                                    StandardExceptions::ExcInvalidIterator());
                                  copy_data.vec_local_neumann_matrix_row_i
                                    .back()(j) +=
                                    ((D * this->quadNormals[cell][*pos]) *
                                     this->quadShapeFunValues[cell][*pos][j] *
                                     this->quadJxW[cell][*pos]);
                                  copy_data.vec_local_dirichlet_matrix_row_i
                                    .back()(j) +=
                                    (s *
                                     this->quadShapeFunValues[cell][*pos][j] *
                                     this->quadJxW[cell][*pos]);
                                }
                            }
                        } // end if
                      else
                        {
                          // after some checks, we have to create the singular
                          // quadrature: here the quadrature points of the cell
                          // will be IGNORED, and the singular quadrature points
                          // are instead used. the 3d and 2d quadrature rules
                          // are different

                          // QUESTO E' IL SOLITO STEP 34, VEDI SE CAMBIARE CON
                          // QUELLO NUOVO PER STOKES
                          Assert(singular_index !=
                                   numbers::invalid_unsigned_int,
                                 ExcInternalError());

                          // TODO: validate
                          /* //reference
                          const Quadrature<dim - 1> *singular_quadrature =
                            (dim == 2 ?
                               dynamic_cast<Quadrature<dim - 1> *>(
                                 &sing_quadratures[singular_index]) :
                               (dim == 3 ?
                                  dynamic_cast<Quadrature<dim - 1> *>(
                                    &sing_quadratures[singular_index]) :
                                  0));
                          */
                          // it's templated anyway, there wouldn't be any need
                          // for ifs
                          const Quadrature<dim - 1> *singular_quadrature =
                            dynamic_cast<Quadrature<dim - 1> *>(
                              &sing_quadratures[singular_index]);
                          Assert(singular_quadrature, ExcInternalError());

                          // once the singular quadrature has been created, we
                          // employ it to create the corresponding fe_values
                          FEValues<dim - 1, dim> fe_v_singular(
                            *(this->fma_mapping),
                            this->fma_dh->get_fe(),
                            *(singular_quadrature),
                            update_jacobians | update_values |
                              update_normal_vectors | update_quadrature_points);
                          fe_v_singular.reinit(cell);

                          // here are the vectors of the quad points and normals
                          // vectors
                          const std::vector<Tensor<1, dim>> &singular_normals =
                            fe_v_singular.get_normal_vectors();
                          const std::vector<Point<dim>> &singular_q_points =
                            fe_v_singular.get_quadrature_points();

                          // and here is the integrals computation: note how in
                          // this case the values for shape functions & co. are
                          // not taken from the precomputed ones in
                          // ComputationalDomain class
                          for (unsigned int q = 0;
                               q < singular_quadrature->size();
                               ++q)
                            {
                              Point<dim> D;
                              double     s = 0.;

                              const Tensor<1, dim> R =
                                singular_q_points[q] -
                                support_points[nodeIndex];
                              LaplaceKernel::kernels(R, D, s);

                              for (unsigned int j = 0;
                                   j < this->fma_dh->get_fe().dofs_per_cell;
                                   ++j)
                                {
                                  copy_data.vec_local_neumann_matrix_row_i
                                    .back()(j) +=
                                    ((D * singular_normals[q]) *
                                     fe_v_singular.shape_value(j, q) *
                                     fe_v_singular.JxW(q));

                                  copy_data.vec_local_dirichlet_matrix_row_i
                                    .back()(j) +=
                                    (s * fe_v_singular.shape_value(j, q) *
                                     fe_v_singular.JxW(q));
                                }
                            }

                          // TODO: validate: it's allocated elsewhere, but dim=2
                          // is invalid anyway
                          if (dim == 2)
                            {
                              delete singular_quadrature;
                            }
                        } // end else

                      helper_index += 1;
                    } // end loop on cells of the intList
                }     // end check on this_cpu_set
            }         // end loop over nodes of block1
        }             // end if (nodes in block > 0)
    };

  // The copier function, it copies the value from the local array to the global
  // matrix. The copier function is composed by three nested for loops and
  // neeeds to copy all the contributions of all the nodes of the block the
  // worker has taken care off.
  auto f_copier_direct = [this](const DirectCopyData &copy_data) {
    // Finally, we need to add
    // the contributions of the
    // current cell to the
    // global matrix.
    if (copy_data.vec_node_index.size() > 0)
      {
        for (types::global_dof_index ii = 0;
             ii < copy_data.vec_node_index.size();
             ++ii)
          {
            // We need the helper function to know the proper contributions
            // associated with each block
            types::global_dof_index foo_start = copy_data.vec_start_helper[ii];
            types::global_dof_index foo_end =
              copy_data.vec_local_dof_indices.size();
            if (ii < copy_data.vec_node_index.size() - 1)
              {
                foo_end = copy_data.vec_start_helper[ii + 1];
              }

            for (types::global_dof_index kk = foo_start; kk < foo_end; ++kk)
              {
                for (unsigned int j = 0;
                     j < this->fma_dh->get_fe().dofs_per_cell;
                     ++j)
                  {
                    this->prec_neumann_matrix.add(
                      copy_data.vec_node_index[ii],
                      copy_data.vec_local_dof_indices[kk][j],
                      copy_data.vec_local_neumann_matrix_row_i[kk](j));
                    this->prec_dirichlet_matrix.add(
                      copy_data.vec_node_index[ii],
                      copy_data.vec_local_dof_indices[kk][j],
                      copy_data.vec_local_dirichlet_matrix_row_i[kk](j));

                    if ((*(this->dirichlet_nodes))(
                          copy_data.vec_local_dof_indices[kk][j]) > 0.8)
                      {
                        this->init_preconditioner.add(
                          copy_data.vec_node_index[ii],
                          copy_data.vec_local_dof_indices[kk][j],
                          -copy_data.vec_local_dirichlet_matrix_row_i[kk](j));
                      }
                    else
                      {
                        this->init_preconditioner.add(
                          copy_data.vec_node_index[ii],
                          copy_data.vec_local_dof_indices[kk][j],
                          copy_data.vec_local_neumann_matrix_row_i[kk](j));
                      }
                  }
              } // end loop on everything in the non int list of the node of the
                // block
          }     // end loop on nodes in block
      }
  };

  DirectScratchData direct_childless_scratch_data;
  DirectCopyData    direct_childless_copy_data;

  // Since our worker function needs an additional paramenter we need to perform
  // a bind to let WorkStream see a function requiring only 3 arguments
  WorkStream::run(childlessList.begin(),
                  childlessList.end(),
                  std::bind(f_worker_direct_childless_non_int_list,
                            std::placeholders::_1,
                            std::placeholders::_2,
                            std::placeholders::_3,
                            support_points,
                            sing_quadratures),
                  f_copier_direct,
                  direct_childless_scratch_data,
                  direct_childless_copy_data);

  // We need a worker function that takes care of
  // all the blocks in the nonIntlist of a block that are bigger than the block
  // itslef.
  auto f_worker_direct_bigger_blocks =
    [this](typename std::vector<types::global_dof_index>::iterator block_it,
           DirectScratchData &,
           DirectCopyData                &copy_data,
           const std::vector<Point<dim>> &support_points,
           types::global_dof_index        startBlockLevel) {
      copy_data.vec_local_dof_indices.clear();
      copy_data.vec_local_neumann_matrix_row_i.clear();
      copy_data.vec_local_dirichlet_matrix_row_i.clear();
      copy_data.vec_node_index.clear();
      copy_data.vec_start_helper.clear();

      types::global_dof_index blockId = *block_it;
      OctreeBlock<dim> *      block1  = this->blocks[blockId];
      // std::vector<types::global_dof_index>
      const auto &            nodesBlk1Ids = block1->GetBlockNodeList();
      types::global_dof_index helper_index = 0;
      for (types::global_dof_index i = 0; i < nodesBlk1Ids.size(); i++)
        {
          // for each block containing nodes, loop over all sublevels in his NN
          // list (this is because if a block remains childless BEFORE the last
          // level, at this point we need to compute all its contributions up to
          // the bottom level)
          types::global_dof_index nodeIndex = nodesBlk1Ids[i];
          copy_data.vec_node_index.push_back(nodeIndex);
          copy_data.vec_start_helper.push_back(helper_index);

          if (this->this_cpu_set.is_element(nodeIndex))
            {
              std::map<cell_it, std::set<types::global_dof_index>>
                directQuadPoints;

              for (unsigned int subLevel = 0;
                   subLevel < block1->NumNearNeighLevels();
                   subLevel++)
                {
                  // std::set<types::global_dof_index>
                  const auto &nonIntList = block1->GetNonIntList(subLevel);

                  // loop over well separated blocks of higher size
                  // (level)-----> in this case
                  // we must use direct evaluation (luckily being childless they
                  // only contain 1 element)
                  for (auto pos = nonIntList.begin();
                       pos != nonIntList.lower_bound(startBlockLevel);
                       pos++)
                    {
                      OctreeBlock<dim> *block2 = this->blocks[*pos];
                      // std::map<cell_it, std::vector<types::global_dof_index>>
                      const auto &blockQuadPointsList =
                        block2->GetBlockQuadPointsList();
                      for (auto it = blockQuadPointsList.begin();
                           it != blockQuadPointsList.end();
                           it++)
                        {
                          for (types::global_dof_index ii = 0;
                               ii < (*it).second.size();
                               ii++)
                            {
                              directQuadPoints[(*it).first].insert(
                                (*it).second[ii]);
                            }
                        }
                    } // end loop over blocks of a sublevel of nonIntList
                }     // end loop over sublevels

              for (auto it = directQuadPoints.begin();
                   it != directQuadPoints.end();
                   it++)
                {
                  // the vectors with the local integrals for the cell must
                  // first be zeroed
                  /*
                  copy_data.vec_local_neumann_matrix_row_i.push_back(
                    Vector<double>(this->fma_dh->get_fe().dofs_per_cell));
                  copy_data.vec_local_dirichlet_matrix_row_i.push_back(
                    Vector<double>(this->fma_dh->get_fe().dofs_per_cell));
                  */
                  copy_data.vec_local_neumann_matrix_row_i.emplace_back(
                    this->fma_dh->get_fe().dofs_per_cell);
                  copy_data.vec_local_dirichlet_matrix_row_i.emplace_back(
                    this->fma_dh->get_fe().dofs_per_cell);

                  // we get the first entry of the map, i.e. the cell pointer
                  // here the quadrature is regular as the cell is well
                  // separated
                  cell_it cell = (*it).first;
                  /*
                  copy_data.vec_local_dof_indices.push_back(
                    std::vector<types::global_dof_index>(
                      this->fma_dh->get_fe().dofs_per_cell));
                  */
                  copy_data.vec_local_dof_indices.emplace_back(
                    this->fma_dh->get_fe().dofs_per_cell);
                  cell->get_dof_indices(copy_data.vec_local_dof_indices.back());

                  // we copy the cell quad points in this set
                  // std::set<types::global_dof_index>
                  const auto &cellQuadPoints = (*it).second;

                  // we start looping on the quad points of the cell: *pos will
                  // be the index of the quad point
                  for (auto pos = cellQuadPoints.begin();
                       pos != cellQuadPoints.end();
                       pos++)
                    {
                      // here we compute the distance R between the node and the
                      // quad point
                      const Tensor<1, dim> R =
                        quadPoints[cell][*pos] - support_points[nodeIndex];
                      Point<dim> D;
                      double     s = 0.;
                      LaplaceKernel::kernels(R, D, s);

                      // and here are the integrals for each of the degrees of
                      // freedom of the cell: note how the quadrature values
                      // (position, normals, jacobianXweight, shape functions)
                      // are taken from the precomputed ones in
                      // ComputationalDomain class

                      for (unsigned int j = 0;
                           j < fma_dh->get_fe().dofs_per_cell;
                           ++j)
                        {
                          copy_data.vec_local_neumann_matrix_row_i.back()(j) +=
                            ((D * this->quadNormals[cell][*pos]) *
                             this->quadShapeFunValues[cell][*pos][j] *
                             this->quadJxW[cell][*pos]);
                          copy_data.vec_local_dirichlet_matrix_row_i.back()(
                            j) += (s * this->quadShapeFunValues[cell][*pos][j] *
                                   this->quadJxW[cell][*pos]);

                        } // end loop over the dofs in the cell
                    }     // end loop over the quad points in a cell

                  helper_index += 1;
                  // Finally, we need to add
                  // the contributions of the
                  // current cell to the
                  // global matrix.

                } // end loop over quad points in the direct quad points list
            }     // end check on proc

        } // end loop over nodes in a block
    };

  for (unsigned int level = 1; level < num_octree_levels + 1;
       level++) // loop over levels
    {
      types::global_dof_index startBlockLevel = startLevel[level];

      // For every level we need to run on all the blocks and check for bigger
      // blocks in their nonIntLists. Once again we need a bind to define the
      // proper 3 arguments function needed by WorkStream.
      DirectScratchData direct_bigger_scratch_data;
      DirectCopyData    direct_bigger_copy_data;
      WorkStream::run(dofs_filled_blocks[level].begin(),
                      dofs_filled_blocks[level].end(),
                      std::bind(f_worker_direct_bigger_blocks,
                                std::placeholders::_1,
                                std::placeholders::_2,
                                std::placeholders::_3,
                                support_points,
                                startBlockLevel),
                      f_copier_direct,
                      direct_bigger_scratch_data,
                      direct_bigger_copy_data);
    } // end loop over octree levels

  // as said, the direct integrals must not be computed only for the
  // quadPoints in the intList: if a bigger block is in the nonIntList of
  // another block, the bound for the multipole expansion application does
  // not hold, and so we must compute direct integrals. Here we scan the
  // nonIntlists of each block al each level to look for bigger blocks and
  // compute the direct integral contribution for the quadNodes in such
  // blocks

  pcout << "...done computing direct integrals" << std::endl;
}

template <int dim>
void
BEMFMA<dim>::direct_integrals_omp()
{
  pcout << "Computing direct integrals (OpenMP)..." << std::endl;
  Teuchos::TimeMonitor LocalTimer(*DirInt);
  // The following function performs the direct integrals for the fast multipole
  // algorithm and saves the results into two sparse matrices which will be also
  // used for precondictioning: the actual preconditioner is a third sparse
  // matrix declaration of the 3d singular quadrature to be used

  // We compute a vector containing all the possible singular quadratures. We
  // need them to properly treat the direct contributions among nodes on the
  // same block.
  // TODO: is this the same as BEMProblem::get_singular_quadrature()?
  // potentially, there's a different dh
  const auto                    dofs_per_cell = fma_dh->get_fe().dofs_per_cell;
  std::vector<QTelles<dim - 1>> sing_quadratures;
  sing_quadratures.reserve(dofs_per_cell);
  for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      sing_quadratures.push_back(
        QTelles<dim - 1>(singular_quadrature_order,
                         fma_dh->get_fe().get_unit_support_points()[i]));
    }

  // vector containing the ids of the dofs
  // of each cell: it will be used to transfer
  // the computed local rows of the matrices
  // into the global matrices
  // std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // vector to store parts of rows of neumann
  // and dirichlet matrix obtained in local
  // operations
  Vector<double> local_neumann_matrix_row_i(dofs_per_cell);
  Vector<double> local_dirichlet_matrix_row_i(dofs_per_cell);

  // Now that we have checked that
  // the number of vertices is equal
  // to the number of degrees of
  // freedom, we construct a vector
  // of support points which will be
  // used in the local integrations:
  std::vector<Point<dim>> support_points(fma_dh->n_dofs());
  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*fma_mapping,
                                                     *fma_dh,
                                                     support_points);

  // After doing so, we can start the
  // integration loop over all cells,
  // where we first initialize the
  // FEValues object and get the
  // values of $\mathbf{\tilde v}$ at
  // the quadrature points (this
  // vector field should be constant,
  // but it doesn't hurt to be more
  // general):

  // first, we (re)initialize the
  // preconditioning matricies by
  // generating the corresponding
  // sparsity pattern, obtained by
  // means of the octree blocking
  // the idea here is that we take
  // each childless block containing
  // at least a dof (such dof index is i)
  // and its interaction list.
  // for each of the dofs i
  // we create a set with all the
  // dofs (this dof index is j) of the
  // elements having
  // at least one quad point in the
  // interaction list blocks: these
  // dofs will determine the non null
  // elements ij of the precondition
  // matrix

  /// TODO: understand the bandwith of the preconditioner
  types::global_dof_index preconditioner_band =
    125 * fma_dh->get_fe().dofs_per_cell;
  init_prec_sparsity_pattern.reinit(this_cpu_set,
                                    mpi_communicator,
                                    preconditioner_band);

  /*
  shared(this,
         dofs_per_cell,
         sing_quadratures,
         local_dof_indices,
         local_neumann_matrix_row_i,
         local_dirichlet_matrix_row_i,
         support_points)
  */

// 1: assemble the preconditioner sparsity pattern
//  1.1: assemble from (the interacting blocks of) the leaf octcells
#pragma omp parallel
  {
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < this->childlessList.size(); ++i)
      {
        const OctreeBlock<dim> *block1 = this->blocks[childlessList[i]];
        // if the block doesn't contain dofs belonging to the current mpi
        // proc, it can be skipped
        bool proceed = false;
        for (const auto idx : block1->GetBlockNodeList())
          {
            if (this->this_cpu_set.is_element(idx))
              {
                proceed = true;
                break;
              }
          }

        if (proceed)
          {
            // if block1 contains nodes,
            // we need to get all the quad points
            // in the intList blocks of block1
            // (such quad points will be used for
            // direct integrals)

            // get the deepest interacting octcells
            const auto &block1IntList =
              block1->GetIntList(block1->GetIntListSize() - 1);

            // in this set we will put all the
            // dofs of the cell to whom
            // the quad points belong
            std::set<types::global_dof_index> directNodes;

            // start looping on the intList
            // blocks (block2 here)
            for (const auto blockId2 : block1IntList)
              {
                // get the list of quad points
                // in block2 and loop on it; but it's only using the
                // cell, it the loop container was a vector it would be
                // better
                for (const auto &pair :
                     blocks[blockId2]->GetBlockQuadPointsList())
                  {
                    // the key of the map (*it.first pointer) is
                    // the cell of the quad point: we will
                    // get its dofs and put them in the set
                    // of direct nodes
                    pair.first->get_dof_indices(local_dof_indices);
                    directNodes.insert(local_dof_indices.cbegin(),
                                       local_dof_indices.cend());
                  }
              }
            // directNodes contains all the dofs that interact with the
            // current block

            // the loop over blocks in intList
            // is over: for all the nodes in
            // block1, we know nodes in directNodes
            // list have direct integrals, so
            // we use them to create the
            // direct contributions matrices
            // sparsity pattern
            boost::container::flat_set<types::global_dof_index> flattened(
              directNodes.cbegin(), directNodes.cend());
            for (const auto idx : block1->GetBlockNodeList())
              {
                if (this->this_cpu_set.is_element(idx))
                  {
                    // unfortunately, the sparsity pattern is not yet
                    // allocated entirely and it requires synch
#pragma omp critical
                    {
                      /*
                      std::cout << "rank " << this_mpi_process << " omp thread "
                                << omp_get_thread_num()
                                << " before calling add_entries" << std::endl;
                      */
                      this->init_prec_sparsity_pattern.add_entries(
                        idx, flattened.cbegin(), flattened.cend(), true);
                    }
                  }
              }
          }
      }
  }

//  1.2: assemble from (the non interacting blocks of) the
#pragma omp parallel
  {
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

// tree traversal is a more interesting scenario than the flat leaf nodes list
// so, we will use tasks
#pragma omp single
    {
      for (unsigned int level = 1; level < num_octree_levels + 1; level++)
        {
          for (const auto &blockId : this->dofs_filled_blocks[level])
            {
              // blockId has at least a dof
              const OctreeBlock<dim> *block1 = this->blocks[blockId];
              // each task processes a block; it will only execute after the
              // parent's; again, there's synch when inserting new positions in
              // the sparsity pattern
#pragma omp task depend(in                                     \
                        : this->blocks[block1->GetParentId()]) \
  depend(out                                                   \
         : this->blocks[blockId]) firstprivate(blockId, block1, level)
              {
                // again, this block is interesting only if it contains dofs of
                // this mpi proc
                // if the block doesn't contain dofs belonging to the current
                // mpi proc, it can be skipped
                bool proceed = false;
                for (const auto idx : block1->GetBlockNodeList())
                  {
                    if (this->this_cpu_set.is_element(idx))
                      {
                        proceed = true;
                        break;
                      }
                  }

                if (proceed)
                  {
                    for (unsigned int subLevel = 0;
                         subLevel < block1->NumNearNeighLevels();
                         subLevel++)
                      {
                        std::set<types::global_dof_index> directNodes;

                        const auto &nonIntList =
                          block1->GetNonIntList(subLevel);

                        // loop over well separated blocks of higher size
                        // (level): in this case we must use direct evaluation:
                        // for each block we get the quad points list
                        for (auto pos = nonIntList.begin();
                             pos !=
                             nonIntList.lower_bound(this->startLevel[level]);
                             pos++)
                          {
                            const OctreeBlock<dim> *block2 = this->blocks[*pos];

                            // get the list of quad points
                            // in block2 and loop on it; but it's only using
                            // the cell, it the loop container was a vector it
                            // would be better
                            for (const auto &pair :
                                 block2->GetBlockQuadPointsList())
                              {
                                // the key of the map (*it.first pointer) is
                                // the cell of the quad point: we will
                                // get its dofs and put them in the set
                                // of direct nodes
                                pair.first->get_dof_indices(local_dof_indices);
                                directNodes.insert(local_dof_indices.cbegin(),
                                                   local_dof_indices.cend());
                              }
                          } // end loop over blocks of a sublevel of nonIntList

                        // we use the nodes in directList, to create the
                        // sparsity pattern
                        boost::container::flat_set<types::global_dof_index>
                          flattened(directNodes.cbegin(), directNodes.cend());
                        for (const auto idx : block1->GetBlockNodeList())
                          {
                            if (this->this_cpu_set.is_element(idx))
                              {
                                // unfortunately, the sparsity pattern is not
                                // yet allocated entirely and it requires synch
#pragma omp critical
                                {
                                  this->init_prec_sparsity_pattern.add_entries(
                                    idx,
                                    flattened.cbegin(),
                                    flattened.cend(),
                                    true);
                                }
                              }
                          }
                      } // end loop over sublevels
                  }
              }
            }
        }
    }
  }

  // sparsity pattern is ready and can be compressed; the direct matrices
  // and the preconditioner one are the initialized with the sparsity
  // pattern just computed
  init_prec_sparsity_pattern.compress();
  double filling_percentage =
    double(init_prec_sparsity_pattern.n_nonzero_elements()) /
    double(fma_dh->n_dofs() * fma_dh->n_dofs()) * 100.;

  pcout << "Preconditioner: " << init_prec_sparsity_pattern.n_nonzero_elements()
        << " Nonzeros out of " << fma_dh->n_dofs() * fma_dh->n_dofs() << ":  "
        << filling_percentage << "%" << std::endl;

  prec_neumann_matrix.reinit(init_prec_sparsity_pattern);
  prec_dirichlet_matrix.reinit(init_prec_sparsity_pattern);
  init_preconditioner.reinit(init_prec_sparsity_pattern);

  // 2: assemble the actual preconditioners
  //  2.1: loop over leaf blocks
#pragma omp parallel
  {
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < this->childlessList.size(); ++i)
      {
        const OctreeBlock<dim> *block1 = this->blocks[childlessList[i]];
        // if the block doesn't contain dofs belonging to the current mpi
        // proc, it can be skipped
        bool proceed = false;
        for (const auto idx : block1->GetBlockNodeList())
          {
            if (this->this_cpu_set.is_element(idx))
              {
                proceed = true;
                break;
              }
          }

        if (proceed)
          {
            std::map<cell_it, std::set<types::global_dof_index>>
                                                 directQuadPoints;
            std::vector<types::global_dof_index> local_dofs(dofs_per_cell);

            // loop over the deepest blocks that are well-separated from
            // block1
            unsigned int intListNumLevs = block1->GetIntListSize();
            for (const auto blockId : block1->GetIntList(intListNumLevs - 1))
              {
                const OctreeBlock<dim> *block2 = this->blocks[blockId];

                // merge their quadrature points, cell by cell
                for (const auto &pair : block2->GetBlockQuadPointsList())
                  {
                    directQuadPoints[pair.first].insert(pair.second.begin(),
                                                        pair.second.end());
                  }
              }

            // actual integration from
            for (const auto idx : block1->GetBlockNodeList())
              {
                if (this->this_cpu_set.is_element(idx))
                  {
                    for (const auto &pair : directQuadPoints)
                      {
                        const auto  cell       = pair.first;
                        const auto &qpointIdxs = pair.second;
                        cell->get_dof_indices(local_dofs);
                        // allocate buffers
                        std::vector<double> vec_local_neumann_matrix_row_i(
                          dofs_per_cell, 0);
                        std::vector<double> vec_local_dirichlet_matrix_row_i(
                          dofs_per_cell, 0);

                        // this is the contribution from cell: pair.first
                        // if dof idx is part (either itself or a double)
                        // of said cell, it requires a special quadrature

                        bool         is_singular = false;
                        unsigned int singular_index =
                          numbers::invalid_unsigned_int;
                        for (unsigned int j = 0; j < dofs_per_cell; ++j)
                          {
                            // recall const
                            // std::vector<std::set<types::global_dof_index>>
                            // *double_nodes_set
                            if ((*double_nodes_set)[idx].count(local_dofs[j]) >
                                0)
                              {
                                singular_index = j;
                                is_singular    = true;
                                break;
                              }
                          }

                        if (!is_singular)
                          {
                            for (const auto &qpointIdx : qpointIdxs)
                              {
                                // here we compute the distance R between
                                // the node and the quad point
                                Point<dim> D;
                                double     s = 0.;

                                Assert(
                                  this->quadPoints.count(cell) > 0,
                                  StandardExceptions::ExcInvalidIterator());
                                const Tensor<1, dim> R =
                                  this->quadPoints[cell][qpointIdx] -
                                  support_points[idx];
                                LaplaceKernel::kernels(R, D, s);

                                Assert(
                                  this->quadNormals.count(cell) > 0,
                                  StandardExceptions::ExcInvalidIterator());
                                Assert(
                                  this->quadShapeFunValues.count(cell) > 0,
                                  StandardExceptions::ExcInvalidIterator());
                                Assert(
                                  this->quadJxW.count(cell) > 0,
                                  StandardExceptions::ExcInvalidIterator());

                                // and here are the integrals for each of
                                // the degrees of freedom of the cell:
                                // note how the quadrature values
                                // (position, normals, jacobianXweight,
                                // shape functions) are taken from the
                                // precomputed ones in
                                // ComputationalDomain class
                                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                                  {
                                    const double m =
                                      this->quadShapeFunValues[cell][qpointIdx]
                                                              [j] *
                                      quadJxW[cell][qpointIdx];

                                    vec_local_neumann_matrix_row_i[j] +=
                                      (D * this->quadNormals[cell][qpointIdx]) *
                                      m;
                                    vec_local_dirichlet_matrix_row_i[j] +=
                                      s * m;
                                  }
                              }
                          }
                        else
                          {
                            Assert(singular_index !=
                                     numbers::invalid_unsigned_int,
                                   ExcInternalError());

                            const Quadrature<dim - 1> *singular_quadrature =
                              dynamic_cast<Quadrature<dim - 1> *>(
                                &sing_quadratures[singular_index]);
                            Assert(singular_quadrature, ExcInternalError());

                            // once the singular quadrature has been
                            // created, we employ it to create the
                            // corresponding fe_values
                            FEValues<dim - 1, dim> fe_v_singular(
                              *fma_mapping,
                              this->fma_dh->get_fe(),
                              *singular_quadrature,
                              update_jacobians | update_values |
                                update_normal_vectors |
                                update_quadrature_points);
                            fe_v_singular.reinit(cell);

                            const auto &singular_normals =
                              fe_v_singular.get_normal_vectors();
                            const auto &singular_q_points =
                              fe_v_singular.get_quadrature_points();

                            for (unsigned int q = 0;
                                 q < singular_quadrature->size();
                                 ++q)
                              {
                                Point<dim> D;
                                double     s = 0.;

                                const Tensor<1, dim> R =
                                  singular_q_points[q] - support_points[idx];
                                LaplaceKernel::kernels(R, D, s);

                                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                                  {
                                    const double m =
                                      fe_v_singular.shape_value(j, q) *
                                      fe_v_singular.JxW(q);

                                    vec_local_neumann_matrix_row_i[j] +=
                                      (D * singular_normals[q]) * m;
                                    vec_local_dirichlet_matrix_row_i[j] +=
                                      s * m;
                                  }
                              }
                          }

                        // the preconditioners have been fully allocated,
                        // there is no conflict
                        this->prec_neumann_matrix.add(
                          idx, local_dofs, vec_local_neumann_matrix_row_i);
                        this->prec_dirichlet_matrix.add(
                          idx, local_dofs, vec_local_dirichlet_matrix_row_i);

                        // TODO: these might change for each dof, cannot
                        // use the compact method
                        for (unsigned int j = 0; j < dofs_per_cell; ++j)
                          {
                            if ((*dirichlet_nodes)(local_dofs[j]) > 0.8)
                              {
                                this->init_preconditioner.add(
                                  idx,
                                  local_dofs[j],
                                  -vec_local_dirichlet_matrix_row_i[j]);
                              }
                            else
                              {
                                this->init_preconditioner.add(
                                  idx,
                                  local_dofs[j],
                                  vec_local_neumann_matrix_row_i[j]);
                              }
                          }
                      }
                  }
              }
          }
      }
  }
//  2.2: loop over the well separated blocks from levels above the current
#pragma omp parallel
  {
    std::vector<types::global_dof_index> local_dofs(dofs_per_cell);
// again, use tasks instead of flat loop;
#pragma omp single
    {
      for (unsigned int level = 1; level < num_octree_levels + 1; level++)
        {
          for (const auto &blockId : this->dofs_filled_blocks[level])
            {
              // blockId has at least a dof
              const OctreeBlock<dim> *block1 = this->blocks[blockId];
              // each task processes a block; it will only execute after the
              // parent's; again, there's synch when inserting new positions in
              // the sparsity pattern
#pragma omp task depend(in                                     \
                        : this->blocks[block1->GetParentId()]) \
  depend(out                                                   \
         : this->blocks[blockId]) firstprivate(blockId, block1, level)
              {
                for (const auto idx : block1->GetBlockNodeList())
                  {
                    if (this->this_cpu_set.is_element(idx))
                      {
                        std::map<cell_it, std::set<types::global_dof_index>>
                          directQuadPoints;

                        for (unsigned int subLevel = 0;
                             subLevel < block1->NumNearNeighLevels();
                             subLevel++)
                          {
                            const auto &nonIntList =
                              block1->GetNonIntList(subLevel);
                            for (auto pos = nonIntList.begin();
                                 pos != nonIntList.lower_bound(
                                          this->startLevel[level]);
                                 pos++)
                              {
                                const OctreeBlock<dim> *block2 =
                                  this->blocks[*pos];

                                // merge their quadrature points, cell by cell
                                for (const auto &pair :
                                     block2->GetBlockQuadPointsList())
                                  {
                                    directQuadPoints[pair.first].insert(
                                      pair.second.begin(), pair.second.end());
                                  }
                              }
                          }

                        for (const auto &pair : directQuadPoints)
                          {
                            const auto  cell       = pair.first;
                            const auto &qpointIdxs = pair.second;
                            cell->get_dof_indices(local_dofs);
                            // allocate buffers
                            std::vector<double> vec_local_neumann_matrix_row_i(
                              dofs_per_cell, 0);
                            std::vector<double>
                              vec_local_dirichlet_matrix_row_i(dofs_per_cell,
                                                               0);

                            for (const auto qpointIdx : qpointIdxs)
                              {
                                const Tensor<1, dim> R =
                                  this->quadPoints[cell][qpointIdx] -
                                  support_points[idx];
                                Point<dim> D;
                                double     s = 0.;
                                LaplaceKernel::kernels(R, D, s);

                                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                                  {
                                    const double m =
                                      this->quadShapeFunValues[cell][qpointIdx]
                                                              [j] *
                                      this->quadJxW[cell][qpointIdx];

                                    vec_local_neumann_matrix_row_i[j] +=
                                      (D * this->quadNormals[cell][qpointIdx]) *
                                      m;
                                    vec_local_dirichlet_matrix_row_i[j] +=
                                      s * m;
                                  }
                              }

                            // the preconditioners have been fully allocated,
                            // there is no conflict
                            this->prec_neumann_matrix.add(
                              idx, local_dofs, vec_local_neumann_matrix_row_i);
                            this->prec_dirichlet_matrix.add(
                              idx,
                              local_dofs,
                              vec_local_dirichlet_matrix_row_i);

                            // TODO: these might change for each dof, cannot
                            // use the compact method
                            for (unsigned int j = 0; j < dofs_per_cell; ++j)
                              {
                                if ((*dirichlet_nodes)(local_dofs[j]) > 0.8)
                                  {
                                    this->init_preconditioner.add(
                                      idx,
                                      local_dofs[j],
                                      -vec_local_dirichlet_matrix_row_i[j]);
                                  }
                                else
                                  {
                                    this->init_preconditioner.add(
                                      idx,
                                      local_dofs[j],
                                      vec_local_neumann_matrix_row_i[j]);
                                  }
                              }
                          }
                      }
                  }
              }
            }
        }
    }
  }

  pcout << "...done computing direct integrals" << std::endl;
}

template <>
void
BEMFMA<2>::multipole_integrals()
{}

// The following function set up the structure needed to generate the multipole
// expansions with a Boundary Element Method.
template <int dim>
void
BEMFMA<dim>::multipole_integrals()
{
  pcout << "Computing multipole integrals..." << std::endl;
  Teuchos::TimeMonitor LocalTimer(*MultInt);
  // we start clearing the two structures in which we will
  // store the integrals. these objects are quite complicated:
  // for each block we are going to get the portion of
  // the cells contained in it, and by means if their
  // quad points in the block, perform the integrals
  // (but remember that there is going to be one integral
  // for each shape function/cell dof, and that the integral
  // will be stored in a multipole expansion).
  // thus, this structure is a map associating the
  // blockId (an unsigned int) to the dofs_per_cell integrals
  // for each cell, i.e. to a map associating cells
  // to vectors of multipole expansions (vectors of size
  // dofs_per_cell)

  elemMultipoleExpansionsKer1.clear();
  elemMultipoleExpansionsKer2.clear();

  // we set a useful integer variable and perform some check

  const unsigned int dofs_per_cell = fma_dh->get_fe().dofs_per_cell;

  // AssertThrow(dofs_per_cell == GeometryInfo<dim-1>::vertices_per_cell,
  //             ExcMessage("The code in this function can only be used for "
  //                        "the usual Q1 elements."));

  // We need to set up elemMultipoleExpansionsKer1 and
  // elemMultipoleExpansionsKer2 these are quite complicated objects so we need
  // great care. Since the creation of new elements inside a map is not thread
  // safe we need to use WorkStream to ensure that everything is set up
  // properly. Basically we just need to set up a structure so we can allow for
  // an empty Scratch. The WorkStream will replace a loop over the childlessList
  // that used to set up all the structure. All the workers need to create their
  // own multipole structures and then let the copiers copy everything in the
  // global memory  using the capture of the lambda functions.
  struct MultipoleScratch
  {};

  struct MultipoleData
  {
    MultipoleData()
    {
      myelemMultipoleExpansionsKer1.clear();
      myelemMultipoleExpansionsKer2.clear();
    };

    MultipoleData(const MultipoleData &in_data)
    {
      myelemMultipoleExpansionsKer1 = in_data.myelemMultipoleExpansionsKer1;
      myelemMultipoleExpansionsKer2 = in_data.myelemMultipoleExpansionsKer2;
    };

    ~MultipoleData()
    {}

    // The local element multipole exapansions that the workers will fill.
    std::map<types::global_dof_index,
             std::map<cell_it, std::vector<MultipoleExpansion>>>
      myelemMultipoleExpansionsKer1;
    std::map<types::global_dof_index,
             std::map<cell_it, std::vector<MultipoleExpansion>>>
      myelemMultipoleExpansionsKer2;
  };

  // The worker lambda function that sets up the local structures.
  auto f_worker_multipole_integral =
    [this, dofs_per_cell](
      std::vector<types::global_dof_index>::iterator blocky,
      MultipoleScratch &,
      MultipoleData &copy_data) //, const unsigned int dofs_per_cell)
  {
    types::global_dof_index blockId = *blocky;
    OctreeBlock<dim>       *block   = this->blocks[blockId];
    double                  delta   = block->GetDelta();
    Point<dim>              deltaHalf;
    for (unsigned int i = 0; i < dim; i++)
      {
        deltaHalf(i) = delta / 2.;
      }
    Point<dim> blockCenter = block->GetPMin() + deltaHalf;

    // at this point, we get the list of quad nodes for the current block,
    // and loop over it
    // std::map<cell_it, std::vector<types::global_dof_index>>
    const auto &blockQuadPointsList = block->GetBlockQuadPointsList();
    std::vector<std::complex<double>> cache;
    for (auto it = blockQuadPointsList.begin(); it != blockQuadPointsList.end();
         it++)
      {
        // for each cell in the list, we get the list of its quad nodes
        // present in the current block
        cell_it cell = (*it).first;
        // std::vector<types::global_dof_index>
        const auto &cellQuadPoints = (*it).second;

        // the vectors in the structures that we have previously cleared
        // neet to be resized
        copy_data.myelemMultipoleExpansionsKer1[blockId][cell].resize(
          dofs_per_cell);
        copy_data.myelemMultipoleExpansionsKer2[blockId][cell].resize(
          dofs_per_cell);

        // std::map <cell_it, std::vector <MultipoleExpansion > > *foo_ker1
        // the vectors are now initialized with an empty multipole expansion
        // centered in the current block center
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {
            copy_data.myelemMultipoleExpansionsKer1[blockId][cell][j] =
              MultipoleExpansion(this->trunc_order,
                                 blockCenter,
                                 &(this->assLegFunction));
            copy_data.myelemMultipoleExpansionsKer2[blockId][cell][j] =
              MultipoleExpansion(this->trunc_order,
                                 blockCenter,
                                 &(this->assLegFunction));
          }

        // the contribution of each quadrature node (which can be seen as a
        // source with a certain strength) is introduced in the
        // multipole expansion with the appropriate methods (AddNormDer
        // for neumann matrix integrals, Add for dirichlet matrix
        // integrals)
        for (types::global_dof_index k = 0; k < cellQuadPoints.size(); ++k)
          {
            types::global_dof_index q = cellQuadPoints[k];
            for (unsigned int j = 0; j < this->fma_dh->get_fe().dofs_per_cell;
                 ++j)
              {
                Assert(this->quadShapeFunValues.count(cell) > 0,
                       StandardExceptions::ExcInvalidIterator());
                Assert(this->quadJxW.count(cell) > 0,
                       StandardExceptions::ExcInvalidIterator());
                Assert(this->quadPoints.count(cell) > 0,
                       StandardExceptions::ExcInvalidIterator());
                Assert(this->quadNormals.count(cell) > 0,
                       StandardExceptions::ExcInvalidIterator());
                copy_data.myelemMultipoleExpansionsKer1[blockId][cell][j]
                  .AddNormDer(this->quadShapeFunValues[cell][q][j] *
                                this->quadJxW[cell][q] / 4 / numbers::PI,
                              this->quadPoints[cell][q],
                              this->quadNormals[cell][q],
                              cache);
                copy_data.myelemMultipoleExpansionsKer2[blockId][cell][j].Add(
                  this->quadShapeFunValues[cell][q][j] *
                    this->quadJxW[cell][q] / 4 / numbers::PI,
                  this->quadPoints[cell][q],
                  cache);
              }
          } // end loop on cell quadrature points in the block
      }
  };

  // The copier function copies the local structures in the global memory.
  auto f_copier_multipole_integral = [this](const MultipoleData &copy_data) {
    this->elemMultipoleExpansionsKer1.insert(
      copy_data.myelemMultipoleExpansionsKer1.cbegin(),
      copy_data.myelemMultipoleExpansionsKer1.cend());
    this->elemMultipoleExpansionsKer2.insert(
      copy_data.myelemMultipoleExpansionsKer2.cbegin(),
      copy_data.myelemMultipoleExpansionsKer2.cend());
  };

  MultipoleScratch foo_scratch;
  MultipoleData    copy_data;

  WorkStream::run(childlessList.begin(),
                  childlessList.end(),
                  f_worker_multipole_integral,
                  f_copier_multipole_integral,
                  foo_scratch,
                  copy_data);

  pcout << "...done computing multipole integrals" << std::endl;
}

template <>
void
BEMFMA<2>::generate_multipole_expansions(
  const TrilinosWrappers::MPI::Vector &,
  const TrilinosWrappers::MPI::Vector &) const
{
  AssertThrow(true, ExcMessage("BEMFMA only 3D"));
}

// The following function performs the ascending phase of the algorithm. We
// need the values of the two traces of the solutions to fill the multipole
// expansions and then to let them be translated along the octree.
template <int dim>
void
BEMFMA<dim>::generate_multipole_expansions(
  const TrilinosWrappers::MPI::Vector &phi_values_in,
  const TrilinosWrappers::MPI::Vector &dphi_dn_values_in) const
{
  pcout << "Generating multipole expansions..." << std::endl;
  Teuchos::TimeMonitor LocalTimer(*MultGen);
  // TODO I DON'T KNOW IF WE CAN SPLIT THIS OPERATION, IN CASE THIS WOULD NOT BE
  // EASY COMMUNICATION, WE WOULD NEED TO COMMUNICATE THE ENTIRE CLASSES THAT
  // HOLD THE MULTIPOLE<->LOCAL EXPANSIONS
  const Vector<double> phi_values(phi_values_in);
  const Vector<double> dphi_dn_values(dphi_dn_values_in);
  // also here we clear the structures storing the multipole
  // expansions for Dirichlet and Neumann matrices

  blockMultipoleExpansionsKer1.clear();
  blockMultipoleExpansionsKer2.clear();

  // we reisze them: there's going to be an expansion per block
  blockMultipoleExpansionsKer1.resize(num_blocks);
  blockMultipoleExpansionsKer2.resize(num_blocks);

  // these two variables will be handy in the following
  std::vector<types::global_dof_index> local_dof_indices(
    fma_dh->get_fe().dofs_per_cell);

  // we loop on blocks and for each of them we create an empty multipole
  // expansion centered in the block center

  // First of all we need to create all the empty expansiones for all the
  // blocks. This is an embarassingly parallel operation that we can perform
  // using the ThreadGroup strategy without requiring any synchronization time.
  auto f_creation_tbb = [this](blocked_range<unsigned int> r) {
    for (unsigned int ii = r.begin(); ii < r.end(); ++ii)
      {
        Point<dim> blockCenter = blocks[ii]->GetCenter();

        blockMultipoleExpansionsKer1[ii] =
          MultipoleExpansion(trunc_order, blockCenter, &(assLegFunction));
        blockMultipoleExpansionsKer2[ii] =
          MultipoleExpansion(trunc_order, blockCenter, &(assLegFunction));
      }
  };

  parallel_for(blocked_range<unsigned int>(0, num_blocks, tbb_granularity),
               f_creation_tbb);

  // we now begin the rising phase of the algorithm: starting from the lowest
  // block levels (childless blocks) we get all the values of the multipole
  // integrals and aggregate them in the multipole expansion for each blocks

  // We need to create the multipole expansions for all the blocks in the
  // childlessList. Once again this is an embarassingly parallel operation. We
  // set up a new ThreadGroup strategy to let all the blocks to run in parallel.
  auto f_childless_creation_tbb = [this, &phi_values, &dphi_dn_values](
                                    blocked_range<types::global_dof_index> r) {
    for (types::global_dof_index kk = r.begin(); kk < r.end(); ++kk)
      {
        // for each block we get the center and the quad points
        types::global_dof_index blockId = childlessList[kk];
        OctreeBlock<dim>       *block   = blocks[blockId];

        double     delta = blocks[blockId]->GetDelta();
        Point<dim> deltaHalf;
        for (unsigned int i = 0; i < dim; i++)
          {
            deltaHalf(i) = delta / 2.;
          }

        // std::map<cell_it, std::vector<types::global_dof_index>>
        const auto &blockQuadPointsList = block->GetBlockQuadPointsList();

        // we loop on the cells of the quad points in the block: remember that
        // for each cell with a node in the block, we had created a set of
        // dofs_per_cell multipole expansion, representing the (partial)
        // integral on each cell
        std::vector<types::global_dof_index> my_local_dof_indices(
          fma_dh->get_fe().dofs_per_cell);

        std::vector<std::complex<double>> cache;
        for (auto it = blockQuadPointsList.begin();
             it != blockQuadPointsList.end();
             it++)
          {
            cell_it cell = (*it).first;
            cell->get_dof_indices(my_local_dof_indices);

            // for each cell we get the dof_indices, and add to the block
            // multipole expansion, the integral previously computed, multiplied
            // by the phi or dphi_dn value at the corresponding dof of the cell.
            // A suitable MultipoleExpansion class method has been created for
            // this purpose
            for (unsigned int jj = 0; jj < fma_dh->get_fe().dofs_per_cell; ++jj)
              {
                AssertIndexRange(blockId, blockMultipoleExpansionsKer1.size());
                blockMultipoleExpansionsKer1[blockId].Add(
                  elemMultipoleExpansionsKer1[blockId][cell][jj],
                  phi_values(my_local_dof_indices[jj]));

                AssertIndexRange(blockId, blockMultipoleExpansionsKer2.size());
                blockMultipoleExpansionsKer2[blockId].Add(
                  elemMultipoleExpansionsKer2[blockId][cell][jj],
                  dphi_dn_values(my_local_dof_indices[jj]));
              }
          } // end loop ond block elements
      }
  };

  parallel_for(blocked_range<types::global_dof_index>(0,
                                                      childlessList.size(),
                                                      tbb_granularity),
               f_childless_creation_tbb);

  // Now all the lower level blocks have a multipole expansion containing the
  // contribution to the integrals of all the quad points in them: we now begin
  // summing the child multipole expansion to the the parents expansions: to do
  // that we need to translate che children expansions to the parent block
  // center: there is a MultipoleExpansion class for this too. In this case we
  // have to set up some syncrhonization since more blocks may have the same
  // parent and this would lead to race conditions among thread. For this reason
  // we need WorkStream to copy all the expansion up along the tree.

  // The scratch data is nothing
  struct AscendScratchData
  {};

  // Basically we are applying a local to global operation so we need only great
  // care in the copy.
  struct AscendCopyData
  {
    // We need an input BEMFMA to set up the local MultipoleExpansion classes.
    AscendCopyData(const BEMFMA<dim> *dummy_fma)
    {
      // each thread will hold a local copy of Multipole expansions. here they
      // are initialized in a very dumb way, but they're always overwritten
      // so...
      Point<dim>         zero;
      MultipoleExpansion dummy(dummy_fma->trunc_order,
                               zero,
                               &(dummy_fma->assLegFunction));
      translatedBlockMultipoleExpansionKer1 = dummy;
      translatedBlockMultipoleExpansionKer2 = dummy;
    };

    // The working copy constructor for the copy structure
    AscendCopyData(const AscendCopyData &in_vec)
    {
      translatedBlockMultipoleExpansionKer1 =
        in_vec.translatedBlockMultipoleExpansionKer1;
      translatedBlockMultipoleExpansionKer2 =
        in_vec.translatedBlockMultipoleExpansionKer2;
      start = in_vec.start;
    };

    // The Destructor needs to make foo_fma to point to NULL (for this reason it
    // is mutable const)
    ~AscendCopyData(){};

    types::global_dof_index start;
    types::global_dof_index parentId;
    MultipoleExpansion      translatedBlockMultipoleExpansionKer1;
    MultipoleExpansion      translatedBlockMultipoleExpansionKer2;
    // The pointer we use to copy everything back.
  };

  // The worker function, it copies the value from the memory at a certain level
  // to local array in copy.data
  auto f_worker_ascend =
    [this](typename std::vector<OctreeBlock<dim> *>::iterator block_it,
           AscendScratchData &,
           AscendCopyData         &copy_data,
           types::global_dof_index start) {
      types::global_dof_index kk =
        std::distance(this->blocks.begin(), block_it);

      Point<dim>         zero;
      MultipoleExpansion dummy(this->trunc_order,
                               zero,
                               &(this->assLegFunction));
      copy_data.translatedBlockMultipoleExpansionKer1 = dummy;
      copy_data.translatedBlockMultipoleExpansionKer2 = dummy;
      copy_data.start                                 = start;
      copy_data.parentId = (*block_it)->GetParentId();

      // We set the center to the parentId center.
      AssertIndexRange(copy_data.parentId, blockMultipoleExpansionsKer1.size());
      copy_data.translatedBlockMultipoleExpansionKer1.SetCenter(
        this->blockMultipoleExpansionsKer1[copy_data.parentId].GetCenter());
      AssertIndexRange(copy_data.parentId, blockMultipoleExpansionsKer2.size());
      copy_data.translatedBlockMultipoleExpansionKer2.SetCenter(
        this->blockMultipoleExpansionsKer2[copy_data.parentId].GetCenter());

      std::vector<std::complex<double>> cache;
      // We translate the children blocks to the local array.
      AssertIndexRange(kk, blockMultipoleExpansionsKer1.size());
      copy_data.translatedBlockMultipoleExpansionKer1.Add(
        this->blockMultipoleExpansionsKer1[kk], cache);
      AssertIndexRange(kk, blockMultipoleExpansionsKer2.size());
      copy_data.translatedBlockMultipoleExpansionKer2.Add(
        this->blockMultipoleExpansionsKer2[kk], cache);
    };

  // The copier function, it copies the value from the local array to the parent
  // blocks
  auto f_copier_ascend = [this](const AscendCopyData &copy_data) {
    // Now we just need to copy on the real global memory. Since the center is
    // the same we just need to add the value to the expansion wothout any
    // further translation. The MultipoleExpansion class takes care of that
    // automatically.
    std::vector<std::complex<double>> cache;
    AssertIndexRange(copy_data.parentId, blockMultipoleExpansionsKer1.size());
    this->blockMultipoleExpansionsKer1[copy_data.parentId].Add(
      copy_data.translatedBlockMultipoleExpansionKer1, cache);
    AssertIndexRange(copy_data.parentId, blockMultipoleExpansionsKer2.size());
    this->blockMultipoleExpansionsKer2[copy_data.parentId].Add(
      copy_data.translatedBlockMultipoleExpansionKer2, cache);
  };

  for (unsigned int level = num_octree_levels; level > 0; level--)
    {
      AscendScratchData sample_scratch;
      AscendCopyData    sample_copy(this);

      // We need to be sure that there actually are some blocks in the current
      // level, then we need to loop over all blocks on the level and perform
      // all the translations. We need a bind for the worker in order to deal
      // with the additional unsigned int startLevel[level].
      if (endLevel[level] >= startLevel[level])
        {
          WorkStream::run(blocks.begin() + startLevel[level],
                          blocks.begin() + endLevel[level] + 1,
                          std::bind(f_worker_ascend,
                                    std::placeholders::_1,
                                    std::placeholders::_2,
                                    std::placeholders::_3,
                                    startLevel[level]),
                          f_copier_ascend,
                          sample_scratch,
                          sample_copy);
        }
    } // end loop over levels

  pcout << "...done generating multipole expansions" << std::endl;
}

template <>
void
BEMFMA<2>::multipole_matr_vect_products(const TrilinosWrappers::MPI::Vector &,
                                        const TrilinosWrappers::MPI::Vector &,
                                        TrilinosWrappers::MPI::Vector &,
                                        TrilinosWrappers::MPI::Vector &) const
{}

template <>
void
BEMFMA<2>::multipole_matr_vect_products_tbb(
  const TrilinosWrappers::MPI::Vector &,
  const TrilinosWrappers::MPI::Vector &,
  TrilinosWrappers::MPI::Vector &,
  TrilinosWrappers::MPI::Vector &) const
{}

template <>
void
BEMFMA<2>::multipole_matr_vect_products_omp(
  const TrilinosWrappers::MPI::Vector &,
  const TrilinosWrappers::MPI::Vector &,
  TrilinosWrappers::MPI::Vector &,
  TrilinosWrappers::MPI::Vector &) const
{}

// The following functions takes care of the descending phase of the FMA.
template <int dim>
void
BEMFMA<dim>::multipole_matr_vect_products(
  const TrilinosWrappers::MPI::Vector &phi_values,
  const TrilinosWrappers::MPI::Vector &dphi_dn_values,
  TrilinosWrappers::MPI::Vector       &matrVectProdN,
  TrilinosWrappers::MPI::Vector       &matrVectProdD) const
{
#ifdef _OPENMP
  multipole_matr_vect_products_omp(phi_values,
                                   dphi_dn_values,
                                   matrVectProdN,
                                   matrVectProdD);
#else
  multipole_matr_vect_products_tbb(phi_values,
                                   dphi_dn_values,
                                   matrVectProdN,
                                   matrVectProdD);
#endif
}

template <int dim>
void
BEMFMA<dim>::multipole_matr_vect_products_tbb(
  const TrilinosWrappers::MPI::Vector &phi_values,
  const TrilinosWrappers::MPI::Vector &dphi_dn_values,
  TrilinosWrappers::MPI::Vector &      matrVectProdN,
  TrilinosWrappers::MPI::Vector &      matrVectProdD) const
{
  pcout << "Computing multipole matrix-vector products (TBB)..." << std::endl;
  Teuchos::TimeMonitor LocalTimer(*MatrVec);

  // and here we compute the direct integral contributions (stored in two
  // sparse matrices)
  prec_neumann_matrix.vmult(matrVectProdN, phi_values);
  prec_dirichlet_matrix.vmult(matrVectProdD, dphi_dn_values);

  // from here on, we compute the multipole expansions contributions
  // we start cleaning past sessions
  // store old values

  // TODO: would not need clear-resize if LocalExpansion's copy assignement is
  // correct
  blockLocalExpansionsKer1.clear();
  blockLocalExpansionsKer2.clear();

  blockLocalExpansionsKer1.resize(num_blocks);
  blockLocalExpansionsKer2.resize(num_blocks);

  // we declare some familiar variables that will be useful in the method
  std::vector<Point<dim>> support_points(fma_dh->n_dofs());
  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*fma_mapping,
                                                     *fma_dh,
                                                     support_points);
  // std::vector<types::global_dof_index> local_dof_indices(
  //  fma_dh->get_fe().dofs_per_cell);

  // First of all we need to create all the empty expansiones for all the
  // blocks. This is an embarassingly parallel operation that we can perform
  // using the ThreadGroup strategy without requiring any synchronization
  // time.
  auto f_local_creation_tbb = [this](blocked_range<types::global_dof_index> r) {
    for (types::global_dof_index ii = r.begin(); ii < r.end(); ++ii)
      {
        const auto blockCenter = blocks[ii]->GetCenter();

        blockLocalExpansionsKer1[ii] =
          LocalExpansion(trunc_order, blockCenter, &(assLegFunction));
        blockLocalExpansionsKer2[ii] =
          LocalExpansion(trunc_order, blockCenter, &(assLegFunction));
      }
  };

  // each instance is independent
  parallel_for(blocked_range<types::global_dof_index>(0,
                                                      num_blocks,
                                                      tbb_granularity),
               f_local_creation_tbb);

  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////

  // In order to perform the descending phase properly we need WorkStream.
  // Inside the worker we check the IndexSet for the MPI parallelisation. The
  // scaracth is empty once again since this class only has to copy from
  // parents to children. As in the ascending phase we use captures in order
  // to pass global information to the worker and the copier
  struct DescendScratchData
  {};

  // Basically we are applying a local to global operation so we need only
  // great care in the copy.
  struct DescendCopyData
  {
    DescendCopyData(const BEMFMA<dim> *dummy_fma)
    {
      // each thread will hold a local copy of Multipole expansions. here they
      // are initialized in a very dumb way, but they're always overwritten
      // so...
      Point<dim>     zero;
      LocalExpansion dummy(dummy_fma->trunc_order,
                           zero,
                           &(dummy_fma->assLegFunction));
      blockLocalExpansionKer1 = dummy;
      blockLocalExpansionKer2 = dummy;

      localTimeEvalNumCalls = 0;
    };

    // The working copy constructor for the copy structure
    DescendCopyData(const DescendCopyData &in_vec)
    {
      blockLocalExpansionKer1 = in_vec.blockLocalExpansionKer1;
      blockLocalExpansionKer2 = in_vec.blockLocalExpansionKer2;
      start                   = in_vec.start;

      localTimeEvalNumCalls = in_vec.localTimeEvalNumCalls;
    };

    // The Destructor needs to make foo_fma to point to NULL (for this reason
    // it is mutable const)
    ~DescendCopyData(){};

    types::global_dof_index              start;
    types::global_dof_index              blockId;
    LocalExpansion                       blockLocalExpansionKer1;
    LocalExpansion                       blockLocalExpansionKer2;
    Vector<double>                       matrVectorProductContributionKer1;
    Vector<double>                       matrVectorProductContributionKer2;
    std::vector<types::global_dof_index> nodesBlk1Ids;

    int localTimeEvalNumCalls;
  };

  // The worker function, it copies the value about the parent from the global
  // memory at a certain level to local array in copy.data

  auto f_worker_Descend =
    [this, &support_points](
      std::vector<types::global_dof_index>::const_iterator block_it_id,
      DescendScratchData &,
      DescendCopyData &             copy_data,
      const types::global_dof_index start) {
      copy_data.start                  = start;
      copy_data.blockId                = *block_it_id;
      types::global_dof_index kk       = *block_it_id;
      OctreeBlock<dim> *      block_it = this->blocks[*block_it_id];

      //*****************definire chi e' on_process qui
      AssertIndexRange(kk, blockLocalExpansionsKer1.size());
      Point<dim>     center = this->blockLocalExpansionsKer1[kk].GetCenter();
      LocalExpansion dummy(this->trunc_order, center, &(this->assLegFunction));
      copy_data.blockLocalExpansionKer1 = dummy;
      copy_data.blockLocalExpansionKer2 = dummy;
      // TODO: could the level be captured by the lambda? as they go in waves
      unsigned int level = 0;
      for (unsigned int lev = 0; lev < this->num_octree_levels + 1; lev++)
        {
          if (kk >= this->startLevel[lev] && kk <= this->endLevel[lev])
            {
              level = lev;
              break;
            }
        }

      types::global_dof_index startBlockLevel = this->startLevel[level];
      types::global_dof_index endBlockLevel   = this->endLevel[level];
      // std::vector<types::global_dof_index>
      const auto &nodesBlk1Ids = block_it->GetBlockNodeList();
      bool        on_process   = false;
      for (auto ind : nodesBlk1Ids)
        {
          if (this->this_cpu_set.is_element(ind))
            {
              on_process = true;
              break;
            }
        }

      copy_data.nodesBlk1Ids = nodesBlk1Ids;
      copy_data.matrVectorProductContributionKer1.reinit(nodesBlk1Ids.size());
      copy_data.matrVectorProductContributionKer2.reinit(nodesBlk1Ids.size());

      if (on_process)
        {
          std::vector<std::complex<double>> cache;
          // the local expansion of the parent must be translated down into
          // the current block
          types::global_dof_index parentId = block_it->GetParentId();
          AssertIndexRange(parentId, blockLocalExpansionsKer1.size());
          copy_data.blockLocalExpansionKer1.Add(
            this->blockLocalExpansionsKer1[parentId], cache);
          AssertIndexRange(parentId, blockLocalExpansionsKer2.size());
          copy_data.blockLocalExpansionKer2.Add(
            this->blockLocalExpansionsKer2[parentId], cache);

          for (unsigned int subLevel = 0;
               subLevel < block_it->NumNearNeighLevels();
               subLevel++)
            {
              // we get the nonIntList of each block
              // std::set<dealii::types::boundary_id>
              const auto &nonIntList = block_it->GetNonIntList(subLevel);

              // we start converting into local expansions, all the multipole
              // expansions of all the blocks in nonIntList, that are of the
              // same size (level) of the current block. To perform the
              // conversion, we use another member of the LocalExpansion
              // class. Note that all the contributions to the integrals of
              // blocks bigger than current block had already been considered
              // in the direct integrals method (the bounds of the local and
              // multipole expansion do not hold in such case)
              for (auto pos1 = nonIntList.lower_bound(startBlockLevel);
                   pos1 != nonIntList.upper_bound(endBlockLevel);
                   pos1++) // loop over NNs of NNs and add them to intList
                {
                  types::global_dof_index block2Id = *pos1;

                  copy_data.blockLocalExpansionKer1.Add(
                    this->blockMultipoleExpansionsKer1[block2Id], cache);
                  copy_data.blockLocalExpansionKer2.Add(
                    this->blockMultipoleExpansionsKer2[block2Id], cache);
                } // end loop over well separated blocks of the same size
                  // (level)

              // loop over well separated blocks of the smaller size
              // (level)-----> use multipoles without local expansions

              // we now have to loop over blocks in the nonIntList that are
              // smaller than the current blocks: in this case the bound for
              // the conversion of a multipole into local expansion does not
              // hold, but the bound for the evaluation of the multipole
              // expansions does hold: thus, we will simply evaluate the
              // multipole expansions of such blocks for each node in the
              // block
              for (auto pos1 = nonIntList.upper_bound(endBlockLevel);
                   pos1 != nonIntList.end();
                   pos1++)
                {
                  types::global_dof_index block2Id = *pos1;

                  // TODO: restore functionality without forcing thread safety
                  // Teuchos::TimeMonitor LocalTimer(*LocEval);
                  // copy_data.localTimeEvalNumCalls++;

                  for (types::global_dof_index ii = 0; ii < nodesBlk1Ids.size();
                       ii++) // loop over each node of (*block_it)
                    {
                      const Point<dim> &nodeBlk1 =
                        support_points[nodesBlk1Ids[ii]];

                      copy_data.matrVectorProductContributionKer1(ii) +=
                        this->blockMultipoleExpansionsKer1[block2Id].Evaluate(
                          nodeBlk1, cache);
                      copy_data.matrVectorProductContributionKer2(ii) +=
                        this->blockMultipoleExpansionsKer2[block2Id].Evaluate(
                          nodeBlk1, cache);
                    }
                } // end loop over well separated blocks of smaller size
                  // (level)
            }     // end loop over all sublevels in  nonIntlist
        }
    };

  int localTimeEvalNumCalls = 0;
  // The copier function, it copies the value from the local array to the
  // parent blocks and it fills the actual parallel vector of the matrix
  // vector products.
  auto f_copier_Descend =
    [this, &matrVectProdD, &matrVectProdN, &localTimeEvalNumCalls](
      const DescendCopyData &copy_data) {
      std::vector<std::complex<double>> cache;
      AssertIndexRange(copy_data.blockId, blockLocalExpansionsKer1.size());
      this->blockLocalExpansionsKer1[copy_data.blockId].Add(
        copy_data.blockLocalExpansionKer1, cache);
      AssertIndexRange(copy_data.blockId, blockLocalExpansionsKer2.size());
      this->blockLocalExpansionsKer2[copy_data.blockId].Add(
        copy_data.blockLocalExpansionKer2, cache);

      for (types::global_dof_index i = 0; i < copy_data.nodesBlk1Ids.size();
           ++i)
        {
          matrVectProdN(copy_data.nodesBlk1Ids[i]) +=
            copy_data.matrVectorProductContributionKer1(i);
          matrVectProdD(copy_data.nodesBlk1Ids[i]) +=
            copy_data.matrVectorProductContributionKer2(i);
        }

      localTimeEvalNumCalls += copy_data.localTimeEvalNumCalls;
    };

  // so, here we loop over levels, starting form lower levels (bigger blocks)
  for (unsigned int level = 1; level < num_octree_levels + 1; level++)
    {
      types::global_dof_index startBlockLevel = startLevel[level];
      types::global_dof_index endBlockLevel   = endLevel[level];

      // to reduce computational cost, we decide to loop on the list of blocks
      // which contain at least one node (the local and multipole expansion
      // will be finally evaluated at the nodes positions)

      // TODO WE COULD SPLIT THIS LOOP OVER THE BLOCKS OVER ALL THE
      // PROCESSORS, THEN DO A TRILINOS.ADD WITH DIFFERENT MAPS. HERE WE NEED
      // A FULL ONE.
      DescendScratchData sample_scratch;
      DescendCopyData    sample_copy(this);

      // The Workstream has to run only if there are blocks in the current
      // level. Then it basically performs a loop over all the blocks in the
      // current level.
      // localTimeEvalNumCalls = 0;
      if (endBlockLevel >= startBlockLevel)
        {
          localTimeEvalNumCalls = 0;
          {
            Teuchos::TimeMonitor LocalTimer(*LocEval);
            WorkStream::run(dofs_filled_blocks[level].begin(),
                            dofs_filled_blocks[level].end(),
                            std::bind(f_worker_Descend,
                                      std::placeholders::_1,
                                      std::placeholders::_2,
                                      std::placeholders::_3,
                                      startBlockLevel),
                            f_copier_Descend,
                            sample_scratch,
                            sample_copy);
          }
          // add the actual number of calls done as if the timer was created
          // inside the workers
          for (auto i = 0; i < localTimeEvalNumCalls; ++i)
            {
              LocEval->incrementNumCalls();
            }
        }
    }

  // finally, when the loop over levels is done, we need to evaluate local
  // expansions of all childless blocks, at each block node(s). This is an
  // embarassingly parallel operation so it can be easily performed using
  // ThreadGroup. We also check the IndexSet to perform the mixed TBB-MPI
  // parallelisation.
  auto f_local_evaluation_tbb =
    [this, &matrVectProdD, &matrVectProdN, &support_points](
      blocked_range<types::global_dof_index> r) {
      std::vector<std::complex<double>> cache;
      for (types::global_dof_index kk = r.begin(); kk < r.end(); ++kk)
        {
          types::global_dof_index              block1Id = childlessList[kk];
          OctreeBlock<dim>                    *block1   = blocks[block1Id];
          std::vector<types::global_dof_index> nodesBlk1Ids =
            block1->GetBlockNodeList();

          // loop over nodes of block
          for (types::global_dof_index ii = 0; ii < nodesBlk1Ids.size();
               ii++) // loop over each node of block1
            {
              if (this_cpu_set.is_element(nodesBlk1Ids[ii]))
                {
                  const Point<dim> &nodeBlk1 = support_points[nodesBlk1Ids[ii]];
                  matrVectProdD(nodesBlk1Ids[ii]) +=
                    (blockLocalExpansionsKer2[block1Id])
                      .Evaluate(nodeBlk1, cache);
                  matrVectProdN(nodesBlk1Ids[ii]) +=
                    (blockLocalExpansionsKer1[block1Id])
                      .Evaluate(nodeBlk1, cache);
                }
            } // end loop over nodes
        }
    };

  // nodes are being used with no sharing; only leaf nodes are expanded
  parallel_for(blocked_range<types::global_dof_index>(0,
                                                      childlessList.size(),
                                                      tbb_granularity),
               f_local_evaluation_tbb);

  pcout << "...done computing multipole matrix-vector products" << std::endl;
}


template <int dim>
void
BEMFMA<dim>::multipole_matr_vect_products_omp(
  const TrilinosWrappers::MPI::Vector &phi_values,
  const TrilinosWrappers::MPI::Vector &dphi_dn_values,
  TrilinosWrappers::MPI::Vector &      matrVectProdN,
  TrilinosWrappers::MPI::Vector &      matrVectProdD) const
{
  pcout << "Computing multipole matrix-vector products (OpenMP)..."
        << std::endl;
  Teuchos::TimeMonitor LocalTimer(*MatrVec);

  // and here we compute the direct integral contributions (stored in two
  // sparse matrices)
  prec_neumann_matrix.vmult(matrVectProdN, phi_values);
  prec_dirichlet_matrix.vmult(matrVectProdD, dphi_dn_values);

  // from here on, we compute the multipole expansions contributions
  // we start cleaning past sessions
  // store old values

  // TODO: would not need clear-resize if LocalExpansion's copy assignement is
  // correct
  blockLocalExpansionsKer1.clear();
  blockLocalExpansionsKer2.clear();

  blockLocalExpansionsKer1.resize(num_blocks);
  blockLocalExpansionsKer2.resize(num_blocks);

  // we declare some familiar variables that will be useful in the method
  std::vector<Point<dim>> support_points(fma_dh->n_dofs());
  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*fma_mapping,
                                                     *fma_dh,
                                                     support_points);

  // reset current expansions
#pragma omp parallel for
  for (unsigned int i = 0; i < this->num_blocks; ++i)
    {
      const auto blockCenter = this->blocks[i]->GetCenter();

      this->blockLocalExpansionsKer1[i] =
        LocalExpansion(this->trunc_order, blockCenter, &(this->assLegFunction));
      this->blockLocalExpansionsKer2[i] =
        LocalExpansion(this->trunc_order, blockCenter, &(this->assLegFunction));
    }

  // this scope is only for the local evaluation timer
  {
    Teuchos::TimeMonitor LocalTimer(*LocEval);
    // first, convert non-well separated blocks' expansions to the interacting
    // blocks
    unsigned int localTimeEvalNumCalls = 0;
#pragma omp parallel reduction(+ : localTimeEvalNumCalls)
    {
      // one cache for each thread, will be reused in the various stages
      std::vector<std::complex<double>> cache;
      // TODO: these can be probably deleted safely, writing directly on the
      // final objects since the work is well partitioned and the tasks have the
      // correct dependencies
      // LocalExpansion thread_blockLocalExpansionKer1;
      // LocalExpansion thread_blockLocalExpansionKer2;
      Vector<double> thread_matrVectorProductContributionKer1;
      Vector<double> thread_matrVectorProductContributionKer2;

#pragma omp single
      {
        // traverse the octree;
        for (unsigned int level = 1; level < this->num_octree_levels + 1;
             level++)
          {
            auto startBlockLevel = this->startLevel[level];
            auto endBlockLevel   = this->endLevel[level];

            // The Workstream has to run only if there are blocks in the current
            // level. Then it basically performs a loop over all the blocks in
            // the current level.
            for (const auto blockId : this->dofs_filled_blocks[level])
              {
                const OctreeBlock<dim> *block1 = this->blocks[blockId];

#pragma omp task depend(in                                     \
                        : this->blocks[block1->GetParentId()]) \
  depend(out                                                   \
         : this->blocks[blockId])                              \
    firstprivate(blockId, block1, level, startBlockLevel, endBlockLevel)
                {
                  const auto &nodesBlk1Ids = block1->GetBlockNodeList();

                  bool proceed = false;
                  for (auto idx : nodesBlk1Ids)
                    {
                      if (this->this_cpu_set.is_element(idx))
                        {
                          proceed = true;
                          break;
                        }
                    }

                  if (proceed)
                    {
                      AssertIndexRange(blockId,
                                       this->blockLocalExpansionsKer1.size());
                      Point<dim> center =
                        this->blockLocalExpansionsKer1[blockId].GetCenter();

                      // reset the local expansions and Vectors with the
                      // contributions
                      /*
                      LocalExpansion dummy(this->trunc_order,
                                           center,
                                           &(this->assLegFunction));
                      */
                      // thread_blockLocalExpansionKer1 = dummy;
                      // thread_blockLocalExpansionKer2 = dummy;

                      thread_matrVectorProductContributionKer1.reinit(
                        nodesBlk1Ids.size());
                      thread_matrVectorProductContributionKer2.reinit(
                        nodesBlk1Ids.size());

                      // init with the expansion from the parent
                      types::global_dof_index parentId = block1->GetParentId();
                      AssertIndexRange(parentId,
                                       this->blockLocalExpansionsKer1.size());
                      // thread_blockLocalExpansionKer1.Add(
                      this->blockLocalExpansionsKer1[blockId].Add(
                        this->blockLocalExpansionsKer1[parentId], cache);
                      AssertIndexRange(parentId,
                                       this->blockLocalExpansionsKer2.size());
                      // thread_blockLocalExpansionKer2.Add(
                      this->blockLocalExpansionsKer2[blockId].Add(
                        this->blockLocalExpansionsKer2[parentId], cache);

                      for (unsigned int subLevel = 0;
                           subLevel < block1->NumNearNeighLevels();
                           subLevel++)
                        {
                          // for each level, for each non interacting block in
                          // that level
                          const auto &nonIntList =
                            block1->GetNonIntList(subLevel);

                          // first, combine the expansions from blocks at the
                          // current level

                          // we start converting into local
                          // expansions, all the multipole expansions of all
                          // the blocks in nonIntList, that are of the same
                          // size (level) of the current block. To perform the
                          // conversion, we use another member of the
                          // LocalExpansion class. Note that all the
                          // contributions to the integrals of blocks bigger
                          // than current block had already been considered in
                          // the direct integrals method (the bounds of the
                          // local and multipole expansion do not hold in such
                          // case)
                          for (auto pos2 =
                                 nonIntList.lower_bound(startBlockLevel);
                               pos2 != nonIntList.upper_bound(endBlockLevel);
                               pos2++) // loop over NNs of NNs and add them to
                                       // intList
                            {
                              types::global_dof_index block2Id = *pos2;

                              AssertIndexRange(
                                block2Id,
                                this->blockMultipoleExpansionsKer1.size());
                              // thread_blockLocalExpansionKer1.Add(
                              this->blockLocalExpansionsKer1[blockId].Add(
                                this->blockMultipoleExpansionsKer1[block2Id],
                                cache);
                              AssertIndexRange(
                                block2Id,
                                this->blockMultipoleExpansionsKer2.size());
                              // thread_blockLocalExpansionKer2.Add(
                              this->blockLocalExpansionsKer2[blockId].Add(
                                this->blockMultipoleExpansionsKer2[block2Id],
                                cache);
                            }

                          // then, expansions of blocks from deeper levels
                          // (that is, smaller blocks)

                          // we now have to loop over blocks in the nonIntList
                          // that are smaller than the current blocks: in this
                          // case the bound for the conversion of a multipole
                          // into local expansion does not hold, but the bound
                          // for the evaluation of the multipole expansions
                          // does hold: thus, we will simply evaluate the
                          // multipole expansions of such blocks for each node
                          // in the block
                          for (auto pos2 =
                                 nonIntList.upper_bound(endBlockLevel);
                               pos2 != nonIntList.end();
                               pos2++)
                            {
                              types::global_dof_index block2Id = *pos2;

                              for (types::global_dof_index ii = 0;
                                   ii < nodesBlk1Ids.size();
                                   ii++)
                                {
                                  const Point<dim> &nodeBlk1 =
                                    support_points[nodesBlk1Ids[ii]];

                                  // TODO: if going directly on the
                                  // this->matrVector* objs, can ditch ii
                                  thread_matrVectorProductContributionKer1(
                                    ii) +=
                                    this->blockMultipoleExpansionsKer1[block2Id]
                                      .Evaluate(nodeBlk1, cache);
                                  thread_matrVectorProductContributionKer2(
                                    ii) +=
                                    this->blockMultipoleExpansionsKer2[block2Id]
                                      .Evaluate(nodeBlk1, cache);
                                }
                            }
                        }

                        // everything has been accumulated in the thread_*
                        // variables
                        /*
                        AssertIndexRange(blockId,
                                         this->blockLocalExpansionsKer1.size());
                        this->blockLocalExpansionsKer1[blockId].Add(
                          thread_blockLocalExpansionKer1, cache);
                        AssertIndexRange(blockId,
                                         this->blockLocalExpansionsKer2.size());
                        this->blockLocalExpansionsKer2[blockId].Add(
                          thread_blockLocalExpansionKer2, cache);
                        */

                        // synch here shouldn't be needed; but on an i7, using
                        // 2mpi 2omp there's an error in the Petra vector add
                        // which does not manifest with Xmpi 1omp (nor with 1mpi
                        // Xomp, but that seems more like a fluke)
                        // is the add method somehow conflicting across threads?
#pragma omp critical(matrVectProdN_add)
                      {
                        matrVectProdN.add(
                          nodesBlk1Ids,
                          thread_matrVectorProductContributionKer1);
                      }
#pragma omp critical(matrVectProdD_add)
                      {
                        matrVectProdD.add(
                          nodesBlk1Ids,
                          thread_matrVectorProductContributionKer2);
                      }
                    }
                }
              }
          }
      }
    }

    // add the actual number of calls done as if the timer was created
    // inside the workers
    for (unsigned int i = 0; i < localTimeEvalNumCalls; ++i)
      {
        LocEval->incrementNumCalls();
      }
  }

// finally, when the loop over levels is done, we need to evaluate local
// expansions of all childless blocks, at each block node(s). This is an
// embarassingly parallel operation so it can be easily performed using
// ThreadGroup. We also check the IndexSet to perform the mixed TBB-MPI
// parallelisation.
#pragma omp parallel
  {
    std::vector<std::complex<double>> cache;

#pragma omp for
    for (unsigned int i = 0; i < this->childlessList.size(); ++i)
      {
        types::global_dof_index blockId = this->childlessList[i];
        const OctreeBlock<dim> *block1  = this->blocks[blockId];

        const auto &nodesBlk1Ids = block1->GetBlockNodeList();
        for (auto idx : nodesBlk1Ids)
          {
            if (this->this_cpu_set.is_element(idx))
              {
                const Point<dim> &nodeBlk1 = support_points[idx];

                matrVectProdN(idx) +=
                  this->blockLocalExpansionsKer1[blockId].Evaluate(nodeBlk1,
                                                                   cache);
                matrVectProdD(idx) +=
                  this->blockLocalExpansionsKer2[blockId].Evaluate(nodeBlk1,
                                                                   cache);
              }
          }
      }
  }

  pcout << "...done computing multipole matrix-vector products" << std::endl;
}

// this method computes the preconditioner needed for the GMRES:
// to do it, it needs to receive the alpha vector from the bem_problem
// class, along with the constraint matrix of the bem problem
template <int dim>
TrilinosWrappers::PreconditionILU &
BEMFMA<dim>::FMA_preconditioner(const TrilinosWrappers::MPI::Vector &alpha,
                                AffineConstraints<double> &          c)
{
#ifdef _OPENMP
  return FMA_preconditioner_omp(alpha, c);
#else
  return FMA_preconditioner_tbb(alpha, c);
#endif
}

template <int dim>
TrilinosWrappers::PreconditionILU &
BEMFMA<dim>::FMA_preconditioner_tbb(
  const TrilinosWrappers::MPI::Vector &alpha,
  AffineConstraints<double>           &c) // TO BE CHANGED!!!
{
  pcout << "Computing FMA preconditioner (TBB)" << std::endl;
  Teuchos::TimeMonitor LocalTimer(*PrecondTime);
  // the final preconditioner (with constraints) has a slightly different
  // sparsity pattern with respect to the non constrained one. we must here
  // initialize such sparsity pattern
  final_prec_sparsity_pattern.reinit(this_cpu_set,
                                     mpi_communicator,
                                     (types::global_dof_index)125 *
                                       fma_dh->get_fe().dofs_per_cell);

  // As before we will use the captures to simplify the calls of the worker
  // copier mechanism. For this reason we build an empty scratch structure.
  struct PrecScratch
  {};

  // We need to fill a vector that memorise the entries of the row associated
  // with the index we are treating.
  struct PrecCopy
  {
    PrecCopy()
      : row(numbers::invalid_unsigned_int)
      , sparsity_row(0){};

    PrecCopy(const PrecCopy &in_copy)
    {
      row          = in_copy.row;
      sparsity_row = in_copy.sparsity_row;
    };

    types::global_dof_index              row;
    std::vector<types::global_dof_index> sparsity_row;
  };

  auto f_worker_prec = [this, &c](IndexSet::ElementIterator iter,
                                  PrecScratch &,
                                  PrecCopy &copy_data) {
    unsigned int i = *iter;
    copy_data.sparsity_row.clear();

    copy_data.row = i;
    if (c.is_constrained(i))
      {
        // constrained nodes entries are taken from the bem problem
        // constraint matrix
        copy_data.sparsity_row.push_back(i);
        const std::vector<std::pair<types::global_dof_index, double>> *entries =
          c.get_constraint_entries(i);
        for (types::global_dof_index j = 0; j < entries->size(); ++j)
          {
            copy_data.sparsity_row.push_back((*entries)[j].first);
          }
      }
    else
      {
        // other nodes entries are taken from the unconstrained
        // preconditioner matrix
        for (unsigned int j = 0; j < fma_dh->n_dofs(); ++j)
          {
            if (this->init_prec_sparsity_pattern.exists(i, j))
              {
                copy_data.sparsity_row.push_back(j);
              }
          }
      }
  };

  // We only need a for cycle to add the indices to the sparsity pattern.
  auto f_copier_prec = [this](const PrecCopy &copy_data) {
    this->final_prec_sparsity_pattern.add_entries(
      copy_data.row,
      copy_data.sparsity_row.begin(),
      copy_data.sparsity_row.end(),
      false);
  };

  PrecCopy    foo_copy;
  PrecScratch foo_scratch;

  // The following Workstream replaces a for cycle on all dofs to check all
  // the constraints. WorkStream::run(0, fma_dh->n_dofs(), f_worker_prec,
  // f_copier_prec, foo_scratch, foo_copy);
  WorkStream::run(this_cpu_set.begin(),
                  this_cpu_set.end(),
                  f_worker_prec,
                  f_copier_prec,
                  foo_scratch,
                  foo_copy);

  final_prec_sparsity_pattern.compress();
  pcout << "Sparsity pattern nonzeros: "
        << final_prec_sparsity_pattern.n_nonzero_elements() << std::endl;
  final_preconditioner.reinit(final_prec_sparsity_pattern);

  // now we assemble the final preconditioner matrix: the loop works
  // exactly like the previous one

  // We need a worker function that fills the final sparisty pattern once its
  // sparsity pattern has been set up. In this case no race condition occurs
  // in the worker so we can let it copy in the global memory.
  auto f_sparsity_filler_tbb = [this, &c](unsigned int pos_begin,
                                          unsigned int pos_end) {
    for (unsigned int iter = pos_begin; iter != pos_end; ++iter)
      {
        unsigned int i = this->this_cpu_set.nth_index_in_set(iter);
        if (c.is_constrained(i))
          {
            final_preconditioner.set(i, i, 1);
            // constrainednodes entries are taken from the bem problem
            // constraint matrix
            for (const auto &entry : *c.get_constraint_entries(i))
              {
                final_preconditioner.set(i, entry.first, entry.second);
              }
          }
        else
          {
            // other nodes entries are taken from the unconstrained
            // preconditioner matrix
            for (unsigned int j = 0; j < fma_dh->n_dofs(); ++j)
              {
                // QUI CHECK SU NEUMANN - DIRICHLET PER METTERE A POSTO,
                // tanto lui già conosce le matrici.
                if (init_prec_sparsity_pattern.exists(i, j))
                  {
                    final_preconditioner.set(i, j, init_preconditioner(i, j));
                  }
              }
          }
      }
  };

  parallel::apply_to_subranges(0,
                               this_cpu_set.n_elements(),
                               f_sparsity_filler_tbb,
                               tbb_granularity);

  // The compress operation makes all the vectors on different processors
  // compliant.
  final_preconditioner.compress(VectorOperation::insert);

  // In order to add alpha we can again use the parallel_for strategy.
  auto f_alpha_adder_tbb = [this, &c, &alpha](unsigned int pos_begin,
                                              unsigned int pos_end) {
    for (unsigned int iter = pos_begin; iter != pos_end; ++iter)
      {
        unsigned int i = this->this_cpu_set.nth_index_in_set(iter);
        if ((*(dirichlet_nodes))(i) == 0 && !(c.is_constrained(i)))
          {
            final_preconditioner.add(i, i, alpha(i));
          }
        else // this is just to avoid a deadlock. we need a better strategy
          {
            final_preconditioner.add(i, i, 0);
          }
      }
  };

  parallel::apply_to_subranges(0,
                               this_cpu_set.n_elements(),
                               f_alpha_adder_tbb,
                               tbb_granularity);

  final_preconditioner.compress(VectorOperation::add);
  final_preconditioner.compress(VectorOperation::insert);

  // Finally we can initialize the ILU final preconditioner.
  preconditioner.initialize(final_preconditioner);
  pcout << "...done computing FMA preconditioner" << std::endl;

  return preconditioner;
}

template <int dim>
TrilinosWrappers::PreconditionILU &
BEMFMA<dim>::FMA_preconditioner_omp(
  const TrilinosWrappers::MPI::Vector &alpha,
  AffineConstraints<double> &          c) // TO BE CHANGED!!!
{
  pcout << "Computing FMA preconditioner (OpenMP)" << std::endl;
  Teuchos::TimeMonitor LocalTimer(*PrecondTime);

  auto max_threads = omp_get_max_threads();

  auto dofs_per_cell = fma_dh->get_fe().dofs_per_cell;
  // the final preconditioner (with constraints) has a slightly different
  // sparsity pattern with respect to the non constrained one. we must here
  // initialize such sparsity pattern
  final_prec_sparsity_pattern.reinit(this_cpu_set,
                                     mpi_communicator,
                                     (types::global_dof_index)125 *
                                       dofs_per_cell);

#pragma omp parallel for schedule(dynamic)
  for (unsigned int i = 0; i < this->this_cpu_set.n_elements(); ++i)
    {
      unsigned int idx = this->this_cpu_set.nth_index_in_set(i);
      if (c.is_constrained(idx))
        {
          this->final_prec_sparsity_pattern.add(idx, idx);
          // constrained nodes entries are taken from the bem problem
          // constraint matrix
          for (const auto &pair : *c.get_constraint_entries(idx))
            {
#pragma omp critical
              {
                this->final_prec_sparsity_pattern.add(idx, pair.first);
              }
            }
        }
      else
        {
          // other nodes entries are taken from the unconstrained
          // preconditioner matrix
          for (unsigned int j = 0; j < this->fma_dh->n_dofs(); ++j)
            {
              if (this->init_prec_sparsity_pattern.exists(idx, j))
                {
#pragma omp critical
                  {
                    this->final_prec_sparsity_pattern.add(idx, j);
                  }
                }
            }
        }
    }

  final_prec_sparsity_pattern.compress();
  pcout << "Sparsity pattern nonzeros: "
        << final_prec_sparsity_pattern.n_nonzero_elements() << std::endl;
  final_preconditioner.reinit(final_prec_sparsity_pattern);

  // now we assemble the final preconditioner matrix: the loop works
  // exactly like the previous one

  // We need a worker function that fills the final sparisty pattern once its
  // sparsity pattern has been set up. In this case no race condition occurs
  // in the worker so we can let it copy in the global memory.

#pragma omp parallel for schedule(dynamic)
  for (unsigned int i = 0; i < this->this_cpu_set.n_elements(); ++i)
    {
      unsigned int idx = this->this_cpu_set.nth_index_in_set(i);
      if (c.is_constrained(idx))
        {
          this->final_preconditioner.set(idx, idx, 1);
          // constrainednodes entries are taken from the bem problem
          // constraint matrix
          for (const auto &pair : *c.get_constraint_entries(idx))
            {
              this->final_preconditioner.set(idx, pair.first, pair.second);
            }
        }
      else
        {
          // other nodes entries are taken from the unconstrained
          // preconditioner matrix
          for (unsigned int j = 0; j < this->fma_dh->n_dofs(); ++j)
            {
              if (this->init_prec_sparsity_pattern.exists(idx, j))
                {
                  this->final_preconditioner.set(
                    idx, j, this->init_preconditioner(idx, j));
                }
            }
        }
    }

  // The compress operation makes all the vectors on different processors
  // compliant.
  final_preconditioner.compress(VectorOperation::insert);

#pragma omp parallel for schedule(dynamic)
  for (unsigned int i = 0; i < this->this_cpu_set.n_elements(); ++i)
    {
      unsigned int idx = this->this_cpu_set.nth_index_in_set(i);
      if ((*(this->dirichlet_nodes))(idx) == 0 && !(c.is_constrained(idx)))
        {
          this->final_preconditioner.add(idx, idx, alpha(idx));
        }
      else // this is just to avoid a deadlock. we need a better strategy
        {
          this->final_preconditioner.add(idx, idx, 0);
        }
    }

  final_preconditioner.compress(VectorOperation::add);
  final_preconditioner.compress(VectorOperation::insert);

  // Finally we can initialize the ILU final preconditioner.
  preconditioner.initialize(final_preconditioner);
  pcout << "...done computing FMA preconditioner" << std::endl;

  omp_set_num_threads(max_threads);

  return preconditioner;
}

template <int dim>
void
BEMFMA<dim>::compute_geometry_cache()
{
  pcout << "Generating geometry cache..." << std::endl;

  FESystem<dim - 1, dim>   gradient_fe(fma_dh->get_fe(), dim);
  DoFHandler<dim - 1, dim> gradient_dh(fma_dh->get_triangulation());

  std::vector<Point<dim>> support_points(fma_dh->n_dofs());

  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*fma_mapping,
                                                     *fma_dh,
                                                     support_points);

  std::vector<types::global_dof_index> dofs(fma_dh->get_fe().dofs_per_cell);
  std::vector<types::global_dof_index> gradient_dofs(
    fma_dh->get_fe().dofs_per_cell);
  // mappa che associa ad ogni dof le celle cui esso appartiene
  dof_to_elems.clear();

  // mappa che associa ad ogni gradient dof le celle cui esso appartiene
  gradient_dof_to_elems.clear();

  // vettore che associa ad ogni gradient dof la sua componente
  gradient_dof_components.clear();
  gradient_dof_components.resize(gradient_dh.n_dofs());

  // mappa che associa ad ogni cella un set contenente le celle circostanti
  // elem_to_surr_elems.clear();

  // for the gradient dofs finding coupled
  // dofs is a little bit difficult, as the
  // gradient is a vectorial function: usually
  // in such case the dofs are numbered
  // so that for each support point dim
  // consecutive dofs represent each component
  // of the vector field: in this case
  // (and only in this case) the following
  // piece of code works

  // the only coupling is the index check
  cell_it gradient_cell = gradient_dh.begin_active(),
          gradient_endc = gradient_dh.end();
  cell_it cell          = fma_dh->begin_active();
  for (; gradient_cell != gradient_endc; ++cell, ++gradient_cell)
    {
      Assert(cell->index() == gradient_cell->index(), ExcInternalError());

      cell->get_dof_indices(dofs);
      for (unsigned int j = 0; j < fma_dh->get_fe().dofs_per_cell; ++j)
        {
          dof_to_elems[dofs[j]].push_back(cell);
        }

      gradient_cell->get_dof_indices(gradient_dofs);
      for (unsigned int j = 0; j < gradient_fe.dofs_per_cell; ++j)
        {
          gradient_dof_to_elems[gradient_dofs[j]].push_back(gradient_cell);
          gradient_dof_components[gradient_dofs[j]] =
            gradient_fe.system_to_component_index(j).first;
        }
    }

  // TODO: deprecated
  // qui viene creata la mappa dei elmenti che circondano ciascun elemento
  // for (cell = fma_dh->begin_active(); cell != endc; ++cell)
  //   {
  //     cell->get_dof_indices(dofs);
  //     for (unsigned int j = 0; j < fma_dh->get_fe().dofs_per_cell; ++j)
  //       {
  //         // std::set<types::global_dof_index>
  //         const auto &duplicates = (*double_nodes_set)[dofs[j]];
  //         for (auto pos = duplicates.begin(); pos != duplicates.end();
  //         pos++)
  //           {
  //             /*
  //             std::vector<cell_it> dof_cell_list = dof_to_elems[*pos];
  //             for (unsigned int k = 0; k < dof_cell_list.size(); ++k)
  //               {
  //                 elem_to_surr_elems[cell].insert(dof_cell_list[k]);
  //               }
  //             */
  //             const auto &dof_cell_list = dof_to_elems[*pos];
  //             elem_to_surr_elems[cell].insert(dof_cell_list.begin(),
  //                                             dof_cell_list.end());
  //           }
  //       }
  //   }

  pcout << "...done" << std::endl;
}

// The following is the function
// which creates the octree blocking
// for the fast multipole algorithm
template <int dim>
void
BEMFMA<dim>::generate_octree_blocking()
{
  pcout << "Generating octree blocking... " << std::endl;
  Teuchos::TimeMonitor LocalTimer(*ListCreat);

  std::vector<Point<dim>> support_points(fma_dh->n_dofs());
  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*fma_mapping,
                                                     *fma_dh,
                                                     support_points);

  // !!!TO BE CHANGED
  quadrature = std::make_shared<QGauss<dim - 1>>(quadrature_order);
  FEValues<dim - 1, dim> fe_v(*fma_mapping,
                              fma_dh->get_fe(),
                              *quadrature,
                              update_values | update_normal_vectors |
                                update_quadrature_points | update_JxW_values);

  double max_coor_value = 0;

  for (types::global_dof_index i = 0; i < fma_dh->n_dofs(); i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          max_coor_value =
            std::max(max_coor_value, std::abs(support_points[i](j)));
        }
    }

  if (blocks.size() > 0)
    {
      for (types::global_dof_index ii = 0; ii < num_blocks; ii++)
        {
          delete blocks[ii];
        }
    }

  types::global_dof_index maxNumBlocks =
    num_octree_levels * fma_dh->get_triangulation().n_active_cells() *
    fe_v.n_quadrature_points;

  blocks.clear();
  blocks.reserve(maxNumBlocks);
  blocks.resize(maxNumBlocks);

  types::global_dof_index blocksCount = 0;
  startLevel.resize(num_octree_levels + 1);
  endLevel.resize(num_octree_levels + 1);

  parentList.clear();
  parentList.resize(num_octree_levels + 1);
  parentList[0].push_back(0);

  childlessList.clear();
  types::global_dof_index numChildless = 0;
  numParent.resize(num_octree_levels + 1);

  // qui di seguito vengono reinizializzate strutture utili al multipolo

  // mappa che associa ad ogni dof un vettore con i blocchi cui essa
  // appartiene per ogni livello
  // TODO: validate
  // dof_to_block.clear();

  // mappa che associa ad ogni quad point un vettore con i blocchi cui essa
  // appartiene per ogni livello
  quad_point_to_block.clear();

  // vettore di vettori contenente per ogni livello, gli ids dei blocchi
  // contenenti almeno un dof
  dofs_filled_blocks.clear();

  // vettore di vettori contenente per ogni livello, gli ids dei blocchi
  // contenenti almeno un quad point
  // TODO: validate
  // quad_points_filled_blocks.clear();

  quadPoints.clear();
  quadNormals.clear();
  quadShapeFunValues.clear();
  quadJxW.clear();

  dofs_filled_blocks.resize(num_octree_levels + 1);

  // TODO: validate
  //  quad_points_filled_blocks.resize(num_octree_levels + 1);

  for (unsigned int ii = 0; ii < num_octree_levels + 1; ii++)
    {
      numParent[ii] = 0;
    }

  Point<dim> pMin;
  for (int i = 0; i < dim; i++)
    {
      pMin(i) = -1.1 * max_coor_value;
    }

  // delta is the edge length of the cube
  double delta = 2.2 * max_coor_value;

  OctreeBlock<dim> *block = new OctreeBlock<dim>(0, 0, pMin, delta);

  std::vector<types::global_dof_index> local_dof_indices(
    fma_dh->get_fe().dofs_per_cell);

  for (const auto &cell : fma_dh->active_cell_iterators())
    {
      fe_v.reinit(cell);
      const unsigned int n_q_points = fe_v.n_quadrature_points;
      quadPoints[cell]              = fe_v.get_quadrature_points();
      quadNormals[cell]             = fe_v.get_normal_vectors();
      quadJxW[cell].resize(n_q_points);
      quadShapeFunValues[cell].resize(n_q_points);
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          quadJxW[cell][q] = fe_v.JxW(q);
          for (unsigned int j = 0; j < fma_dh->get_fe().dofs_per_cell; ++j)
            {
              quadShapeFunValues[cell][q].push_back(fe_v.shape_value(j, q));
            }
        }

      quad_point_to_block[cell].resize(n_q_points);
      for (unsigned int j = 0; j < n_q_points; ++j)
        {
          block->AddQuadPoint(cell, j);
          quad_point_to_block[cell][j].push_back(0);
        }

      cell->get_dof_indices(local_dof_indices);
      for (unsigned int j = 0; j < fma_dh->get_fe().dofs_per_cell; ++j)
        {
          dof_to_elems[local_dof_indices[j]].push_back(cell);
        }
    }

  for (types::global_dof_index ii = 0; ii < fma_dh->n_dofs(); ii++)
    {
      block->AddNode(ii);
      // TODO: validate
      // dof_to_block[ii].push_back(0);
    }

  blocks[0]    = block;
  numParent[0] = 1;

  types::global_dof_index quadPointsInChildless = 0;
  types::global_dof_index nodesInChildless      = 0;
  for (unsigned int level = 1; level < num_octree_levels + 1; level++)
    {
      types::global_dof_index quadPointsCheck = quadPointsInChildless;
      types::global_dof_index nodesCheck      = nodesInChildless;
      delta /= 2.;

      for (types::global_dof_index kk = 0; kk < numParent[level - 1]; kk++)
        {
          types::global_dof_index jj     = parentList[level - 1][kk];
          OctreeBlock<dim> *      parent = blocks[jj];

          pMin = parent->GetPMin();
          unsigned int num_children_per_block =
            int(pow((double)2, (double)dim));
          std::vector<OctreeBlock<dim> *> children(num_children_per_block);

          if (dim == 3)
            {
              children[0] = new OctreeBlock<dim>(level, jj, pMin, delta);
              children[1] = new OctreeBlock<dim>(
                level, jj, pMin + Point<dim>(delta, 0., 0.), delta);
              children[2] = new OctreeBlock<dim>(
                level, jj, pMin + Point<dim>(delta, delta, 0.), delta);
              children[3] = new OctreeBlock<dim>(
                level, jj, pMin + Point<dim>(0., delta, 0.), delta);
              children[4] = new OctreeBlock<dim>(
                level, jj, pMin + Point<dim>(0., 0., delta), delta);
              children[5] = new OctreeBlock<dim>(
                level, jj, pMin + Point<dim>(delta, 0., delta), delta);
              children[6] = new OctreeBlock<dim>(
                level, jj, pMin + Point<dim>(delta, delta, delta), delta);
              children[7] = new OctreeBlock<dim>(
                level, jj, pMin + Point<dim>(0., delta, delta), delta);
            }

          if (dim == 2)
            {
              children[0] = new OctreeBlock<dim>(level, jj, pMin, delta);
              children[1] = new OctreeBlock<dim>(level,
                                                 jj,
                                                 pMin + Point<dim>(delta, 0.),
                                                 delta);
              children[2] = new OctreeBlock<dim>(
                level, jj, pMin + Point<dim>(delta, delta), delta);
              children[3] = new OctreeBlock<dim>(level,
                                                 jj,
                                                 pMin + Point<dim>(0., delta),
                                                 delta);
            }

          // std::map<cell_it, std::vector<types::global_dof_index>>
          const auto &blockQuadPointsList = parent->GetBlockQuadPointsList();

          // std::vector<types::global_dof_index>
          const auto &blockNodeList = parent->GetBlockNodeList();

          if (dim == 3)
            {
              for (types::global_dof_index i = 0; i < blockNodeList.size(); i++)
                {
                  Point<dim> node = support_points[blockNodeList[i]];
                  // assegnamento nodi del blocco padre ai blocchi figli

                  if (node(2) <= parent->GetPMin()(2) + delta)
                    {
                      if (node(1) <= parent->GetPMin()(1) + delta)
                        {
                          if (node(0) <= parent->GetPMin()(0) + delta)
                            {
                              children[0]->AddNode(blockNodeList[i]);
                            }
                          else
                            {
                              children[1]->AddNode(blockNodeList[i]);
                            }
                        }
                      else
                        {
                          if (node(0) <= parent->GetPMin()(0) + delta)
                            {
                              children[3]->AddNode(blockNodeList[i]);
                            }
                          else
                            {
                              children[2]->AddNode(blockNodeList[i]);
                            }
                        }
                    }
                  else
                    {
                      if (node(1) <= parent->GetPMin()(1) + delta)
                        {
                          if (node(0) <= parent->GetPMin()(0) + delta)
                            {
                              children[4]->AddNode(blockNodeList[i]);
                            }
                          else
                            {
                              children[5]->AddNode(blockNodeList[i]);
                            }
                        }
                      else
                        {
                          if (node(0) <= parent->GetPMin()(0) + delta)
                            {
                              children[7]->AddNode(blockNodeList[i]);
                            }
                          else
                            {
                              children[6]->AddNode(blockNodeList[i]);
                            }
                        }
                    } // fine assegnazione nodi del padre ai blocchi figli
                }     // fine loop nodi del blocco

              for (auto it = blockQuadPointsList.begin();
                   it != blockQuadPointsList.end();
                   it++)
                {
                  for (types::global_dof_index pp = 0; pp < (*it).second.size();
                       pp++)
                    {
                      Point<dim> quadPoint =
                        quadPoints[(*it).first][(*it).second[pp]];
                      // assegnamento punti quadratura del blocco padre ai
                      // blocchi figli
                      if (quadPoint(2) <= parent->GetPMin()(2) + delta)
                        {
                          if (quadPoint(1) <= parent->GetPMin()(1) + delta)
                            {
                              if (quadPoint(0) <= parent->GetPMin()(0) + delta)
                                {
                                  children[0]->AddQuadPoint((*it).first,
                                                            (*it).second[pp]);
                                }
                              else
                                {
                                  children[1]->AddQuadPoint((*it).first,
                                                            (*it).second[pp]);
                                }
                            }
                          else
                            {
                              if (quadPoint(0) <= parent->GetPMin()(0) + delta)
                                {
                                  children[3]->AddQuadPoint((*it).first,
                                                            (*it).second[pp]);
                                }
                              else
                                {
                                  children[2]->AddQuadPoint((*it).first,
                                                            (*it).second[pp]);
                                }
                            }
                        }
                      else
                        {
                          if (quadPoint(1) <= parent->GetPMin()(1) + delta)
                            {
                              if (quadPoint(0) <= parent->GetPMin()(0) + delta)
                                {
                                  children[4]->AddQuadPoint((*it).first,
                                                            (*it).second[pp]);
                                }
                              else
                                {
                                  children[5]->AddQuadPoint((*it).first,
                                                            (*it).second[pp]);
                                }
                            }
                          else
                            {
                              if (quadPoint(0) <= parent->GetPMin()(0) + delta)
                                {
                                  children[7]->AddQuadPoint((*it).first,
                                                            (*it).second[pp]);
                                }
                              else
                                {
                                  children[6]->AddQuadPoint((*it).first,
                                                            (*it).second[pp]);
                                }
                            }
                        } // fine assegnazione punti quadratura del padre ai
                          // blocchi figli
                    }
                } // fine loop punti quadratura del blocco

              for (unsigned int j = 0; j < num_children_per_block; j++)
                {
                  if (children[j]->GetBlockNodeList().size() +
                        children[j]->GetBlockQuadPointsList().size() >
                      0)
                    {
                      blocksCount += 1;
                      // TODO: is there a more compact pattern?
                      /* //reference
                      blocks[blocksCount] = new OctreeBlock<dim>();
                      blocks[blocksCount]->CopyContent(children[j]);
                      delete children[j];
                      */
                      blocks[blocksCount] =
                        new OctreeBlock<dim>(std::move(*children[j]));
                      delete children[j];

                      parent->AddChild(blocksCount);
                      // std::map<cell_it,
                      // std::vector<types::global_dof_index>>
                      const auto &blockQuadPointsList =
                        blocks[blocksCount]->GetBlockQuadPointsList();
                      for (auto it = blockQuadPointsList.begin();
                           it != blockQuadPointsList.end();
                           it++)
                        {
                          // cell_it cell = (*it).first;
                          for (types::global_dof_index kk = 0;
                               kk < (*it).second.size();
                               kk++)
                            {
                              quad_point_to_block[(*it).first][(*it).second[kk]]
                                .push_back(blocksCount);
                            }
                        }
                      // TODO: deprecated
                      /*
                      // std::vector<types::global_dof_index>
                      const auto &blockNodesList =
                        blocks[jj]->GetBlockNodeList();
                      for (types::global_dof_index k = 0;
                          k < blockNodesList.size();
                          k++)
                        {
                          dof_to_block[blockNodesList[k]].push_back(jj);
                        }
                      */
                    }
                  else
                    {
                      delete children[j];
                    }
                } // fine loop sui blocchi figlio appena creati
            }     // fine ramo dim = 3 dell'if
          else
            {
              for (types::global_dof_index i = 0; i < blockNodeList.size(); i++)
                {
                  // assegnamento nodi del blocco padre ai blocchi figli
                  Point<dim> node = support_points[blockNodeList[i]];

                  if (node(1) <= parent->GetPMin()(1) + delta)
                    {
                      if (node(0) <= parent->GetPMin()(0) + delta)
                        {
                          children[0]->AddNode(blockNodeList[i]);
                        }
                      else
                        {
                          children[1]->AddNode(blockNodeList[i]);
                        }
                    }
                  else
                    {
                      if (node(0) <= parent->GetPMin()(0) + delta)
                        {
                          children[3]->AddNode(blockNodeList[i]);
                        }
                      else
                        {
                          children[2]->AddNode(blockNodeList[i]);
                        }
                    } // fine assegnazione blocchi del padre ai blocchi figli
                }     // fine loop nodi del blocco

              for (auto it = blockQuadPointsList.begin();
                   it != blockQuadPointsList.end();
                   it++)
                {
                  for (types::global_dof_index pp = 0; pp < (*it).second.size();
                       pp++)
                    {
                      // assegnamento quad points del blocco padre ai blocchi
                      // figli
                      Point<dim> quadPoint =
                        quadPoints[(*it).first][(*it).second[pp]];
                      if (quadPoint(1) <= parent->GetPMin()(1) + delta)
                        {
                          if (quadPoint(0) <= parent->GetPMin()(0) + delta)
                            {
                              children[0]->AddQuadPoint((*it).first,
                                                        (*it).second[pp]);
                            }
                          else
                            {
                              children[1]->AddQuadPoint((*it).first,
                                                        (*it).second[pp]);
                            }
                        }
                      else
                        {
                          if (quadPoint(0) <= parent->GetPMin()(0) + delta)
                            {
                              children[3]->AddQuadPoint((*it).first,
                                                        (*it).second[pp]);
                            }
                          else
                            {
                              children[2]->AddQuadPoint((*it).first,
                                                        (*it).second[pp]);
                            }
                        } // fine assegnazione blocchi del padre ai blocchi
                          // figli
                    }
                }

              for (unsigned int j = 0; j < num_children_per_block; j++)
                {
                  if (children[j]->GetBlockNodeList().size() +
                        children[j]->GetBlockQuadPointsList().size() >
                      0)
                    {
                      blocksCount += 1;
                      // TODO: is there a more compact pattern?
                      /* //reference
                      blocks[blocksCount] = new OctreeBlock<dim>();
                      blocks[blocksCount]->CopyContent(children[j]);
                      delete children[j];
                      */
                      blocks[blocksCount] =
                        new OctreeBlock<dim>(std::move(*children[j]));
                      delete children[j];

                      parent->AddChild(blocksCount);
                      // std::map<cell_it,
                      // std::vector<types::global_dof_index>>
                      const auto &blockQuadPointsList =
                        blocks[blocksCount]->GetBlockQuadPointsList();
                      for (auto it = blockQuadPointsList.begin();
                           it != blockQuadPointsList.end();
                           it++)
                        {
                          for (types::global_dof_index kk = 0;
                               kk < (*it).second.size();
                               kk++)
                            {
                              quad_point_to_block[(*it).first][(*it).second[kk]]
                                .push_back(blocksCount);
                            }
                        }

                      // TODO: validate
                      /*
                      // std::vector<types::global_dof_index>
                      const auto &blockNodesList =
                        blocks[jj]->GetBlockNodeList();
                      for (types::global_dof_index k = 0;
                            k < blockNodesList.size();
                            k++)
                        {
                          dof_to_block[blockNodesList[k]].push_back(jj);
                        }
                      */
                    }
                  else
                    {
                      delete children[j];
                    }
                } // fine loop sui blocchi figlio appena creati
            }     // fine ramo dim == 2 dell'if
        }         // fine loop blocchi livello precedente

      startLevel[level] = endLevel[level - 1] + 1;
      endLevel[level]   = blocksCount;

      // here we loop over the blocks of the newly created level and
      // we decide if each block is to be split again in the next level:
      // if it contains more
      // than a node or quad point, it will be placed in the parent list.
      // Instead, if it only contains a node or quad point, it goes in the
      // childless list, and not be refined any more. It is important to
      // account for the presence of double nodes: if not, blocks will be
      // always refined
      for (types::global_dof_index jj = startLevel[level];
           jj < endLevel[level] + 1;
           jj++)
        {
          // here we get the number of nodes in the block
          // std::vector<types::global_dof_index>
          const auto &nodesId       = blocks[jj]->GetBlockNodeList();
          double      blockNumNodes = 0.0;

          // now we compute the number of the nodes that are double of others
          for (types::global_dof_index kk = 0; kk < nodesId.size(); kk++)
            {
              blockNumNodes +=
                1.0 / (double((*double_nodes_set)[nodesId[kk]].size()));
            }

          // here we compute the number of quad points in the block
          int blockNumQuadPoints = 0;
          std::map<cell_it, std::vector<types::global_dof_index>>
            blockQuadPointsList = blocks[jj]->GetBlockQuadPointsList();
          typename std::map<cell_it,
                            std::vector<types::global_dof_index>>::iterator it;
          for (it = blockQuadPointsList.begin();
               it != blockQuadPointsList.end();
               it++)
            {
              blockNumQuadPoints += (int)(*it).second.size();
            }

          quadPointsCheck += blockNumQuadPoints;
          nodesCheck += (types::global_dof_index)blockNumNodes;
          // here we decide if a block is to be placed in the parent
          // or childless listì
          if (round(blockNumNodes) <= max_num_nodes_per_block)
            {
              numChildless += 1;
              childlessList.push_back(jj);
              quadPointsInChildless += blockNumQuadPoints;
              nodesInChildless += (types::global_dof_index)blockNumNodes;

              // TODO: deprecated
              /*
                      // if a block is childless, we must assign now the nodes
                 and quad
                      // points that belong to it for all the next levels
                      for (types::global_dof_index kk = 0; kk <
                 nodesId.size(); kk++)
                        {
                          for (unsigned int j = level + 1; j <
                 num_octree_levels
                 + 1; j++)
                            {
                              dof_to_block[nodesId[kk]].push_back(jj);
                            }
                        }
              */

              for (auto it = blockQuadPointsList.begin();
                   it != blockQuadPointsList.end();
                   it++)
                {
                  for (types::global_dof_index i = 0; i < (*it).second.size();
                       i++)
                    {
                      for (unsigned int j = level + 1;
                           j < num_octree_levels + 1;
                           j++)
                        {
                          quad_point_to_block[(*it).first][(*it).second[i]]
                            .push_back(jj);
                        }
                    }
                }
            }
          else
            {
              numParent[level] += 1;
              parentList[level].push_back(jj);
            }

          // let's update the list of node filled block
          if (blockNumNodes > 0)
            {
              dofs_filled_blocks[level].push_back(jj);
            }

          // let's update the list of quad point filled block
          // TODO: deprecated
          /*
          if (blockNumQuadPoints > 0)
            {
              quad_points_filled_blocks[level].push_back(jj);
            }
          */
        }

      pcout << " Total nodes at level " << level << " of " << num_octree_levels
            << " are " << nodesCheck << std::endl;
      pcout << " Total quad points at level " << level << " of "
            << num_octree_levels << " are " << quadPointsCheck << std::endl;
      pcout << " Blocks at level " << level << " of " << num_octree_levels
            << " are " << endLevel[level] - endLevel[level - 1] << " out of "
            << std::pow(8, level) << std::endl;
      pcout << " Total blocks up to level " << level << " of "
            << num_octree_levels << " are " << endLevel[level] + 1 << std::endl;
      pcout << std::endl;
    } // fine loop livelli

  childlessList.resize(childlessList.size() +
                       parentList[num_octree_levels].size());

  for (types::global_dof_index jj = 0;
       jj < parentList[num_octree_levels].size();
       jj++)
    {
      childlessList[numChildless + jj] = parentList[num_octree_levels][jj];
    }

  num_blocks = blocksCount + 1;

  pcout << "...done generating octree blocking" << std::endl;
  pcout << "Computing proximity lists for blocks" << std::endl;

  // ricerca blocchi nearest neighbors
  for (types::global_dof_index ii = startLevel[1]; ii < endLevel[1] + 1; ii++)
    {
      for (types::global_dof_index jj = startLevel[1]; jj < endLevel[1] + 1;
           jj++)
        {
          blocks[ii]->AddNearNeigh(0, jj);
        }
    }

  for (types::global_dof_index level = 2; level < num_octree_levels + 1;
       level++)
    {
      for (types::global_dof_index kk = startLevel[level];
           kk < endLevel[level] + 1;
           kk++)
        {
          OctreeBlock<dim> *block1 = blocks[kk];
          block1->AddNearNeigh(0, kk); // a block is NearNeigh of itself

          double     delta1 = block1->GetDelta();
          Point<dim> PMin1  = block1->GetPMin();
          Point<dim> Center1;
          for (unsigned int iii = 0; iii < dim; iii++)
            {
              Center1(iii) = delta1 / 2.;
            }
          Point<dim> PMax1 = 2. * Center1;
          PMax1 += PMin1;
          Center1 += PMin1;
          types::global_dof_index parentId = block1->GetParentId();
          // std::set<types::global_dof_index>
          const auto &parentNNeighs = blocks[parentId]->GetNearNeighs(0);

          // the nearest neighbors are searched among the father's nearest
          // neighbors children
          for (auto pos = parentNNeighs.begin(); pos != parentNNeighs.end();
               pos++)
            {
              if (blocks[*pos]->GetBlockChildrenNum() ==
                  0) // if a parent's near neigh is childless, he can be a
                     // near neigh: let's check
                {
                  types::global_dof_index block2Id = *pos;
                  OctreeBlock<dim>       *block2   = blocks[block2Id];
                  double                  delta2   = block2->GetDelta();
                  Point<dim>              PMin2    = block2->GetPMin();
                  Point<dim>              Center2;
                  for (unsigned int iii = 0; iii < dim; iii++)
                    {
                      Center2(iii) = delta2 / 2.;
                    }
                  Point<dim> PMax2 = 2. * Center2;
                  PMax2 += PMin2;
                  Center2 += PMin2;

                  if (dim == 3)
                    {
                      if ((fabs(PMin1(0) - PMax2(0)) <= TOLL) ||
                          (fabs(PMax1(0) - PMin2(0)) <= TOLL))
                        {
                          if ((PMin1(1) - TOLL <= PMax2(1)) &&
                              (PMax1(1) + TOLL >= PMin2(1)))
                            {
                              if ((PMin1(2) - TOLL <= PMax2(2)) &&
                                  (PMax1(2) + TOLL >= PMin2(2)))
                                {
                                  block1->AddNearNeigh(0, block2Id);
                                }
                            }
                        }

                      if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) ||
                          (fabs(PMax1(1) - PMin2(1)) <= TOLL))
                        {
                          if ((PMin1(0) - TOLL <= PMax2(0)) &&
                              (PMax1(0) + TOLL >= PMin2(0)))
                            {
                              if ((PMin1(2) - TOLL <= PMax2(2)) &&
                                  (PMax1(2) + TOLL >= PMin2(2)))
                                {
                                  block1->AddNearNeigh(0, block2Id);
                                }
                            }
                        }

                      if ((fabs(PMin1(2) - PMax2(2)) <= TOLL) ||
                          (fabs(PMax1(2) - PMin2(2)) <= TOLL))
                        {
                          if ((PMin1(1) - TOLL <= PMax2(1)) &&
                              (PMax1(1) + TOLL >= PMin2(1)))
                            {
                              if ((PMin1(0) - TOLL <= PMax2(0)) &&
                                  (PMax1(0) + TOLL >= PMin2(0)))
                                {
                                  block1->AddNearNeigh(0, block2Id);
                                }
                            }
                        }
                    } // fine caso dim ==3
                  else if (dim == 2)
                    {
                      if ((fabs(PMin1(0) - PMax2(0)) <= TOLL) ||
                          (fabs(PMax1(0) - PMin2(0)) <= TOLL))
                        {
                          if ((PMin1(1) - TOLL <= PMax2(1)) &&
                              (PMax1(1) + TOLL >= PMin2(1)))
                            {
                              block1->AddNearNeigh(0, block2Id);
                            }
                        }

                      if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) ||
                          (fabs(PMax1(1) - PMin2(1)) <= TOLL))
                        {
                          if ((PMin1(0) - TOLL <= PMax2(0)) &&
                              (PMax1(0) + TOLL >= PMin2(0)))
                            {
                              block1->AddNearNeigh(0, block2Id);
                            }
                        }
                    } // fine caso dim == 2
                }

              for (unsigned int ii = 0;
                   ii < blocks[*pos]->GetBlockChildrenNum();
                   ii++)
                {
                  types::global_dof_index block2Id =
                    blocks[*pos]->GetChildId(ii);
                  OctreeBlock<dim> *block2 = blocks[block2Id];
                  double            delta2 = block2->GetDelta();
                  Point<dim>        PMin2  = block2->GetPMin();
                  Point<dim>        Center2;
                  for (unsigned int iii = 0; iii < dim; iii++)
                    {
                      Center2(iii) = delta2 / 2.;
                    }
                  Point<dim> PMax2 = 2. * Center2;
                  PMax2 += PMin2;
                  Center2 += PMin2;

                  if (dim == 3)
                    {
                      if ((fabs(PMin1(0) - PMax2(0)) <= TOLL) ||
                          (fabs(PMax1(0) - PMin2(0)) <= TOLL))
                        {
                          if ((PMin1(1) - TOLL <= PMax2(1)) &&
                              (PMax1(1) + TOLL >= PMin2(1)))
                            {
                              if ((PMin1(2) - TOLL <= PMax2(2)) &&
                                  (PMax1(2) + TOLL >= PMin2(2)))
                                {
                                  block1->AddNearNeigh(0, block2Id);
                                }
                            }
                        }

                      if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) ||
                          (fabs(PMax1(1) - PMin2(1)) <= TOLL))
                        {
                          if ((PMin1(0) - TOLL <= PMax2(0)) &&
                              (PMax1(0) + TOLL >= PMin2(0)))
                            {
                              if ((PMin1(2) - TOLL <= PMax2(2)) &&
                                  (PMax1(2) + TOLL >= PMin2(2)))
                                {
                                  block1->AddNearNeigh(0, block2Id);
                                }
                            }
                        }

                      if ((fabs(PMin1(2) - PMax2(2)) <= TOLL) ||
                          (fabs(PMax1(2) - PMin2(2)) <= TOLL))
                        {
                          if ((PMin1(1) - TOLL <= PMax2(1)) &&
                              (PMax1(1) + TOLL >= PMin2(1)))
                            {
                              if ((PMin1(0) - TOLL <= PMax2(0)) &&
                                  (PMax1(0) + TOLL >= PMin2(0)))
                                {
                                  block1->AddNearNeigh(0, block2Id);
                                }
                            }
                        }
                    } // fine caso dim ==3
                  else if (dim == 2)
                    {
                      if ((fabs(PMin1(0) - PMax2(0)) <= TOLL) ||
                          (fabs(PMax1(0) - PMin2(0)) <= TOLL))
                        {
                          if ((PMin1(1) - TOLL <= PMax2(1)) &&
                              (PMax1(1) + TOLL >= PMin2(1)))
                            {
                              block1->AddNearNeigh(0, block2Id);
                            }
                        }

                      if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) ||
                          (fabs(PMax1(1) - PMin2(1)) <= TOLL))
                        {
                          if ((PMin1(0) - TOLL <= PMax2(0)) &&
                              (PMax1(0) + TOLL >= PMin2(0)))
                            {
                              block1->AddNearNeigh(0, block2Id);
                            }
                        }
                    } // fine caso dim == 2
                }     // fine loop sui figli di un nearest neighbor del padre
            }         // fine loop sui nearest neighbors del padre

          if ((block1->GetBlockChildrenNum() ==
               0)) // if the block is childless we must compute now its
                   // nearneigh at all residual levels
            {
              block1->SetNearNeighSize(num_octree_levels - level + 1);
              block1->SetIntListSize(num_octree_levels - level +
                                     1); // intList is a vector of sets with the
                                         // same number of members of nearNeigh
              block1->SetNonIntListSize(
                num_octree_levels - level +
                1); // nonIntList is a vector of sets with the same number of
                    // members of nearNeigh

              for (unsigned int subLevel = 1;
                   subLevel < num_octree_levels - level + 1;
                   subLevel++)
                {
                  // std::set<types::global_dof_index>
                  const auto &upperLevelNNeighs =
                    block1->GetNearNeighs(subLevel - 1);
                  for (auto pos = upperLevelNNeighs.begin();
                       pos != upperLevelNNeighs.end();
                       pos++)
                    {
                      if (blocks[*pos]->GetBlockChildrenNum() == 0)
                        { // if nearneigh is childless, it will stay a near
                          // neigh
                          block1->AddNearNeigh(subLevel, *pos);
                        }

                      for (unsigned int ii = 0;
                           ii < blocks[*pos]->GetBlockChildrenNum();
                           ii++)
                        {
                          types::global_dof_index block2Id =
                            blocks[*pos]->GetChildId(ii);
                          OctreeBlock<dim> *block2 = blocks[block2Id];
                          double            delta2 = block2->GetDelta();
                          Point<dim>        PMin2  = block2->GetPMin();
                          Point<dim>        Center2;
                          for (unsigned int iii = 0; iii < dim; iii++)
                            {
                              Center2(iii) = delta2 / 2.;
                            }
                          Point<dim> PMax2 = 2. * Center2;
                          PMax2 += PMin2;
                          Center2 += PMin2;

                          if (dim == 3)
                            {
                              if ((fabs(PMin1(0) - PMax2(0)) <= TOLL) ||
                                  (fabs(PMax1(0) - PMin2(0)) <= TOLL))
                                {
                                  if ((PMin1(1) - TOLL <= PMax2(1)) &&
                                      (PMax1(1) + TOLL >= PMin2(1)))
                                    {
                                      if ((PMin1(2) - TOLL <= PMax2(2)) &&
                                          (PMax1(2) + TOLL >= PMin2(2)))
                                        {
                                          block1->AddNearNeigh(subLevel,
                                                               block2Id);
                                        }
                                    }
                                }

                              if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) ||
                                  (fabs(PMax1(1) - PMin2(1)) <= TOLL))
                                {
                                  if ((PMin1(0) - TOLL <= PMax2(0)) &&
                                      (PMax1(0) + TOLL >= PMin2(0)))
                                    {
                                      if ((PMin1(2) - TOLL <= PMax2(2)) &&
                                          (PMax1(2) + TOLL >= PMin2(2)))
                                        {
                                          block1->AddNearNeigh(subLevel,
                                                               block2Id);
                                        }
                                    }
                                }

                              if ((fabs(PMin1(2) - PMax2(2)) <= TOLL) ||
                                  (fabs(PMax1(2) - PMin2(2)) <= TOLL))
                                {
                                  if ((PMin1(1) - TOLL <= PMax2(1)) &&
                                      (PMax1(1) + TOLL >= PMin2(1)))
                                    {
                                      if ((PMin1(0) - TOLL <= PMax2(0)) &&
                                          (PMax1(0) + TOLL >= PMin2(0)))
                                        {
                                          block1->AddNearNeigh(subLevel,
                                                               block2Id);
                                        }
                                    }
                                }
                            } // fine caso dim ==3
                          else if (dim == 2)
                            {
                              if ((fabs(PMin1(0) - PMax2(0)) <= TOLL) ||
                                  (fabs(PMax1(0) - PMin2(0)) <= TOLL))
                                {
                                  if ((PMin1(1) - TOLL <= PMax2(1)) &&
                                      (PMax1(1) + TOLL >= PMin2(1)))
                                    {
                                      block1->AddNearNeigh(subLevel, block2Id);
                                    }
                                }

                              if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) ||
                                  (fabs(PMax1(1) - PMin2(1)) <= TOLL))
                                {
                                  if ((PMin1(0) - TOLL <= PMax2(0)) &&
                                      (PMax1(0) + TOLL >= PMin2(0)))
                                    {
                                      block1->AddNearNeigh(subLevel, block2Id);
                                    }
                                }
                            } // fine caso dim == 2
                        }     // fine loop sui figli di ciascun nearest neighbor
                              // del blocco childless
                    } // fine loop sui nearest neighbors del blocco childless
                }     // fine loop sui subLevels (da quello del blocco childless
                      // all'ultimo)
            }         // fine if (il blocco e' childless?)
        }             // fine loop sui blocchi di un livello
    }                 // fine loop sui livelli

  // search for interaction list blocks (NearNeigh + NearNeighOfNearNeigh)
  // and for non interaction list blocks (nonIntList for a block B is composed
  // by blocks that are children of blocks being in intList of B's parent, but
  // are not in intList of B)
  for (types::global_dof_index ii = startLevel[1]; ii < endLevel[1] + 1;
       ii++) // at level 1, add all blocks to intList
    {
      for (types::global_dof_index jj = startLevel[1]; jj < endLevel[1] + 1;
           jj++)
        {
          blocks[ii]->AddBlockToIntList(0, jj);
        }
    }

  for (unsigned int level = 2; level < num_octree_levels + 1;
       level++) // loop over levels
    {
      for (unsigned int jj = startLevel[level]; jj < endLevel[level] + 1;
           jj++) // loop over blocks of each level
        {
          OctreeBlock<dim> *block1 = blocks[jj];
          for (unsigned int subLevel = 0;
               subLevel < block1->NumNearNeighLevels();
               subLevel++)
            {
              // std::set<types::global_dof_index>
              const auto &NNList = block1->GetNearNeighs(subLevel);

              for (auto pos1 = NNList.begin(); pos1 != NNList.end();
                   pos1++) // loop over blocks in NN list and get their NNs
                {
                  block1->AddBlockToIntList(subLevel, *pos1);
                }

              // std::vector<types::global_dof_index>
              const auto &nodeIds = block1->GetBlockNodeList();
              for (types::global_dof_index pp = 0; pp < nodeIds.size(); pp++)
                {
                  // std::set<types::global_dof_index>
                  const auto &doubleNodes = (*double_nodes_set)[nodeIds[pp]];
                  for (auto pos = doubleNodes.begin(); pos != doubleNodes.end();
                       pos++)
                    {
                      for (types::global_dof_index k = 0;
                           k < dof_to_elems[*pos].size();
                           k++)
                        {
                          cell_it cell = dof_to_elems[*pos][k];
                          for (unsigned int j = 0; j < quadPoints[cell].size();
                               j++)
                            {
                              block1->AddBlockToIntList(
                                subLevel,
                                quad_point_to_block[cell][j][level + subLevel]);
                            }
                        }
                    }
                }
            }

          for (unsigned int subLevel = 0; subLevel < block1->GetNearNeighSize();
               subLevel++) // for each block, loop over all sublevels in his
                           // NN list (to account for childless blocks)
            {
              // now use intList to compute nonIntList
              // std::set<types::global_dof_index>
              const auto &intList = block1->GetIntList(subLevel);
              // TODO: this is a less than desirable pattern, there's no need
              // to duplicate this set
              const typename OctreeBlock<dim>::small_set *parentIntList =
                nullptr;
              if (subLevel == 0)
                { // if a block is childless we get its intList at the
                  // previous level, otherwise we get its parent's intList
                  parentIntList =
                    &(blocks[block1->GetParentId()]->GetIntList(0));
                }
              else
                {
                  parentIntList = &(block1->GetIntList(subLevel - 1));
                }

              for (auto pos1 = parentIntList->cbegin();
                   pos1 != parentIntList->cend();
                   pos1++) // loop over blocks in parentIntList
                {
                  OctreeBlock<dim> *block2 = blocks[*pos1];
                  if (block2->GetBlockChildrenNum() == 0)
                    { // if blocks in parentIntList are childless, don't look
                      // for their children, but see if they are in
                      // nonIntList
                      if (intList.count(*pos1) == 0)
                        { // if these blocks are not in intList
                          block1->AddBlockToNonIntList(
                            subLevel, *pos1); // then they go in nonIntList
                        }
                    }
                  else
                    { // if blocks in parentIntList are not childless, do
                      // the same test on all their children
                      for (types::global_dof_index kk = 0;
                           kk < block2->GetBlockChildrenNum();
                           kk++) // loop over children of blocks in
                                 // parentIntList
                        {
                          if (intList.count(block2->GetChildId(kk)) == 0)
                            { // if these blocks are not in intList
                              block1->AddBlockToNonIntList(
                                subLevel,
                                block2->GetChildId(
                                  kk)); // then they go in nonIntList
                            }
                        } // end loop over children of blocks in
                          // parentIntList
                    }
                } // loop over blocks in parentIntList
            }     // end loop over subLevels of each block's intList
        }         // end loop over blocks of a level
    }             // end loop over levels

  pcout << "Done computing proximity lists for blocks" << std::endl;
} // end method for octree blocking generation

template class BEMFMA<2>;
template class BEMFMA<3>;
