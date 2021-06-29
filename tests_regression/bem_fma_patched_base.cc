#define TOLL 0.001
// #define MAXELEMENTSPERBLOCK 1

#include "../include/bem_fma.h"

#include "../include/laplace_kernel.h"
#include "Teuchos_TimeMonitor.hpp"

using Teuchos::RCP;
using Teuchos::Time;
using Teuchos::TimeMonitor;
using namespace tbb;

RCP<Time> MatrVec =
  Teuchos::TimeMonitor::getNewTimer("Multipole MatrVec Products Time");
RCP<Time> MultGen =
  Teuchos::TimeMonitor::getNewTimer("Multipole Generation Time");
RCP<Time> MultInt =
  Teuchos::TimeMonitor::getNewTimer("Multipole Integral Time");
RCP<Time> ListCreat =
  Teuchos::TimeMonitor::getNewTimer("Octree Generation Time");
RCP<Time> DirInt = Teuchos::TimeMonitor::getNewTimer("Direct Integral Time");
RCP<Time> PrecondTime =
  Teuchos::TimeMonitor::getNewTimer("FMA_preconditioner Time");
RCP<Time> LocEval = Teuchos::TimeMonitor::getNewTimer("Local Evaluation Time");

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
        delete blocks[ii];
    }
}

template <int dim>
void
BEMFMA<dim>::init_fma(
  const DoFHandler<dim - 1, dim> &                      input_dh,
  const std::vector<std::set<types::global_dof_index>> &db_in,
  const TrilinosWrappers::MPI::Vector &                 input_sn,
  const Mapping<dim - 1, dim> &                         input_mapping,
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
  this_cpu_set.add_indices(
    input_sn
      .locally_owned_elements()); // = new const
                                  // IndexSet(input_sn.locally_owned_elements());
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

  add_parameter(prm,
                &tbb_granularity,
                "Granularity for TBB simple for cycles",
                "10",
                Patterns::Integer());
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
}



template <int dim>
void
BEMFMA<dim>::direct_integrals()
{
  pcout << "Computing direct integrals..." << std::endl;
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
  std::vector<QTelles<dim - 1>> sing_quadratures;
  for (unsigned int i = 0; i < fma_dh->get_fe().dofs_per_cell; ++i)
    sing_quadratures.push_back(
      QTelles<dim - 1>(singular_quadrature_order,
                       fma_dh->get_fe().get_unit_support_points()[i]));
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
  // pcout<<this_cpu_set.size()<<" "<<this_cpu_set.n_elements()<<std::endl;
  types::global_dof_index preconditioner_band =
    125 * fma_dh->get_fe().dofs_per_cell;
  // preconditioner_sparsity_pattern.reinit(sol.vector_partitioner(), (unsigned
  // int) preconditioner_band);
  TrilinosWrappers::MPI::Vector helper(this_cpu_set, mpi_communicator);
  // TODO WHY IT DOES NOT WORK????
  // init_prec_sparsity_pattern.reinit(this_cpu_set.make_trilinos_map(mpi_communicator),preconditioner_band);//,125*fma_dh->get_fe().dofs_per_cell);
  // init_prec_sparsity_pattern.reinit(helper.vector_partitioner(),
  // preconditioner_band);//,125*fma_dh->get_fe().dofs_per_cell);
  init_prec_sparsity_pattern.reinit(
    this_cpu_set,
    mpi_communicator,
    preconditioner_band); //,125*fma_dh->get_fe().dofs_per_cell);

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
    copy_data.block_indices.resize(0);
    copy_data.col_indices.resize(0);
    // for each block in the childless
    // list we get the list of nodes and
    // we check if it contains nodes:
    // if no nodes are contained there is
    // nothing to do

    types::global_dof_index blockId = this->childlessList[kk];

    OctreeBlock<dim> *block1 = this->blocks[blockId];

    std::vector<types::global_dof_index> block1Nodes =
      block1->GetBlockNodeList();

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
        const std::set<types::global_dof_index> &block1IntList =
          block1->GetIntList(intListSubLevs - 1);

        // in this set we will put all the
        // dofs of the cell to whom
        // the quad points belong

        std::set<types::global_dof_index> directNodes;

        // start looping on the intList
        // blocks (block2 here)

        for (std::set<types::global_dof_index>::iterator pos =
               block1IntList.begin();
             pos != block1IntList.end();
             pos++)
          {
            OctreeBlock<dim> *block2 = this->blocks[*pos];
            std::map<cell_it, std::vector<types::global_dof_index>>
              blockQuadPointsList = block2->GetBlockQuadPointsList();

            // get the list of quad points
            // in block2 and loop on it

            typename std::map<cell_it,
                              std::vector<types::global_dof_index>>::iterator
              it;
            for (it = blockQuadPointsList.begin();
                 it != blockQuadPointsList.end();
                 it++)
              {
                // the key of the map (*it.first pointer) is
                // the cell of the quad point: we will
                // get its dofs and put them in the set
                // of direct nodes

                cell_it cell =
                  (*it)
                    .first; // pcout<<cell<<"  end
                            // "<<(*blockQuadPointsList.end()).first<<std::endl;
                cell->get_dof_indices(local_dof_indices);
                for (unsigned int j = 0;
                     j < this->fma_dh->get_fe().dofs_per_cell;
                     j++)
                  directNodes.insert(local_dof_indices[j]);
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
          if (this_cpu_set.is_element(block1Nodes[i]))
            {
              copy_data.block_indices.push_back(block1Nodes[i]);
              copy_data.col_indices.push_back(
                std::vector<types::global_dof_index>());
              for (std::set<types::global_dof_index>::iterator pos =
                     directNodes.begin();
                   pos != directNodes.end();
                   pos++)
                {
                  copy_data.col_indices.back().push_back(*pos);
                  // init_prec_sparsity_pattern.add(block1Nodes[i],*pos);
                }
            }
        // std::cout<<copy_data.col_indices.size()<<"
        // "<<copy_data.block_indices.size()<<std::endl;
      }
  };

  // The copier function uses the InitPrecCopy structure to know the global
  // indices to add to the global initial sparsity pattern. We use once again
  // the capture to access the global memory.
  auto f_init_prec_copier = [this](const InitPrecCopy &copy_data) {
    for (types::global_dof_index i = 0; i < copy_data.col_indices.size(); ++i)
      {
        for (types::global_dof_index j = 0; j < copy_data.col_indices[i].size();
             ++j)
          this->init_prec_sparsity_pattern.add(copy_data.block_indices[i],
                                               copy_data.col_indices[i][j]);
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

  // for (unsigned int kk = 0; kk < childlessList.size(); kk++)
  //   {
  //     // for each block in the childless
  //     // list we get the list of nodes and
  //     // we check if it contains nodes:
  //     // if no nodes are contained there is
  //     // nothing to do
  //
  //     unsigned int blockId = childlessList[kk];
  //
  //     OctreeBlock<dim> *block1 =  blocks[blockId];
  //
  //     std::vector <unsigned int> block1Nodes = block1->GetBlockNodeList();
  //
  //     if  (block1Nodes.size() > 0)
  //       {
  //
  //         // if block1 contains nodes,
  //         // we need to get all the quad points
  //         // in the intList blocks of block1
  //         // (such quad points will be used for
  //         // direct integrals)
  //
  //         unsigned int intListSubLevs = block1->GetIntListSize();
  //         const std::set<unsigned int> &block1IntList =
  //         block1->GetIntList(intListSubLevs-1);
  //
  //         // in this set we will put all the
  //         // dofs of the cell to whom
  //         // the quad points belong
  //
  //         std::set<unsigned int> directNodes;
  //
  //         // start looping on the intList
  //         // blocks (block2 here)
  //
  //         for (std::set<unsigned int>::iterator pos = block1IntList.begin();
  //         pos != block1IntList.end(); pos++)
  //           {
  //             OctreeBlock<dim> *block2 =  blocks[*pos];
  //             std::map <cell_it, std::vector<types::global_dof_index> >
  //             blockQuadPointsList = block2->GetBlockQuadPointsList();
  //
  //             // get the list of quad points
  //             // in block2 and loop on it
  //
  //             typename std::map <cell_it,
  //             std::vector<types::global_dof_index> >::iterator it; for (it =
  //             blockQuadPointsList.begin(); it != blockQuadPointsList.end();
  //             it++)
  //               {
  //                 // the key of the map (*it.first pointer) is
  //                 // the cell of the quad point: we will
  //                 // get its dofs and put them in the set
  //                 // of direct nodes
  //
  //                 cell_it cell = (*it).first;//pcout<<cell<<"  end
  //                 "<<(*blockQuadPointsList.end()).first<<std::endl;
  //                 cell->get_dof_indices(local_dof_indices);
  //                 for (unsigned int j = 0; j < dofs_per_cell; j++)
  //                   directNodes.insert(local_dof_indices[j]);
  //               }
  //           }
  //         // the loop over blocks in intList
  //         // is over: for all the nodes in
  //         // block1, we know nodes in directNodes
  //         // list have direct integrals, so
  //         // we use them to create the
  //         // direct contributions matrices
  //         // sparsity pattern
  //
  //         for (unsigned int i = 0; i < block1Nodes.size(); i++)
  //           if (this_cpu_set.is_element(block1Nodes[i]))
  //             for (std::set<unsigned int>::iterator pos =
  //             directNodes.begin(); pos != directNodes.end(); pos++)
  //               {
  //                 init_prec_sparsity_pattern.add(block1Nodes[i],*pos);
  //               }
  //       }
  //
  //   }

  // unfortunately, the direct integrals must not be computed only for the
  // quadPoints in the intList: if a bigger block is in the nonIntList of
  // another block, the bound for the multipole expansion application does
  // not hold, and so we must compute direct integrals. Here we scan the
  // nonIntlists of each block at each level to look for bigger blocks and
  // initialize the prec matrices sparsity pattern with the corresponding nodes

  for (unsigned int level = 1; level < num_octree_levels + 1;
       level++) // loop over levels

    {
      // std::vector<types::global_dof_index>
      // dofs_filled_blocks =  dofs_filled_blocks[level];
      types::global_dof_index startBlockLevel = startLevel[level];

      // For each level we need again WorkStream to compute the entries in the
      // sparisity pattern that belong to block that are in the nonIntList of
      // the current block but that are of greater size. Once again we use
      // IndexSet to know if we the node in the block belong to the current
      // proccesor.
      auto f_init_prec_level_worker =
        [this, &startBlockLevel, &level](types::global_dof_index jj,
                                         InitPrecScratch &,
                                         InitPrecCopy &copy_data) {
          copy_data.block_indices.resize(0);
          copy_data.col_indices.resize(0);

          OctreeBlock<dim> *block1 =
            this->blocks[this->dofs_filled_blocks[level][jj]];
          const std::vector<types::global_dof_index> &nodesBlk1Ids =
            block1->GetBlockNodeList();
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

                  std::set<types::global_dof_index>        directNodes;
                  const std::set<types::global_dof_index> &nonIntList =
                    block1->GetNonIntList(subLevel);

                  // loop over well separated blocks of higher size (level): in
                  // this case we must use direct evaluation: for each block we
                  // get the quad points list
                  for (std::set<types::global_dof_index>::iterator pos =
                         nonIntList.begin();
                       pos != nonIntList.lower_bound(startBlockLevel);
                       pos++)
                    {
                      OctreeBlock<dim> *block2 = this->blocks[*pos];
                      std::map<cell_it, std::vector<types::global_dof_index>>
                        blockQuadPointsList = block2->GetBlockQuadPointsList();

                      // we loop on the cells of the quad blocks (*it.first
                      // pointer) and put their dofs in the direct list

                      typename std::map<
                        cell_it,
                        std::vector<types::global_dof_index>>::iterator it;
                      for (it = blockQuadPointsList.begin();
                           it != blockQuadPointsList.end();
                           it++)
                        {
                          cell_it cell = (*it).first;
                          cell->get_dof_indices(local_dof_indices);
                          for (unsigned int j = 0;
                               j < this->fma_dh->get_fe().dofs_per_cell;
                               j++)
                            directNodes.insert(local_dof_indices[j]);
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
                          copy_data.col_indices.push_back(
                            std::vector<types::global_dof_index>());
                          for (std::set<types::global_dof_index>::iterator pos =
                                 directNodes.begin();
                               pos != directNodes.end();
                               pos++)
                            copy_data.col_indices.back().push_back(*pos);
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

      // for (unsigned int jj = 0; jj < dofs_filled_blocks[level].size();  jj++)
      //   {
      //     OctreeBlock<dim> *block1 =  blocks[dofs_filled_blocks[level][jj]];
      //     const std::vector <unsigned int> &nodesBlk1Ids =
      //     block1->GetBlockNodeList();
      //
      //     // again, no need to perform next operations if block has no nodes
      //
      //     if  (nodesBlk1Ids.size() > 0)// !!!CHECK, IT SEEMS TO BE USELESS
      //       {
      //         // for each block containing nodes, loop over all sublevels in
      //         his NN list (this is because if a
      //         // block remains childless BEFORE the last level, at this point
      //         we need to compute
      //         // all its contributions up to the bottom level)
      //
      //
      //         for (unsigned int subLevel = 0; subLevel <
      //         block1->NumNearNeighLevels();  subLevel++)
      //           {
      //
      //             // in this vectors we are saving the nodes needing direct
      //             integrals
      //
      //             std::set <unsigned int> directNodes;
      //             const std::set <unsigned int> &nonIntList =
      //             block1->GetNonIntList(subLevel);
      //
      //             // loop over well separated blocks of higher size (level):
      //             in this case
      //             // we must use direct evaluation: for each block we get the
      //             quad points
      //             // list
      //             for (std::set<unsigned int>::iterator pos =
      //             nonIntList.begin(); pos
      //             !=nonIntList.lower_bound(startBlockLevel); pos++)
      //               {
      //                 OctreeBlock<dim> *block2 =  blocks[*pos];
      //                 std::map <cell_it, std::vector<types::global_dof_index>
      //                 > blockQuadPointsList =
      //                 block2->GetBlockQuadPointsList();
      //
      //                 // we loop on the cells of the quad blocks (*it.first
      //                 pointer)
      //                 // and put their dofs in the direct list
      //
      //                 typename std::map <cell_it,
      //                 std::vector<types::global_dof_index> >::iterator it;
      //                 for (it = blockQuadPointsList.begin(); it !=
      //                 blockQuadPointsList.end(); it++)
      //                   {
      //                     cell_it cell = (*it).first;
      //                     cell->get_dof_indices(local_dof_indices);
      //                     for (unsigned int j = 0; j < dofs_per_cell; j++)
      //                       directNodes.insert(local_dof_indices[j]);
      //                   }
      //               } // end loop over blocks of a sublevel of nonIntList
      //
      //             // we use the nodes in directList, to create the sparsity
      //             pattern
      //
      //             for (unsigned int i = 0; i < nodesBlk1Ids.size(); i++)
      //               {
      //                 if (this_cpu_set.is_element(nodesBlk1Ids[i]))
      //                   for (std::set<unsigned int>::iterator pos =
      //                   directNodes.begin(); pos != directNodes.end(); pos++)
      //                     init_prec_sparsity_pattern.add(nodesBlk1Ids[i],*pos);
      //               }
      //
      //           } // end loop over sublevels
      //       } // end if: is there any node in the block?
      //   }// end loop over block of a level
    } // end loop over octree levels


  // sparsity pattern is ready and can be compressed; the direct matrices
  // and the preconditioner one are the initialized with the sparsity
  // pattern just computed

  init_prec_sparsity_pattern.compress();
  double filling_percentage =
    double(init_prec_sparsity_pattern.n_nonzero_elements()) /
    double(fma_dh->n_dofs() * fma_dh->n_dofs()) * 100.;
  pcout << init_prec_sparsity_pattern.n_nonzero_elements()
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
           DirectCopyData &               copy_data,
           const std::vector<Point<dim>> &support_points,
           std::vector<QTelles<dim - 1>> &sing_quadratures) {
      // pcout<<"processing block "<<kk <<"  of
      // "<<cMesh->GetNumChildlessBlocks()<<std::endl; pcout<<"block
      // "<<cMesh->GetChildlessBlockId(kk) <<"  of  "<<cMesh->GetNumBlocks()<<"
      // in block list"<<std::endl;

      // this is the Id of the block
      copy_data.vec_local_dof_indices.resize(0);
      copy_data.vec_local_neumann_matrix_row_i.resize(0);
      copy_data.vec_local_dirichlet_matrix_row_i.resize(0);
      copy_data.vec_node_index.resize(0);
      copy_data.vec_start_helper.resize(0);
      types::global_dof_index blockId = *block_it;
      // unsigned int blockId =  this->childlessList[block_it];
      // and this is the block pointer
      OctreeBlock<dim> *block1 = this->blocks[blockId];
      // we get the block node list
      const std::vector<types::global_dof_index> &block1Nodes =
        block1->GetBlockNodeList();

      // if a block has no nodes (if it only contains quad points), there is
      // nothing to do if instead there are nodes, we start integrating
      if (block1Nodes.size() > 0)
        {
          // std::cout<<"Nodes in childless block :
          // "<<block1Nodes.size()<<std::endl; we first get all the blocks in
          // the intList of the current block (block1) and loop over these
          // blocks, to create a list of ALL the quadrature points that lie in
          // the interaction list blocks: these quad points have to be
          // integrated directly. the list of direct quad points has to be a
          // std::map of std::set of integers, meaning that to each cell, we
          // associate a std::set containing all the direct quad point ids
          unsigned int intListNumLevs = block1->GetIntListSize();
          std::set<types::global_dof_index> block1IntList =
            block1->GetIntList(intListNumLevs - 1);

          std::map<cell_it, std::set<types::global_dof_index>> directQuadPoints;
          for (std::set<types::global_dof_index>::iterator pos =
                 block1IntList.begin();
               pos != block1IntList.end();
               pos++)
            {
              // now for each block block2 we get the list of quad points
              OctreeBlock<dim> *block2 = this->blocks[*pos];
              std::map<cell_it, std::vector<types::global_dof_index>>
                blockQuadPointsList = block2->GetBlockQuadPointsList();

              // we now loop among the cells of the list and for each cell we
              // loop among its quad points, to copy them into the direct quad
              // points list
              typename std::map<cell_it,
                                std::vector<types::global_dof_index>>::iterator
                it;
              for (it = blockQuadPointsList.begin();
                   it != blockQuadPointsList.end();
                   it++)
                {
                  for (types::global_dof_index i = 0; i < (*it).second.size();
                       i++)
                    {
                      directQuadPoints[(*it).first].insert((*it).second[i]);

                      /*//////////this is for a check///////////////////
                      for (types::global_dof_index kk=0; kk<block1Nodes.size();
                      kk++) integralCheck[block1Nodes[kk]][(*it).first] += 1;
                      ///////////////////////////*/
                    }
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

                  typename std::map<cell_it,
                                    std::set<types::global_dof_index>>::iterator
                    it;
                  // we loop on the list of quad points to be treated directly
                  for (it = directQuadPoints.begin();
                       it != directQuadPoints.end();
                       it++)
                    {
                      // the vectors with the local integrals for the cell must
                      // first be zeroed
                      copy_data.vec_local_neumann_matrix_row_i.push_back(
                        Vector<double>(this->fma_dh->get_fe().dofs_per_cell));
                      copy_data.vec_local_dirichlet_matrix_row_i.push_back(
                        Vector<double>(this->fma_dh->get_fe().dofs_per_cell));

                      // we get the first entry of the map, i.e. the cell
                      // pointer and we check if the cell contains the current
                      // node, to decide if singular of regular quadrature is to
                      // be used
                      cell_it cell = (*it).first;
                      copy_data.vec_local_dof_indices.push_back(
                        std::vector<types::global_dof_index>(
                          this->fma_dh->get_fe().dofs_per_cell));
                      cell->get_dof_indices(
                        copy_data.vec_local_dof_indices.back());

                      // we copy the cell quad points in this set
                      std::set<types::global_dof_index> &cellQuadPoints =
                        (*it).second;
                      bool         is_singular = false;
                      unsigned int singular_index =
                        numbers::invalid_unsigned_int;

                      for (unsigned int j = 0;
                           j < this->fma_dh->get_fe().dofs_per_cell;
                           ++j)
                        if ((*(this->double_nodes_set))[nodeIndex].count(
                              copy_data.vec_local_dof_indices.back()[j]) > 0)
                          {
                            singular_index = j;
                            is_singular    = true;
                            break;
                          }
                      // first case: the current node does not belong to the
                      // current cell: we use regular quadrature
                      if (is_singular == false)
                        {
                          // pcout<<"Node "<<i<<"  Elem "<<cell<<" (Direct)
                          // Nodes: "; for(unsigned int j=0; j<fe.dofs_per_cell;
                          // ++j) pcout<<" "<<local_dof_indices[j];
                          // pcout<<std::endl;

                          // we start looping on the quad points of the cell:
                          // *pos will be the index of the quad point
                          for (std::set<types::global_dof_index>::iterator pos =
                                 cellQuadPoints.begin();
                               pos != cellQuadPoints.end();
                               pos++)
                            {
                              // here we compute the distance R between the node
                              // and the quad point

                              // MAGARI USARE FEVALUES CON IL DOFHANDLER CRETINO
                              // DISCONTINUO E IL MAPPING bem_fma
                              Point<dim> D;
                              double     s = 0.;

                              const Tensor<1, dim> R =
                                this->quadPoints.at(cell)[*pos] -
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
                                  copy_data.vec_local_neumann_matrix_row_i
                                    .back()(j) +=
                                    ((D * this->quadNormals.at(cell)[*pos]) *
                                     this->quadShapeFunValues.at(
                                       cell)[*pos][j] *
                                     this->quadJxW.at(cell)[*pos]);
                                  copy_data.vec_local_dirichlet_matrix_row_i
                                    .back()(j) +=
                                    (s *
                                     this->quadShapeFunValues.at(
                                       cell)[*pos][j] *
                                     this->quadJxW.at(cell)[*pos]);
                                  // if(std::abs(copy_data.vec_local_neumann_matrix_row_i.back()(j))<1e-12)
                                  //   std::cout<<D<<" "<<
                                  //   this->quadNormals.at(cell)[*pos]<<"
                                  //   "<<copy_data.vec_local_neumann_matrix_row_i.back()(j)<<std::endl;
                                  // pcout<<
                                  // quadShapeFunValues[cell][*pos][j]<<" ";
                                  // pcout<< quadJxW[cell][*pos]<<std::endl;
                                  // std::cout<<D<<std::endl<<"
                                  // "<<this->quadNormals.at(cell)[*pos]<<std::endl;
                                  // std::cout<<copy_data.vec_local_neumann_matrix_row_i.back()(j)<<"
                                  // "<<copy_data.vec_local_dirichlet_matrix_row_i.back()(j)<<std::endl;
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

                          const Quadrature<dim - 1> *singular_quadrature =
                            (dim == 2 ?
                               dynamic_cast<Quadrature<dim - 1> *>(
                                 &sing_quadratures[singular_index]) :
                               (dim == 3 ?
                                  dynamic_cast<Quadrature<dim - 1> *>(
                                    &sing_quadratures[singular_index]) :
                                  0));
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
                          if (dim == 2)
                            delete singular_quadrature;

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
        // std::cout<<"Sizes of vec_node_index and vec_start_helper:
        // "<<copy_data.vec_node_index.size()<<"
        // "<<copy_data.vec_start_helper.size()<<std::endl;
        // for(auto fdj : copy_data.vec_start_helper)
        //   std::cout<<fdj<<" ";


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
              foo_end = copy_data.vec_start_helper[ii + 1];

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
                        //  std::cout<<"DIANE"<<std::endl;
                        this->init_preconditioner.add(
                          copy_data.vec_node_index[ii],
                          copy_data.vec_local_dof_indices[kk][j],
                          -copy_data.vec_local_dirichlet_matrix_row_i[kk](j));
                      }
                    else
                      this->init_preconditioner.add(
                        copy_data.vec_node_index[ii],
                        copy_data.vec_local_dof_indices[kk][j],
                        copy_data.vec_local_neumann_matrix_row_i[kk](j));
                    //  std::cout<<this->init_preconditioner(copy_data.vec_node_index[ii],copy_data.vec_local_dof_indices[kk][j])<<"
                    //  ";
                    // std::cout<<copy_data.vec_local_dirichlet_matrix_row_i[kk][j]<<"
                    // "<<std::endl;
                  }


                // std::cout<<std::endl;
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
           DirectCopyData &               copy_data,
           const std::vector<Point<dim>> &support_points,
           types::global_dof_index        startBlockLevel) {
      copy_data.vec_local_dof_indices.resize(0);
      copy_data.vec_local_neumann_matrix_row_i.resize(0);
      copy_data.vec_local_dirichlet_matrix_row_i.resize(0);
      copy_data.vec_node_index.resize(0);
      copy_data.vec_start_helper.resize(0);

      types::global_dof_index blockId = *block_it;
      OctreeBlock<dim> *      block1  = this->blocks[blockId];
      const std::vector<types::global_dof_index> &nodesBlk1Ids =
        block1->GetBlockNodeList();
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

          if (this->this_cpu_set.is_element(
                nodeIndex)) //(m2l_flags[level][jj]==this_mpi_process)
            {
              std::map<cell_it, std::set<types::global_dof_index>>
                directQuadPoints;

              for (unsigned int subLevel = 0;
                   subLevel < block1->NumNearNeighLevels();
                   subLevel++)
                {
                  const std::set<types::global_dof_index> &nonIntList =
                    block1->GetNonIntList(subLevel);

                  // loop over well separated blocks of higher size
                  // (level)-----> in this case
                  // we must use direct evaluation (luckily being childless they
                  // only contain 1 element)
                  for (std::set<types::global_dof_index>::iterator pos =
                         nonIntList.begin();
                       pos != nonIntList.lower_bound(startBlockLevel);
                       pos++)
                    {
                      OctreeBlock<dim> *block2 = this->blocks[*pos];
                      std::map<cell_it, std::vector<types::global_dof_index>>
                        blockQuadPointsList = block2->GetBlockQuadPointsList();
                      typename std::map<
                        cell_it,
                        std::vector<types::global_dof_index>>::iterator it;
                      for (it = blockQuadPointsList.begin();
                           it != blockQuadPointsList.end();
                           it++)
                        {
                          for (types::global_dof_index ii = 0;
                               ii < (*it).second.size();
                               ii++)
                            {
                              directQuadPoints[(*it).first].insert(
                                (*it).second[ii]);

                              /*////////this is for a check/////////////////////
                                         integralCheck[nodesBlk1Ids[i]][(*it).first]
                                 += 1;
                                          ////////////////////////////*/
                            }
                        }
                    } // end loop over blocks of a sublevel of nonIntList
                }     // end loop over sublevels

              typename std::map<cell_it,
                                std::set<types::global_dof_index>>::iterator it;
              for (it = directQuadPoints.begin(); it != directQuadPoints.end();
                   it++)
                {
                  // the vectors with the local integrals for the cell must
                  // first be zeroed
                  copy_data.vec_local_neumann_matrix_row_i.push_back(
                    Vector<double>(this->fma_dh->get_fe().dofs_per_cell));
                  copy_data.vec_local_dirichlet_matrix_row_i.push_back(
                    Vector<double>(this->fma_dh->get_fe().dofs_per_cell));

                  // we get the first entry of the map, i.e. the cell pointer
                  // here the quadrature is regular as the cell is well
                  // separated
                  cell_it cell = (*it).first;
                  copy_data.vec_local_dof_indices.push_back(
                    std::vector<types::global_dof_index>(
                      this->fma_dh->get_fe().dofs_per_cell));
                  cell->get_dof_indices(copy_data.vec_local_dof_indices.back());

                  // we copy the cell quad points in this set
                  std::set<types::global_dof_index> &cellQuadPoints =
                    (*it).second;

                  // pcout<<"Node "<<i<<"  Elem "<<cell<<" (Direct) Nodes: ";
                  // for(unsigned int j=0; j<fe.dofs_per_cell; ++j) pcout<<"
                  // "<<local_dof_indices[j]; pcout<<std::endl;

                  // we start looping on the quad points of the cell: *pos will
                  // be the index of the quad point
                  for (std::set<types::global_dof_index>::iterator pos =
                         cellQuadPoints.begin();
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

  // here we finally start computing the direct integrals: we
  // first loop among the childless blocks

  // for (unsigned int kk = 0; kk <  childlessList.size(); kk++)
  //
  //   {
  //     //pcout<<"processing block "<<kk <<"  of
  //     "<<cMesh->GetNumChildlessBlocks()<<std::endl;
  //     //pcout<<"block "<<cMesh->GetChildlessBlockId(kk) <<"  of
  //     "<<cMesh->GetNumBlocks()<<"  in block list"<<std::endl;
  //
  //     // this is the Id of the block
  //     unsigned int blockId =  childlessList[kk];
  //     // and this is the block pointer
  //     OctreeBlock<dim> *block1 =  blocks[blockId];
  //     // we get the block node list
  //     const std::vector <unsigned int> &block1Nodes =
  //     block1->GetBlockNodeList();
  //
  //     // if a block has no nodes (if it only contains quad points), there is
  //     nothing to do
  //     // if instead there are nodes, we start integrating
  //     if  (block1Nodes.size() > 0)
  //       {
  //         // we first get all the blocks in the intList of the current block
  //         (block1)
  //         // and loop over these blocks, to create a list of ALL the
  //         quadrature points that
  //         // lie in the interaction list blocks: these quad points have to be
  //         integrated
  //         // directly. the list of direct quad points has to be a std::map of
  //         std::set of
  //         // integers, meaning that to each cell, we associate a std::set
  //         containing all
  //         // the direct quad point ids
  //         unsigned int intListNumLevs = block1->GetIntListSize();
  //         std::set <unsigned int> block1IntList =
  //         block1->GetIntList(intListNumLevs-1);
  //
  //         std::map <cell_it,std::set<unsigned int> > directQuadPoints;
  //         for (std::set<unsigned int>::iterator pos = block1IntList.begin();
  //         pos != block1IntList.end(); pos++)
  //           {
  //             // now for each block block2 we get the list of quad points
  //             OctreeBlock<dim> *block2 =  blocks[*pos];
  //             std::map <cell_it, std::vector<types::global_dof_index> >
  //             blockQuadPointsList = block2->GetBlockQuadPointsList();
  //
  //             // we now loop among the cells of the list and for each cell we
  //             loop
  //             // among its quad points, to copy them into the direct quad
  //             points list typename std::map <cell_it,
  //             std::vector<types::global_dof_index> >::iterator it; for (it =
  //             blockQuadPointsList.begin(); it != blockQuadPointsList.end();
  //             it++)
  //               {
  //                 for (unsigned int i=0; i<(*it).second.size(); i++)
  //                   {
  //                     directQuadPoints[(*it).first].insert((*it).second[i]);
  //
  //                     /*//////////this is for a check///////////////////
  //                     for (unsigned int kk=0; kk<block1Nodes.size(); kk++)
  //                            integralCheck[block1Nodes[kk]][(*it).first] +=
  //                            1;
  //                     ///////////////////////////*/
  //                   }
  //               }
  //           }
  //         // we are now ready to go: for each node, we know which quad points
  //         are to be
  //         // treated directly, and for each node, we will now perform the
  //         integral.
  //         // we then start looping on the nodes of the block
  //         for (unsigned int i=0; i<block1Nodes.size(); i++)
  //           {
  //             unsigned int nodeIndex = block1Nodes[i];
  //             if(this_cpu_set.is_element(nodeIndex))
  //             {
  //               typename std::map <cell_it, std::set<unsigned int>
  //               >::iterator it;
  //               // we loop on the list of quad points to be treated directly
  //               for (it = directQuadPoints.begin(); it !=
  //               directQuadPoints.end(); it++)
  //                 {
  //                   // the vectors with the local integrals for the cell must
  //                   first
  //                   // be zeroed
  //                   local_neumann_matrix_row_i = 0;
  //                   local_dirichlet_matrix_row_i = 0;
  //
  //                   // we get the first entry of the map, i.e. the cell
  //                   pointer
  //                   // and we check if the cell contains the current node, to
  //                   // decide if singular of regular quadrature is to be used
  //                   cell_it cell = (*it).first;
  //                   cell->get_dof_indices(local_dof_indices);
  //
  //                   // we copy the cell quad points in this set
  //                   std::set<unsigned int> &cellQuadPoints = (*it).second;
  //                   bool is_singular = false;
  //                   unsigned int singular_index =
  //                   numbers::invalid_unsigned_int;
  //
  //                   for (unsigned int j=0; j<fma_dh->get_fe().dofs_per_cell;
  //                   ++j)
  //                     if (
  //                     (*double_nodes_set)[nodeIndex].count(local_dof_indices[j])
  //                     > 0)
  //                       {
  //                         singular_index = j;
  //                         is_singular = true;
  //                         break;
  //                       }
  //                   // first case: the current node does not belong to the
  //                   current cell:
  //                   // we use regular quadrature
  //                   if (is_singular == false)
  //                     {
  //                       //pcout<<"Node "<<i<<"  Elem "<<cell<<" (Direct)
  //                       Nodes: ";
  //                       //for(unsigned int j=0; j<fe.dofs_per_cell; ++j)
  //                       pcout<<" "<<local_dof_indices[j];
  //                       //pcout<<std::endl;
  //
  //                       // we start looping on the quad points of the cell:
  //                       *pos will be the
  //                       // index of the quad point
  //                       for (std::set<unsigned int>::iterator
  //                       pos=cellQuadPoints.begin();
  //                       pos!=cellQuadPoints.end(); pos++)
  //                         {
  //                           // here we compute the distance R between the
  //                           node and the quad point
  //
  //                           //MAGARI USARE FEVALUES CON IL DOFHANDLER CRETINO
  //                           DISCONTINUO E IL MAPPING bem_fma const Tensor<1,
  //                           dim> R =  quadPoints[cell][*pos] -
  //                           support_points[nodeIndex];
  //                           LaplaceKernel::kernels(R, D, s);
  //
  //                           // and here are the integrals for each of the
  //                           degrees of freedom of the cell: note
  //                           // how the quadrature values (position, normals,
  //                           jacobianXweight, shape functions)
  //                           // are taken from the precomputed ones in
  //                           ComputationalDomain class for (unsigned int j=0;
  //                           j<fma_dh->get_fe().dofs_per_cell; ++j)
  //                             {
  //                               local_neumann_matrix_row_i(j) += ( ( D *
  //                                                                    quadNormals[cell][*pos]
  //                                                                    ) *
  //                                                                  quadShapeFunValues[cell][*pos][j]
  //                                                                  *
  //                                                                  quadJxW[cell][*pos]
  //                                                                  );
  //                               local_dirichlet_matrix_row_i(j) += ( s *
  //                                                                    quadShapeFunValues[cell][*pos][j]
  //                                                                    *
  //                                                                    quadJxW[cell][*pos]
  //                                                                    );
  //                               //pcout<<D<<" "<< quadNormals[cell][*pos]<<"
  //                               ";
  //                               //pcout<<
  //                               quadShapeFunValues[cell][*pos][j]<<" ";
  //                               //pcout<< quadJxW[cell][*pos]<<std::endl;
  //                             }
  //                         }
  //                     } // end if
  //                   else
  //                     {
  //                       // after some checks, we have to create the singular
  //                       quadrature:
  //                       // here the quadrature points of the cell will be
  //                       IGNORED,
  //                       // and the singular quadrature points are instead
  //                       used.
  //                       // the 3d and 2d quadrature rules are different
  //
  //                       // QUESTO E' IL SOLITO STEP 34, VEDI SE CAMBIARE CON
  //                       QUELLO NUOVO PER STOKES Assert(singular_index !=
  //                       numbers::invalid_unsigned_int,
  //                              ExcInternalError());
  //
  //                       const Quadrature<dim-1> *
  //                       singular_quadrature
  //                         = (dim == 2
  //                            ?
  //                            dynamic_cast<Quadrature<dim-1>*>(
  //                              &sing_quadratures[singular_index])
  //                            :
  //                            (dim == 3
  //                             ?
  //                             dynamic_cast<Quadrature<dim-1>*>(
  //                               &sing_quadratures[singular_index])
  //                             :
  //                             0));
  //                       Assert(singular_quadrature, ExcInternalError());
  //
  //                       // once the singular quadrature has been created, we
  //                       employ it
  //                       // to create the corresponding fe_values
  //
  //                       FEValues<dim-1,dim> fe_v_singular (*fma_mapping,
  //                       fma_dh->get_fe(), *singular_quadrature,
  //                                                          update_jacobians |
  //                                                          update_values |
  //                                                          update_cell_normal_vectors
  //                                                          |
  //                                                          update_quadrature_points
  //                                                          );
  //
  //                       fe_v_singular.reinit(cell);
  //
  //                       // here are the vectors of the quad points and
  //                       normals vectors
  //
  //                       const std::vector<Point<dim> > &singular_normals =
  //                       fe_v_singular.get_normal_vectors(); const
  //                       std::vector<Point<dim> > &singular_q_points =
  //                       fe_v_singular.get_quadrature_points();
  //
  //
  //                       // and here is the integrals computation: note how in
  //                       this case the
  //                       // values for shape functions & co. are not taken
  //                       from the precomputed
  //                       // ones in ComputationalDomain class
  //
  //                       for (unsigned int q=0; q<singular_quadrature->size();
  //                       ++q)
  //                         {
  //                           const Tensor<1, dim> R = singular_q_points[q] -
  //                           support_points[nodeIndex];
  //                           LaplaceKernel::kernels(R, D, s);
  //                           for (unsigned int j=0;
  //                           j<fma_dh->get_fe().dofs_per_cell; ++j)
  //                             {
  //                               local_neumann_matrix_row_i(j) += (( D *
  //                                                                   singular_normals[q])
  //                                                                   *
  //                                                                 fe_v_singular.shape_value(j,q)
  //                                                                 *
  //                                                                 fe_v_singular.JxW(q)
  //                                                                 );
  //
  //                               local_dirichlet_matrix_row_i(j) += ( s   *
  //                                                                    fe_v_singular.shape_value(j,q)
  //                                                                    *
  //                                                                    fe_v_singular.JxW(q)
  //                                                                    );
  //                             }
  //                         }
  //                       if (dim==2)
  //                         delete singular_quadrature;
  //
  //                     } // end else
  //
  //                   // Finally, we need to add
  //                   // the contributions of the
  //                   // current cell to the
  //                   // global matrix.
  //
  //                   for (unsigned int j=0; j<fma_dh->get_fe().dofs_per_cell;
  //                   ++j)
  //                     {
  //                         prec_neumann_matrix.add(nodeIndex,local_dof_indices[j],local_neumann_matrix_row_i(j));
  //                         prec_dirichlet_matrix.add(nodeIndex,local_dof_indices[j],local_dirichlet_matrix_row_i(j));
  //                         if ((*dirichlet_nodes)(local_dof_indices[j]) > 0.8)
  //                           init_preconditioner.add(nodeIndex,local_dof_indices[j],-local_dirichlet_matrix_row_i(j));
  //                         else
  //                           init_preconditioner.add(nodeIndex,local_dof_indices[j],
  //                           local_neumann_matrix_row_i(j));
  //                     }
  //
  //                 } // end loop on cells of the intList
  //             }
  //           } // end loop over nodes of block1
  //       } // end if (nodes in block > 0)
  //   } // end loop over childless blocks


  // as said, the direct integrals must not be computed only for the
  // quadPoints in the intList: if a bigger block is in the nonIntList of
  // another block, the bound for the multipole expansion application does
  // not hold, and so we must compute direct integrals. Here we scan the
  // nonIntlists of each block al each level to look for bigger blocks and
  // compute the direct integral contribution for the quadNodes in such
  // blocks

  // for (unsigned int level = 1; level <  num_octree_levels + 1;  level++) //
  // loop over levels
  //
  //   {
  //     unsigned int startBlockLevel =  startLevel[level];
  //     // !!! Io spezzerei qui per poi comunicare alla fine (se vogliamo, ma
  //     questo viene chiamato poche volte). for (unsigned int jj = 0; jj <
  //     dofs_filled_blocks[level].size();  jj++) // loop over blocks of each
  //     level
  //       {
  //         OctreeBlock<dim> *block1 =  blocks[ dofs_filled_blocks[level][jj]];
  //         const std::vector <unsigned int> &nodesBlk1Ids =
  //         block1->GetBlockNodeList();
  //
  //         for (unsigned int i = 0; i < nodesBlk1Ids.size(); i++)
  //           {
  //             // for each block containing nodes, loop over all sublevels in
  //             his NN list (this is because if a
  //             // block remains childless BEFORE the last level, at this point
  //             we need to compute
  //             // all its contributions up to the bottom level)
  //             unsigned int nodeIndex = nodesBlk1Ids[i];
  //             if(this_cpu_set.is_element(nodeIndex))//(m2l_flags[level][jj]==this_mpi_process)
  //             {
  //
  //             std::map <cell_it,std::set<unsigned int> > directQuadPoints;
  //
  //             for (unsigned int subLevel = 0; subLevel <
  //             block1->NumNearNeighLevels();  subLevel++)
  //               {
  //                 const std::set <unsigned int> &nonIntList =
  //                 block1->GetNonIntList(subLevel);
  //
  //                 // loop over well separated blocks of higher size
  //                 (level)-----> in this case
  //                 //we must use direct evaluation (luckily being childless
  //                 they only contain 1 element) for (std::set<unsigned
  //                 int>::iterator pos = nonIntList.begin(); pos
  //                 !=nonIntList.lower_bound(startBlockLevel); pos++)
  //                   {
  //                     OctreeBlock<dim> *block2 =  blocks[*pos];
  //                     std::map <cell_it, std::vector<types::global_dof_index>
  //                     > blockQuadPointsList =
  //                     block2->GetBlockQuadPointsList(); typename std::map
  //                     <cell_it, std::vector<types::global_dof_index>
  //                     >::iterator it; for (it = blockQuadPointsList.begin();
  //                     it != blockQuadPointsList.end(); it++)
  //                       {
  //                         for (unsigned int ii=0; ii<(*it).second.size();
  //                         ii++)
  //                           {
  //                             directQuadPoints[(*it).first].insert((*it).second[ii]);
  //
  //                             /*////////this is for a
  //                             check/////////////////////
  //                                        integralCheck[nodesBlk1Ids[i]][(*it).first]
  //                                        += 1;
  //                                         ////////////////////////////*/
  //                           }
  //                       }
  //                   } // end loop over blocks of a sublevel of nonIntList
  //               } // end loop over sublevels
  //
  //             typename std::map <cell_it, std::set<unsigned int> >::iterator
  //             it; for (it = directQuadPoints.begin(); it !=
  //             directQuadPoints.end(); it++)
  //               {
  //                 // the vectors with the local integrals for the cell must
  //                 first
  //                 // be zeroed
  //                 local_neumann_matrix_row_i = 0;
  //                 local_dirichlet_matrix_row_i = 0;
  //
  //                 // we get the first entry of the map, i.e. the cell pointer
  //                 // here the quadrature is regular as the cell is well
  //                 // separated
  //                 cell_it cell = (*it).first;
  //                 cell->get_dof_indices(local_dof_indices);
  //                 // we copy the cell quad points in this set
  //                 std::set<unsigned int> &cellQuadPoints = (*it).second;
  //
  //                 //pcout<<"Node "<<i<<"  Elem "<<cell<<" (Direct) Nodes: ";
  //                 //for(unsigned int j=0; j<fe.dofs_per_cell; ++j) pcout<<"
  //                 "<<local_dof_indices[j];
  //                 //pcout<<std::endl;
  //
  //                 // we start looping on the quad points of the cell: *pos
  //                 will be the
  //                 // index of the quad point
  //                 for (std::set<unsigned int>::iterator
  //                 pos=cellQuadPoints.begin(); pos!=cellQuadPoints.end();
  //                 pos++)
  //                   {
  //                     // here we compute the distance R between the node and
  //                     the quad point const Tensor<1,dim> R =
  //                     quadPoints[cell][*pos] - support_points[nodeIndex];
  //                     LaplaceKernel::kernels(R, D, s);
  //
  //                     // and here are the integrals for each of the degrees
  //                     of freedom of the cell: note
  //                     // how the quadrature values (position, normals,
  //                     jacobianXweight, shape functions)
  //                     // are taken from the precomputed ones in
  //                     ComputationalDomain class
  //
  //                     for (unsigned int j=0;
  //                     j<fma_dh->get_fe().dofs_per_cell; ++j)
  //                       {
  //                         local_neumann_matrix_row_i(j) += ( ( D *
  //                                                              quadNormals[cell][*pos]
  //                                                              ) *
  //                                                            quadShapeFunValues[cell][*pos][j]
  //                                                            *
  //                                                            quadJxW[cell][*pos]
  //                                                            );
  //                         local_dirichlet_matrix_row_i(j) += ( s *
  //                                                              quadShapeFunValues[cell][*pos][j]
  //                                                              *
  //                                                              quadJxW[cell][*pos]
  //                                                              );
  //
  //                       } // end loop over the dofs in the cell
  //                   } // end loop over the quad points in a cell
  //
  //                 // Finally, we need to add
  //                 // the contributions of the
  //                 // current cell to the
  //                 // global matrix.
  //
  //                 for (unsigned int j=0; j<fma_dh->get_fe().dofs_per_cell;
  //                 ++j)
  //                   {
  //                     // if(this_cpu_set.is_element(local_dof_indices[j]))
  //                     // {
  //                       prec_neumann_matrix.add(nodeIndex,local_dof_indices[j],local_neumann_matrix_row_i(j));
  //                       prec_dirichlet_matrix.add(nodeIndex,local_dof_indices[j],local_dirichlet_matrix_row_i(j));
  //
  //                       if ((*dirichlet_nodes)(local_dof_indices[j]) > 0.8)
  //                         init_preconditioner.add(nodeIndex,local_dof_indices[j],-local_dirichlet_matrix_row_i(j));
  //                       else
  //                         init_preconditioner.add(nodeIndex,local_dof_indices[j],
  //                         local_neumann_matrix_row_i(j));
  //                     // }
  //                   }
  //
  //
  //               } // end loop over quad points in the direct quad points list
  //             }// end check on proc
  //
  //           } // end loop over nodes in a block
  //       }// end loop over block of a level
  //   }//end loop over octree levels



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
    OctreeBlock<dim> *      block   = this->blocks[blockId];
    double                  delta   = block->GetDelta();
    Point<dim>              deltaHalf;
    for (unsigned int i = 0; i < dim; i++)
      deltaHalf(i) = delta / 2.;
    Point<dim> blockCenter = block->GetPMin() + deltaHalf;

    // at this point, we get the list of quad nodes for the current block,
    // and loop over it
    std::map<cell_it, std::vector<types::global_dof_index>>
      blockQuadPointsList = block->GetBlockQuadPointsList();

    typename std::map<cell_it, std::vector<types::global_dof_index>>::iterator
      it;
    for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end();
         it++)
      {
        // for each cell in the list, we get the list of its quad nodes
        // present in the current block
        cell_it                               cell           = (*it).first;
        std::vector<types::global_dof_index> &cellQuadPoints = (*it).second;

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
                copy_data.myelemMultipoleExpansionsKer1[blockId][cell][j]
                  .AddNormDer(this->quadShapeFunValues.at(cell)[q][j] *
                                this->quadJxW.at(cell)[q] / 4 / numbers::PI,
                              this->quadPoints.at(cell)[q],
                              this->quadNormals.at(cell)[q]);
                copy_data.myelemMultipoleExpansionsKer2[blockId][cell][j].Add(
                  this->quadShapeFunValues.at(cell)[q][j] *
                    this->quadJxW.at(cell)[q] / 4 / numbers::PI,
                  this->quadPoints.at(cell)[q]);
              }
          } // end loop on cell quadrature points in the block
      }
  };

  // The copier function copies the local structures in the global memory.
  auto f_copier_multipole_integral = [this](const MultipoleData &copy_data) {
    this->elemMultipoleExpansionsKer1.insert(
      copy_data.myelemMultipoleExpansionsKer1.begin(),
      copy_data.myelemMultipoleExpansionsKer1.end());
    this->elemMultipoleExpansionsKer2.insert(
      copy_data.myelemMultipoleExpansionsKer2.begin(),
      copy_data.myelemMultipoleExpansionsKer2.end());
  };

  MultipoleScratch foo_scratch;
  MultipoleData    copy_data;

  WorkStream::run(childlessList.begin(),
                  childlessList.end(),
                  f_worker_multipole_integral,
                  f_copier_multipole_integral,
                  foo_scratch,
                  copy_data);



  // for (unsigned int kk = 0; kk <  childlessList.size(); kk++)
  //   f_generate_multipole_integral(kk, dofs_per_cell, this);
  // for (unsigned int kk = 0; kk <  childlessList.size(); kk++)
  //
  //   {
  //     //pcout<<"processing block "<<kk <<"  of
  //     "<<cMesh->GetNumChildlessBlocks()<<std::endl;
  //     //pcout<<"block "<<cMesh->GetChildlessBlockId(kk) <<"  of
  //     "<<cMesh->GetNumBlocks()<<"  in block list"<<std::endl;
  //
  //     // we get the current block and its Id, and then we
  //     // compute its center, which is needed to construct the
  //     // multipole expansion in which we store the integrals
  //
  //     unsigned int blockId =  childlessList[kk];
  //     OctreeBlock<dim> *block =  blocks[blockId];
  //     double delta = block->GetDelta();
  //     Point<dim> deltaHalf;
  //     for (unsigned int i=0; i<dim; i++)
  //       deltaHalf(i) = delta/2.;
  //     Point<dim> blockCenter = block->GetPMin()+deltaHalf;
  //
  //     // at this point, we get the list of quad nodes for the current block,
  //     // and loop over it
  //     std::map <cell_it, std::vector <unsigned int> > blockQuadPointsList =
  //     block->GetBlockQuadPointsList();
  //
  //     typename std::map <cell_it, std::vector<types::global_dof_index>
  //     >::iterator it; for (it = blockQuadPointsList.begin(); it !=
  //     blockQuadPointsList.end(); it++)
  //       {
  //         // for each cell in the list, we get the list of its quad nodes
  //         // present in the current block
  //         cell_it cell = (*it).first;
  //         std::vector <unsigned int> &cellQuadPoints = (*it).second;
  //
  //         // the vectors in the structures that we have previously cleared
  //         // neet to be resized
  //         elemMultipoleExpansionsKer1[blockId][cell].resize(dofs_per_cell);
  //         elemMultipoleExpansionsKer2[blockId][cell].resize(dofs_per_cell);
  //
  //         // the vectors are now initialized with an empty multipole
  //         expansion
  //         // centered in the current block center
  //         for (unsigned int j=0; j<fma_dh->get_fe().dofs_per_cell; ++j)
  //           {
  //             elemMultipoleExpansionsKer1[blockId][cell][j] =
  //               MultipoleExpansion(trunc_order, blockCenter,
  //               &assLegFunction);
  //             elemMultipoleExpansionsKer2[blockId][cell][j] =
  //               MultipoleExpansion(trunc_order, blockCenter,
  //               &assLegFunction);
  //           }
  //
  //         // the contribution of each quadrature node (which can be seen as a
  //         // source with a certain strength) is introduced in the
  //         // multipole expansion with the appropriate methods (AddNormDer
  //         // for neumann matrix integrals, Add for dirichlet matrix
  //         // integrals)
  //         for (unsigned int k=0; k<cellQuadPoints.size(); ++k)
  //           {
  //             unsigned int q = cellQuadPoints[k];
  //             for (unsigned int j=0; j<fma_dh->get_fe().dofs_per_cell; ++j)
  //               {
  //                 elemMultipoleExpansionsKer1[blockId][cell][j].AddNormDer(
  //                 quadShapeFunValues[cell][q][j]*
  //                 quadJxW[cell][q]/4/numbers::PI, quadPoints[cell][q],
  //                 quadNormals[cell][q]);
  //                 elemMultipoleExpansionsKer2[blockId][cell][j].Add(
  //                 quadShapeFunValues[cell][q][j]*
  //                 quadJxW[cell][q]/4/numbers::PI, quadPoints[cell][q]);
  //               }
  //           } // end loop on cell quadrature points in the block
  //
  //       } // end loop over cells in the block
  //
  //   } // fine loop sui blocchi childless



  pcout << "...done computing multipole integrals" << std::endl;
}

// // A helper function that we may use to know the overall nemuber of M2L
// operations
// // on each processor.
// template <int dim>
// void BEMFMA<dim>::compute_m2l_flags()
// {
//   pcout<<"Partitioning the descending phase Multipole2Local"<<std::endl;
//   m2l_flags.resize(num_octree_levels+1);
//   std::vector<types::global_dof_index>
//   m2l_operations_per_level(num_octree_levels+1);
//   std::vector<std::vector<types::global_dof_index> >
//   m2l_operations_per_block(num_octree_levels+1); types::global_dof_index
//   my_total_operations=0;
//
//   for (unsigned int level = 1; level <  num_octree_levels + 1;  level++)
//     {
//       m2l_flags[level].resize(dofs_filled_blocks[level].size());
//       m2l_operations_per_block[level].resize(dofs_filled_blocks[level].size());
//       m2l_operations_per_level[level] = 0;
//       for (types::global_dof_index k = 0; k <
//       dofs_filled_blocks[level].size();  k++) // loop over blocks of each
//       level
//         {
//           m2l_operations_per_block[level][k] = 0;
//           types::global_dof_index jj =  dofs_filled_blocks[level][k];
//
//           OctreeBlock<dim> *block1 = blocks[jj];
//           for (unsigned int subLevel = 0; subLevel <
//           block1->NumNearNeighLevels();  subLevel++)
//             {
//               m2l_operations_per_block[level][k] +=
//               block1->GetNonIntList(subLevel).size();
//               m2l_operations_per_level[level] +=
//               block1->GetNonIntList(subLevel).size();
//             }
//         }
//       unsigned int test = n_mpi_processes;
//       types::global_dof_index operations_per_proc =
//       m2l_operations_per_level[level]/test;     //(int)
//       ceil(m2l_operations_per_level[level]/test); int rest_op =
//       m2l_operations_per_level[level]%test; types::global_dof_index
//       my_operations=0; types::global_dof_index cumulative_check=0; unsigned
//       int proc=0; types::global_dof_index k=0;
//       std::vector<types::global_dof_index> m2l_operations_per_proc(test);
//       std::vector<types::global_dof_index> blocks_per_proc(test);
//       for (types::global_dof_index k = 0; k <
//       dofs_filled_blocks[level].size();  k++) // loop over blocks of each
//       level
//         {
//           // m2l_operations_per_block[level][k] = 0;
//           types::global_dof_index jj =  dofs_filled_blocks[level][k];
//           OctreeBlock<dim> *block1 = blocks[jj];
//           std::vector <types::global_dof_index> nodesBlk1Ids =
//           block1->GetBlockNodeList(); bool on_process = false; for (auto ind
//           : nodesBlk1Ids)
//             {
//               if (this_cpu_set.is_element(ind))
//                 {
//                   on_process = true;
//                   break;
//                 }
//             }
//           if (on_process)
//             {
//               for (unsigned int subLevel = 0; subLevel <
//               block1->NumNearNeighLevels();  subLevel++)
//                 {
//                   my_operations += m2l_operations_per_block[level][k];
//                   my_total_operations += m2l_operations_per_block[level][k];
//                 }
//             }
//         }
//       std::cout<<"level --- mpi_proc --- ops"<<std::endl;
//       std::cout<<level<<" --- "<<this_mpi_process<<" ---
//       "<<my_operations<<std::endl;
//       // while (k <  dofs_filled_blocks[level].size() && proc<test) // loop
//       over blocks of each level
//       // {
//       //   if(my_operations >= operations_per_proc)
//       //   {
//       //
//       //     rest_op -= my_operations - operations_per_proc;
//       //     // pcout<<"On processor "<< proc<<", we have
//       "<<my_operations<<"operations, cumulatively we are taking
//       "<<cumulative_check<<" m2l ops over a total of
//       "<<m2l_operations_per_level[level]<<std::endl;
//       //     proc += 1;
//       //     my_operations = 0;
//       //   }
//       //   m2l_flags[level][k]=proc;
//       //   my_operations += m2l_operations_per_block[level][k];
//       //   cumulative_check += m2l_operations_per_block[level][k];
//       //   m2l_operations_per_proc[proc] +=
//       m2l_operations_per_block[level][k];
//       //   blocks_per_proc[proc] += 1;
//       //   k+=1;
//       // }
//       // pcout<<"LEVEL "<<level<<std::endl;
//       // pcout<<"Rest is "<<rest_op<<" last block is "<<k<<" , total of
//       "<<dofs_filled_blocks[level].size()<<
//       //        " operations taken "<<cumulative_check<<" over a total of
//       "<<m2l_operations_per_level[level]<<std::endl;
//       // for(proc = 0 ; proc < test; ++proc)
//       //   pcout<<"On processor "<< proc<<", we have
//       "<<m2l_operations_per_proc[proc]<<" m2l operations, and "
//       //        <<blocks_per_proc[proc]<<" blocks over a total of
//       "<<dofs_filled_blocks[level].size()<<std::endl;
//     }
//
//   std::cout<<"finalcountmpi_proc --- ops"<<std::endl;
//   std::cout<<this_mpi_process<<" --- "<<my_total_operations<<std::endl;
//
//
//
// }

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
        double     delta = blocks[ii]->GetDelta();
        Point<dim> deltaHalf;
        for (unsigned int i = 0; i < dim; i++)
          deltaHalf(i) = delta / 2.;

        Point<dim> blockCenter = blocks[ii]->GetPMin() + deltaHalf;
        blockMultipoleExpansionsKer1[ii] =
          MultipoleExpansion(trunc_order, blockCenter, &(assLegFunction));
        blockMultipoleExpansionsKer2[ii] =
          MultipoleExpansion(trunc_order, blockCenter, &(assLegFunction));
      }
  };


  parallel_for(blocked_range<unsigned int>(0, num_blocks, tbb_granularity),
               f_creation_tbb);

  // auto f_creation = [] (types::global_dof_index ii, const BEMFMA<dim>
  // *foo_fma)
  // {
  //
  //   double delta = foo_fma->blocks[ii]->GetDelta();
  //   Point<dim> deltaHalf;
  //   for (unsigned int i=0; i<dim; i++)
  //     deltaHalf(i) = delta/2.;
  //
  //   Point<dim> blockCenter =  foo_fma->blocks[ii]->GetPMin()+deltaHalf;
  //   foo_fma->blockMultipoleExpansionsKer1[ii] =
  //   MultipoleExpansion(foo_fma->trunc_order, blockCenter,
  //   &(foo_fma->assLegFunction)); foo_fma->blockMultipoleExpansionsKer2[ii] =
  //   MultipoleExpansion(foo_fma->trunc_order, blockCenter,
  //   &(foo_fma->assLegFunction));
  // };
  //
  // Threads::TaskGroup<> group_creation;
  // for (types::global_dof_index ii = 0; ii <  num_blocks ; ii++)
  //   group_creation += Threads::new_task ( static_cast<void
  //   (*)(types::global_dof_index, const BEMFMA<dim> *)> (f_creation), ii,
  //   this);
  // group_creation.join_all();
  // After a discussion with Wolfgang it should be like the following
  // group_childless_creation += Threads::new_task ( std::function<void ()>
  // (std::bind( f_childless_creation, kk, this, phi_values,
  // dphi_dn_values)));//( static_cast<void (*)(types::global_dof_index, const
  // BEMFMA<dim> *, const Vector<double> &, const Vector<double> &)>
  // (f_childless_creation), kk, this, phi_values, dphi_dn_values); However if
  // we do so it does not work...

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
        OctreeBlock<dim> *      block   = blocks[blockId];

        double     delta = blocks[blockId]->GetDelta();
        Point<dim> deltaHalf;
        for (unsigned int i = 0; i < dim; i++)
          deltaHalf(i) = delta / 2.;
        // Point<dim> blockCenter =  blocks[blockId]->GetPMin()+deltaHalf;

        std::map<cell_it, std::vector<types::global_dof_index>>
          blockQuadPointsList = block->GetBlockQuadPointsList();

        // we loop on the cells of the quad points in the block: remember that
        // for each cell with a node in the block, we had created a set of
        // dofs_per_cell multipole expansion, representing the (partial)
        // integral on each cell

        std::vector<types::global_dof_index> my_local_dof_indices(
          fma_dh->get_fe().dofs_per_cell);

        typename std::map<cell_it,
                          std::vector<types::global_dof_index>>::iterator it;
        for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end();
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
                blockMultipoleExpansionsKer2.at(blockId).Add(
                  elemMultipoleExpansionsKer2[blockId][cell][jj],
                  dphi_dn_values(my_local_dof_indices[jj]));
                blockMultipoleExpansionsKer1.at(blockId).Add(
                  elemMultipoleExpansionsKer1[blockId][cell][jj],
                  phi_values(my_local_dof_indices[jj]));
              }
          } // end loop ond block elements
      }
  };

  parallel_for(blocked_range<types::global_dof_index>(0,
                                                      childlessList.size(),
                                                      tbb_granularity),
               f_childless_creation_tbb);

  // auto f_childless_creation = [] (types::global_dof_index kk, const
  // BEMFMA<dim> *foo_fma, const Vector<double> &phi_values, const
  // Vector<double> &dphi_dn_values)
  // {
  //
  //   // for each block we get the center and the quad points
  //
  //   types::global_dof_index blockId =  foo_fma->childlessList[kk];
  //   OctreeBlock<dim> *block =  foo_fma->blocks[blockId];
  //
  //   double delta =  foo_fma->blocks[blockId]->GetDelta();
  //   Point<dim> deltaHalf;
  //   for (unsigned int i=0; i<dim; i++)
  //     deltaHalf(i) = delta/2.;
  //   //Point<dim> blockCenter =  blocks[blockId]->GetPMin()+deltaHalf;
  //
  //   std::map <cell_it, std::vector <types::global_dof_index> >
  //   blockQuadPointsList = block->GetBlockQuadPointsList();
  //
  //   // we loop on the cells of the quad points in the block: remember that
  //   for each cell with a node in the
  //   // block, we had created a set of dofs_per_cell multipole expansion,
  //   representing
  //   // the (partial) integral on each cell
  //
  //   std::vector<types::global_dof_index>
  //   my_local_dof_indices(foo_fma->fma_dh->get_fe().dofs_per_cell);
  //
  //   typename std::map <cell_it, std::vector<types::global_dof_index>
  //   >::iterator it; for (it = blockQuadPointsList.begin(); it !=
  //   blockQuadPointsList.end(); it++)
  //     {
  //       cell_it cell = (*it).first;
  //       cell->get_dof_indices(my_local_dof_indices);
  //
  //       // for each cell we get the dof_indices, and add to the block
  //       multipole expansion,
  //       // the integral previously computed, multiplied by the phi or dphi_dn
  //       value at the
  //       // corresponding dof of the cell. A suitable MultipoleExpansion class
  //       method has been
  //       // created for this purpose
  //
  //       for (unsigned int jj=0; jj < foo_fma->fma_dh->get_fe().dofs_per_cell;
  //       ++jj)
  //         {
  //           foo_fma->blockMultipoleExpansionsKer2.at(blockId).Add(foo_fma->elemMultipoleExpansionsKer2[blockId][cell][jj],dphi_dn_values(my_local_dof_indices[jj]));
  //           foo_fma->blockMultipoleExpansionsKer1.at(blockId).Add(foo_fma->elemMultipoleExpansionsKer1[blockId][cell][jj],phi_values(my_local_dof_indices[jj]));
  //         }
  //     } //end loop ond block elements
  // };
  //
  // Threads::TaskGroup<> group_childless_creation;
  // for (types::global_dof_index kk = 0; kk <  childlessList.size(); kk++)
  //   group_childless_creation += Threads::new_task ( static_cast<void
  //   (*)(types::global_dof_index, const BEMFMA<dim> *, const Vector<double> &,
  //   const Vector<double> &)> (f_childless_creation), kk, this, phi_values,
  //   dphi_dn_values);
  // group_childless_creation.join_all();

  // {

  // // for each block we get the center and the quad points
  //
  // unsigned int blockId =  childlessList[kk];
  // OctreeBlock<dim> *block =  blocks[blockId];
  //
  // delta =  blocks[blockId]->GetDelta();
  // Point<dim> deltaHalf;
  // for (unsigned int i=0; i<dim; i++)
  //   deltaHalf(i) = delta/2.;
  // //Point<dim> blockCenter =  blocks[blockId]->GetPMin()+deltaHalf;
  //
  // std::map <cell_it, std::vector <unsigned int> > blockQuadPointsList =
  // block->GetBlockQuadPointsList();
  //
  // // we loop on the cells of the quad points in the block: remember that for
  // each cell with a node in the
  // // block, we had created a set of dofs_per_cell multipole expansion,
  // representing
  // // the (partial) integral on each cell
  //
  // typename std::map <cell_it, std::vector<types::global_dof_index>
  // >::iterator it; for (it = blockQuadPointsList.begin(); it !=
  // blockQuadPointsList.end(); it++)
  //   {
  //     cell_it cell = (*it).first;
  //     cell->get_dof_indices(local_dof_indices);
  //
  //     // for each cell we get the dof_indices, and add to the block multipole
  //     expansion,
  //     // the integral previously computed, multiplied by the phi or dphi_dn
  //     value at the
  //     // corresponding dof of the cell. A suitable MultipoleExpansion class
  //     method has been
  //     // created for this purpose
  //
  //     for (unsigned int jj=0; jj < fma_dh->get_fe().dofs_per_cell; ++jj)
  //       {
  //         blockMultipoleExpansionsKer2.at(blockId).Add(elemMultipoleExpansionsKer2[blockId][cell][jj],dphi_dn_values(local_dof_indices[jj]));
  //         blockMultipoleExpansionsKer1.at(blockId).Add(elemMultipoleExpansionsKer1[blockId][cell][jj],phi_values(local_dof_indices[jj]));
  //       }
  //   } //end loop ond block elements
  // } // end loop on childless blocks


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
           AscendCopyData &        copy_data,
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
      // copy_data.local_level_indices[kk] = start + kk;
      copy_data.parentId = (*block_it)->GetParentId();
      // We set the center to the parentId center.
      copy_data.translatedBlockMultipoleExpansionKer1.SetCenter(
        this->blockMultipoleExpansionsKer1.at(copy_data.parentId).GetCenter());
      copy_data.translatedBlockMultipoleExpansionKer2.SetCenter(
        this->blockMultipoleExpansionsKer2.at(copy_data.parentId).GetCenter());
      // std::cout<<"Before Add:
      // "<<copy_data.translatedBlockMultipoleExpansionKer1.GetCoeff(5, 1)<<"
      // "<<copy_data.translatedBlockMultipoleExpansionKer2.GetCoeff(5,
      // 1)<<std::endl;

      // We translate the children blocks to the local array.
      copy_data.translatedBlockMultipoleExpansionKer1.Add(
        this->blockMultipoleExpansionsKer1.at(kk));
      copy_data.translatedBlockMultipoleExpansionKer2.Add(
        this->blockMultipoleExpansionsKer2.at(kk));
      // std::cout<<"After Add:
      // "<<copy_data.translatedBlockMultipoleExpansionKer1.GetCoeff(5, 1)<<"
      // "<<copy_data.translatedBlockMultipoleExpansionKer2.GetCoeff(5,
      // 1)<<std::endl;
    };

  // The copier function, it copies the value from the local array to the parent
  // blocks
  auto f_copier_ascend = [this](const AscendCopyData &copy_data) {
    // Now we just need to copy on the real global memory. Since the center is
    // the same we just need to add the value to the expansion wothout any
    // further translation. The MultipoleExpansion class takes care of that
    // automatically.
    this->blockMultipoleExpansionsKer1.at(copy_data.parentId)
      .Add(copy_data.translatedBlockMultipoleExpansionKer1);
    this->blockMultipoleExpansionsKer2.at(copy_data.parentId)
      .Add(copy_data.translatedBlockMultipoleExpansionKer2);
  };
  for (unsigned int level = num_octree_levels; level > 0; level--)

    {
      AscendScratchData sample_scratch;

      AscendCopyData sample_copy(this);

      // We need to be sure that there actually are some blocks in the current
      // level, then we need to loop over all blocks on the level and perform
      // all the translations. We need a bind for the worker in order to deal
      // with the additional unsigned int startLevel[level].
      if (endLevel[level] >= startLevel[level])
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

      // for (typename std::vector<OctreeBlock<dim> *>::iterator
      // blocky=blocks.begin()+startLevel[level];
      // blocky<blocks.begin()+endLevel[level]+1; ++blocky)
      //
      //   {
      //     f_worker_ascend(blocky, sample_scratch, sample_copy,
      //     startLevel[level]);
      //   }
      //   f_copier_ascend(sample_copy);

      // for (unsigned int kk =  startLevel[level]; kk <  endLevel[level]+1;
      // kk++)
      //
      //   {
      //     unsigned int parentId =  blocks[kk]->GetParentId();
      //     std::cout<<kk<<std::endl;
      //
      //     blockMultipoleExpansionsKer1.at(parentId).Add(sample_copy.local_level_indicesblockMultipoleExpansionKer1.at(kk-startLevel[level]));
      //     blockMultipoleExpansionsKer2.at(parentId).Add(sample_copy.local_level_indicesblockMultipoleExpansionKer2.at(kk-startLevel[level]));
      //
      //   } // end loop over blocks of a level


      // for (unsigned int kk =  startLevel[level]; kk <  endLevel[level]+1;
      // kk++)
      //
      //   {
      //     unsigned int parentId =  blocks[kk]->GetParentId();
      //     std::cout<<kk<<std::endl;
      //
      //     blockMultipoleExpansionsKer1.at(parentId).Add(blockMultipoleExpansionsKer1.at(kk));
      //     blockMultipoleExpansionsKer2.at(parentId).Add(blockMultipoleExpansionsKer2.at(kk));
      //
      //   } // end loop over blocks of a level

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
// The following functions takes care of the descending phase of the FMA.
template <int dim>
void
BEMFMA<dim>::multipole_matr_vect_products(
  const TrilinosWrappers::MPI::Vector &phi_values,
  const TrilinosWrappers::MPI::Vector &dphi_dn_values,
  TrilinosWrappers::MPI::Vector &      matrVectProdN,
  TrilinosWrappers::MPI::Vector &      matrVectProdD) const
{
  pcout << "Computing multipole matrix-vector products... " << std::endl;
  Teuchos::TimeMonitor LocalTimer(*MatrVec);
  // we start re-initializing matrix-vector-product vectors
  // matrVectProdN = 0;
  // matrVectProdD = 0;

  // and here we compute the direct integral contributions (stored in two sparse
  // matrices)
  // const Vector<double> phi_values_loc(phi_values);
  // const Vector<double> dphi_dn_values_loc(dphi_dn_values);
  // Vector<double> matrVectProdD(fma_dh->n_dofs());
  // Vector<double> matrVectProdN(fma_dh->n_dofs());
  prec_neumann_matrix.vmult(matrVectProdN, phi_values);
  prec_dirichlet_matrix.vmult(matrVectProdD, dphi_dn_values);

  // if(this_cpu_set.is_element(55))
  //   std::cout<<matrVectProdN[55]<<" "<<phi_values[55]<<"
  //   "<<matrVectProdD[55]<<" "<<dphi_dn_values[55]<<std::endl;

  // from here on, we compute the multipole expansions contributions

  // we start cleaning past sessions
  blockLocalExpansionsKer1.clear();
  blockLocalExpansionsKer2.clear();

  blockLocalExpansionsKer1.resize(num_blocks);
  blockLocalExpansionsKer2.resize(num_blocks);

  // we declare some familiar variables that will be useful in the method
  std::vector<Point<dim>> support_points(fma_dh->n_dofs());
  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*fma_mapping,
                                                     *fma_dh,
                                                     support_points);
  std::vector<types::global_dof_index> local_dof_indices(
    fma_dh->get_fe().dofs_per_cell);



  // First of all we need to create all the empty expansiones for all the
  // blocks. This is an embarassingly parallel operation that we can perform
  // using the ThreadGroup strategy without requiring any synchronization time.

  auto f_local_creation_tbb = [this](blocked_range<types::global_dof_index> r) {
    for (types::global_dof_index ii = r.begin(); ii < r.end(); ++ii)
      {
        double     delta = blocks[ii]->GetDelta();
        Point<dim> deltaHalf;
        for (unsigned int i = 0; i < dim; i++)
          deltaHalf(i) = delta / 2.;

        Point<dim> blockCenter = blocks[ii]->GetPMin() + deltaHalf;
        blockLocalExpansionsKer1[ii] =
          LocalExpansion(trunc_order, blockCenter, &(assLegFunction));
        blockLocalExpansionsKer2[ii] =
          LocalExpansion(trunc_order, blockCenter, &(assLegFunction));
      }
  };

  parallel_for(blocked_range<types::global_dof_index>(0,
                                                      num_blocks,
                                                      tbb_granularity),
               f_local_creation_tbb);

  // auto f_local_creation = [] (types::global_dof_index ii, const BEMFMA<dim>
  // *foo_fma)
  // {
  //
  //   double delta = foo_fma->blocks[ii]->GetDelta();
  //   Point<dim> deltaHalf;
  //   for (unsigned int i=0; i<dim; i++)
  //     deltaHalf(i) = delta/2.;
  //
  //   Point<dim> blockCenter =  foo_fma->blocks[ii]->GetPMin()+deltaHalf;
  //   foo_fma->blockLocalExpansionsKer1[ii] =
  //   LocalExpansion(foo_fma->trunc_order, blockCenter,
  //   &(foo_fma->assLegFunction)); foo_fma->blockLocalExpansionsKer2[ii] =
  //   LocalExpansion(foo_fma->trunc_order, blockCenter,
  //   &(foo_fma->assLegFunction));
  // };
  // Threads::TaskGroup<> group_local_creation;
  // for (types::global_dof_index ii = 0; ii <  num_blocks ; ii++)
  //   group_local_creation += Threads::new_task ( static_cast<void
  //   (*)(types::global_dof_index, const BEMFMA<dim> *)> (f_local_creation),
  //   ii, this);
  // group_local_creation.join_all();

  // // here we loop on all the blocks and build an empty local expansion for
  // // each of them
  // for (unsigned int ii = 0; ii <  num_blocks ; ii++)
  //
  //   {
  //     delta =  blocks[ii]->GetDelta();
  //     Point<dim> deltaHalf;
  //     for (unsigned int i=0; i<dim; i++)
  //       deltaHalf(i) = delta/2.;
  //     Point<dim> blockCenter =  blocks[ii]->GetPMin()+deltaHalf;
  //
  //     blockLocalExpansionsKer1[ii] =  LocalExpansion(trunc_order,
  //     blockCenter, &assLegFunction); blockLocalExpansionsKer2[ii] =
  //     LocalExpansion(trunc_order, blockCenter, &assLegFunction);
  //   }



  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////

  // In order to perform the descending phase properly we need WorkStream.
  // Inside the worker we check the IndexSet for the MPI parallelisation. The
  // scaracth is empty once again since this class only has to copy from parents
  // to children. As in the ascending phase we use captures in order to pass
  // global information to the worker and the copier
  struct DescendScratchData
  {};

  // Basically we are applying a local to global operation so we need only great
  // care in the copy.
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
    };

    // The working copy constructor for the copy structure
    DescendCopyData(const DescendCopyData &in_vec)
    {
      blockLocalExpansionKer1 = in_vec.blockLocalExpansionKer1;
      blockLocalExpansionKer2 = in_vec.blockLocalExpansionKer2;
      start                   = in_vec.start;
    };

    // The Destructor needs to make foo_fma to point to NULL (for this reason it
    // is mutable const)
    ~DescendCopyData(){};

    types::global_dof_index              start;
    types::global_dof_index              blockId;
    LocalExpansion                       blockLocalExpansionKer1;
    LocalExpansion                       blockLocalExpansionKer2;
    Vector<double>                       matrVectorProductContributionKer1;
    Vector<double>                       matrVectorProductContributionKer2;
    std::vector<types::global_dof_index> nodesBlk1Ids;
  };


  // The worker function, it copies the value about the parent from the global
  // memory at a certain level to local array in copy.data
  auto f_worker_Descend = [this, &support_points](
                            std::vector<types::global_dof_index>::const_iterator
                              block_it_id,
                            DescendScratchData &,
                            DescendCopyData &             copy_data,
                            const types::global_dof_index start) {
    // types::global_dof_index kk = std::distance(this->blocks.begin(),
    // block_it);


    copy_data.start            = start;
    copy_data.blockId          = *block_it_id;
    types::global_dof_index kk = *block_it_id;
    OctreeBlock<dim> *      block_it;
    block_it = this->blocks[*block_it_id];
    // copy_data.local_level_indices[kk] = start + kk;
    //*****************definire chi e' on_process qui
    Point<dim>     center = this->blockLocalExpansionsKer1.at(kk).GetCenter();
    LocalExpansion dummy(this->trunc_order, center, &(this->assLegFunction));
    copy_data.blockLocalExpansionKer1 = dummy;
    copy_data.blockLocalExpansionKer2 = dummy;
    // std::cout<<copy_data.blockLocalExpansionKer1.GetCoeff(5,1)<<"
    // "<<copy_data.blockLocalExpansionKer2.GetCoeff(5,1)<<std::endl;
    // pcout<<"Block "<<jj<<std::endl;
    unsigned int level = 0;
    for (unsigned int lev = 0; lev < this->num_octree_levels + 1; lev++)
      {
        // std::cout<<this->startLevel[lev]<<"  "<<kk<<"
        // "<<this->endLevel[lev]<<std::endl;
        if (kk >= this->startLevel[lev] && kk <= this->endLevel[lev])
          {
            level = lev;
            break;
          }
      }

    types::global_dof_index startBlockLevel = this->startLevel[level];
    types::global_dof_index endBlockLevel   = this->endLevel[level];
    std::vector<types::global_dof_index> nodesBlk1Ids =
      block_it->GetBlockNodeList();
    bool on_process = false;
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
        // std::cout<<"Level: "<<level<<" "<<kk<<"("<<jj<<")
        // "<<nodesBlk1Ids.size()<<"  "<<startBlockLevel<<std::endl;
        // the local expansion of the parent must be translated down into the
        // current block
        types::global_dof_index parentId = block_it->GetParentId();
        copy_data.blockLocalExpansionKer1.Add(
          this->blockLocalExpansionsKer1.at(parentId));
        copy_data.blockLocalExpansionKer2.Add(
          this->blockLocalExpansionsKer2.at(parentId));


        for (unsigned int subLevel = 0;
             subLevel < block_it->NumNearNeighLevels();
             subLevel++)
          {
            // we get the nonIntList of each block

            std::set<types::global_dof_index> nonIntList =
              block_it->GetNonIntList(subLevel);
            std::set<cell_it> nonIntListElemsBlk1;

            // we start converting into local expansions, all the multipole
            // expansions of all the blocks in nonIntList, that are of the same
            // size (level) of the current block. To perform the conversion, we
            // use another member of the LocalExpansion class. Note that all the
            // contributions to the integrals of blocks bigger than current
            // block had already been considered in the direct integrals method
            // (the bounds of the local and multipole expansion do not hold in
            // such case)


            for (std::set<types::global_dof_index>::iterator pos1 =
                   nonIntList.lower_bound(startBlockLevel);
                 pos1 != nonIntList.upper_bound(endBlockLevel);
                 pos1++) // loop over NNs of NNs and add them to intList
              {
                // std::cout<<"NonIntListPart2 Blocks: "<<*pos1<<" "<<std::endl;
                types::global_dof_index block2Id = *pos1;
                copy_data.blockLocalExpansionKer1.Add(
                  this->blockMultipoleExpansionsKer1[block2Id]);
                copy_data.blockLocalExpansionKer2.Add(
                  this->blockMultipoleExpansionsKer2[block2Id]);

                /*////////this is for a check///////////////////////
                  OctreeBlock<dim>* block2 =  blocks[block2Id];
                std::vector<types::global_dof_index> nodesBlk1Ids =
                (*block_it)->GetBlockNodeList(); std::map <cell_it,
                std::vector<types::global_dof_index> > blockQuadPointsList =
                block2->GetBlockQuadPointsList(); typename std::map <cell_it,
                std::vector<types::global_dof_index> >::iterator it; for (it =
                blockQuadPointsList.begin(); it != blockQuadPointsList.end();
                it++)
                                      {
                                for (unsigned int i=0; i<(*it).second.size();
                i++)
                                    {
                                    for (unsigned int kk=0;
                kk<nodesBlk1Ids.size(); kk++)
                                       integralCheck[nodesBlk1Ids[kk]][(*it).first]
                += 1;
                        }
                                            }
                                        ////////////////////////////*/
              } // end loop over well separated blocks of the same size (level)


            // loop over well separated blocks of the smaller size (level)----->
            // use multipoles without local expansions

            // we now have to loop over blocks in the nonIntList that are
            // smaller than the current blocks: in this case the bound for the
            // conversion of a multipole into local expansion does not hold, but
            // the bound for the evaluation of the multipole expansions does
            // hold: thus, we will simply evaluate the multipole expansions of
            // such blocks for each node in the block
            for (std::set<types::global_dof_index>::iterator pos1 =
                   nonIntList.upper_bound(endBlockLevel);
                 pos1 != nonIntList.end();
                 pos1++)
              {
                // pcout<<"NonIntListPart3 Blocks: "<<*pos1<<" ";
                types::global_dof_index block2Id = *pos1;
                // std::vector <cell_it> elemBlk2Ids =
                // block2.GetBlockElementsList();
                //Teuchos::TimeMonitor LocalTimer(*LocEval);

                for (types::global_dof_index ii = 0; ii < nodesBlk1Ids.size();
                     ii++) // loop over each node of (*block_it)
                  {
                    // std::cout<<nodesBlk1Ids.at(ii)<<std::endl;
                    const Point<dim> &nodeBlk1 =
                      support_points[nodesBlk1Ids.at(ii)];
                    copy_data.matrVectorProductContributionKer1(ii) +=
                      this->blockMultipoleExpansionsKer1[block2Id].Evaluate(
                        nodeBlk1);
                    copy_data.matrVectorProductContributionKer2(ii) +=
                      this->blockMultipoleExpansionsKer2[block2Id].Evaluate(
                        nodeBlk1);

                    /*//////this is for a check/////////////////////////
                          OctreeBlock<dim>* block2 =  blocks[block2Id];
                          std::map <cell_it,
                       std::vector<types::global_dof_index> >
                                            blockQuadPointsList =
                       block2->GetBlockQuadPointsList(); typename std::map
                       <cell_it, std::vector<types::global_dof_index>
                       >::iterator it; for (it = blockQuadPointsList.begin(); it
                       != blockQuadPointsList.end(); it++)
                                                {
                                          for (unsigned int i=0;
                       i<(*it).second.size(); i++)
                                              {
                                           integralCheck[nodesBlk1Ids[ii]][(*it).first]
                       += 1;
                                  }
                                                      }
                                                  ////////////////////////////*/
                  }
              } // end loop over well separated blocks of smaller size (level)


          } // end loop over all sublevels in  nonIntlist
      }

    // std::cout<<"After Add: "<<copy_data.blockLocalExpansionKer1.GetCoeff(5,
    // 1)<<" "<<copy_data.blockLocalExpansionKer2.GetCoeff(5, 1)<<std::endl;
  };

  // The copier function, it copies the value from the local array to the parent
  // blocks and it fills the actual parallel vector of the matrix vector
  // products.
  auto f_copier_Descend = [this, &matrVectProdD, &matrVectProdN](
                            const DescendCopyData &copy_data) {
    this->blockLocalExpansionsKer1.at(copy_data.blockId)
      .Add(copy_data.blockLocalExpansionKer1);
    this->blockLocalExpansionsKer2.at(copy_data.blockId)
      .Add(copy_data.blockLocalExpansionKer2);

    for (types::global_dof_index i = 0; i < copy_data.nodesBlk1Ids.size(); ++i)
      {
        matrVectProdD(copy_data.nodesBlk1Ids[i]) +=
          copy_data.matrVectorProductContributionKer2(i);
        matrVectProdN(copy_data.nodesBlk1Ids[i]) +=
          copy_data.matrVectorProductContributionKer1(i);
      }
  };
  // so, here we loop over levels, starting form lower levels (bigger blocks)
  for (unsigned int level = 1; level < num_octree_levels + 1; level++)
    {
      Teuchos::TimeMonitor LocalTimer(*LocEval);
      types::global_dof_index startBlockLevel = startLevel[level];
      types::global_dof_index endBlockLevel   = endLevel[level];

      // to reduce computational cost, we decide to loop on the list of blocks
      // which contain at least one node (the local and multipole expansion will
      // be finally evaluated at the nodes positions)

      // TODO WE COULD SPLIT THIS LOOP OVER THE BLOCKS OVER ALL THE PROCESSORS,
      // THEN DO A TRILINOS.ADD WITH DIFFERENT MAPS. HERE WE NEED A FULL ONE.


      DescendScratchData sample_scratch;

      DescendCopyData sample_copy(this);

      // The Workstream has to run only if there are blocks in the current
      // level. Then it basically performs a loop over all the blocks in the
      // current level.
      if (endBlockLevel >= startBlockLevel)
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
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  /*
    // here the descending phase of the algorithm begins: we will start from the
    top level, and loop
    // on each block of each level: all the needed multipole integrals
    contributions must be
    // introduced in each specific block local expansion. to do this, for each
    block, we first need to
    // translate into the block center the parent block local expansion (it
    accounts for all far blocks
    // contribution of the parent, and all blocks that are far from the parent,
    are far from its
    // children); then, we have to consider all the blocks in the current
    block's nonIntList,
    // and convert  their multipole expansions to the local expansion of the
    current block
    // (there is only one exception to this procedure, and we'll comment it in
    the following)

    // so, here we loop over levels, starting form lower levels (bigger blocks)
    for (unsigned int level = 1; level <  num_octree_levels + 1;  level++)

      {
        // pcout<<"processing level "<<level <<"  of  "<<
    num_octree_levels<<std::endl;

        // we get the ids of the first and last block of the level
        unsigned int startBlockLevel =  startLevel[level];
        unsigned int endBlockLevel =  endLevel[level];

        // to reduce computational cost, we decide to loop on the list of blocks
    which
        // contain at least one node (the local and multipole expansion will be
    finally evaluated
        // at the nodes positions)

        // TODO WE COULD SPLIT THIS LOOP OVER THE BLOCKS OVER ALL THE
    PROCESSORS, THEN DO A TRILINOS.ADD WITH
        // DIFFERENT MAPS. HERE WE NEED A FULL ONE.
        for (unsigned int k = 0; k <  dofs_filled_blocks[level].size();  k++) //
    loop over blocks of each level
          {

            //pcout<<"Block "<<jj<<std::endl;
            unsigned int jj =  dofs_filled_blocks[level][k];
            OctreeBlock<dim> *block1 =  blocks[jj];
            unsigned int block1Parent = block1->GetParentId();
            std::vector <unsigned int> nodesBlk1Ids =
    block1->GetBlockNodeList(); bool on_process = false; for(auto ind :
    nodesBlk1Ids)
            {
              if(this_cpu_set.is_element(ind))
              {
                on_process = true;
                break;
              }
            }
            if(on_process)
            {
            // here we get parent contribution to local expansion and transfer
    it in current
            // block: this operation requires a local expansion translation,
    implemented
            // in a specific LocalExpansion class member

            blockLocalExpansionsKer1[jj].Add(blockLocalExpansionsKer1[block1Parent]);
            blockLocalExpansionsKer2[jj].Add(blockLocalExpansionsKer2[block1Parent]);

            // for each block, loop over all sublevels in his NN list (this is
    because if a
            // block remains childless BEFORE the last level, at this point we
    need to compute
            // all its contributions up to the bottom level)


            for (unsigned int subLevel = 0; subLevel <
    block1->NumNearNeighLevels();  subLevel++)
              {
                // we get the nonIntList of each block

                std::set <unsigned int> nonIntList =
    block1->GetNonIntList(subLevel); std::set <cell_it> nonIntListElemsBlk1;

                // we start converting into local expansions, all the multipole
    expansions
                // of all the blocks in nonIntList, that are of the same size
    (level) of the
                // current block. To perform the conversion, we use another
    member of the
                // LocalExpansion class. Note that all the contributions to the
    integrals
                // of blocks bigger than current block had already been
    considered in
                // the direct integrals method (the bounds of the local and
    multipole
                // expansion do not hold in such case)


                for (std::set <unsigned int>::iterator pos1 =
    nonIntList.lower_bound(startBlockLevel); pos1 !=
    nonIntList.upper_bound(endBlockLevel);  pos1++) // loop over NNs of NNs and
    add them to intList
                  {
                    //pcout<<"NonIntListPart2 Blocks: "<<*pos1<<" ";
                    unsigned int block2Id = *pos1;
                    blockLocalExpansionsKer1[jj].Add(blockMultipoleExpansionsKer1[block2Id]);
                    blockLocalExpansionsKer2[jj].Add(blockMultipoleExpansionsKer2[block2Id]);

                   ////////this is for a check///////////////////////
                   //   OctreeBlock<dim>* block2 =  blocks[block2Id];
                   // std::vector<types::global_dof_index> nodesBlk1Ids =
    block1->GetBlockNodeList();
                   // std::map <cell_it, std::vector<types::global_dof_index> >
                    //                  blockQuadPointsList =
    block2->GetBlockQuadPointsList();
                   //                   typename std::map <cell_it,
    std::vector<types::global_dof_index> >::iterator it;
                   //                         for (it =
    blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
                    //                      {
                   //                 for (unsigned int i=0;
    i<(*it).second.size(); i++)
                   //                     {
                   //                     for (unsigned int kk=0;
    kk<nodesBlk1Ids.size(); kk++)
                   // integralCheck[nodesBlk1Ids[kk]][(*it).first] += 1;
                   //         }
                   //                             }
                                            ////////////////////////////
                  } // end loop over well separated blocks of the same size
    (level)


                // loop over well separated blocks of the smaller size
    (level)-----> use multipoles without local expansions

                // we now have to loop over blocks in the nonIntList that are
    smaller than the current
                // blocks: in this case the bound for the conversion of a
    multipole into local expansion
                // do not hold, but the bound for the evaluation of the
    multipole expansions does hold:
                // thus, we will simply evaluate the multipole expansions of
    such blocks for each node in
                // the clock

                for (std::set <unsigned int>::iterator pos1 =
    nonIntList.upper_bound(endBlockLevel); pos1 != nonIntList.end();  pos1++)
                  {
                    //pcout<<"NonIntListPart3 Blocks: "<<*pos1<<" ";
                    unsigned int block2Id = *pos1;
                    //std::vector <cell_it> elemBlk2Ids =
    block2.GetBlockElementsList(); Teuchos::TimeMonitor LocalTimer(*LocEval);

                    for (unsigned int ii = 0; ii < nodesBlk1Ids.size(); ii++)
    //loop over each node of block1
                      {
                        Point<dim> &nodeBlk1 =
    support_points[nodesBlk1Ids.at(ii)]; matrVectProdD(nodesBlk1Ids[ii]) +=
    blockMultipoleExpansionsKer2[block2Id].Evaluate(nodeBlk1);
                        matrVectProdN(nodesBlk1Ids[ii]) +=
    blockMultipoleExpansionsKer1[block2Id].Evaluate(nodeBlk1);

                        //////this is for a check/////////////////////////
                         //     OctreeBlock<dim>* block2 =  blocks[block2Id];
                         //     std::map <cell_it,
    std::vector<types::global_dof_index> >
                         //                       blockQuadPointsList =
    block2->GetBlockQuadPointsList();
                         //                       typename std::map <cell_it,
    std::vector<types::global_dof_index> >::iterator it;
                         //                             for (it =
    blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
                         //                           {
                         //                     for (unsigned int i=0;
    i<(*it).second.size(); i++)
                         //                         {
                         // integralCheck[nodesBlk1Ids[ii]][(*it).first] += 1;
                         //             }
                         //                                 }
                                                      ////////////////////////////
                      }
                  } // end loop over well separated blocks of smaller size
    (level)


              } // end loop over all sublevels in  nonIntlist
            }//end if for m2l flags
          } // end loop over blocks of each level

      } // end loop over levels

  */
  // finally, when the loop over levels is done, we need to evaluate local
  // expansions of all childless blocks, at each block node(s). This is an
  // embarassingly parallel operation so it can be easily performed using
  // ThreadGroup. We also check the IndexSet to perform the mixed TBB-MPI
  // parallelisation.
  auto f_local_evaluation_tbb =
    [this, &matrVectProdD, &matrVectProdN, &support_points](
      blocked_range<types::global_dof_index> r) {
      for (types::global_dof_index kk = r.begin(); kk < r.end(); ++kk)
        {
          types::global_dof_index              block1Id = childlessList[kk];
          OctreeBlock<dim> *                   block1   = blocks[block1Id];
          std::vector<types::global_dof_index> nodesBlk1Ids =
            block1->GetBlockNodeList();

          // loop over nodes of block
          for (types::global_dof_index ii = 0; ii < nodesBlk1Ids.size();
               ii++) // loop over each node of block1
            {
              if (this_cpu_set.is_element(nodesBlk1Ids[ii]))
                {
                  // Teuchos::TimeMonitor LocalTimer(*LocEval);

                  const Point<dim> &nodeBlk1 =
                    support_points[nodesBlk1Ids.at(ii)];
                  matrVectProdD(nodesBlk1Ids[ii]) +=
                    (blockLocalExpansionsKer2[block1Id]).Evaluate(nodeBlk1);
                  matrVectProdN(nodesBlk1Ids[ii]) +=
                    (blockLocalExpansionsKer1[block1Id]).Evaluate(nodeBlk1);
                }
            } // end loop over nodes
        }
    };

  parallel_for(blocked_range<types::global_dof_index>(0,
                                                      childlessList.size(),
                                                      tbb_granularity),
               f_local_evaluation_tbb);

  // auto f_local_evaluation = [] (types::global_dof_index kk,
  // TrilinosWrappers::MPI::Vector & matrVectProdD,
  // TrilinosWrappers::MPI::Vector & matrVectProdN, const BEMFMA<dim> *foo_fma,
  // const std::vector<Point<dim> > &support_points)
  // {
  //
  //   types::global_dof_index block1Id =  foo_fma->childlessList[kk];
  //   OctreeBlock<dim> *block1 =  foo_fma->blocks[block1Id];
  //   std::vector <types::global_dof_index> nodesBlk1Ids =
  //   block1->GetBlockNodeList();
  //
  //   // loop over nodes of block
  //   for (types::global_dof_index ii = 0; ii < nodesBlk1Ids.size(); ii++)
  //   //loop over each node of block1
  //     {
  //       if (foo_fma->this_cpu_set.is_element(nodesBlk1Ids[ii]))
  //         {
  //           //Teuchos::TimeMonitor LocalTimer(*LocEval);
  //
  //           const Point<dim> &nodeBlk1 = support_points[nodesBlk1Ids.at(ii)];
  //           matrVectProdD(nodesBlk1Ids[ii]) +=
  //           (foo_fma->blockLocalExpansionsKer2[block1Id]).Evaluate(nodeBlk1);
  //           matrVectProdN(nodesBlk1Ids[ii]) +=
  //           (foo_fma->blockLocalExpansionsKer1[block1Id]).Evaluate(nodeBlk1);
  //         }
  //     } // end loop over nodes
  // };
  //
  // Threads::TaskGroup<> group_local_evaluation;
  // for (types::global_dof_index kk = 0; kk <  childlessList.size(); kk++)
  //   group_local_evaluation += Threads::new_task ( static_cast<void
  //   (*)(types::global_dof_index, TrilinosWrappers::MPI::Vector &,
  //   TrilinosWrappers::MPI::Vector &, const BEMFMA<dim> *, const
  //   std::vector<Point<dim> > &)> (f_local_evaluation), kk, matrVectProdD,
  //   matrVectProdN, this, support_points);
  // group_local_evaluation.join_all();

  // for (unsigned int kk = 0; kk <  childlessList.size(); kk++)
  //
  //   {
  //     unsigned int block1Id =  childlessList[kk];
  //     OctreeBlock<dim> *block1 =  blocks[block1Id];
  //     std::vector <unsigned int> nodesBlk1Ids = block1->GetBlockNodeList();
  //
  //     // loop over nodes of block
  //     for (unsigned int ii = 0; ii < nodesBlk1Ids.size(); ii++) //loop over
  //     each node of block1
  //       {
  //         if(this_cpu_set.is_element(nodesBlk1Ids[ii]))
  //         {
  //           Teuchos::TimeMonitor LocalTimer(*LocEval);
  //
  //           Point<dim> &nodeBlk1 = support_points[nodesBlk1Ids.at(ii)];
  //           matrVectProdD(nodesBlk1Ids[ii]) +=
  //           (blockLocalExpansionsKer2[block1Id]).Evaluate(nodeBlk1);
  //           matrVectProdN(nodesBlk1Ids[ii]) +=
  //           (blockLocalExpansionsKer1[block1Id]).Evaluate(nodeBlk1);
  //         }
  //       } // end loop over nodes
  //
  //   } // end loop over childless blocks

  /*////////this is for a check//////////////////////
  for (unsigned int i = 0; i < fma_dh->n_dofs(); i++)
      {
    for (cell_it cell = fma_dh->begin_active(); cell != fma_dh->end(); ++cell)
          {
      pcout<<i<<" "<<cell<<" "<< integralCheck[i][cell]<<std::endl;
       integralCheck[i][cell] = 0;
    }
      }
  //////////////////////////////*/

  // std::cout<<matrVectProdN[23]<<" "<<matrVectProdD[23]<<std::endl;

  // TODO I WOULD COMMUNICATE HERE, AT THE END OF ALL THINGS
  // TrilinosWrappers::MPI::Vector
  // Ndummy(matrVectProdN_in.locally_owned_elements(), mpi_communicator);
  // TrilinosWrappers::MPI::Vector
  // Ddummy(matrVectProdD_in.locally_owned_elements(), mpi_communicator); Ndummy
  // = matrVectProdN; Ddummy = matrVectProdD; matrVectProdN_in.add( Ndummy);
  // matrVectProdD_in.add( Ddummy);

  pcout << "...done computing multipole matrix-vector products" << std::endl;
}

// this method computes the preconditioner needed for the GMRES:
// to do it, it needs to receive the alpha vector from the bem_problem
// class, along with the constraint matrix of the bem problem

template <int dim>
TrilinosWrappers::PreconditionILU &
BEMFMA<dim>::FMA_preconditioner(
  const TrilinosWrappers::MPI::Vector &alpha,
  AffineConstraints<double> &          c) // TO BE CHANGED!!!
{
  Teuchos::TimeMonitor LocalTimer(*PrecondTime);
  // the final preconditioner (with constraints) has a slightly different
  // sparsity pattern with respect to the non constrained one. we must here
  // initialize such sparsity pattern
  // final_prec_sparsity_pattern.reinit(alpha.vector_partitioner(),(types::global_dof_index)
  // 125*fma_dh->get_fe().dofs_per_cell);
  final_prec_sparsity_pattern.reinit(this_cpu_set,
                                     mpi_communicator,
                                     (types::global_dof_index)125 *
                                       fma_dh->get_fe().dofs_per_cell);

  // final_prec_sparsity_pattern.reinit(fma_dh->n_dofs(),fma_dh->n_dofs(),125*fma_dh->get_fe().dofs_per_cell);

  // As before we will use the captures to simplify the calls of the worker
  // copier mechanism. For this reason we build an empty scratch structure.
  struct PrecScratch
  {};

  // We need to fill a vector that memorise the entries of the row associated
  // with the index we are treating.
  struct PrecCopy
  {
    PrecCopy()
    {
      row = numbers::invalid_unsigned_int;
      sparsity_row.resize(0);
    };
    PrecCopy(const PrecCopy &in_copy)
    {
      row          = in_copy.row;
      sparsity_row = in_copy.sparsity_row;
    };
    types::global_dof_index              row;
    std::vector<types::global_dof_index> sparsity_row;
  };
  auto f_worker_prec =
    [this, &c](types::global_dof_index i, PrecScratch &, PrecCopy &copy_data) {
      copy_data.sparsity_row.resize(0);
      // unsigned int i = *index_it;
      if (this->this_cpu_set.is_element(i))
        {
          copy_data.row = i;
          if (c.is_constrained(i))
            {
              // cout<<i<<"  (c):"<<endl;
              // constrained nodes entries are taken from the bem problem
              // constraint matrix
              copy_data.sparsity_row.push_back(i);
              const std::vector<std::pair<types::global_dof_index, double>>
                *entries = c.get_constraint_entries(i);
              for (types::global_dof_index j = 0; j < entries->size(); ++j)
                copy_data.sparsity_row.push_back((*entries)[j].first);
              // std::cout<<"demonio proco";
            }
          else
            {
              // cout<<i<<"  (nc): ";
              // other nodes entries are taken from the unconstrained
              // preconditioner matrix
              for (unsigned int j = 0; j < fma_dh->n_dofs(); ++j)
                {
                  if (this->init_prec_sparsity_pattern.exists(i, j))
                    {
                      copy_data.sparsity_row.push_back(j);
                      // std::cout<<"boia"<<j<<" ";
                    }
                }
              // cout<<endl;
            }
          // std::cout<<copy_data.sparsity_row.size()<<std::endl;
        }
    };

  // We only need a for cycle to add the indices to the sparsity pattern.
  auto f_copier_prec = [this](const PrecCopy &copy_data) {
    // std::cout<<"COPY"<<std::endl;
    if (this->this_cpu_set.is_element(copy_data.row))
      {
        // std::cout<<copy_data.row<<"
        // "<<copy_data.sparsity_row.size()<<std::endl;
        for (types::global_dof_index i = 0; i < copy_data.sparsity_row.size();
             ++i)
          {
            this->final_prec_sparsity_pattern.add(copy_data.row,
                                                  copy_data.sparsity_row[i]);
            // std::cout<<copy_data.row<<" "<<i<<std::endl;
          }
      }
  };

  PrecCopy    foo_copy;
  PrecScratch foo_scratch;

  // The following Workstream replaces a for cycle on all dofs to check all the
  // constraints.
  WorkStream::run(
    0, fma_dh->n_dofs(), f_worker_prec, f_copier_prec, foo_scratch, foo_copy);

  final_prec_sparsity_pattern.compress();
  std::cout << final_prec_sparsity_pattern.n_nonzero_elements() << std::endl;
  final_preconditioner.reinit(final_prec_sparsity_pattern);

  // std::cout<<"ok prec sparsity pattern"<<std::endl;
  // now we assemble the final preconditioner matrix: the loop works
  // exactly like the previous one

  // We need a worker function that fills the final sparisty pattern once its
  // saprsity pattern has been set up. In this case no race condition occurs in
  // the worker so we can let it copy in the global memory.
  auto f_sparsity_filler_tbb = [this,
                                &c](blocked_range<types::global_dof_index> r) {
    for (types::global_dof_index i = r.begin(); i < r.end(); ++i)
      {
        if (this_cpu_set.is_element(i))
          {
            if (c.is_constrained(i))
              {
                final_preconditioner.set(i, i, 1);
                // pcout<<i<<" "<<i<<"  **
                // "<<final_preconditioner(i,i)<<std::endl;
                // constrainednodes entries are taken from the bem problem
                // constraint matrix
                const std::vector<std::pair<types::global_dof_index, double>>
                  *entries = c.get_constraint_entries(i);
                for (types::global_dof_index j = 0; j < entries->size(); ++j)
                  {
                    final_preconditioner.set(i,
                                             (*entries)[j].first,
                                             (*entries)[j].second);
                    // pcout<<i<<" "<<(*entries)[j].first<<"  *
                    // "<<(*entries)[j].second<<std::endl;
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
                        final_preconditioner.set(i,
                                                 j,
                                                 init_preconditioner(i, j));
                        // foo_fma->pcout<<i<<" "<<j<<"
                        // "<<foo_fma->init_preconditioner(i,j)<<std::endl;
                      }
                  }
              }
          }
      }
  };

  parallel_for(blocked_range<types::global_dof_index>(0,
                                                      fma_dh->n_dofs(),
                                                      tbb_granularity),
               f_sparsity_filler_tbb);
  // auto f_sparsity_filler = [] (types::global_dof_index i,
  // TrilinosWrappers::SparseMatrix &final_preconditioner, const
  // AffineConstraints<double> &c, const BEMFMA<dim> *foo_fma)
  // {
  //   if (foo_fma->this_cpu_set.is_element(i))
  //     {
  //       if (c.is_constrained(i))
  //         {
  //           final_preconditioner.set(i,i,1);
  //           //pcout<<i<<" "<<i<<"  **
  //           "<<final_preconditioner(i,i)<<std::endl;
  //           // constrainednodes entries are taken from the bem problem
  //           constraint matrix const std::vector< std::pair <
  //           types::global_dof_index, double > > *entries =
  //           c.get_constraint_entries (i); for (types::global_dof_index j=0;
  //           j< entries->size(); ++j)
  //             {
  //               final_preconditioner.set(i,(*entries)[j].first,(*entries)[j].second);
  //               //pcout<<i<<" "<<(*entries)[j].first<<"  *
  //               "<<(*entries)[j].second<<std::endl;
  //             }
  //         }
  //       else
  //         {
  //           // other nodes entries are taken from the unconstrained
  //           preconditioner matrix for (unsigned int j=0;
  //           j<foo_fma->fma_dh->n_dofs(); ++j)
  //             {
  //               // QUI CHECK SU NEUMANN - DIRICHLET PER METTERE A POSTO,
  //               tanto lui già conosce le matrici. if
  //               (foo_fma->init_prec_sparsity_pattern.exists(i,j))
  //                 {
  //                   final_preconditioner.set(i,j,foo_fma->init_preconditioner(i,j));
  //                   // foo_fma->pcout<<i<<" "<<j<<"
  //                   "<<foo_fma->init_preconditioner(i,j)<<std::endl;
  //                 }
  //             }
  //
  //         }
  //
  //     }
  // };
  //
  // // We use the ThreadGroup function to let different parallel thread fill
  // the preconditioner. In this case
  // // we allow for a greater parallelism ensuring that no race condition
  // occurs. Since such a strategy does
  // // not allow for the use of lambda function we must cast them as standard
  // function. For this reason the capture
  // // must be empty.
  // Threads::TaskGroup<> prec_filler;
  // for (types::global_dof_index ii = 0; ii <  fma_dh->n_dofs() ; ii++)
  //   prec_filler += Threads::new_task ( static_cast<void
  //   (*)(types::global_dof_index, TrilinosWrappers::SparseMatrix &, const
  //   AffineConstraints<double> &, const BEMFMA<dim> *)> (f_sparsity_filler),
  //   ii, final_preconditioner, c, this);
  // prec_filler.join_all();

  // The compress operation makes all the vectors on different processors
  // compliant.
  final_preconditioner.compress(VectorOperation::insert);

  // In order to add alpha we can again use the parallel_for strategy.
  auto f_alpha_adder_tbb = [this, &c, &alpha](
                             blocked_range<types::global_dof_index> r) {
    for (types::global_dof_index i = r.begin(); i < r.end(); ++i)
      {
        if (this_cpu_set.is_element(i))
          {
            if ((*(dirichlet_nodes))(i) == 0 && !(c.is_constrained(i)))
              {
                final_preconditioner.add(i, i, alpha(i));
                // pcout<<i<<" "<<i<<" "<<final_preconditioner(i,i)<<std::endl;
              }
            else // this is just to avoid a deadlock. we need a better strategy
              {
                final_preconditioner.add(i, i, 0);
              }
          }
      }
  };

  parallel_for(blocked_range<types::global_dof_index>(0,
                                                      fma_dh->n_dofs(),
                                                      tbb_granularity),
               f_alpha_adder_tbb);
  // auto f_alpha_adder = [] (types::global_dof_index i,
  // TrilinosWrappers::SparseMatrix &final_preconditioner, const
  // AffineConstraints<double> &c, const TrilinosWrappers::MPI::Vector &alpha,
  // const BEMFMA<dim> *foo_fma)
  // {
  //   if (foo_fma->this_cpu_set.is_element(i))
  //     {
  //       if ( (*(foo_fma->dirichlet_nodes))(i) == 0 && !(c.is_constrained(i)))
  //         {
  //           final_preconditioner.add(i,i,alpha(i));
  //           //pcout<<i<<" "<<i<<" "<<final_preconditioner(i,i)<<std::endl;
  //         }
  //       else // this is just to avoid a deadlock. we need a better strategy
  //         {
  //           final_preconditioner.add(i,i,0);
  //         }
  //     }
  //
  // };
  //
  // Threads::TaskGroup<> alpha_adder;
  // for (types::global_dof_index ii = 0; ii <  fma_dh->n_dofs() ; ii++)
  //   alpha_adder += Threads::new_task ( static_cast<void
  //   (*)(types::global_dof_index, TrilinosWrappers::SparseMatrix &, const
  //   AffineConstraints<double> &, const TrilinosWrappers::MPI::Vector &, const
  //   BEMFMA<dim> *)> (f_alpha_adder), ii, final_preconditioner, c, alpha,
  //   this);
  // alpha_adder.join_all();
  final_preconditioner.compress(VectorOperation::add);
  final_preconditioner.compress(VectorOperation::insert);

  // for (unsigned int i=0; i < fma_dh->n_dofs(); i++)
  //   {
  //     if (this_cpu_set.is_element(i))
  //       {
  //         if (c.is_constrained(i))
  //           {
  //             final_preconditioner.set(i,i,1);
  //             //pcout<<i<<" "<<i<<"  **
  //             "<<final_preconditioner(i,i)<<std::endl;
  //             // constrainednodes entries are taken from the bem problem
  //             constraint matrix const std::vector< std::pair <
  //             types::global_dof_index, double > > *entries =
  //             c.get_constraint_entries (i); for (unsigned int j=0; j<
  //             entries->size(); ++j)
  //               {
  //                 final_preconditioner.set(i,(*entries)[j].first,(*entries)[j].second);
  //                 //pcout<<i<<" "<<(*entries)[j].first<<"  *
  //                 "<<(*entries)[j].second<<std::endl;
  //               }
  //           }
  //         else
  //           {
  //             // other nodes entries are taken from the unconstrained
  //             preconditioner matrix for (unsigned int j=0;
  //             j<fma_dh->n_dofs(); ++j)
  //               {
  //                 // QUI CHECK SU NEUMANN - DIRICHLET PER METTERE A POSTO,
  //                 tanto lui già conosce le matrici. if
  //                 (init_prec_sparsity_pattern.exists(i,j))
  //                   {
  //                     final_preconditioner.set(i,j,init_preconditioner(i,j));
  //                     //pcout<<i<<" "<<j<<"
  //                     "<<init_preconditioner(i,j)<<std::endl;
  //                   }
  //               }
  //
  //           }
  //       }
  //   }
  // // std::cout<<"now alpha"<<std::endl;
  // // finally, we have to add the alpha values on the diagonal, whenever
  // dealing with a
  // // neumann (in such nodes the potential phi is an unknown) and non
  // constrained node
  //
  // for (unsigned int i=0; i < fma_dh->n_dofs(); i++)
  //   {
  //     if (this_cpu_set.is_element(i))
  //       {
  //         if ( (*dirichlet_nodes)(i) == 0 && !(c.is_constrained(i)))
  //           {
  //             final_preconditioner.add(i,i,alpha(i));
  //             //pcout<<i<<" "<<i<<" "<<final_preconditioner(i,i)<<std::endl;
  //           }
  //         else // this is just to avoid a deadlock. we need a better strategy
  //           {
  //             final_preconditioner.add(i,i,0);
  //           }
  //       }
  //   }

  // preconditioner.print_formatted(pcout,4,true,0," 0 ",1.);

  // Finally we can initialize the ILU final preconditioner.
  preconditioner.initialize(final_preconditioner);

  return preconditioner;
}

template <int dim>
void
BEMFMA<dim>::compute_geometry_cache()
{
  pcout << "Generating geometry cache..." << std::endl;


  FESystem<dim - 1, dim>   gradient_fe(fma_dh->get_fe(), dim);
  DoFHandler<dim - 1, dim> gradient_dh(fma_dh->get_triangulation());

  // double tol = 1e-8;
  std::vector<Point<dim>> support_points(fma_dh->n_dofs());

  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*fma_mapping,
                                                     *fma_dh,
                                                     support_points);

  // for (unsigned int i=0; i<fma_dh->n_dofs(); ++i)
  //   {
  //     for (unsigned int j=0; j<fma_dh->n_dofs(); ++j)
  //       {
  //         //pcout<<"i "<<i<<" ("<<support_points[i]<<")  j "<<j<<"
  //         ("<<support_points[j]<<")  distance
  //         "<<support_points[i].distance(support_points[j])<<std::endl;
  //       }
  //
  //   }

  // for the gradient dofs finding coupled
  // dofs is a little bit difficult, as the
  // gradient is a vectorial function: usually
  // in such case the dofs are numbered
  // so that for each support point dim
  // consecutive dofs represent each component
  // of the vector field: in this case
  // (and only in this case) the following
  // piece of code works


  cell_it gradient_cell = gradient_dh.begin_active(),
          gradient_endc = gradient_dh.end();

  cell_it cell = fma_dh->begin_active(), endc = fma_dh->end();

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
  elem_to_surr_elems.clear();



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

  // qui viene creata la mappa dei elmenti che circondano ciascun elemento
  for (cell = fma_dh->begin_active(); cell != endc; ++cell)
    {
      cell->get_dof_indices(dofs);
      for (unsigned int j = 0; j < fma_dh->get_fe().dofs_per_cell; ++j)
        {
          std::set<types::global_dof_index> duplicates =
            (*double_nodes_set)[dofs[j]];
          for (std::set<types::global_dof_index>::iterator pos =
                 duplicates.begin();
               pos != duplicates.end();
               pos++)
            {
              std::vector<cell_it> dof_cell_list = dof_to_elems[*pos];
              for (unsigned int k = 0; k < dof_cell_list.size(); ++k)
                elem_to_surr_elems[cell].insert(dof_cell_list[k]);
            }
        }
    }

  // gradient_cell = gradient_dh.begin_active();
  // cell = fma_dh->begin_active();
  //
  // for (; gradient_cell!=gradient_endc; ++cell,++gradient_cell)
  //     {
  //     Assert(cell->index() == gradient_cell->index(), ExcInternalError());
  //
  //     if (cell->material_id() == wall_sur_ID1 ||
  //         cell->material_id() == wall_sur_ID2 ||
  //         cell->material_id() == wall_sur_ID3)
  //  {
  //        // This is a free surface node.
  //        gradient_cell->get_dof_indices(gradient_dofs);
  //        for (unsigned int j=0; j<gradient_fe.dofs_per_cell; ++j)
  //            {
  //      std::set <unsigned int> duplicates =
  //      gradient_double_nodes_set[gradient_dofs[j]]; for (std::set<unsigned
  //      int>::iterator pos = duplicates.begin(); pos !=duplicates.end();
  //      pos++)
  //          {
  //          cell_it duplicate_cell =  gradient_dof_to_elems[*pos][0];
  //   if (duplicate_cell->material_id() == free_sur_ID1 ||
  //                    duplicate_cell->material_id() == free_sur_ID2 ||
  //                    duplicate_cell->material_id() == free_sur_ID3)
  //      {
  //      free_surf_and_boat_nodes.insert(gradient_dofs[j]);
  //      }
  //
  //          }
  //      }
  //            //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"
  //            surfNodes: "<<dirichlet_nodes(dofs[i])<<"  otherNodes:
  //            "<<other_nodes(dofs[i])<<std::endl;
  //        }
  //     }
  //
  // gradient_cell = gradient_dh.begin_active();
  // cell = fma_dh->begin_active();
  //
  // for (; gradient_cell!=gradient_endc; ++cell,++gradient_cell)
  //     {
  //     Assert(cell->index() == gradient_cell->index(), ExcInternalError());
  //
  //     if (gradient_cell->material_id() == wall_sur_ID1)
  //  {
  //        gradient_cell->get_dof_indices(gradient_dofs);
  //        for (unsigned int j=0; j<gradient_fe.dofs_per_cell; ++j)
  //            {
  //      std::set <unsigned int> duplicates =
  //      gradient_double_nodes_set[gradient_dofs[j]]; for (std::set<unsigned
  //      int>::iterator pos = duplicates.begin(); pos !=duplicates.end();
  //      pos++)
  //          {
  //          cell_it duplicate_cell =  gradient_dof_to_elems[*pos][0];
  //   if (
  //                    duplicate_cell->material_id() == wall_sur_ID2 ||
  //                    duplicate_cell->material_id() == wall_sur_ID3)
  //      {
  //      boat_keel_nodes.insert(gradient_dofs[j]);
  //      }
  //
  //          }
  //      }
  //            //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"
  //            surfNodes: "<<dirichlet_nodes(dofs[i])<<"  otherNodes:
  //            "<<other_nodes(dofs[i])<<std::endl;
  //        }
  //     }
  // gradient_cell = gradient_dh.begin_active();
  // cell = fma_dh->begin_active();
  //
  // for (; gradient_cell!=gradient_endc; ++cell,++gradient_cell)
  //     {
  //     Assert(cell->index() == gradient_cell->index(), ExcInternalError());
  //
  //     if (gradient_cell->material_id() == wall_sur_ID2)
  //  {
  //        // This is a free surface node.
  //        gradient_cell->get_dof_indices(gradient_dofs);
  //        for (unsigned int j=0; j<gradient_fe.dofs_per_cell; ++j)
  //            {
  //      std::set <unsigned int> duplicates =
  //      gradient_double_nodes_set[gradient_dofs[j]]; for (std::set<unsigned
  //      int>::iterator pos = duplicates.begin(); pos !=duplicates.end();
  //      pos++)
  //          {
  //          cell_it duplicate_cell =  gradient_dof_to_elems[*pos][0];
  //   if (duplicate_cell->material_id() == wall_sur_ID3)
  //      {
  //      boat_keel_nodes.insert(gradient_dofs[j]);
  //      }
  //
  //          }
  //      }
  //            //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"
  //            surfNodes: "<<dirichlet_nodes(dofs[i])<<"  otherNodes:
  //            "<<other_nodes(dofs[i])<<std::endl;
  //        }
  //     }


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
      // for printout
      // pcout<<"Node "<<i<< "["<<support_points[i]<<"] "<<std::endl;
      for (unsigned int j = 0; j < dim; j++)
        {
          max_coor_value =
            std::max(max_coor_value, std::abs(support_points[i](j)));
        }
    }

  if (blocks.size() > 0)
    {
      for (types::global_dof_index ii = 0; ii < num_blocks; ii++)
        delete blocks[ii];
    }

  types::global_dof_index maxNumBlocks =
    num_octree_levels * fma_dh->get_triangulation().n_active_cells() *
    fe_v.n_quadrature_points;
  // unsigned int maxNumBlocks = 0;
  // for (unsigned int ii = 0; ii < num_octree_levels + 1;  ii++)
  //  {
  //   maxNumBlocks += int(pow(8.,double(ii)));
  //  }

  blocks.clear();
  blocks.reserve(maxNumBlocks);
  blocks.resize(maxNumBlocks);

  types::global_dof_index blocksCount = 0;
  startLevel.resize(num_octree_levels + 1);
  endLevel.resize(num_octree_levels + 1);

  // for (unsigned int j=0; j < num_octree_levels + 1; j++)
  //     parentList[j].clear();
  parentList.clear();
  parentList.resize(num_octree_levels + 1);
  parentList[0].push_back(0);


  childlessList.clear();
  types::global_dof_index numChildless = 0;
  numParent.resize(num_octree_levels + 1);

  // qui di seguito vengono reinizializzate strutture utili al multipolo

  // mappa che associa ad ogni dof un vettore con i blocchi cui essa appartiene
  // per ogni livello
  dof_to_block.clear();

  // mappa che associa ad ogni quad point un vettore con i blocchi cui essa
  // appartiene per ogni livello
  quad_point_to_block.clear();

  // vettore di vettori contenente per ogni livello, gli ids dei blocchi
  // contenenti almeno un dof
  dofs_filled_blocks.clear();

  // vettore di vettori contenente per ogni livello, gli ids dei blocchi
  // contenenti almeno un quad point
  quad_points_filled_blocks.clear();

  quadPoints.clear();
  quadNormals.clear();
  quadShapeFunValues.clear();
  quadJxW.clear();

  dofs_filled_blocks.resize(num_octree_levels + 1);

  quad_points_filled_blocks.resize(num_octree_levels + 1);



  for (unsigned int ii = 0; ii < num_octree_levels + 1; ii++)
    {
      numParent[ii] = 0;
    }



  Point<dim> pMin;
  for (int i = 0; i < dim; i++)
    pMin(i) = -1.1 * max_coor_value;

  // delta e' il lato del kazzo di kubo...
  double delta = 2.2 * max_coor_value;

  OctreeBlock<dim> *block = new OctreeBlock<dim>(0, 0, pMin, delta);

  std::vector<types::global_dof_index> local_dof_indices(
    fma_dh->get_fe().dofs_per_cell);
  cell_it cell = fma_dh->begin_active(), endc = fma_dh->end();
  for (cell = fma_dh->begin_active(); cell != endc; ++cell)
    {
      fe_v.reinit(cell);
      const unsigned int n_q_points = fe_v.n_quadrature_points;
      quadPoints[cell]              = fe_v.get_quadrature_points();
      // quadNormals[cell] = fe_v.get_normal_vectors();
      quadNormals[cell] = fe_v.get_normal_vectors();
      quadJxW[cell].resize(n_q_points);
      quadShapeFunValues[cell].resize(n_q_points);
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          quadJxW[cell][q] = fe_v.JxW(q);
          for (unsigned int j = 0; j < fma_dh->get_fe().dofs_per_cell; ++j)
            quadShapeFunValues[cell][q].push_back(fe_v.shape_value(j, q));
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
      dof_to_block[ii].push_back(0);
    }


  // just for output
  /*for (cell = fma_dh->begin_active(); cell != endc; ++cell)
      {
      std::set<cell_it> surr_elems = elem_to_surr_elems[cell];
      pcout<<std::endl<<"cell "<<cell<<"  surrounded by: ";
      for (typename std::set<cell_it>::iterator pos = surr_elems.begin(); pos
     !=surr_elems.end(); pos++) pcout<<" "<<*pos;
      }*/

  blocks[0]    = block;
  numParent[0] = 1;

  // pcout<<"blocks[0].GetBlockChildrenNum()
  // "<<blocks[0].GetBlockChildrenNum()<<std::endl;

  /*pcout<<std::endl;
  pcout<<blocks[0].GetPMin()<<std::endl;
  pcout<<pMin<<std::endl;
  pcout<<block.GetDelta()<<std::endl;
  pcout<<block.GetBlockNodeList()[0]<<std::endl;
  pcout<<block.GetBlockElementsList()[1]<<std::endl;
  pcout<<delta<<std::endl;
  pcout<<std::endl;//*/

  types::global_dof_index quadPointsInChildless = 0;
  types::global_dof_index nodesInChildless      = 0;

  for (unsigned int level = 1; level < num_octree_levels + 1; level++)

    {
      types::global_dof_index quadPointsCheck = quadPointsInChildless;
      types::global_dof_index nodesCheck      = nodesInChildless;
      delta /= 2.;

      for (types::global_dof_index kk = 0; kk < numParent[level - 1]; kk++)

        {
          types::global_dof_index jj = parentList[level - 1][kk];
          // pcout<<" level "<<level<<"     block "<<jj<<std::endl;
          OctreeBlock<dim> *parent = blocks[jj];
          // pcout<<"parent.GetBlockChildrenNum()
          // "<<parent.GetBlockChildrenNum()<<std::endl; pcout<<" Pmin
          // "<<parent.GetPMin()(0)<<", "<<parent.GetPMin()(1)<<",
          // "<<parent.GetPMin()(2)<<" "<<std::endl; pcout<<" delta
          // "<<parent.GetDelta()<<" "<<std::endl;

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

          std::map<cell_it, std::vector<types::global_dof_index>>
            blockQuadPointsList = parent->GetBlockQuadPointsList();

          std::vector<types::global_dof_index> blockNodeList =
            parent->GetBlockNodeList();

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
                              // pcout<<" Sono in 1 "<<std::endl;
                              children[0]->AddNode(blockNodeList[i]);
                            }
                          else
                            {
                              // pcout<<" Sono in 2 "<<std::endl;
                              children[1]->AddNode(blockNodeList[i]);
                            }
                        }
                      else
                        {
                          if (node(0) <= parent->GetPMin()(0) + delta)
                            {
                              // pcout<<" Sono in 4 "<<std::endl;
                              children[3]->AddNode(blockNodeList[i]);
                            }
                          else
                            {
                              // pcout<<" Sono in 3 "<<std::endl;
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
                              // pcout<<" Sono in 5 "<<std::endl;
                              children[4]->AddNode(blockNodeList[i]);
                            }
                          else
                            {
                              // pcout<<" Sono in 6 "<<std::endl;
                              children[5]->AddNode(blockNodeList[i]);
                            }
                        }
                      else
                        {
                          if (node(0) <= parent->GetPMin()(0) + delta)
                            {
                              // pcout<<" Sono in 8 "<<std::endl;
                              children[7]->AddNode(blockNodeList[i]);
                            }
                          else
                            {
                              // pcout<<" Sono in 7 "<<std::endl;
                              children[6]->AddNode(blockNodeList[i]);
                            }
                        }
                    } // fine assegnazione nodi del padre ai blocchi figli

                } // fine loop nodi del blocco

              typename std::map<cell_it,
                                std::vector<types::global_dof_index>>::iterator
                it;
              for (it = blockQuadPointsList.begin();
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
                                  // pcout<<" Sono in 1 "<<std::endl;
                                  children[0]->AddQuadPoint((*it).first,
                                                            (*it).second[pp]);
                                }
                              else
                                {
                                  // pcout<<" Sono in 2 "<<std::endl;
                                  children[1]->AddQuadPoint((*it).first,
                                                            (*it).second[pp]);
                                }
                            }
                          else
                            {
                              if (quadPoint(0) <= parent->GetPMin()(0) + delta)
                                {
                                  // pcout<<" Sono in 4 "<<std::endl;
                                  children[3]->AddQuadPoint((*it).first,
                                                            (*it).second[pp]);
                                }
                              else
                                {
                                  // pcout<<" Sono in 3 "<<std::endl;
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
                                  // pcout<<" Sono in 5 "<<std::endl;
                                  children[4]->AddQuadPoint((*it).first,
                                                            (*it).second[pp]);
                                }
                              else
                                {
                                  // pcout<<" Sono in 6 "<<std::endl;
                                  children[5]->AddQuadPoint((*it).first,
                                                            (*it).second[pp]);
                                }
                            }
                          else
                            {
                              if (quadPoint(0) <= parent->GetPMin()(0) + delta)
                                {
                                  // pcout<<" Sono in 8 "<<std::endl;
                                  children[7]->AddQuadPoint((*it).first,
                                                            (*it).second[pp]);
                                }
                              else
                                {
                                  // pcout<<" Sono in 7 "<<std::endl;
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
                      blocks[blocksCount] = new OctreeBlock<dim>();
                      blocks[blocksCount]->CopyContent(children[j]);
                      delete children[j];

                      parent->AddChild(blocksCount);
                      std::map<cell_it, std::vector<types::global_dof_index>>
                        blockQuadPointsList =
                          blocks[blocksCount]->GetBlockQuadPointsList();
                      typename std::map<
                        cell_it,
                        std::vector<types::global_dof_index>>::iterator it;
                      for (it = blockQuadPointsList.begin();
                           it != blockQuadPointsList.end();
                           it++)
                        {
                          cell_it cell = (*it).first;
                          for (types::global_dof_index kk = 0;
                               kk < (*it).second.size();
                               kk++)
                            {
                              quad_point_to_block[(*it).first][(*it).second[kk]]
                                .push_back(blocksCount);
                            }
                        }
                      std::vector<types::global_dof_index> blockNodesList =
                        blocks[jj]->GetBlockNodeList();
                      for (types::global_dof_index k = 0;
                           k < blockNodesList.size();
                           k++)
                        dof_to_block[blockNodesList[k]].push_back(jj);
                    }
                  else
                    {
                      delete children[j];
                    }

                } // fine loop sui blocchi figlio appena creati

            } // fine ramo dim = 3 dell'if
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
                          // pcout<<" Sono in 1 "<<std::endl;
                          children[0]->AddNode(blockNodeList[i]);
                        }
                      else
                        {
                          // pcout<<" Sono in 2 "<<std::endl;
                          children[1]->AddNode(blockNodeList[i]);
                        }
                    }
                  else
                    {
                      if (node(0) <= parent->GetPMin()(0) + delta)
                        {
                          // pcout<<" Sono in 4 "<<std::endl;
                          children[3]->AddNode(blockNodeList[i]);
                        }
                      else
                        {
                          // pcout<<" Sono in 3 "<<std::endl;
                          children[2]->AddNode(blockNodeList[i]);
                        }
                    } // fine assegnazione blocchi del padre ai blocchi figli

                } // fine loop nodi del blocco

              typename std::map<cell_it,
                                std::vector<types::global_dof_index>>::iterator
                it;
              for (it = blockQuadPointsList.begin();
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
                              // pcout<<" Sono in 1 "<<std::endl;
                              children[0]->AddQuadPoint((*it).first,
                                                        (*it).second[pp]);
                            }
                          else
                            {
                              // pcout<<" Sono in 2 "<<std::endl;
                              children[1]->AddQuadPoint((*it).first,
                                                        (*it).second[pp]);
                            }
                        }
                      else
                        {
                          if (quadPoint(0) <= parent->GetPMin()(0) + delta)
                            {
                              // pcout<<" Sono in 4 "<<std::endl;
                              children[3]->AddQuadPoint((*it).first,
                                                        (*it).second[pp]);
                            }
                          else
                            {
                              // pcout<<" Sono in 3 "<<std::endl;
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
                      blocks[blocksCount] = new OctreeBlock<dim>();
                      blocks[blocksCount]->CopyContent(children[j]);
                      delete children[j];

                      parent->AddChild(blocksCount);
                      std::map<cell_it, std::vector<types::global_dof_index>>
                        blockQuadPointsList =
                          blocks[blocksCount]->GetBlockQuadPointsList();
                      typename std::map<
                        cell_it,
                        std::vector<types::global_dof_index>>::iterator it;
                      for (it = blockQuadPointsList.begin();
                           it != blockQuadPointsList.end();
                           it++)
                        {
                          cell_it cell = (*it).first;
                          for (types::global_dof_index kk = 0;
                               kk < (*it).second.size();
                               kk++)
                            {
                              quad_point_to_block[(*it).first][(*it).second[kk]]
                                .push_back(blocksCount);
                            }
                        }
                      std::vector<types::global_dof_index> blockNodesList =
                        blocks[jj]->GetBlockNodeList();
                      for (types::global_dof_index k = 0;
                           k < blockNodesList.size();
                           k++)
                        dof_to_block[blockNodesList[k]].push_back(jj);
                    }
                  else
                    {
                      delete children[j];
                    }

                } // fine loop sui blocchi figlio appena creati


            } // fine ramo dim == 2 dell'if

        } // fine loop blocchi livello precedente


      // double elemCheck = numChildless;

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
          std::vector<types::global_dof_index> nodesId =
            blocks[jj]->GetBlockNodeList();
          double blockNumNodes = 0.0;

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
            blockNumQuadPoints += (int)(*it).second.size();
          // pcout<<"Level "<<level<<" Block "<<jj<<"  nodes "<<blockNumNodes<<"
          // + quadPoints "<<blockNumQuadPoints<<std::endl;

          quadPointsCheck += blockNumQuadPoints;
          nodesCheck += blockNumNodes;
          // here we decide if a block is to be placed in the parent
          // or childless list
          // if (blockNumNodes + blockNumQuadPoints - numDoubleNodes < 2)
          if (round(blockNumNodes) <= max_num_nodes_per_block)
            {
              numChildless += 1;
              childlessList.push_back(jj);
              quadPointsInChildless += blockNumQuadPoints;
              nodesInChildless += blockNumNodes;


              // if a block is childless, we must assign now the nodes and quad
              // points that belong to it for all the next levels

              for (types::global_dof_index kk = 0; kk < nodesId.size(); kk++)
                for (unsigned int j = level + 1; j < num_octree_levels + 1; j++)
                  dof_to_block[nodesId[kk]].push_back(jj);

              for (it = blockQuadPointsList.begin();
                   it != blockQuadPointsList.end();
                   it++)
                for (types::global_dof_index i = 0; i < (*it).second.size();
                     i++)
                  for (unsigned int j = level + 1; j < num_octree_levels + 1;
                       j++)
                    quad_point_to_block[(*it).first][(*it).second[i]].push_back(
                      jj);
            }
          else
            {
              numParent[level] += 1;
              parentList[level].push_back(jj);
            }

          // let's update the list of node filled block
          if (blockNumNodes > 0)
            dofs_filled_blocks[level].push_back(jj);

          // let's update the list of quad point filled block
          if (blockNumQuadPoints > 0)
            quad_points_filled_blocks[level].push_back(jj);
          // elemCheck += blockNumNodes + blockNumQuadPoints;
        }


      pcout << " Total nodes at level " << level << " of " << num_octree_levels
            << " are " << nodesCheck << std::endl;
      pcout << " Total quad points at level " << level << " of "
            << num_octree_levels << " are " << quadPointsCheck << std::endl;
      pcout << " Blocks at level " << level << " of " << num_octree_levels
            << " are " << endLevel[level] - endLevel[level - 1] << std::endl;
      pcout << " Total blocks at level " << level << " of " << num_octree_levels
            << " are " << endLevel[level] + 1 << std::endl;
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

  // just for output
  /*for (cell = fma_dh->begin_active(); cell != endc; ++cell)
      {
      unsigned int levelCheck = elem_to_blocks[cell].size();
      pcout<<std::endl<<"Elem "<<cell<<" is in the "<<levelCheck<<" blocks: ";
      for (unsigned int zz = 0; zz < levelCheck; zz++)
           pcout<<elem_to_blocks[cell][zz]<<" ";
      }*/

  /*for (cell_it cell = fma_dh->begin_active(); cell != endc; ++cell)
      for (unsigned int j=0; j < quadPoints[cell].size(); j++)
          pcout<<"Cell "<<cell<<"  QP "<<j<<"  of "<<quadPoints[cell].size()<<":
     "<<quadPoints[cell][j]<<std::endl;//*/

  /*for (cell_it cell = fma_dh->begin_active(); cell != endc; ++cell)
      for (unsigned int j=0; j < quad_point_to_block[cell].size(); j++)
          {
    pcout<<"Cell "<<cell<<"  QP "<<j<<"  of
    "<<quad_point_to_block[cell].size()<<": "; for (unsigned int i=0; i <
    quad_point_to_block[cell][j].size(); i++)
              pcout<<quad_point_to_block[cell][j][i]<<" ";
    pcout<<std::endl;
    }  //*/


  /*for (unsigned int i=0; i < fma_dh->n_dofs(); i++)
      {
      pcout<<"Node "<<i<<"  doubles: ";
      std::set <unsigned int> doubleNodes = double_nodes_set[i];
      for (std::set<unsigned int>::iterator pos = doubleNodes.begin(); pos !=
    doubleNodes.end(); pos++)
          {
    pcout<<*pos<<"( ";
    for (unsigned int j=0; j < dof_to_elems[*pos].size(); j++)
        pcout<<" "<<dof_to_elems[*pos][j];
    pcout<<") ";
          }
      pcout<<std::endl;
      } //*/



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
            Center1(iii) = delta1 / 2.;
          Point<dim> PMax1 = 2. * Center1;
          PMax1 += PMin1;
          Center1 += PMin1;
          types::global_dof_index           parentId = block1->GetParentId();
          std::set<types::global_dof_index> parentNNeighs =
            blocks[parentId]->GetNearNeighs(0);

          // the nearest neighbors are searched among the father's nearest
          // neighbors children
          for (std::set<types::global_dof_index>::iterator pos =
                 parentNNeighs.begin();
               pos != parentNNeighs.end();
               pos++)

            {
              if (blocks[*pos]->GetBlockChildrenNum() ==
                  0) // if a parent's near neigh is childless, he can be a near
                     // neigh: let's check
                {
                  types::global_dof_index block2Id = *pos;
                  OctreeBlock<dim> *      block2   = blocks[block2Id];
                  double                  delta2   = block2->GetDelta();
                  Point<dim>              PMin2    = block2->GetPMin();
                  Point<dim>              Center2;
                  for (unsigned int iii = 0; iii < dim; iii++)
                    Center2(iii) = delta2 / 2.;
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
                                  // pcout<<" *"<<block2Id;
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
                                  // pcout<<" *"<<block2Id;
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
                                  // pcout<<" *"<<block2Id;
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
                              // pcout<<block2Id<<" ";
                            }
                        }

                      if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) ||
                          (fabs(PMax1(1) - PMin2(1)) <= TOLL))
                        {
                          if ((PMin1(0) - TOLL <= PMax2(0)) &&
                              (PMax1(0) + TOLL >= PMin2(0)))
                            {
                              block1->AddNearNeigh(0, block2Id);
                              // pcout<<block2Id<<" ";
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
                    Center2(iii) = delta2 / 2.;
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
                                  // pcout<<" "<<block2Id;
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
                                  // pcout<<" "<<block2Id;
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
                                  // pcout<<" "<<block2Id;
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
                              // pcout<<block2Id<<" ";
                            }
                        }

                      if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) ||
                          (fabs(PMax1(1) - PMin2(1)) <= TOLL))
                        {
                          if ((PMin1(0) - TOLL <= PMax2(0)) &&
                              (PMax1(0) + TOLL >= PMin2(0)))
                            {
                              block1->AddNearNeigh(0, block2Id);
                              // pcout<<block2Id<<" ";
                            }
                        }

                    } // fine caso dim == 2

                } // fine loop sui figli di un nearest neighbor del padre


            } // fine loop sui nearest neighbors del padre



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
                  std::set<types::global_dof_index> upperLevelNNeighs =
                    block1->GetNearNeighs(subLevel - 1);
                  for (std::set<types::global_dof_index>::iterator pos =
                         upperLevelNNeighs.begin();
                       pos != upperLevelNNeighs.end();
                       pos++)

                    {
                      if (blocks[*pos]->GetBlockChildrenNum() ==
                          0) // if nearneigh is childless, it will stay a near
                             // neigh
                        block1->AddNearNeigh(subLevel, *pos);

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
                            Center2(iii) = delta2 / 2.;
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
                                          // pcout<<block2Id<<" ";
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
                                          // pcout<<block2Id<<" ";
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
                                          // pcout<<block2Id<<" ";
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
                                      // pcout<<block2Id<<" ";
                                    }
                                }

                              if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) ||
                                  (fabs(PMax1(1) - PMin2(1)) <= TOLL))
                                {
                                  if ((PMin1(0) - TOLL <= PMax2(0)) &&
                                      (PMax1(0) + TOLL >= PMin2(0)))
                                    {
                                      block1->AddNearNeigh(subLevel, block2Id);
                                      // pcout<<block2Id<<" ";
                                    }
                                }

                            } // fine caso dim == 2

                        } // fine loop sui figli di ciascun nearest neighbor del
                          // blocco childless


                    } // fine loop sui nearest neighbors del blocco childless


                } // fine loop sui subLevels (da quello del blocco childless
                  // all'ultimo)


            } // fine if (il blocco e' childless?)



        } // fine loop sui blocchi di un livello

    } // fine loop sui livelli

  // for printout

  /*pcout<<std::endl;
  pcout<<"-------------------------------- "<<std::endl;
  pcout<<"-------------------------------- "<<std::endl;
  pcout<<std::endl;

  //for printout
  for (cell=fma_dh->begin_active();cell!=endc; cell++)
      {
      pcout<<std::endl;
      pcout<<"-------------------------------- "<<std::endl;
      pcout<<"Cell "<<cell<<"  elementPlot(";
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int j = 0; j<local_dof_indices.size(); j++)
          pcout<<"["<<support_points[local_dof_indices[j]]<<"],";
          pcout<<"'r')";
      }

  pcout<<std::endl;
  pcout<<"-------------------------------- "<<std::endl;
  pcout<<"-------------------------------- "<<std::endl;
  pcout<<std::endl;
  pcout<<std::endl;//*/


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
          // pcout<<"??Out "<<jj<<" "<<block1.GetNearNeighSize()<<std::endl;
          for (unsigned int subLevel = 0;
               subLevel < block1->NumNearNeighLevels();
               subLevel++)
            {
              std::set<types::global_dof_index> NNList =
                block1->GetNearNeighs(subLevel);

              for (std::set<types::global_dof_index>::iterator pos1 =
                     NNList.begin();
                   pos1 != NNList.end();
                   pos1++) // loop over blocks in NN list and get their NNs
                {
                  block1->AddBlockToIntList(subLevel, *pos1);
                }
              // pcout<<std::endl<<"Sublevel "<<subLevel<<"
              // elem("<<block1.GetBlockElementsList()[0]<<") NearNeighs: ";
              std::vector<types::global_dof_index> nodeIds =
                block1->GetBlockNodeList();
              // pcout<<std::endl<<"Level "<<level<<"  Block1: "<<kk<<"
              // NumElements: "<<block1.GetBlockElementsNum()<<"
              // "<<elemIds.size()<<std::endl; pcout<<"Nearest Neighbors
              // Found:"<<std::endl;
              for (types::global_dof_index pp = 0; pp < nodeIds.size(); pp++)
                {
                  // pcout<<"Node "<<nodeIds[pp]<<std::endl;
                  std::set<types::global_dof_index> doubleNodes =
                    (*double_nodes_set)[nodeIds[pp]];
                  for (std::set<types::global_dof_index>::iterator pos =
                         doubleNodes.begin();
                       pos != doubleNodes.end();
                       pos++)
                    {
                      // pcout<<"Node Double"<<*pos<<std::endl;
                      // std::vector<cell_it > surrCellIds = dof_to_elems[*pos];
                      for (types::global_dof_index k = 0;
                           k < dof_to_elems[*pos].size();
                           k++)
                        {
                          cell_it cell = dof_to_elems[*pos][k];
                          // pcout<<cell<<std::endl;
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
               subLevel++) // for each block, loop over all sublevels in his NN
                           // list (to account for childless blocks)
            {
              // now use intList to compute nonIntList
              std::set<types::global_dof_index> intList =
                block1->GetIntList(subLevel);
              std::set<types::global_dof_index>
                parentIntList; // intList at the  previous level
              if (subLevel ==
                  0) // if a block is childless we get its intList at the
                     // previous level, otherwise we get its parent's intList
                parentIntList = blocks[block1->GetParentId()]->GetIntList(0);
              else
                parentIntList = block1->GetIntList(subLevel - 1);

              for (std::set<types::global_dof_index>::iterator pos1 =
                     parentIntList.begin();
                   pos1 != parentIntList.end();
                   pos1++) // loop over blocks in parentIntList
                {
                  OctreeBlock<dim> *block2 = blocks[*pos1];
                  if (block2->GetBlockChildrenNum() ==
                      0) // if blocks in parentIntList are childless, don't look
                         // for their children, but see if they are in
                         // nonIntList
                    {
                      if (intList.count(*pos1) ==
                          0) // if these blocks are not in intList
                        block1->AddBlockToNonIntList(
                          subLevel, *pos1); // then they go in nonIntList
                    }
                  else // if blocks in parentIntList are not childless, do the
                       // same test on all their children
                    {
                      for (types::global_dof_index kk = 0;
                           kk < block2->GetBlockChildrenNum();
                           kk++) // loop over children of blocks in
                                 // parentIntList
                        {
                          // pcout<<"Sublevel "<<subLevel<<" Block1 "<<jj<<"
                          // Block2 "<<*pos1<<" child(kk)
                          // "<<block2.GetChildId(kk)<<std::endl;
                          if (intList.count(block2->GetChildId(kk)) ==
                              0) // if these blocks are not in intList
                            block1->AddBlockToNonIntList(
                              subLevel,
                              block2->GetChildId(
                                kk)); // then they go in nonIntList
                        } // end loop over children of blocks in parentIntList
                    }
                } // loop over blocks in parentIntList

            } // end loop over subLevels of each block's intList

          //      for printout

          /*
                pcout<<std::endl;
                                  pcout<<"--------------------------------
          "<<std::endl; pcout<<"Block "<<jj<<": Parent "<<block1->GetParentId();
                pcout<<"
          cubePlot(["<<block1->GetPMin()<<"],"<<block1->GetDelta()<<",'m')"<<std::endl;
                pcout<<"Quad points of block "<<jj<<": "<<std::endl;
                std::map <cell_it,std::vector<types::global_dof_index> >
          quadPointsList = block1->GetBlockQuadPointsList(); typename std::map
          <cell_it, std::vector<types::global_dof_index> >::iterator it; for (it
          = quadPointsList.begin(); it != quadPointsList.end(); it++)
                                {
                    pcout<<(*it).first<<"( ";
                    for (unsigned int zz = 0; zz < (*it).second.size();  zz++)
                        pcout<<(*it).second[zz]<<" ";
                    pcout<<")  ";
                    }
                pcout<<std::endl;

                pcout<<"Nodes of block "<<jj<<" :"<<std::endl;
                std::vector <unsigned int> nodeList =
          block1->GetBlockNodeList(); for (unsigned int zz = 0; zz <
          block1->GetBlockNodesNum();  zz++)
                  {
                  pcout<<nodeList.at(zz)<<" ";
                  }
                pcout<<std::endl;

                for (unsigned int subLevel = 0; subLevel <
          block1->GetNearNeighSize();  subLevel++)
                  {
                  std::set <unsigned int> NNList =
          block1->GetNearNeighs(subLevel); pcout<<"NearNeigh for block "<<jj<<"
          at level "<<level+subLevel<<":"<<std::endl; for (std::set <unsigned
          int>::iterator pos1 = NNList.begin(); pos1 != NNList.end();  pos1++)
                    {
                    pcout<<*pos1<<" ";
                    }
                  pcout<<std::endl;

                  std::set <unsigned int> intList =
          block1->GetIntList(subLevel); pcout<<"IntList for block "<<jj<<" at
          level "<<level+subLevel<<":"<<std::endl; for (std::set <unsigned
          int>::iterator pos1 = intList.begin(); pos1 != intList.end();  pos1++)
                    {
                    pcout<<*pos1<<" ";
                    }
                  pcout<<std::endl;

                  std::set <unsigned int> nonIntList =
          block1->GetNonIntList(subLevel); pcout<<"NonIntList for block "<<jj<<"
          at level "<<level+subLevel<<":"<<std::endl; for (std::set <unsigned
          int>::iterator pos1 = nonIntList.begin(); pos1 != nonIntList.end();
          pos1++)
                    {
                    pcout<<*pos1<<" ";
                    }
                  pcout<<std::endl;
                  }
          //*/



        } // end loop over blocks of a level

    } // end loop over levels

  // if (blocks.size() > 0)
  //   {
  //   for (unsigned int ii = 0; ii < num_blocks;  ii++)
  //        delete blocks[ii];
  //   }

  //    integralCheck.clear();
  // for (unsigned int i = 0; i < fma_dh->n_dofs(); i++)
  // {
  //     for (cell_it cell = fma_dh->begin_active(); cell != endc; ++cell)
  //     {
  //         //pcout<<i<<" "<<cell<<" "<<integralCheck[i][cell]<<std::endl;
  //         integralCheck[i][cell] = 0;
  //     }
  // }



  pcout << "Done computing proximity lists for blocks" << std::endl;
} // end method for octree blocking generation

template class BEMFMA<2>;
template class BEMFMA<3>;
