//----------------------------  step-34.cc  ---------------------------
//    $Id: step-34.cc 18734 2009-04-25 13:36:48Z heltai $
//    Version: $Name$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//    Authors: Luca Heltai, Cataldo Manigrasso, Andrea Mola
//
//----------------------------  step-34.cc  ---------------------------
#define TOLL 0.001
#define MAXELEMENTSPERBLOCK 1

#include "../include/bem_fma.h"
#include "../include/laplace_kernel.h"
#include "utilities.h"


#include "Teuchos_TimeMonitor.hpp"

using Teuchos::Time;
using Teuchos::TimeMonitor;
using Teuchos::RCP;

RCP<Time> MatrVec = TimeMonitor::getNewTimer("Multipole MatrVec Products Time");
RCP<Time> MultGen = TimeMonitor::getNewTimer("Multipole Generation Time");
RCP<Time> MultInt = TimeMonitor::getNewTimer("Multipole Integral Time");
RCP<Time> ListCreat = TimeMonitor::getNewTimer("Octree Generation Time");
RCP<Time> DirInt = TimeMonitor::getNewTimer("Direct Integral Time");
RCP<Time> PrecondTime = TimeMonitor::getNewTimer("FMA_preconditioner Time");

template <int dim>
BEMFMA<dim>::BEMFMA(MPI_Comm mpi_commy)
  :
  singular_quadrature_order(5),//TO BE CHANGED WITH A PARSER
  mpi_communicator (mpi_commy),
  n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
  this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator)),
  pcout(std::cout,
        (this_mpi_process
         == 0))
{}


template <int dim>
BEMFMA<dim>::~BEMFMA()
{

  if (blocks.size() > 0)
    {
      for (unsigned int ii = 0; ii < num_blocks;  ii++)
        delete blocks[ii];
    }

  // fma_dh = SP();
  // fma_fe = SP();
  // dirichlet_nodes = SP();//per quadratura singolare e octree generator
  // double_nodes_set = SP();//da passare al metodo che fa il precondizionatore
  // fma_mapping = SP();


}

template <int dim>
void BEMFMA<dim>::init_fma(const DoFHandler<dim-1,dim> &input_dh,
                           const std::vector<std::set<unsigned int> > &db_in,
                           const TrilinosWrappers::MPI::Vector &input_sn,
                           const Mapping<dim-1,dim> &input_mapping)
{
  fma_dh = &input_dh;
  fma_fe = &(input_dh.get_fe());
  dirichlet_nodes = new const Vector<double>(input_sn);//per quadratura singolare e octree generator
  this_cpu_set.clear();
  this_cpu_set.set_size(fma_dh->n_dofs());
  this_cpu_set.add_indices(input_sn.locally_owned_elements());// = new const IndexSet(input_sn.locally_owned_elements());
  this_cpu_set.compress();
  double_nodes_set = &db_in;//da passare al metodo che fa il precondizionatore
  fma_mapping = &input_mapping;

}

template <int dim>
void BEMFMA<dim>::declare_parameters (ParameterHandler &prm)
{
  prm.enter_subsection("Octree Params");
  {
    prm.declare_entry("Number of Octree Levels", "10", Patterns::Integer());

    prm.declare_entry("Maximum Number of Collocation Points per Childless Block", "20", Patterns::Integer());
  }
  prm.leave_subsection();

  prm.enter_subsection("FMA Params");
  {
    prm.declare_entry("FMA Truncation Order", "6", Patterns::Integer());
  }
  prm.leave_subsection();

}

template <int dim>
void BEMFMA<dim>::parse_parameters (ParameterHandler &prm)
{
  prm.enter_subsection("Octree Params");
  {
    num_octree_levels = prm.get_integer("Number of Octree Levels");
    max_num_nodes_per_block = prm.get_integer("Maximum Number of Collocation Points per Childless Block");
  }
  prm.leave_subsection();

  prm.enter_subsection("FMA Params");
  {
    trunc_order = prm.get_integer("FMA Truncation Order");
  }
  prm.leave_subsection();

}



template <int dim>
void BEMFMA<dim>::direct_integrals()
{
  pcout<<"Computing direct integrals..."<<std::endl;
  TimeMonitor LocalTimer(*DirInt);
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


  std::vector<QTelles<dim-1> > sing_quadratures;
  for (unsigned int i=0; i<fma_fe->dofs_per_cell; ++i)
    sing_quadratures.push_back
    (QTelles<dim-1>(singular_quadrature_order,
                    fma_fe->get_unit_support_points()[i]));
  const unsigned int dofs_per_cell = fma_fe->dofs_per_cell;

  // vector containing the ids of the dofs
  // of each cell: it will be used to transfer
  // the computed local rows of the matrices
  // into the global matrices

  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // vector to store parts of rows of neumann
  // and dirichlet matrix obtained in local
  // operations

  Vector<double>      local_neumann_matrix_row_i(fma_fe->dofs_per_cell);
  Vector<double>      local_dirichlet_matrix_row_i(fma_fe->dofs_per_cell);


  // Now that we have checked that
  // the number of vertices is equal
  // to the number of degrees of
  // freedom, we construct a vector
  // of support points which will be
  // used in the local integrations:

  std::vector<Point<dim> > support_points(fma_dh->n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>(*fma_mapping, *fma_dh, support_points);


  // After doing so, we can start the
  // integration loop over all cells,
  // where we first initialize the
  // FEValues object and get the
  // values of $\mathbf{\tilde v}$ at
  // the quadrature points (this
  // vector field should be constant,
  // but it doesn't hurt to be more
  // general):


  cell_it
  cell = fma_dh->begin_active(),
  endc = fma_dh->end();

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
  unsigned int preconditioner_band = 125*fma_fe->dofs_per_cell;
  // preconditioner_sparsity_pattern.reinit(sol.vector_partitioner(), (unsigned int) preconditioner_band);
  TrilinosWrappers::MPI::Vector helper(this_cpu_set, mpi_communicator);
  // TODO WHY IT DOES NOT WORK????
  // init_prec_sparsity_pattern.reinit(this_cpu_set.make_trilinos_map(mpi_communicator),preconditioner_band);//,125*fma_fe->dofs_per_cell);
  init_prec_sparsity_pattern.reinit(helper.vector_partitioner(),preconditioner_band);//,125*fma_fe->dofs_per_cell);

  for (unsigned int kk = 0; kk < childlessList.size(); kk++)
    {
      // for each block in the childless
      // list we get the list of nodes and
      // we check if it contains nodes:
      // if no nodes are contained there is
      // nothing to do

      unsigned int blockId = childlessList[kk];

      OctreeBlock<dim> *block1 =  blocks[blockId];

      std::vector <unsigned int> block1Nodes = block1->GetBlockNodeList();

      if  (block1Nodes.size() > 0)
        {

          // if block1 contains nodes,
          // we need to get all the quad points
          // in the intList blocks of block1
          // (such quad points will be used for
          // direct integrals)

          unsigned int intListSubLevs = block1->GetIntListSize();
          const std::set<unsigned int> &block1IntList = block1->GetIntList(intListSubLevs-1);

          // in this set we will put all the
          // dofs of the cell to whom
          // the quad points belong

          std::set<unsigned int> directNodes;

          // start looping on the intList
          // blocks (block2 here)

          for (std::set<unsigned int>::iterator pos = block1IntList.begin(); pos != block1IntList.end(); pos++)
            {
              OctreeBlock<dim> *block2 =  blocks[*pos];
              std::map <cell_it, std::vector<unsigned int> >
              blockQuadPointsList = block2->GetBlockQuadPointsList();

              // get the list of quad points
              // in block2 and loop on it

              typename std::map <cell_it, std::vector<unsigned int> >::iterator it;
              for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
                {
                  // the key of the map (*it.first pointer) is
                  // the cell of the quad point: we will
                  // get its dofs and put them in the set
                  // of direct nodes

                  cell_it cell = (*it).first;//pcout<<cell<<"  end "<<(*blockQuadPointsList.end()).first<<std::endl;
                  cell->get_dof_indices(local_dof_indices);
                  for (unsigned int j = 0; j < dofs_per_cell; j++)
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

          for (unsigned int i = 0; i < block1Nodes.size(); i++)
            if(this_cpu_set.is_element(block1Nodes[i]))
              for (std::set<unsigned int>::iterator pos = directNodes.begin(); pos != directNodes.end(); pos++)
              {
                init_prec_sparsity_pattern.add(block1Nodes[i],*pos);
              }
        }

    }

  // unfortunately, the direct integrals must not be computed only for the
  // quadPoints in the intList: if a bigger block is in the nonIntList of
  // another block, the bound for the multipole expansion application does
  // not hold, and so we must compute direct integrals. Here we scan the
  // nonIntlists of each block at each level to look for bigger blocks and
  // initialize the prec matrices sparsity pattern with the corresponding nodes

  for (unsigned int level = 1; level <  num_octree_levels + 1;  level++) // loop over levels

    {
      // std::vector<unsigned int>
      // dofs_filled_blocks =  dofs_filled_blocks[level];
      unsigned int startBlockLevel =  startLevel[level];

      // we loop over blocks of each level

      for (unsigned int jj = 0; jj < dofs_filled_blocks[level].size();  jj++)
        {
          OctreeBlock<dim> *block1 =  blocks[dofs_filled_blocks[level][jj]];
          const std::vector <unsigned int> &nodesBlk1Ids = block1->GetBlockNodeList();

          // again, no need to perform next operations if block has no nodes

          if  (nodesBlk1Ids.size() > 0)// !!!CHECK, IT SEEMS TO BE USELESS
            {
              // for each block containing nodes, loop over all sublevels in his NN list (this is because if a
              // block remains childless BEFORE the last level, at this point we need to compute
              // all its contributions up to the bottom level)


              for (unsigned int subLevel = 0; subLevel < block1->NumNearNeighLevels();  subLevel++)
                {

                  // in this vectors we are saving the nodes needing direct integrals

                  std::set <unsigned int> directNodes;
                  const std::set <unsigned int> &nonIntList = block1->GetNonIntList(subLevel);

                  // loop over well separated blocks of higher size (level): in this case
                  // we must use direct evaluation: for each block we get the quad points
                  // list
                  for (std::set<unsigned int>::iterator pos = nonIntList.begin(); pos !=nonIntList.lower_bound(startBlockLevel); pos++)
                    {
                      OctreeBlock<dim> *block2 =  blocks[*pos];
                      std::map <cell_it, std::vector<unsigned int> >
                      blockQuadPointsList = block2->GetBlockQuadPointsList();

                      // we loop on the cells of the quad blocks (*it.first pointer)
                      // and put their dofs in the direct list

                      typename std::map <cell_it, std::vector<unsigned int> >::iterator it;
                      for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
                        {
                          cell_it cell = (*it).first;
                          cell->get_dof_indices(local_dof_indices);
                          for (unsigned int j = 0; j < dofs_per_cell; j++)
                            directNodes.insert(local_dof_indices[j]);
                        }
                    } // end loop over blocks of a sublevel of nonIntList

                  // we use the nodes in directList, to create the sparsity pattern

                  for (unsigned int i = 0; i < nodesBlk1Ids.size(); i++)
                  {
                    if(this_cpu_set.is_element(nodesBlk1Ids[i]))
                      for (std::set<unsigned int>::iterator pos = directNodes.begin(); pos != directNodes.end(); pos++)
                        init_prec_sparsity_pattern.add(nodesBlk1Ids[i],*pos);
                  }

                } // end loop over sublevels
            } // end if: is there any node in the block?
        }// end loop over block of a level
    }//end loop over octree levels


  // sparsity pattern is ready and can be compressed; the direct matrices
  // and the preconditioner one are the initialized with the sparsity
  // pattern just computed

  init_prec_sparsity_pattern.compress();
  double filling_percentage = double(init_prec_sparsity_pattern.n_nonzero_elements())/double(fma_dh->n_dofs()*fma_dh->n_dofs())*100.;
  pcout<<init_prec_sparsity_pattern.n_nonzero_elements()<<" Nonzeros out of "<<fma_dh->n_dofs()*fma_dh->n_dofs()<<":  "<<filling_percentage<<"%"<<std::endl;

  prec_neumann_matrix.reinit(init_prec_sparsity_pattern);
  prec_dirichlet_matrix.reinit(init_prec_sparsity_pattern);
  init_preconditioner.reinit(init_prec_sparsity_pattern);

  Point<dim> D;
  double s;

  // here we finally start computing the direct integrals: we
  // first loop among the childless blocks

  for (unsigned int kk = 0; kk <  childlessList.size(); kk++)

    {
      //pcout<<"processing block "<<kk <<"  of  "<<cMesh->GetNumChildlessBlocks()<<std::endl;
      //pcout<<"block "<<cMesh->GetChildlessBlockId(kk) <<"  of  "<<cMesh->GetNumBlocks()<<"  in block list"<<std::endl;

      // this is the Id of the block
      unsigned int blockId =  childlessList[kk];
      // and this is the block pointer
      OctreeBlock<dim> *block1 =  blocks[blockId];
      // we get the block node list
      const std::vector <unsigned int> &block1Nodes = block1->GetBlockNodeList();

      // if a block has no nodes (if it only contains quad points), there is nothing to do
      // if instead there are nodes, we start integrating
      if  (block1Nodes.size() > 0)
        {
          // we first get all the blocks in the intList of the current block (block1)
          // and loop over these blocks, to create a list of ALL the quadrature points that
          // lie in the interaction list blocks: these quad points have to be integrated
          // directly. the list of direct quad points has to be a std::map of std::set of
          // integers, meaning that to each cell, we associate a std::set containing all
          // the direct quad point ids
          unsigned int intListNumLevs = block1->GetIntListSize();
          std::set <unsigned int> block1IntList = block1->GetIntList(intListNumLevs-1);

          std::map <cell_it,std::set<unsigned int> > directQuadPoints;
          for (std::set<unsigned int>::iterator pos = block1IntList.begin(); pos != block1IntList.end(); pos++)
            {
              // now for each block block2 we get the list of quad points
              OctreeBlock<dim> *block2 =  blocks[*pos];
              std::map <cell_it, std::vector<unsigned int> >
              blockQuadPointsList = block2->GetBlockQuadPointsList();

              // we now loop among the cells of the list and for each cell we loop
              // among its quad points, to copy them into the direct quad points list
              typename std::map <cell_it, std::vector<unsigned int> >::iterator it;
              for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
                {
                  for (unsigned int i=0; i<(*it).second.size(); i++)
                    {
                      directQuadPoints[(*it).first].insert((*it).second[i]);

                      /*//////////this is for a check///////////////////
                      for (unsigned int kk=0; kk<block1Nodes.size(); kk++)
                             integralCheck[block1Nodes[kk]][(*it).first] += 1;
                      ///////////////////////////*/
                    }
                }
            }
          // we are now ready to go: for each node, we know which quad points are to be
          // treated directly, and for each node, we will now perform the integral.
          // we then start looping on the nodes of the block
          for (unsigned int i=0; i<block1Nodes.size(); i++)
            {
              unsigned int nodeIndex = block1Nodes[i];
              if(this_cpu_set.is_element(nodeIndex))
              {
                typename std::map <cell_it, std::set<unsigned int> >::iterator it;
                // we loop on the list of quad points to be treated directly
                for (it = directQuadPoints.begin(); it != directQuadPoints.end(); it++)
                  {
                    // the vectors with the local integrals for the cell must first
                    // be zeroed
                    local_neumann_matrix_row_i = 0;
                    local_dirichlet_matrix_row_i = 0;

                    // we get the first entry of the map, i.e. the cell pointer
                    // and we check if the cell contains the current node, to
                    // decide if singular of regular quadrature is to be used
                    cell_it cell = (*it).first;
                    cell->get_dof_indices(local_dof_indices);

                    // we copy the cell quad points in this set
                    std::set<unsigned int> &cellQuadPoints = (*it).second;
                    bool is_singular = false;
                    unsigned int singular_index = numbers::invalid_unsigned_int;

                    for (unsigned int j=0; j<fma_fe->dofs_per_cell; ++j)
                      if ( (*double_nodes_set)[nodeIndex].count(local_dof_indices[j]) > 0)
                        {
                          singular_index = j;
                          is_singular = true;
                          break;
                        }
                    // first case: the current node does not belong to the current cell:
                    // we use regular quadrature
                    if (is_singular == false)
                      {
                        //pcout<<"Node "<<i<<"  Elem "<<cell<<" (Direct) Nodes: ";
                        //for(unsigned int j=0; j<fe.dofs_per_cell; ++j) pcout<<" "<<local_dof_indices[j];
                        //pcout<<std::endl;

                        // we start looping on the quad points of the cell: *pos will be the
                        // index of the quad point
                        for (std::set<unsigned int>::iterator pos=cellQuadPoints.begin(); pos!=cellQuadPoints.end(); pos++)
                          {
                            // here we compute the distance R between the node and the quad point

                            //MAGARI USARE FEVALUES CON IL DOFHANDLER CRETINO DISCONTINUO E IL MAPPING bem_fma
                            const Tensor<1, dim> R =  quadPoints[cell][*pos] - support_points[nodeIndex];
                            LaplaceKernel::kernels(R, D, s);

                            // and here are the integrals for each of the degrees of freedom of the cell: note
                            // how the quadrature values (position, normals, jacobianXweight, shape functions)
                            // are taken from the precomputed ones in ComputationalDomain class
                            for (unsigned int j=0; j<fma_fe->dofs_per_cell; ++j)
                              {
                                local_neumann_matrix_row_i(j) += ( ( D *
                                                                     quadNormals[cell][*pos] ) *
                                                                   quadShapeFunValues[cell][*pos][j] *
                                                                   quadJxW[cell][*pos] );
                                local_dirichlet_matrix_row_i(j) += ( s *
                                                                     quadShapeFunValues[cell][*pos][j] *
                                                                     quadJxW[cell][*pos] );
                                //pcout<<D<<" "<< quadNormals[cell][*pos]<<" ";
                                //pcout<< quadShapeFunValues[cell][*pos][j]<<" ";
                                //pcout<< quadJxW[cell][*pos]<<std::endl;
                              }
                          }
                      } // end if
                    else
                      {
                        // after some checks, we have to create the singular quadrature:
                        // here the quadrature points of the cell will be IGNORED,
                        // and the singular quadrature points are instead used.
                        // the 3d and 2d quadrature rules are different

                        // QUESTO E' IL SOLITO STEP 34, VEDI SE CAMBIARE CON QUELLO NUOVO PER STOKES
                        Assert(singular_index != numbers::invalid_unsigned_int,
                               ExcInternalError());

                        const Quadrature<dim-1> *
                        singular_quadrature
                          = (dim == 2
                             ?
                             dynamic_cast<Quadrature<dim-1>*>(
                               &sing_quadratures[singular_index])
                             :
                             (dim == 3
                              ?
                              dynamic_cast<Quadrature<dim-1>*>(
                                &sing_quadratures[singular_index])
                              :
                              0));
                        Assert(singular_quadrature, ExcInternalError());

                        // once the singular quadrature has been created, we employ it
                        // to create the corresponding fe_values

                        FEValues<dim-1,dim> fe_v_singular (*fma_mapping, *fma_fe, *singular_quadrature,
                                                           update_jacobians |
                                                           update_values |
                                                           update_cell_normal_vectors |
                                                           update_quadrature_points );

                        fe_v_singular.reinit(cell);

                        // here are the vectors of the quad points and normals vectors

                        const std::vector<Point<dim> > &singular_normals = fe_v_singular.get_normal_vectors();
                        const std::vector<Point<dim> > &singular_q_points = fe_v_singular.get_quadrature_points();


                        // and here is the integrals computation: note how in this case the
                        // values for shape functions & co. are not taken from the precomputed
                        // ones in ComputationalDomain class

                        for (unsigned int q=0; q<singular_quadrature->size(); ++q)
                          {
                            const Tensor<1, dim> R = singular_q_points[q] - support_points[nodeIndex];
                            LaplaceKernel::kernels(R, D, s);
                            for (unsigned int j=0; j<fma_fe->dofs_per_cell; ++j)
                              {
                                local_neumann_matrix_row_i(j) += (( D *
                                                                    singular_normals[q]) *
                                                                  fe_v_singular.shape_value(j,q) *
                                                                  fe_v_singular.JxW(q) );

                                local_dirichlet_matrix_row_i(j) += ( s   *
                                                                     fe_v_singular.shape_value(j,q) *
                                                                     fe_v_singular.JxW(q) );
                              }
                          }
                        if (dim==2)
                          delete singular_quadrature;

                      } // end else

                    // Finally, we need to add
                    // the contributions of the
                    // current cell to the
                    // global matrix.

                    for (unsigned int j=0; j<fma_fe->dofs_per_cell; ++j)
                      {
                          prec_neumann_matrix.add(nodeIndex,local_dof_indices[j],local_neumann_matrix_row_i(j));
                          prec_dirichlet_matrix.add(nodeIndex,local_dof_indices[j],local_dirichlet_matrix_row_i(j));
                          if ((*dirichlet_nodes)(local_dof_indices[j]) > 0.8)
                            init_preconditioner.add(nodeIndex,local_dof_indices[j],-local_dirichlet_matrix_row_i(j));
                          else
                            init_preconditioner.add(nodeIndex,local_dof_indices[j], local_neumann_matrix_row_i(j));
                      }

                  } // end loop on cells of the intList
              }
            } // end loop over nodes of block1
        } // end if (nodes in block > 0)
    } // end loop over childless blocks


  // as said, the direct integrals must not be computed only for the
  // quadPoints in the intList: if a bigger block is in the nonIntList of
  // another block, the bound for the multipole expansion application does
  // not hold, and so we must compute direct integrals. Here we scan the
  // nonIntlists of each block al each level to look for bigger blocks and
  // compute the direct integral contribution for the quadNodes in such
  // blocks

  for (unsigned int level = 1; level <  num_octree_levels + 1;  level++) // loop over levels

    {
      unsigned int startBlockLevel =  startLevel[level];
      // !!! Io spezzerei qui per poi comunicare alla fine (se vogliamo, ma questo viene chiamato poche volte).
      for (unsigned int jj = 0; jj <  dofs_filled_blocks[level].size();  jj++) // loop over blocks of each level
        {
          OctreeBlock<dim> *block1 =  blocks[ dofs_filled_blocks[level][jj]];
          const std::vector <unsigned int> &nodesBlk1Ids = block1->GetBlockNodeList();

          for (unsigned int i = 0; i < nodesBlk1Ids.size(); i++)
            {
              // for each block containing nodes, loop over all sublevels in his NN list (this is because if a
              // block remains childless BEFORE the last level, at this point we need to compute
              // all its contributions up to the bottom level)
              unsigned int nodeIndex = nodesBlk1Ids[i];
              if(this_cpu_set.is_element(nodeIndex))//(m2l_flags[level][jj]==this_mpi_process)
              {

              std::map <cell_it,std::set<unsigned int> > directQuadPoints;

              for (unsigned int subLevel = 0; subLevel < block1->NumNearNeighLevels();  subLevel++)
                {
                  const std::set <unsigned int> &nonIntList = block1->GetNonIntList(subLevel);

                  // loop over well separated blocks of higher size (level)-----> in this case
                  //we must use direct evaluation (luckily being childless they only contain 1 element)
                  for (std::set<unsigned int>::iterator pos = nonIntList.begin(); pos !=nonIntList.lower_bound(startBlockLevel); pos++)
                    {
                      OctreeBlock<dim> *block2 =  blocks[*pos];
                      std::map <cell_it, std::vector<unsigned int> >
                      blockQuadPointsList = block2->GetBlockQuadPointsList();
                      typename std::map <cell_it, std::vector<unsigned int> >::iterator it;
                      for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
                        {
                          for (unsigned int ii=0; ii<(*it).second.size(); ii++)
                            {
                              directQuadPoints[(*it).first].insert((*it).second[ii]);

                              /*////////this is for a check/////////////////////
                                         integralCheck[nodesBlk1Ids[i]][(*it).first] += 1;
                                          ////////////////////////////*/
                            }
                        }
                    } // end loop over blocks of a sublevel of nonIntList
                } // end loop over sublevels

              typename std::map <cell_it, std::set<unsigned int> >::iterator it;
              for (it = directQuadPoints.begin(); it != directQuadPoints.end(); it++)
                {
                  // the vectors with the local integrals for the cell must first
                  // be zeroed
                  local_neumann_matrix_row_i = 0;
                  local_dirichlet_matrix_row_i = 0;

                  // we get the first entry of the map, i.e. the cell pointer
                  // here the quadrature is regular as the cell is well
                  // separated
                  cell_it cell = (*it).first;
                  cell->get_dof_indices(local_dof_indices);
                  // we copy the cell quad points in this set
                  std::set<unsigned int> &cellQuadPoints = (*it).second;

                  //pcout<<"Node "<<i<<"  Elem "<<cell<<" (Direct) Nodes: ";
                  //for(unsigned int j=0; j<fe.dofs_per_cell; ++j) pcout<<" "<<local_dof_indices[j];
                  //pcout<<std::endl;

                  // we start looping on the quad points of the cell: *pos will be the
                  // index of the quad point
                  for (std::set<unsigned int>::iterator pos=cellQuadPoints.begin(); pos!=cellQuadPoints.end(); pos++)
                    {
                      // here we compute the distance R between the node and the quad point
                      const Tensor<1,dim> R =  quadPoints[cell][*pos] - support_points[nodeIndex];
                      LaplaceKernel::kernels(R, D, s);

                      // and here are the integrals for each of the degrees of freedom of the cell: note
                      // how the quadrature values (position, normals, jacobianXweight, shape functions)
                      // are taken from the precomputed ones in ComputationalDomain class

                      for (unsigned int j=0; j<fma_fe->dofs_per_cell; ++j)
                        {
                          local_neumann_matrix_row_i(j) += ( ( D *
                                                               quadNormals[cell][*pos] ) *
                                                             quadShapeFunValues[cell][*pos][j] *
                                                             quadJxW[cell][*pos] );
                          local_dirichlet_matrix_row_i(j) += ( s *
                                                               quadShapeFunValues[cell][*pos][j] *
                                                               quadJxW[cell][*pos] );

                        } // end loop over the dofs in the cell
                    } // end loop over the quad points in a cell

                  // Finally, we need to add
                  // the contributions of the
                  // current cell to the
                  // global matrix.

                  for (unsigned int j=0; j<fma_fe->dofs_per_cell; ++j)
                    {
                      // if(this_cpu_set.is_element(local_dof_indices[j]))
                      // {
                        prec_neumann_matrix.add(nodeIndex,local_dof_indices[j],local_neumann_matrix_row_i(j));
                        prec_dirichlet_matrix.add(nodeIndex,local_dof_indices[j],local_dirichlet_matrix_row_i(j));

                        if ((*dirichlet_nodes)(local_dof_indices[j]) > 0.8)
                          init_preconditioner.add(nodeIndex,local_dof_indices[j],-local_dirichlet_matrix_row_i(j));
                        else
                          init_preconditioner.add(nodeIndex,local_dof_indices[j], local_neumann_matrix_row_i(j));
                      // }
                    }


                } // end loop over quad points in the direct quad points list
              }// end check on proc

            } // end loop over nodes in a block
        }// end loop over block of a level
    }//end loop over octree levels



  pcout<<"...done computing direct integrals"<<std::endl;
}



template <int dim>
void BEMFMA<dim>::multipole_integrals()
{

  pcout<<"Computing multipole integrals..."<<std::endl;
  TimeMonitor LocalTimer(*MultInt);
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

  const unsigned int dofs_per_cell = fma_fe->dofs_per_cell;


  AssertThrow(dofs_per_cell == GeometryInfo<dim-1>::vertices_per_cell,
              ExcMessage("The code in this function can only be used for "
                         "the usual Q1 elements."));

  // now we start looping on the childless blocks to perform the integrals
  for (unsigned int kk = 0; kk <  childlessList.size(); kk++)

    {
      //pcout<<"processing block "<<kk <<"  of  "<<cMesh->GetNumChildlessBlocks()<<std::endl;
      //pcout<<"block "<<cMesh->GetChildlessBlockId(kk) <<"  of  "<<cMesh->GetNumBlocks()<<"  in block list"<<std::endl;

      // we get the current block and its Id, and then we
      // compute its center, which is needed to construct the
      // multipole expansion in which we store the integrals

      unsigned int blockId =  childlessList[kk];
      OctreeBlock<dim> *block =  blocks[blockId];
      double delta = block->GetDelta();
      Point<dim> deltaHalf;
      for (unsigned int i=0; i<dim; i++)
        deltaHalf(i) = delta/2.;
      Point<dim> blockCenter = block->GetPMin()+deltaHalf;

      // at this point, we get the list of quad nodes for the current block,
      // and loop over it
      std::map <cell_it, std::vector <unsigned int> > blockQuadPointsList = block->GetBlockQuadPointsList();

      typename std::map <cell_it, std::vector<unsigned int> >::iterator it;
      for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
        {
          // for each cell in the list, we get the list of its quad nodes
          // present in the current block
          cell_it cell = (*it).first;
          std::vector <unsigned int> &cellQuadPoints = (*it).second;

          // the vectors in the structures that we have previously cleared
          // neet to be resized
          elemMultipoleExpansionsKer1[blockId][cell].resize(dofs_per_cell);
          elemMultipoleExpansionsKer2[blockId][cell].resize(dofs_per_cell);

          // the vectors are now initialized with an empty multipole expansion
          // centered in the current block center
          for (unsigned int j=0; j<fma_fe->dofs_per_cell; ++j)
            {
              elemMultipoleExpansionsKer1[blockId][cell][j] =
                MultipoleExpansion(trunc_order, blockCenter, &assLegFunction);
              elemMultipoleExpansionsKer2[blockId][cell][j] =
                MultipoleExpansion(trunc_order, blockCenter, &assLegFunction);
            }

          // the contribution of each quadrature node (which can be seen as a
          // source with a certain strength) is introduced in the
          // multipole expansion with the appropriate methods (AddNormDer
          // for neumann matrix integrals, Add for dirichlet matrix
          // integrals)
          for (unsigned int k=0; k<cellQuadPoints.size(); ++k)
            {
              unsigned int q = cellQuadPoints[k];
              for (unsigned int j=0; j<fma_fe->dofs_per_cell; ++j)
                {
                  elemMultipoleExpansionsKer1[blockId][cell][j].AddNormDer( quadShapeFunValues[cell][q][j]* quadJxW[cell][q]/4/numbers::PI, quadPoints[cell][q], quadNormals[cell][q]);
                  elemMultipoleExpansionsKer2[blockId][cell][j].Add( quadShapeFunValues[cell][q][j]* quadJxW[cell][q]/4/numbers::PI, quadPoints[cell][q]);
                }
            } // end loop on cell quadrature points in the block

        } // end loop over cells in the block

    } // fine loop sui blocchi childless



  pcout<<"...done computing multipole integrals"<<std::endl;
}

template <int dim>
void BEMFMA<dim>::compute_m2l_flags()
{
  pcout<<"Partitioning the descending phase Multipole2Local"<<std::endl;
  m2l_flags.resize(num_octree_levels+1);
  std::vector<unsigned int> m2l_operations_per_level(num_octree_levels+1);
  std::vector<std::vector<unsigned int> > m2l_operations_per_block(num_octree_levels+1);
  unsigned int my_total_operations=0;

  for(unsigned int level = 1; level <  num_octree_levels + 1;  level++)
  {
    m2l_flags[level].resize(dofs_filled_blocks[level].size());
    m2l_operations_per_block[level].resize(dofs_filled_blocks[level].size());
    m2l_operations_per_level[level] = 0;
    for (unsigned int k = 0; k <  dofs_filled_blocks[level].size();  k++) // loop over blocks of each level
    {
      m2l_operations_per_block[level][k] = 0;
      unsigned int jj =  dofs_filled_blocks[level][k];

      OctreeBlock<dim> *block1 = blocks[jj];
      for (unsigned int subLevel = 0; subLevel < block1->NumNearNeighLevels();  subLevel++)
      {
          m2l_operations_per_block[level][k] += block1->GetNonIntList(subLevel).size();
          m2l_operations_per_level[level] += block1->GetNonIntList(subLevel).size();
      }
    }
    unsigned int test = n_mpi_processes;
    unsigned int operations_per_proc = m2l_operations_per_level[level]/test;     //(int) ceil(m2l_operations_per_level[level]/test);
    int rest_op = m2l_operations_per_level[level]%test;
    unsigned int my_operations=0;
    unsigned int cumulative_check=0;
    unsigned int proc=0;
    unsigned int k=0;
    std::vector<unsigned int> m2l_operations_per_proc(test);
    std::vector<unsigned int> blocks_per_proc(test);
    for (unsigned int k = 0; k <  dofs_filled_blocks[level].size();  k++) // loop over blocks of each level
    {
      // m2l_operations_per_block[level][k] = 0;
      unsigned int jj =  dofs_filled_blocks[level][k];
      OctreeBlock<dim> *block1 = blocks[jj];
      std::vector <unsigned int> nodesBlk1Ids = block1->GetBlockNodeList();
      bool on_process = false;
      for(auto ind : nodesBlk1Ids)
      {
        if(this_cpu_set.is_element(ind))
        {
          on_process = true;
          break;
        }
      }
      if(on_process)
      {
      for (unsigned int subLevel = 0; subLevel < block1->NumNearNeighLevels();  subLevel++)
      {
         my_operations += m2l_operations_per_block[level][k];
         my_total_operations += m2l_operations_per_block[level][k];
      }
      }
    }
    std::cout<<"level --- mpi_proc --- ops"<<std::endl;
    std::cout<<level<<" --- "<<this_mpi_process<<" --- "<<my_operations<<std::endl;
    // while (k <  dofs_filled_blocks[level].size() && proc<test) // loop over blocks of each level
    // {
    //   if(my_operations >= operations_per_proc)
    //   {
    //
    //     rest_op -= my_operations - operations_per_proc;
    //     // pcout<<"On processor "<< proc<<", we have "<<my_operations<<"operations, cumulatively we are taking "<<cumulative_check<<" m2l ops over a total of "<<m2l_operations_per_level[level]<<std::endl;
    //     proc += 1;
    //     my_operations = 0;
    //   }
    //   m2l_flags[level][k]=proc;
    //   my_operations += m2l_operations_per_block[level][k];
    //   cumulative_check += m2l_operations_per_block[level][k];
    //   m2l_operations_per_proc[proc] += m2l_operations_per_block[level][k];
    //   blocks_per_proc[proc] += 1;
    //   k+=1;
    // }
    // pcout<<"LEVEL "<<level<<std::endl;
    // pcout<<"Rest is "<<rest_op<<" last block is "<<k<<" , total of "<<dofs_filled_blocks[level].size()<<
    //        " operations taken "<<cumulative_check<<" over a total of "<<m2l_operations_per_level[level]<<std::endl;
    // for(proc = 0 ; proc < test; ++proc)
    //   pcout<<"On processor "<< proc<<", we have "<<m2l_operations_per_proc[proc]<<" m2l operations, and "
    //        <<blocks_per_proc[proc]<<" blocks over a total of "<<dofs_filled_blocks[level].size()<<std::endl;
  }

  std::cout<<"finalcountmpi_proc --- ops"<<std::endl;
  std::cout<<this_mpi_process<<" --- "<<my_total_operations<<std::endl;



}


template <int dim>
void BEMFMA<dim>::generate_multipole_expansions(const TrilinosWrappers::MPI::Vector &phi_values_in, const TrilinosWrappers::MPI::Vector &dphi_dn_values_in) const
{
  pcout<<"Generating multipole expansions..."<<std::endl;
  TimeMonitor LocalTimer(*MultGen);
  // TODO I DON'T KNOW IF WE CAN SPLIT THIS OPERATION, IN CASE THIS WOULD NOT BE EASY COMMUNICATION, WE WOULD NEED TO
  // COMMUNICATE THE ENTIRE CLASSES THAT HOLD THE MULTIPOLE<->LOCAL EXPANSIONS
  const Vector<double> phi_values(phi_values_in);
  const Vector<double> dphi_dn_values(dphi_dn_values_in);
// also here we clear the structures storing the multipole
// expansions for Dirichlet and Neumann matrices


  blockMultipoleExpansionsKer1.clear();
  blockMultipoleExpansionsKer2.clear();

// we reisze them: there's going to be an expansion per block

  blockMultipoleExpansionsKer1.resize( num_blocks);
  blockMultipoleExpansionsKer2.resize( num_blocks);

// these two variables will be handy in the following

  std::vector<unsigned int> local_dof_indices(fma_fe->dofs_per_cell);
  double delta;

// we loop on blocks and for each of them we create an empty multipole expansion
// centered in the block center

  for (unsigned int ii = 0; ii <  num_blocks ; ii++)

    {
      delta =  blocks[ii]->GetDelta();
      Point<dim> deltaHalf;
      for (unsigned int i=0; i<dim; i++)
        deltaHalf(i) = delta/2.;

      Point<dim> blockCenter =  blocks[ii]->GetPMin()+deltaHalf;

      blockMultipoleExpansionsKer1[ii] = MultipoleExpansion(trunc_order, blockCenter, &assLegFunction);
      blockMultipoleExpansionsKer2[ii] = MultipoleExpansion(trunc_order, blockCenter, &assLegFunction);
    }

// we now begin the rising phase of the algorithm: starting from the lowest block levels (childless blocks)
// we get all the values of the multipole integrals and aggregate them in the multipole expansion for
// each blocks

  for (unsigned int kk = 0; kk <  childlessList.size(); kk++)

    {

      // for each block we get the center and the quad points

      unsigned int blockId =  childlessList[kk];
      OctreeBlock<dim> *block =  blocks[blockId];

      delta =  blocks[blockId]->GetDelta();
      Point<dim> deltaHalf;
      for (unsigned int i=0; i<dim; i++)
        deltaHalf(i) = delta/2.;
      //Point<dim> blockCenter =  blocks[blockId]->GetPMin()+deltaHalf;

      std::map <cell_it, std::vector <unsigned int> > blockQuadPointsList = block->GetBlockQuadPointsList();

      // we loop on the cells of the quad points in the block: remember that for each cell with a node in the
      // block, we had created a set of dofs_per_cell multipole expansion, representing
      // the (partial) integral on each cell

      typename std::map <cell_it, std::vector<unsigned int> >::iterator it;
      for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
        {
          cell_it cell = (*it).first;
          cell->get_dof_indices(local_dof_indices);

          // for each cell we get the dof_indices, and add to the block multipole expansion,
          // the integral previously computed, multiplied by the phi or dphi_dn value at the
          // corresponding dof of the cell. A suitable MultipoleExpansion class method has been
          // created for this purpose

          for (unsigned int jj=0; jj < fma_fe->dofs_per_cell; ++jj)
            {
              blockMultipoleExpansionsKer2.at(blockId).Add(elemMultipoleExpansionsKer2[blockId][cell][jj],dphi_dn_values(local_dof_indices[jj]));
              blockMultipoleExpansionsKer1.at(blockId).Add(elemMultipoleExpansionsKer1[blockId][cell][jj],phi_values(local_dof_indices[jj]));
            }
        } //end loop ond block elements
    } // end loop on childless blocks


// now all the lower level blocks have a multipole expansion containing the contribution to the integrals
// of all the quad points in them: we now begin summing the child multipole expansion to the the parents
// expansions: to do that we need to translate che children expansions to the parent block center: there
// is a MultipoleExpansion clas for this too

// we loop the levels starting from the bottom one

  for (unsigned int level =  num_octree_levels; level > 0; level--)

    {

      // for each block we add the (translated) multipole expansion to the the parent expansion

      // pcout<<"processing level "<<level <<"  of  "<< num_octree_levels<<std::endl;

      for (unsigned int kk =  startLevel[level]; kk <  endLevel[level]+1; kk++)

        {
          unsigned int parentId =  blocks[kk]->GetParentId();

          blockMultipoleExpansionsKer1.at(parentId).Add(blockMultipoleExpansionsKer1.at(kk));
          blockMultipoleExpansionsKer2.at(parentId).Add(blockMultipoleExpansionsKer2.at(kk));

        } // end loop over blocks of a level

    } // end loop over levels


  pcout<<"...done generating multipole expansions"<<std::endl;

}



template <int dim>
void BEMFMA<dim>::multipole_matr_vect_products(const TrilinosWrappers::MPI::Vector &phi_values, const TrilinosWrappers::MPI::Vector &dphi_dn_values,
                                               TrilinosWrappers::MPI::Vector &matrVectProdN,    TrilinosWrappers::MPI::Vector &matrVectProdD) const
{
  pcout<<"Computing multipole matrix-vector products... "<<std::endl;
  TimeMonitor LocalTimer(*MatrVec);
  // we start re-initializing matrix-vector-product vectors
  // matrVectProdN = 0;
  // matrVectProdD = 0;

  //and here we compute the direct integral contributions (stored in two sparse matrices)
  // const Vector<double> phi_values_loc(phi_values);
  // const Vector<double> dphi_dn_values_loc(dphi_dn_values);
  // Vector<double> matrVectProdD(fma_dh->n_dofs());
  // Vector<double> matrVectProdN(fma_dh->n_dofs());
  prec_neumann_matrix.vmult(matrVectProdN, phi_values);
  prec_dirichlet_matrix.vmult(matrVectProdD, dphi_dn_values);

  // if(this_cpu_set.is_element(55))
  //   std::cout<<matrVectProdN[55]<<" "<<phi_values[55]<<" "<<matrVectProdD[55]<<" "<<dphi_dn_values[55]<<std::endl;

  //from here on, we compute the multipole expansions contributions

  // we start cleaning past sessions
  blockLocalExpansionsKer1.clear();
  blockLocalExpansionsKer2.clear();

  blockLocalExpansionsKer1.resize( num_blocks);
  blockLocalExpansionsKer2.resize( num_blocks);

  // we declare some familiar variables that will be useful in the method
  std::vector<Point<dim> > support_points(fma_dh->n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *fma_mapping,
                                                    *fma_dh, support_points);
  std::vector<unsigned int> local_dof_indices(fma_fe->dofs_per_cell);
  double delta;


  // here we loop on all the blocks and build an empty local expansion for
  // each of them
  for (unsigned int ii = 0; ii <  num_blocks ; ii++)

    {
      delta =  blocks[ii]->GetDelta();
      Point<dim> deltaHalf;
      for (unsigned int i=0; i<dim; i++)
        deltaHalf(i) = delta/2.;
      Point<dim> blockCenter =  blocks[ii]->GetPMin()+deltaHalf;

      blockLocalExpansionsKer1[ii] =  LocalExpansion(trunc_order, blockCenter, &assLegFunction);
      blockLocalExpansionsKer2[ii] =  LocalExpansion(trunc_order, blockCenter, &assLegFunction);
    }


  // here the descending phase of the algorithm begins: we will start from the top level, and loop
  // on each block of each level: all the needed multipole integrals contributions must be
  // introduced in each specific block local expansion. to do this, for each block, we first need to
  // translate into the block center the parent block local expansion (it accounts for all far blocks
  // contribution of the parent, and all blocks that are far from the parent, are far from its
  // children); then, we have to consider all the blocks in the current block's nonIntList,
  // and convert  their multipole expansions to the local expansion of the current block
  // (there is only one exception to this procedure, and we'll comment it in the following)

  // so, here we loop over levels, starting form lower levels (bigger blocks)
  for (unsigned int level = 1; level <  num_octree_levels + 1;  level++)

    {
      // pcout<<"processing level "<<level <<"  of  "<< num_octree_levels<<std::endl;

      // we get the ids of the first and last block of the level
      unsigned int startBlockLevel =  startLevel[level];
      unsigned int endBlockLevel =  endLevel[level];

      // to reduce computational cost, we decide to loop on the list of blocks which
      // contain at least one node (the local and multipole expansion will be finally evaluated
      // at the nodes positions)

      // TODO WE COULD SPLIT THIS LOOP OVER THE BLOCKS OVER ALL THE PROCESSORS, THEN DO A TRILINOS.ADD WITH
      // DIFFERENT MAPS. HERE WE NEED A FULL ONE.
      for (unsigned int k = 0; k <  dofs_filled_blocks[level].size();  k++) // loop over blocks of each level
        {

          //pcout<<"Block "<<jj<<std::endl;
          unsigned int jj =  dofs_filled_blocks[level][k];
          OctreeBlock<dim> *block1 =  blocks[jj];
          unsigned int block1Parent = block1->GetParentId();
          std::vector <unsigned int> nodesBlk1Ids = block1->GetBlockNodeList();
          bool on_process = false;
          for(auto ind : nodesBlk1Ids)
          {
            if(this_cpu_set.is_element(ind))
            {
              on_process = true;
              break;
            }
          }
          if(on_process)
          {

          // here we get parent contribution to local expansion and transfer it in current
          // block: this operation requires a local expansion translation, implemented
          // in a specific LocalExpansion class member

          blockLocalExpansionsKer1[jj].Add(blockLocalExpansionsKer1[block1Parent]);
          blockLocalExpansionsKer2[jj].Add(blockLocalExpansionsKer2[block1Parent]);

          // for each block, loop over all sublevels in his NN list (this is because if a
          // block remains childless BEFORE the last level, at this point we need to compute
          // all its contributions up to the bottom level)


          for (unsigned int subLevel = 0; subLevel < block1->NumNearNeighLevels();  subLevel++)
            {
              // we get the nonIntList of each block

              std::set <unsigned int> nonIntList = block1->GetNonIntList(subLevel);
              std::set <cell_it> nonIntListElemsBlk1;

              // we start converting into local expansions, all the multipole expansions
              // of all the blocks in nonIntList, that are of the same size (level) of the
              // current block. To perform the conversion, we use another member of the
              // LocalExpansion class. Note that all the contributions to the integrals
              // of blocks bigger than current block had already been considered in
              // the direct integrals method (the bounds of the local and multipole
              // expansion do not hold in such case)


              for (std::set <unsigned int>::iterator pos1 = nonIntList.lower_bound(startBlockLevel); pos1 != nonIntList.upper_bound(endBlockLevel);  pos1++) // loop over NNs of NNs and add them to intList
                {
                  //pcout<<"NonIntListPart2 Blocks: "<<*pos1<<" ";
                  unsigned int block2Id = *pos1;
                  blockLocalExpansionsKer1[jj].Add(blockMultipoleExpansionsKer1[block2Id]);
                  blockLocalExpansionsKer2[jj].Add(blockMultipoleExpansionsKer2[block2Id]);

                  /*////////this is for a check///////////////////////
                    OctreeBlock<dim>* block2 =  blocks[block2Id];
                  std::vector<unsigned int> nodesBlk1Ids = block1->GetBlockNodeList();
                  std::map <cell_it, std::vector<unsigned int> >
                                    blockQuadPointsList = block2->GetBlockQuadPointsList();
                                    typename std::map <cell_it, std::vector<unsigned int> >::iterator it;
                                          for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
                                        {
                                  for (unsigned int i=0; i<(*it).second.size(); i++)
                                      {
                                      for (unsigned int kk=0; kk<nodesBlk1Ids.size(); kk++)
                                         integralCheck[nodesBlk1Ids[kk]][(*it).first] += 1;
                          }
                                              }
                                          ////////////////////////////*/
                } // end loop over well separated blocks of the same size (level)


              // loop over well separated blocks of the smaller size (level)-----> use multipoles without local expansions

              // we now have to loop over blocks in the nonIntList that are smaller than the current
              // blocks: in this case the bound for the conversion of a multipole into local expansion
              // do not hold, but the bound for the evaluation of the multipole expansions does hold:
              // thus, we will simply evaluate the multipole expansions of such blocks for each node in
              // the clock

              for (std::set <unsigned int>::iterator pos1 = nonIntList.upper_bound(endBlockLevel); pos1 != nonIntList.end();  pos1++)
                {
                  //pcout<<"NonIntListPart3 Blocks: "<<*pos1<<" ";
                  unsigned int block2Id = *pos1;
                  //std::vector <cell_it> elemBlk2Ids = block2.GetBlockElementsList();
                  for (unsigned int ii = 0; ii < nodesBlk1Ids.size(); ii++) //loop over each node of block1
                    {
                      Point<dim> &nodeBlk1 = support_points[nodesBlk1Ids.at(ii)];
                      matrVectProdD(nodesBlk1Ids[ii]) += blockMultipoleExpansionsKer2[block2Id].Evaluate(nodeBlk1);
                      matrVectProdN(nodesBlk1Ids[ii]) += blockMultipoleExpansionsKer1[block2Id].Evaluate(nodeBlk1);

                      /*//////this is for a check/////////////////////////
                            OctreeBlock<dim>* block2 =  blocks[block2Id];
                            std::map <cell_it, std::vector<unsigned int> >
                                              blockQuadPointsList = block2->GetBlockQuadPointsList();
                                              typename std::map <cell_it, std::vector<unsigned int> >::iterator it;
                                                    for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
                                                  {
                                            for (unsigned int i=0; i<(*it).second.size(); i++)
                                                {
                                             integralCheck[nodesBlk1Ids[ii]][(*it).first] += 1;
                                    }
                                                        }
                                                    ////////////////////////////*/
                    }
                } // end loop over well separated blocks of smaller size (level)


            } // end loop over all sublevels in  nonIntlist
          }//end if for m2l flags
        } // end loop over blocks of each level

    } // end loop over levels


  // finally, when the loop over levels is done, we need to evaluate local expansions of all
  // childless blocks, at each block node(s)

  // TODO SPLIT THIS LOOP OVER ALL PROCESSORS THEN COMMUNICATE.
  // if(this_mpi_process==0)
  for (unsigned int kk = 0; kk <  childlessList.size(); kk++)

    {
      unsigned int block1Id =  childlessList[kk];
      OctreeBlock<dim> *block1 =  blocks[block1Id];
      std::vector <unsigned int> nodesBlk1Ids = block1->GetBlockNodeList();

      // loop over nodes of block
      for (unsigned int ii = 0; ii < nodesBlk1Ids.size(); ii++) //loop over each node of block1
        {
          if(this_cpu_set.is_element(nodesBlk1Ids[ii]))
          {
            Point<dim> &nodeBlk1 = support_points[nodesBlk1Ids.at(ii)];
            matrVectProdD(nodesBlk1Ids[ii]) += (blockLocalExpansionsKer2[block1Id]).Evaluate(nodeBlk1);
            matrVectProdN(nodesBlk1Ids[ii]) += (blockLocalExpansionsKer1[block1Id]).Evaluate(nodeBlk1);
          }
        } // end loop over nodes

    } // end loop over childless blocks

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
  // TrilinosWrappers::MPI::Vector Ndummy(matrVectProdN_in.locally_owned_elements(), mpi_communicator);
  // TrilinosWrappers::MPI::Vector Ddummy(matrVectProdD_in.locally_owned_elements(), mpi_communicator);
  // Ndummy = matrVectProdN;
  // Ddummy = matrVectProdD;
  // matrVectProdN_in.add( Ndummy);
  // matrVectProdD_in.add( Ddummy);

  pcout<<"...done computing multipole matrix-vector products"<<std::endl;

}

// this method computes the preconditioner needed for the GMRES:
// to do it, it needs to receive the alpha vector from the bem_problem
// class, along with the constraint matrix of the bem problem

template <int dim>
TrilinosWrappers::PreconditionILU &BEMFMA<dim>::FMA_preconditioner(const TrilinosWrappers::MPI::Vector &alpha, ConstraintMatrix &c )//TO BE CHANGED!!!
{
  TimeMonitor LocalTimer(*PrecondTime);
  // the final preconditioner (with constraints) has a slightly different sparsity pattern with respect
  // to the non constrained one. we must here initialize such sparsity pattern
  final_prec_sparsity_pattern.reinit(alpha.vector_partitioner(),125*fma_fe->dofs_per_cell);
  //final_prec_sparsity_pattern.reinit(fma_dh->n_dofs(),fma_dh->n_dofs(),125*fma_fe->dofs_per_cell);

  // IndexSet this_cpu_set(alpha.locally_owned_elements());
  for (unsigned int i=0; i < fma_dh->n_dofs(); i++)
    {
      if(this_cpu_set.is_element(i))
      {
        if (c.is_constrained(i))
          {
            //cout<<i<<"  (c):"<<endl;
            // constrained nodes entries are taken from the bem problem constraint matrix
            final_prec_sparsity_pattern.add(i,i);
            const std::vector< std::pair < unsigned int, double > >
            *entries = c.get_constraint_entries (i);
            for (unsigned int j=0; j< entries->size(); ++j)
              final_prec_sparsity_pattern.add(i,(*entries)[j].first);
          }
        else
          {
            //cout<<i<<"  (nc): ";
            // other nodes entries are taken from the unconstrained preconditioner matrix
            for (unsigned int j=0; j<fma_dh->n_dofs(); ++j)
              {
                if (init_prec_sparsity_pattern.exists(i,j))
                  {
                    final_prec_sparsity_pattern.add(i,j);
                    //cout<<j<<" ";
                  }
              }
            //cout<<endl;
          }
      }
    }


  final_prec_sparsity_pattern.compress();
  final_preconditioner.reinit(final_prec_sparsity_pattern);

  // std::cout<<"ok prec sparsity pattern"<<std::endl;
  // now we assemble the final preconditioner matrix: the loop works
  // exactly like the previous one
  for (unsigned int i=0; i < fma_dh->n_dofs(); i++)
  {
    if(this_cpu_set.is_element(i))
    {
      if (c.is_constrained(i))
        {
          final_preconditioner.set(i,i,1);
          //pcout<<i<<" "<<i<<"  ** "<<final_preconditioner(i,i)<<std::endl;
          // constrainednodes entries are taken from the bem problem constraint matrix
          const std::vector< std::pair < unsigned int, double > >
          *entries = c.get_constraint_entries (i);
          for (unsigned int j=0; j< entries->size(); ++j)
            {
              final_preconditioner.set(i,(*entries)[j].first,(*entries)[j].second);
              //pcout<<i<<" "<<(*entries)[j].first<<"  * "<<(*entries)[j].second<<std::endl;
            }
        }
      else
        {
          // other nodes entries are taken from the unconstrained preconditioner matrix
          for (unsigned int j=0; j<fma_dh->n_dofs(); ++j)
            {
              // QUI CHECK SU NEUMANN - DIRICHLET PER METTERE A POSTO, tanto lui già conosce le matrici.
              if (init_prec_sparsity_pattern.exists(i,j))
                {
                  final_preconditioner.set(i,j,init_preconditioner(i,j));
                  //pcout<<i<<" "<<j<<" "<<init_preconditioner(i,j)<<std::endl;
                }
            }

        }
    }
  }
  // std::cout<<"now alpha"<<std::endl;
  // finally, we have to add the alpha values on the diagonal, whenever dealing with a
  // neumann (in such nodes the potential phi is an unknown) and non constrained node

  for (unsigned int i=0; i < fma_dh->n_dofs(); i++)
  {
    if(this_cpu_set.is_element(i))
    {
      if ( (*dirichlet_nodes)(i) == 0 && !(c.is_constrained(i)))
        {
          final_preconditioner.add(i,i,alpha(i));
          //pcout<<i<<" "<<i<<" "<<final_preconditioner(i,i)<<std::endl;
        }
      else // this is just to avoid a deadlock. we need a better strategy
      {
        final_preconditioner.add(i,i,0);
      }
    }
  }

  //preconditioner.print_formatted(pcout,4,true,0," 0 ",1.);
  preconditioner.initialize(final_preconditioner);

  return preconditioner;
}

template <int dim>
void BEMFMA<dim>::compute_geometry_cache()
{
  pcout<<"Generating geometry cache..."<<std::endl;

  // @sect5{ComputationalDomain::generate_double_nodes_set}

  // The following is the function
  // which creates a set containing
  // the double nodes.


  FESystem<dim-1,dim> gradient_fe(*fma_fe, dim);
  DoFHandler<dim-1, dim> gradient_dh(fma_dh->get_tria());

  //double tol = 1e-8;
  std::vector<Point<dim> > support_points(fma_dh->n_dofs());

  DoFTools::map_dofs_to_support_points<dim-1, dim>( *fma_mapping,
                                                    *fma_dh, support_points);

  for (unsigned int i=0; i<fma_dh->n_dofs(); ++i)
    {
      for (unsigned int j=0; j<fma_dh->n_dofs(); ++j)
        {
          //pcout<<"i "<<i<<" ("<<support_points[i]<<")  j "<<j<<" ("<<support_points[j]<<")  distance "<<support_points[i].distance(support_points[j])<<std::endl;
        }

    }
  // for the gradient dofs finding coupled
  // dofs is a little bit difficult, as the
  // gradient is a vectorial function: usually
  // in such case the dofs are numbered
  // so that for each support point dim
  // consecutive dofs represent each component
  // of the vector field: in this case
  // (and only in this case) the following
  // piece of code works


  cell_it
  gradient_cell = gradient_dh.begin_active(),
  gradient_endc = gradient_dh.end();

  cell_it
  cell = fma_dh->begin_active(),
  endc = fma_dh->end();

  std::vector<unsigned int> dofs(fma_fe->dofs_per_cell);
  std::vector<unsigned int> gradient_dofs(fma_fe->dofs_per_cell);
  // mappa che associa ad ogni dof le celle cui esso appartiene
  dof_to_elems.clear();

  // mappa che associa ad ogni gradient dof le celle cui esso appartiene
  gradient_dof_to_elems.clear();

  // vettore che associa ad ogni gradient dof la sua componente
  gradient_dof_components.clear();
  gradient_dof_components.resize(gradient_dh.n_dofs());

  // mappa che associa ad ogni cella un set contenente le celle circostanti
  elem_to_surr_elems.clear();



  for (; cell!=endc,gradient_cell!=gradient_endc; ++cell,++gradient_cell)
    {
      Assert(cell->index() == gradient_cell->index(), ExcInternalError());

      cell->get_dof_indices(dofs);
      for (unsigned int j=0; j<fma_fe->dofs_per_cell; ++j)
        {
          dof_to_elems[dofs[j]].push_back(cell);
        }
      gradient_cell->get_dof_indices(gradient_dofs);
      for (unsigned int j=0; j<gradient_fe.dofs_per_cell; ++j)
        {
          gradient_dof_to_elems[gradient_dofs[j]].push_back(gradient_cell);
          gradient_dof_components[gradient_dofs[j]] = gradient_fe.system_to_component_index(j).first;
        }
    }

  // qui viene creata la mappa dei elmenti che circondano ciascun elemento
  for (cell = fma_dh->begin_active(); cell != endc; ++cell)
    {
      cell->get_dof_indices(dofs);
      for (unsigned int j=0; j<fma_fe->dofs_per_cell; ++j)
        {
          std::set <unsigned int> duplicates = (*double_nodes_set)[dofs[j]];
          for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
            {
              std::vector<cell_it>
              dof_cell_list =  dof_to_elems[*pos];
              for (unsigned int k=0; k<dof_cell_list.size(); ++k)
                elem_to_surr_elems[cell].insert(dof_cell_list[k]);
            }
        }
    }

  // gradient_cell = gradient_dh.begin_active();
  // cell = fma_dh->begin_active();
  //
  // for (; cell!=endc,gradient_cell!=gradient_endc; ++cell,++gradient_cell)
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
  //      std::set <unsigned int> duplicates = gradient_double_nodes_set[gradient_dofs[j]];
  //      for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
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
  //            //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"  surfNodes: "<<dirichlet_nodes(dofs[i])<<"  otherNodes: "<<other_nodes(dofs[i])<<std::endl;
  //        }
  //     }
  //
  // gradient_cell = gradient_dh.begin_active();
  // cell = fma_dh->begin_active();
  //
  // for (; cell!=endc,gradient_cell!=gradient_endc; ++cell,++gradient_cell)
  //     {
  //     Assert(cell->index() == gradient_cell->index(), ExcInternalError());
  //
  //     if (gradient_cell->material_id() == wall_sur_ID1)
  //  {
  //        gradient_cell->get_dof_indices(gradient_dofs);
  //        for (unsigned int j=0; j<gradient_fe.dofs_per_cell; ++j)
  //            {
  //      std::set <unsigned int> duplicates = gradient_double_nodes_set[gradient_dofs[j]];
  //      for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
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
  //            //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"  surfNodes: "<<dirichlet_nodes(dofs[i])<<"  otherNodes: "<<other_nodes(dofs[i])<<std::endl;
  //        }
  //     }
  // gradient_cell = gradient_dh.begin_active();
  // cell = fma_dh->begin_active();
  //
  // for (; cell!=endc,gradient_cell!=gradient_endc; ++cell,++gradient_cell)
  //     {
  //     Assert(cell->index() == gradient_cell->index(), ExcInternalError());
  //
  //     if (gradient_cell->material_id() == wall_sur_ID2)
  //  {
  //        // This is a free surface node.
  //        gradient_cell->get_dof_indices(gradient_dofs);
  //        for (unsigned int j=0; j<gradient_fe.dofs_per_cell; ++j)
  //            {
  //      std::set <unsigned int> duplicates = gradient_double_nodes_set[gradient_dofs[j]];
  //      for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
  //          {
  //          cell_it duplicate_cell =  gradient_dof_to_elems[*pos][0];
  //   if (duplicate_cell->material_id() == wall_sur_ID3)
  //      {
  //      boat_keel_nodes.insert(gradient_dofs[j]);
  //      }
  //
  //          }
  //      }
  //            //pcout<<dofs[i]<<"  cellMatId "<<cell->material_id()<<"  surfNodes: "<<dirichlet_nodes(dofs[i])<<"  otherNodes: "<<other_nodes(dofs[i])<<std::endl;
  //        }
  //     }


  pcout<<"...done"<<std::endl;

}
template <int dim>
void BEMFMA<dim>::generate_octree_blocking()
{

  pcout<<"Generating octree blocking... "<<std::endl;
  TimeMonitor LocalTimer(*ListCreat);
  // @sect5{BEMProblem::generate_double_nodes_set}

  // The following is the function
  // which creates the octree blocking
  // for the fast multipole algorithm


  std::vector<Point<dim> > support_points(fma_dh->n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *fma_mapping,
                                                    *fma_dh, support_points);

  // !!!TO BE CHANGED
  quadrature = SP(new QGauss<dim-1>(4));
  FEValues<dim-1,dim> fe_v(*fma_mapping,*fma_fe, *quadrature,
                           update_values |
                           update_cell_normal_vectors |
                           update_quadrature_points |
                           update_JxW_values);

  double max_coor_value = 0;

  for (unsigned int i=0; i < fma_dh->n_dofs(); i++)
    {
      //for printout
      //pcout<<"Node "<<i<< "["<<support_points[i]<<"] "<<std::endl;
      for (unsigned int j=0; j < dim; j++)
        {
          max_coor_value = std::max(max_coor_value,std::abs(support_points[i](j)));
        }
    }

  if (blocks.size() > 0)
    {
      for (unsigned int ii = 0; ii < num_blocks;  ii++)
        delete blocks[ii];
    }

  unsigned int maxNumBlocks = num_octree_levels*fma_dh->get_tria().n_active_cells()*fe_v.n_quadrature_points;
//unsigned int maxNumBlocks = 0;
//for (unsigned int ii = 0; ii < num_octree_levels + 1;  ii++)
//  {
//   maxNumBlocks += int(pow(8.,double(ii)));
//  }

  blocks.clear();
  blocks.reserve(maxNumBlocks);
  blocks.resize(maxNumBlocks);

  unsigned int blocksCount = 0;
  startLevel.resize(num_octree_levels+1);
  endLevel.resize(num_octree_levels+1);

//for (unsigned int j=0; j < num_octree_levels + 1; j++)
//     parentList[j].clear();
  parentList.clear();
  parentList.resize(num_octree_levels+1);
  parentList[0].push_back(0);


  childlessList.clear();
  unsigned int numChildless = 0;
  numParent.resize(num_octree_levels+1);

//qui di seguito vengono reinizializzate strutture utili al multipolo

// mappa che associa ad ogni dof un vettore con i blocchi cui essa appartiene per ogni livello
  dof_to_block.clear();

// mappa che associa ad ogni quad point un vettore con i blocchi cui essa appartiene per ogni livello
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

  dofs_filled_blocks.resize(num_octree_levels+1);

  quad_points_filled_blocks.resize(num_octree_levels+1);



  for (unsigned int ii = 0; ii < num_octree_levels + 1 ;  ii++)
    {
      numParent[ii] = 0;
    }



  Point<dim> pMin;
  for (int i=0; i<dim; i++)
    pMin(i) = -1.1*max_coor_value;

  // delta e' il lato del kazzo di kubo...
  double delta = 2.2*max_coor_value;

  OctreeBlock<dim> *block = new OctreeBlock<dim>(0, 0, pMin, delta);

  std::vector<unsigned int> local_dof_indices(fma_fe->dofs_per_cell);
  cell_it
  cell = fma_dh->begin_active(),
  endc = fma_dh->end();
  for (cell = fma_dh->begin_active(); cell != endc; ++cell)
    {
      fe_v.reinit(cell);
      const unsigned int n_q_points = fe_v.n_quadrature_points;
      quadPoints[cell] = fe_v.get_quadrature_points();
      quadNormals[cell] = fe_v.get_normal_vectors();
      quadJxW[cell].resize(n_q_points);
      quadShapeFunValues[cell].resize(n_q_points);
      for (unsigned int q=0; q<n_q_points; ++q)
        {
          quadJxW[cell][q] = fe_v.JxW(q);
          for (unsigned int j=0; j<fma_fe->dofs_per_cell; ++j)
            quadShapeFunValues[cell][q].push_back(fe_v.shape_value(j,q));
        }

      quad_point_to_block[cell].resize(n_q_points);
      for (unsigned int j=0; j<n_q_points; ++j)
        {
          block->AddQuadPoint(cell,j);
          quad_point_to_block[cell][j].push_back(0);
        }

      cell->get_dof_indices(local_dof_indices);
      for (unsigned int j=0; j<fma_fe->dofs_per_cell; ++j)
        {
          dof_to_elems[local_dof_indices[j]].push_back(cell);
        }
    }

  for (unsigned int ii = 0; ii < fma_dh->n_dofs(); ii++)
    {
      block->AddNode(ii);
      dof_to_block[ii].push_back(0);
    }


// just for output
  /*for (cell = fma_dh->begin_active(); cell != endc; ++cell)
      {
      std::set<cell_it> surr_elems = elem_to_surr_elems[cell];
      pcout<<std::endl<<"cell "<<cell<<"  surrounded by: ";
      for (typename std::set<cell_it>::iterator pos = surr_elems.begin(); pos !=surr_elems.end(); pos++)
           pcout<<" "<<*pos;
      }*/

  blocks[0] = block;
  numParent[0] = 1;

//pcout<<"blocks[0].GetBlockChildrenNum() "<<blocks[0].GetBlockChildrenNum()<<std::endl;

  /*pcout<<std::endl;
  pcout<<blocks[0].GetPMin()<<std::endl;
  pcout<<pMin<<std::endl;
  pcout<<block.GetDelta()<<std::endl;
  pcout<<block.GetBlockNodeList()[0]<<std::endl;
  pcout<<block.GetBlockElementsList()[1]<<std::endl;
  pcout<<delta<<std::endl;
  pcout<<std::endl;//*/

  unsigned int quadPointsInChildless = 0;
  unsigned int nodesInChildless = 0;

  for (unsigned int level = 1; level < num_octree_levels + 1;  level++)

    {
      unsigned int quadPointsCheck = quadPointsInChildless;
      unsigned int nodesCheck = nodesInChildless;
      delta /= 2.;

      for (unsigned int kk = 0; kk < numParent[level-1];  kk++)

        {
          unsigned int jj = parentList[level-1][kk];
          //pcout<<" level "<<level<<"     block "<<jj<<std::endl;
          OctreeBlock<dim> *parent = blocks[jj];
          //pcout<<"parent.GetBlockChildrenNum() "<<parent.GetBlockChildrenNum()<<std::endl;
          //pcout<<" Pmin "<<parent.GetPMin()(0)<<", "<<parent.GetPMin()(1)<<", "<<parent.GetPMin()(2)<<" "<<std::endl;
          //pcout<<" delta "<<parent.GetDelta()<<" "<<std::endl;

          pMin = parent->GetPMin();
          unsigned int num_children_per_block = int(pow((double)2,(double)dim));
          std::vector<OctreeBlock<dim> *> children(num_children_per_block);

          if (dim == 3)
            {
              children[0] = new OctreeBlock<dim>(level, jj, pMin                              , delta);
              children[1] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta   ,0.,   0.), delta);
              children[2] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta,delta,   0.), delta);
              children[3] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(0.   ,delta,   0.), delta);
              children[4] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(0.   ,   0.,delta), delta);
              children[5] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta,   0.,delta), delta);
              children[6] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta,delta,delta), delta);
              children[7] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(   0.,delta,delta), delta);
            }

          if (dim == 2)
            {
              children[0] = new OctreeBlock<dim>(level, jj, pMin                        , delta);
              children[1] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta   ,0.), delta);
              children[2] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta,delta), delta);
              children[3] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(0.   ,delta), delta);

            }

          std::map <cell_it, std::vector <unsigned int> > blockQuadPointsList =
            parent->GetBlockQuadPointsList();

          std::vector <unsigned int> blockNodeList = parent->GetBlockNodeList();

          if (dim == 3)
            {
              for (unsigned int i = 0; i < blockNodeList.size(); i++)
                {
                  Point <dim> node = support_points[blockNodeList[i]];
                  // assegnamento nodi del blocco padre ai blocchi figli

                  if (node(2) <= parent->GetPMin()(2)+delta)
                    {
                      if (node(1) <= parent->GetPMin()(1)+delta)
                        {
                          if (node(0) <= parent->GetPMin()(0)+delta)
                            {
                              //pcout<<" Sono in 1 "<<std::endl;
                              children[0]->AddNode(blockNodeList[i]);
                            }
                          else
                            {
                              //pcout<<" Sono in 2 "<<std::endl;
                              children[1]->AddNode(blockNodeList[i]);
                            }
                        }
                      else
                        {
                          if (node(0) <= parent->GetPMin()(0)+delta)
                            {
                              //pcout<<" Sono in 4 "<<std::endl;
                              children[3]->AddNode(blockNodeList[i]);
                            }
                          else
                            {
                              //pcout<<" Sono in 3 "<<std::endl;
                              children[2]->AddNode(blockNodeList[i]);
                            }
                        }
                    }
                  else
                    {
                      if (node(1) <= parent->GetPMin()(1)+delta)
                        {
                          if (node(0) <= parent->GetPMin()(0)+delta)
                            {
                              //pcout<<" Sono in 5 "<<std::endl;
                              children[4]->AddNode(blockNodeList[i]);
                            }
                          else
                            {
                              //pcout<<" Sono in 6 "<<std::endl;
                              children[5]->AddNode(blockNodeList[i]);
                            }
                        }
                      else
                        {
                          if (node(0) <= parent->GetPMin()(0)+delta)
                            {
                              //pcout<<" Sono in 8 "<<std::endl;
                              children[7]->AddNode(blockNodeList[i]);
                            }
                          else
                            {
                              //pcout<<" Sono in 7 "<<std::endl;
                              children[6]->AddNode(blockNodeList[i]);
                            }
                        }
                    } //fine assegnazione nodi del padre ai blocchi figli

                } //fine loop nodi del blocco

              typename std::map <cell_it, std::vector<unsigned int> >::iterator it;
              for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
                {
                  for (unsigned int pp = 0; pp < (*it).second.size(); pp++)
                    {
                      Point<dim> quadPoint = quadPoints[(*it).first][(*it).second[pp]];
                      // assegnamento punti quadratura del blocco padre ai blocchi figli
                      if (quadPoint(2) <= parent->GetPMin()(2)+delta)
                        {
                          if (quadPoint(1) <= parent->GetPMin()(1)+delta)
                            {
                              if (quadPoint(0) <= parent->GetPMin()(0)+delta)
                                {
                                  //pcout<<" Sono in 1 "<<std::endl;
                                  children[0]->AddQuadPoint((*it).first,(*it).second[pp]);
                                }
                              else
                                {
                                  //pcout<<" Sono in 2 "<<std::endl;
                                  children[1]->AddQuadPoint((*it).first,(*it).second[pp]);
                                }
                            }
                          else
                            {
                              if (quadPoint(0) <= parent->GetPMin()(0)+delta)
                                {
                                  //pcout<<" Sono in 4 "<<std::endl;
                                  children[3]->AddQuadPoint((*it).first,(*it).second[pp]);
                                }
                              else
                                {
                                  //pcout<<" Sono in 3 "<<std::endl;
                                  children[2]->AddQuadPoint((*it).first,(*it).second[pp]);
                                }
                            }
                        }
                      else
                        {
                          if (quadPoint(1) <= parent->GetPMin()(1)+delta)
                            {
                              if (quadPoint(0) <= parent->GetPMin()(0)+delta)
                                {
                                  //pcout<<" Sono in 5 "<<std::endl;
                                  children[4]->AddQuadPoint((*it).first,(*it).second[pp]);
                                }
                              else
                                {
                                  //pcout<<" Sono in 6 "<<std::endl;
                                  children[5]->AddQuadPoint((*it).first,(*it).second[pp]);
                                }
                            }
                          else
                            {
                              if (quadPoint(0) <= parent->GetPMin()(0)+delta)
                                {
                                  //pcout<<" Sono in 8 "<<std::endl;
                                  children[7]->AddQuadPoint((*it).first,(*it).second[pp]);
                                }
                              else
                                {
                                  //pcout<<" Sono in 7 "<<std::endl;
                                  children[6]->AddQuadPoint((*it).first,(*it).second[pp]);
                                }
                            }
                        } //fine assegnazione punti quadratura del padre ai blocchi figli
                    }
                } //fine loop punti quadratura del blocco

              for (unsigned int j=0; j < num_children_per_block; j++ )
                {
                  if (children[j]->GetBlockNodeList().size() +
                      children[j]->GetBlockQuadPointsList().size()> 0)
                    {
                      blocksCount += 1;
                      blocks[blocksCount] = new OctreeBlock<dim>();
                      blocks[blocksCount]->CopyContent(children[j]);
                      delete children[j];

                      parent->AddChild(blocksCount);
                      std::map <cell_it, std::vector<unsigned int> >
                      blockQuadPointsList = blocks[blocksCount]->GetBlockQuadPointsList();
                      typename std::map <cell_it, std::vector<unsigned int> >::iterator it;
                      for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
                        {
                          cell_it cell = (*it).first;
                          for (unsigned int kk = 0; kk < (*it).second.size(); kk++)
                            {
                              quad_point_to_block[(*it).first][(*it).second[kk]].push_back(blocksCount);
                            }
                        }
                      std::vector<unsigned int> blockNodesList = blocks[jj]->GetBlockNodeList();
                      for (unsigned int k = 0; k < blockNodesList.size(); k++)
                        dof_to_block[blockNodesList[k]].push_back(jj);

                    }
                  else
                    {
                      delete children[j];
                    }

                } // fine loop sui blocchi figlio appena creati

            } //fine ramo dim = 3 dell'if
          else
            {
              for (unsigned int i = 0; i < blockNodeList.size(); i++)
                {

                  // assegnamento nodi del blocco padre ai blocchi figli
                  Point <dim> node = support_points[blockNodeList[i]];

                  if (node(1) <= parent->GetPMin()(1)+delta)
                    {
                      if (node(0) <= parent->GetPMin()(0)+delta)
                        {
                          //pcout<<" Sono in 1 "<<std::endl;
                          children[0]->AddNode(blockNodeList[i]);
                        }
                      else
                        {
                          //pcout<<" Sono in 2 "<<std::endl;
                          children[1]->AddNode(blockNodeList[i]);
                        }
                    }
                  else
                    {
                      if (node(0) <= parent->GetPMin()(0)+delta)
                        {
                          //pcout<<" Sono in 4 "<<std::endl;
                          children[3]->AddNode(blockNodeList[i]);
                        }
                      else
                        {
                          //pcout<<" Sono in 3 "<<std::endl;
                          children[2]->AddNode(blockNodeList[i]);
                        }
                    }//fine assegnazione blocchi del padre ai blocchi figli

                } //fine loop nodi del blocco

              typename std::map <cell_it, std::vector<unsigned int> >::iterator it;
              for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
                {
                  for (unsigned int pp = 0; pp < (*it).second.size(); pp++)
                    {
                      // assegnamento quad points del blocco padre ai blocchi figli
                      Point<dim> quadPoint = quadPoints[(*it).first][(*it).second[pp]];
                      if (quadPoint(1) <= parent->GetPMin()(1)+delta)
                        {
                          if (quadPoint(0) <= parent->GetPMin()(0)+delta)
                            {
                              //pcout<<" Sono in 1 "<<std::endl;
                              children[0]->AddQuadPoint((*it).first,(*it).second[pp]);
                            }
                          else
                            {
                              //pcout<<" Sono in 2 "<<std::endl;
                              children[1]->AddQuadPoint((*it).first,(*it).second[pp]);
                            }
                        }
                      else
                        {
                          if (quadPoint(0) <= parent->GetPMin()(0)+delta)
                            {
                              //pcout<<" Sono in 4 "<<std::endl;
                              children[3]->AddQuadPoint((*it).first,(*it).second[pp]);
                            }
                          else
                            {
                              //pcout<<" Sono in 3 "<<std::endl;
                              children[2]->AddQuadPoint((*it).first,(*it).second[pp]);
                            }
                        }//fine assegnazione blocchi del padre ai blocchi figli

                    }
                }

              for (unsigned int j=0; j < num_children_per_block; j++ )
                {
                  if (children[j]->GetBlockNodeList().size() +
                      children[j]->GetBlockQuadPointsList().size()> 0)
                    {
                      blocksCount += 1;
                      blocks[blocksCount] = new OctreeBlock<dim>();
                      blocks[blocksCount]->CopyContent(children[j]);
                      delete children[j];

                      parent->AddChild(blocksCount);
                      std::map <cell_it, std::vector<unsigned int> >
                      blockQuadPointsList = blocks[blocksCount]->GetBlockQuadPointsList();
                      typename std::map <cell_it, std::vector<unsigned int> >::iterator it;
                      for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
                        {
                          cell_it cell = (*it).first;
                          for (unsigned int kk = 0; kk < (*it).second.size(); kk++)
                            {
                              quad_point_to_block[(*it).first][(*it).second[kk]].push_back(blocksCount);
                            }
                        }
                      std::vector<unsigned int> blockNodesList = blocks[jj]->GetBlockNodeList();
                      for (unsigned int k = 0; k < blockNodesList.size(); k++)
                        dof_to_block[blockNodesList[k]].push_back(jj);

                    }
                  else
                    {
                      delete children[j];
                    }

                } // fine loop sui blocchi figlio appena creati


            } // fine ramo dim == 2 dell'if

        } //fine loop blocchi livello precedente


      //double elemCheck = numChildless;

      startLevel[level] = endLevel[level-1] + 1;
      endLevel[level] = blocksCount;

      // here we loop over the blocks of the newly created level and
      // we decide if each block is to be split again in the next level:
      // if it contains more
      // than a node or quad point, it will be placed in the parent list.
      // Instead, if it only contains a node or quad point, it goes in the
      // childless list, and not be refined any more. It is important to
      // account for the presence of double nodes: if not, blocks will be
      // always refined
      for (unsigned int jj = startLevel[level]; jj < endLevel[level]+1;  jj++)
        {

          // here we get the number of nodes in the block
          std::vector<unsigned int> nodesId = blocks[jj]->GetBlockNodeList();
          double blockNumNodes = 0.0;

          // now we compute the number of the nodes that are double of others
          for (unsigned int kk = 0; kk < nodesId.size();  kk++)
            {
              blockNumNodes += 1.0/(double((*double_nodes_set)[nodesId[kk]].size()));
            }

          // here we compute the number of quad points in the block
          int blockNumQuadPoints = 0;
          std::map <cell_it, std::vector<unsigned int> >
          blockQuadPointsList = blocks[jj]->GetBlockQuadPointsList();
          typename std::map <cell_it, std::vector<unsigned int> >::iterator it;
          for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
            blockNumQuadPoints += (int) (*it).second.size();
          //pcout<<"Level "<<level<<" Block "<<jj<<"  nodes "<<blockNumNodes<<" + quadPoints "<<blockNumQuadPoints<<std::endl;

          quadPointsCheck += blockNumQuadPoints;
          nodesCheck += blockNumNodes;
          // here we decide if a block is to be placed in the parent
          // or childless list
          //if (blockNumNodes + blockNumQuadPoints - numDoubleNodes < 2)
          if (round(blockNumNodes) <= max_num_nodes_per_block)
            {
              numChildless += 1;
              childlessList.push_back(jj);
              quadPointsInChildless += blockNumQuadPoints;
              nodesInChildless += blockNumNodes;


              // if a block is childless, we must assign now the nodes and quad points
              // that belong to it for all the next levels

              for (unsigned int kk = 0; kk < nodesId.size();  kk++)
                for (unsigned int j = level+1; j < num_octree_levels+1; j++)
                  dof_to_block[nodesId[kk]].push_back(jj);

              for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
                for (unsigned int i = 0; i < (*it).second.size(); i++)
                  for (unsigned int j = level+1; j < num_octree_levels+1; j++)
                    quad_point_to_block[(*it).first][(*it).second[i]].push_back(jj);

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
          //elemCheck += blockNumNodes + blockNumQuadPoints;
        }


      pcout<<" Total nodes at level "<<level<<" of "<<num_octree_levels<<" are "<<nodesCheck<<std::endl;
      pcout<<" Total quad points at level "<<level<<" of "<<num_octree_levels<<" are "<<quadPointsCheck<<std::endl;
      pcout<<" Blocks at level "<<level<<" of "<<num_octree_levels<<" are "<<endLevel[level]-endLevel[level-1]<<std::endl;
      pcout<<" Total blocks at level "<<level<<" of "<<num_octree_levels<<" are "<<endLevel[level] + 1<<std::endl;
      pcout<<std::endl;

    } //fine loop livelli

  childlessList.resize(childlessList.size()+parentList[num_octree_levels].size());

  for (unsigned int jj = 0; jj < parentList[num_octree_levels].size();  jj++)
    {
      childlessList[numChildless + jj] = parentList[num_octree_levels][jj];
    }




  num_blocks = blocksCount+1;


  pcout<<"...done generating octree blocking"<<std::endl;

  pcout<<"Computing proximity lists for blocks"<<std::endl;

//just for output
  /*for (cell = fma_dh->begin_active(); cell != endc; ++cell)
      {
      unsigned int levelCheck = elem_to_blocks[cell].size();
      pcout<<std::endl<<"Elem "<<cell<<" is in the "<<levelCheck<<" blocks: ";
      for (unsigned int zz = 0; zz < levelCheck; zz++)
           pcout<<elem_to_blocks[cell][zz]<<" ";
      }*/

  /*for (cell_it cell = fma_dh->begin_active(); cell != endc; ++cell)
      for (unsigned int j=0; j < quadPoints[cell].size(); j++)
          pcout<<"Cell "<<cell<<"  QP "<<j<<"  of "<<quadPoints[cell].size()<<": "<<quadPoints[cell][j]<<std::endl;//*/

  /*for (cell_it cell = fma_dh->begin_active(); cell != endc; ++cell)
      for (unsigned int j=0; j < quad_point_to_block[cell].size(); j++)
          {
    pcout<<"Cell "<<cell<<"  QP "<<j<<"  of "<<quad_point_to_block[cell].size()<<": ";
          for (unsigned int i=0; i < quad_point_to_block[cell][j].size(); i++)
              pcout<<quad_point_to_block[cell][j][i]<<" ";
    pcout<<std::endl;
    }  //*/


  /*for (unsigned int i=0; i < fma_dh->n_dofs(); i++)
      {
      pcout<<"Node "<<i<<"  doubles: ";
      std::set <unsigned int> doubleNodes = double_nodes_set[i];
      for (std::set<unsigned int>::iterator pos = doubleNodes.begin(); pos != doubleNodes.end(); pos++)
          {
    pcout<<*pos<<"( ";
    for (unsigned int j=0; j < dof_to_elems[*pos].size(); j++)
        pcout<<" "<<dof_to_elems[*pos][j];
    pcout<<") ";
          }
      pcout<<std::endl;
      } //*/



// ricerca blocchi nearest neighbors

  for (unsigned int ii = startLevel[1]; ii < endLevel[1] + 1;  ii++)
    {
      for (unsigned int jj = startLevel[1]; jj < endLevel[1] + 1;  jj++)
        {
          blocks[ii]->AddNearNeigh(0,jj);
        }
    }


  for (unsigned int level = 2; level < num_octree_levels + 1;  level++)

    {
      for (unsigned int kk = startLevel[level]; kk < endLevel[level]+1;  kk++)

        {
          OctreeBlock<dim> *block1 = blocks[kk];
          block1->AddNearNeigh(0,kk); // a block is NearNeigh of itself

          double delta1 = block1->GetDelta();
          Point<dim> PMin1 = block1->GetPMin();
          Point<dim> Center1;
          for (unsigned int iii = 0; iii < dim; iii++)
            Center1(iii) = delta1/2.;
          Point<dim> PMax1 = 2.*Center1;
          PMax1 += PMin1;
          Center1 += PMin1;
          unsigned int parentId = block1->GetParentId();
          std::set <unsigned int> parentNNeighs = blocks[parentId]->GetNearNeighs(0);

          // the nearest neighbors are searched among the father's nearest neighbors children
          for (std::set <unsigned int>::iterator pos = parentNNeighs.begin(); pos != parentNNeighs.end();  pos++)

            {
              if (blocks[*pos]->GetBlockChildrenNum() == 0) // if a parent's near neigh is childless, he can be a near neigh: let's check
                {
                  unsigned int block2Id = *pos;
                  OctreeBlock<dim> *block2 = blocks[block2Id];
                  double delta2 = block2->GetDelta();
                  Point<dim> PMin2 = block2->GetPMin();
                  Point<dim> Center2;
                  for (unsigned int iii = 0; iii < dim; iii++)
                    Center2(iii) = delta2/2.;
                  Point<dim> PMax2 = 2.*Center2;
                  PMax2 += PMin2;
                  Center2 += PMin2;

                  if (dim == 3)
                    {
                      if  ((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
                        {
                          if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
                            {
                              if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
                                {
                                  block1->AddNearNeigh(0,block2Id);
                                  //pcout<<" *"<<block2Id;
                                }
                            }
                        }

                      if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
                        {
                          if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
                            {
                              if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
                                {
                                  block1->AddNearNeigh(0,block2Id);
                                  //pcout<<" *"<<block2Id;
                                }
                            }
                        }

                      if ((fabs(PMin1(2) - PMax2(2)) <= TOLL) || (fabs(PMax1(2) - PMin2(2)) <= TOLL))
                        {
                          if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
                            {
                              if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
                                {
                                  block1->AddNearNeigh(0,block2Id);
                                  //pcout<<" *"<<block2Id;
                                }
                            }
                        }
                    } //fine caso dim ==3

                  else if (dim == 2)
                    {
                      if  ((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
                        {
                          if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
                            {
                              block1->AddNearNeigh(0,block2Id);
                              //pcout<<block2Id<<" ";
                            }
                        }

                      if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
                        {
                          if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
                            {
                              block1->AddNearNeigh(0,block2Id);
                              //pcout<<block2Id<<" ";
                            }
                        }

                    } // fine caso dim == 2

                }

              for (unsigned int ii = 0; ii < blocks[*pos]->GetBlockChildrenNum();  ii++)
                {
                  unsigned int block2Id = blocks[*pos]->GetChildId(ii);
                  OctreeBlock<dim> *block2 = blocks[block2Id];
                  double delta2 = block2->GetDelta();
                  Point<dim> PMin2 = block2->GetPMin();
                  Point<dim> Center2;
                  for (unsigned int iii = 0; iii < dim; iii++)
                    Center2(iii) = delta2/2.;
                  Point<dim> PMax2 = 2.*Center2;
                  PMax2 += PMin2;
                  Center2 += PMin2;

                  if (dim == 3)
                    {
                      if  ((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
                        {
                          if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
                            {
                              if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
                                {
                                  block1->AddNearNeigh(0,block2Id);
                                  //pcout<<" "<<block2Id;
                                }
                            }
                        }

                      if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
                        {
                          if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
                            {
                              if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
                                {
                                  block1->AddNearNeigh(0,block2Id);
                                  //pcout<<" "<<block2Id;
                                }
                            }
                        }

                      if ((fabs(PMin1(2) - PMax2(2)) <= TOLL) || (fabs(PMax1(2) - PMin2(2)) <= TOLL))
                        {
                          if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
                            {
                              if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
                                {
                                  block1->AddNearNeigh(0,block2Id);
                                  //pcout<<" "<<block2Id;
                                }
                            }
                        }
                    } //fine caso dim ==3

                  else if (dim == 2)
                    {
                      if  ((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
                        {
                          if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
                            {
                              block1->AddNearNeigh(0,block2Id);
                              //pcout<<block2Id<<" ";
                            }
                        }

                      if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
                        {
                          if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
                            {
                              block1->AddNearNeigh(0,block2Id);
                              //pcout<<block2Id<<" ";
                            }
                        }

                    } // fine caso dim == 2

                } // fine loop sui figli di un nearest neighbor del padre


            } // fine loop sui nearest neighbors del padre



          if ((block1->GetBlockChildrenNum() == 0))  // if the block is childless we must compute now its nearneigh at all residual levels
            {
              block1->SetNearNeighSize(num_octree_levels-level+1);
              block1->SetIntListSize(num_octree_levels-level+1); // intList is a vector of sets with the same number of members of nearNeigh
              block1->SetNonIntListSize(num_octree_levels-level+1); // nonIntList is a vector of sets with the same number of members of nearNeigh

              for (unsigned int subLevel = 1; subLevel < num_octree_levels - level + 1; subLevel++)

                {

                  std::set <unsigned int> upperLevelNNeighs = block1->GetNearNeighs(subLevel-1);
                  for (std::set <unsigned int>::iterator pos = upperLevelNNeighs.begin(); pos != upperLevelNNeighs.end();  pos++)

                    {
                      if (blocks[*pos]->GetBlockChildrenNum() == 0) // if nearneigh is childless, it will stay a near neigh
                        block1->AddNearNeigh(subLevel,*pos);

                      for (unsigned int ii = 0; ii < blocks[*pos]->GetBlockChildrenNum();  ii++)
                        {
                          unsigned int block2Id = blocks[*pos]->GetChildId(ii);
                          OctreeBlock<dim> *block2 = blocks[block2Id];
                          double delta2 = block2->GetDelta();
                          Point<dim> PMin2 = block2->GetPMin();
                          Point<dim> Center2;
                          for (unsigned int iii = 0; iii < dim; iii++)
                            Center2(iii) = delta2/2.;
                          Point<dim> PMax2 = 2.*Center2;
                          PMax2 += PMin2;
                          Center2 += PMin2;

                          if (dim == 3)
                            {
                              if  ((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
                                {
                                  if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
                                    {
                                      if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
                                        {
                                          block1->AddNearNeigh(subLevel,block2Id);
                                          //pcout<<block2Id<<" ";
                                        }
                                    }
                                }

                              if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
                                {
                                  if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
                                    {
                                      if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
                                        {
                                          block1->AddNearNeigh(subLevel,block2Id);
                                          //pcout<<block2Id<<" ";
                                        }
                                    }
                                }

                              if ((fabs(PMin1(2) - PMax2(2)) <= TOLL) || (fabs(PMax1(2) - PMin2(2)) <= TOLL))
                                {
                                  if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
                                    {
                                      if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
                                        {
                                          block1->AddNearNeigh(subLevel,block2Id);
                                          //pcout<<block2Id<<" ";
                                        }
                                    }
                                }
                            } //fine caso dim ==3

                          else if (dim == 2)
                            {
                              if  ((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
                                {
                                  if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
                                    {
                                      block1->AddNearNeigh(subLevel,block2Id);
                                      //pcout<<block2Id<<" ";
                                    }
                                }

                              if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
                                {
                                  if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
                                    {
                                      block1->AddNearNeigh(subLevel,block2Id);
                                      //pcout<<block2Id<<" ";
                                    }
                                }

                            } // fine caso dim == 2

                        } // fine loop sui figli di ciascun nearest neighbor del blocco childless


                    } // fine loop sui nearest neighbors del blocco childless


                } // fine loop sui subLevels (da quello del blocco childless all'ultimo)


            } // fine if (il blocco e' childless?)



        } // fine loop sui blocchi di un livello

    } // fine loop sui livelli

//for printout

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
// and for non interaction list blocks (nonIntList for a block B is composed by blocks that are children
// of blocks being in intList of B's parent, but are not in intList of B)

  for (unsigned int ii = startLevel[1]; ii < endLevel[1] + 1;  ii++) // at level 1, add all blocks to intList
    {
      for (unsigned int jj = startLevel[1]; jj < endLevel[1] + 1;  jj++)
        {
          blocks[ii]->AddBlockToIntList(0,jj);
        }
    }


  for (unsigned int level = 2; level < num_octree_levels + 1;  level++) // loop over levels

    {
      for (unsigned int jj = startLevel[level]; jj < endLevel[level] + 1;  jj++) // loop over blocks of each level
        {
          OctreeBlock<dim> *block1 = blocks[jj];
          //pcout<<"??Out "<<jj<<" "<<block1.GetNearNeighSize()<<std::endl;
          for (unsigned int subLevel = 0; subLevel < block1->NumNearNeighLevels(); subLevel++)
            {
              std::set <unsigned int> NNList = block1->GetNearNeighs(subLevel);

              for (std::set <unsigned int>::iterator pos1 = NNList.begin(); pos1 != NNList.end();  pos1++) //loop over blocks in NN list and get their NNs
                {
                  block1->AddBlockToIntList(subLevel,*pos1);
                }
              //pcout<<std::endl<<"Sublevel "<<subLevel<<" elem("<<block1.GetBlockElementsList()[0]<<") NearNeighs: ";
              std::vector <unsigned int> nodeIds = block1->GetBlockNodeList();
              //pcout<<std::endl<<"Level "<<level<<"  Block1: "<<kk<<"  NumElements: "<<block1.GetBlockElementsNum()<<" "<<elemIds.size()<<std::endl;
              //pcout<<"Nearest Neighbors Found:"<<std::endl;
              for (unsigned int pp = 0; pp < nodeIds.size();  pp++)
                {
                  //pcout<<"Node "<<nodeIds[pp]<<std::endl;
                  std::set <unsigned int> doubleNodes = (*double_nodes_set)[nodeIds[pp]];
                  for (std::set<unsigned int>::iterator pos = doubleNodes.begin();
                       pos != doubleNodes.end(); pos++)
                    {
                      //pcout<<"Node Double"<<*pos<<std::endl;
                      //std::vector<cell_it > surrCellIds = dof_to_elems[*pos];
                      for (unsigned int k=0; k < dof_to_elems[*pos].size(); k++)
                        {
                          cell_it cell = dof_to_elems[*pos][k];
                          //pcout<<cell<<std::endl;
                          for (unsigned int j=0; j < quadPoints[cell].size(); j++)
                            {
                              block1->AddBlockToIntList(subLevel,quad_point_to_block[cell][j][level+subLevel]);
                            }
                        }
                    }
                }
            }
          for (unsigned int subLevel = 0; subLevel < block1->GetNearNeighSize();  subLevel++) // for each block, loop over all sublevels in his NN list (to account for childless blocks)
            {
              // now use intList to compute nonIntList
              std::set <unsigned int> intList = block1->GetIntList(subLevel);
              std::set <unsigned int> parentIntList; // intList at the  previous level
              if (subLevel == 0) // if a block is childless we get its intList at the previous level, otherwise we get its parent's intList
                parentIntList = blocks[block1->GetParentId()]->GetIntList(0);
              else
                parentIntList = block1->GetIntList(subLevel-1);

              for (std::set <unsigned int>::iterator pos1 = parentIntList.begin(); pos1 != parentIntList.end();  pos1++) // loop over blocks in parentIntList
                {
                  OctreeBlock<dim> *block2 = blocks[*pos1];
                  if (block2->GetBlockChildrenNum() == 0) // if blocks in parentIntList are childless, don't look for their children, but see if they are in nonIntList
                    {
                      if (intList.count(*pos1) == 0) // if these blocks are not in intList
                        block1->AddBlockToNonIntList(subLevel,*pos1); // then they go in nonIntList
                    }
                  else // if blocks in parentIntList are not childless, do the same test on all their children
                    {
                      for (unsigned int kk = 0; kk < block2->GetBlockChildrenNum();  kk++) // loop over children of blocks in parentIntList
                        {
                          //pcout<<"Sublevel "<<subLevel<<" Block1 "<<jj<<"  Block2 "<<*pos1<<" child(kk) "<<block2.GetChildId(kk)<<std::endl;
                          if (intList.count(block2->GetChildId(kk)) == 0)   // if these blocks are not in intList
                            block1->AddBlockToNonIntList(subLevel,block2->GetChildId(kk)); // then they go in nonIntList
                        } // end loop over children of blocks in parentIntList
                    }
                } // loop over blocks in parentIntList

            } // end loop over subLevels of each block's intList

//      for printout

          /*
                pcout<<std::endl;
                                  pcout<<"-------------------------------- "<<std::endl;
                pcout<<"Block "<<jj<<": Parent "<<block1->GetParentId();
                pcout<<"  cubePlot(["<<block1->GetPMin()<<"],"<<block1->GetDelta()<<",'m')"<<std::endl;
                pcout<<"Quad points of block "<<jj<<": "<<std::endl;
                std::map <cell_it,std::vector<unsigned int> > quadPointsList = block1->GetBlockQuadPointsList();
                typename std::map <cell_it, std::vector<unsigned int> >::iterator it;
                                  for (it = quadPointsList.begin(); it != quadPointsList.end(); it++)
                                {
                    pcout<<(*it).first<<"( ";
                    for (unsigned int zz = 0; zz < (*it).second.size();  zz++)
                        pcout<<(*it).second[zz]<<" ";
                    pcout<<")  ";
                    }
                pcout<<std::endl;

                pcout<<"Nodes of block "<<jj<<" :"<<std::endl;
                std::vector <unsigned int> nodeList = block1->GetBlockNodeList();
                for (unsigned int zz = 0; zz < block1->GetBlockNodesNum();  zz++)
                  {
                  pcout<<nodeList.at(zz)<<" ";
                  }
                pcout<<std::endl;

                for (unsigned int subLevel = 0; subLevel < block1->GetNearNeighSize();  subLevel++)
                  {
                  std::set <unsigned int> NNList = block1->GetNearNeighs(subLevel);
                  pcout<<"NearNeigh for block "<<jj<<" at level "<<level+subLevel<<":"<<std::endl;
                  for (std::set <unsigned int>::iterator pos1 = NNList.begin(); pos1 != NNList.end();  pos1++)
                    {
                    pcout<<*pos1<<" ";
                    }
                  pcout<<std::endl;

                  std::set <unsigned int> intList = block1->GetIntList(subLevel);
                  pcout<<"IntList for block "<<jj<<" at level "<<level+subLevel<<":"<<std::endl;
                  for (std::set <unsigned int>::iterator pos1 = intList.begin(); pos1 != intList.end();  pos1++)
                    {
                    pcout<<*pos1<<" ";
                    }
                  pcout<<std::endl;

                  std::set <unsigned int> nonIntList = block1->GetNonIntList(subLevel);
                  pcout<<"NonIntList for block "<<jj<<" at level "<<level+subLevel<<":"<<std::endl;
                  for (std::set <unsigned int>::iterator pos1 = nonIntList.begin(); pos1 != nonIntList.end();  pos1++)
                    {
                    pcout<<*pos1<<" ";
                    }
                  pcout<<std::endl;
                  }
          //*/



        } // end loop over blocks of a level

    } // end loop over levels

//if (blocks.size() > 0)
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



  pcout<<"Done computing proximity lists for blocks"<<std::endl;
} //end method for octree blocking generation

template class BEMFMA<3>;
