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


template <int dim>
BEMFMA<dim>::BEMFMA(ComputationalDomain<dim> &comp_dom)
		:
		comp_dom(comp_dom),
                mpi_communicator (MPI_COMM_WORLD),
                n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
                this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator))
{}


template <int dim> 
void BEMFMA<dim>::declare_parameters (ParameterHandler &prm)
{
		    
  prm.enter_subsection("FMA Params");
  {    
   prm.declare_entry("FMA Truncation Order", "3", Patterns::Integer());
  }
  prm.leave_subsection();
  
}

template <int dim> 
void BEMFMA<dim>::parse_parameters (ParameterHandler &prm)
{

  prm.enter_subsection("FMA Params");
  {
   trunc_order = prm.get_integer("FMA Truncation Order");
  }
  prm.leave_subsection();
  
}



template <int dim>
void BEMFMA<dim>::direct_integrals()
{
std::cout<<"Computing direct integrals..."<<std::endl;

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

  
  std::vector<QGaussOneOverR<2> > sing_quadratures_3d; 
  for(unsigned int i=0; i<comp_dom.fe.dofs_per_cell; ++i)
    sing_quadratures_3d.push_back
      (QGaussOneOverR<2>(comp_dom.singular_quadrature_order,
			 comp_dom.fe.get_unit_support_points()[i], true));
                                   // number of dofs per cell
				   
  const unsigned int dofs_per_cell = comp_dom.fe.dofs_per_cell;
  
                                   // vector containing the ids of the dofs
				   // of each cell: it will be used to transfer
				   // the computed local rows of the matrices 
				   // into the global matrices

  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

                                   // vector to store parts of rows of neumann
				   // and dirichlet matrix obtained in local
				   // operations

  Vector<double>      local_neumann_matrix_row_i(comp_dom.fe.dofs_per_cell);
  Vector<double>      local_dirichlet_matrix_row_i(comp_dom.fe.dofs_per_cell);
    
    
				   // Now that we have checked that
				   // the number of vertices is equal
				   // to the number of degrees of
				   // freedom, we construct a vector
				   // of support points which will be
				   // used in the local integrations:
				   
  std::vector<Point<dim> > support_points(comp_dom.dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>(*comp_dom.mapping, comp_dom.dh, support_points);

  
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
  cell = comp_dom.dh.begin_active(),
  endc = comp_dom.dh.end();
    
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

  init_prec_sparsity_pattern.reinit(comp_dom.dh.n_dofs(),comp_dom.dh.n_dofs(),125*comp_dom.fe.dofs_per_cell);
  
  for (unsigned int kk = 0; kk < comp_dom.childlessList.size(); kk++)

	{
		 	            // for each block in the childless
				    // list we get the list of nodes and
				    // we check if it contains nodes:
				    // if no nodes are contained there is
				    // nothing to do
				     
	 unsigned int blockId = comp_dom.childlessList[kk];
	 
	 OctreeBlock<dim>* block1 = comp_dom.blocks[blockId];
	 
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
		 OctreeBlock<dim>* block2 = comp_dom.blocks[*pos];
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
		   
		   cell_it cell = (*it).first;//std::cout<<cell<<"  end "<<(*blockQuadPointsList.end()).first<<std::endl;
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
	         for (std::set<unsigned int>::iterator pos = directNodes.begin(); pos != directNodes.end(); pos++)
	             init_prec_sparsity_pattern.add(block1Nodes[i],*pos);
             }
	 
	}

	// unfortunately, the direct integrals must not be computed only for the
	// quadPoints in the intList: if a bigger block is in the nonIntList of
	// another block, the bound for the multipole expansion application does
	// not hold, and so we must compute direct integrals. Here we scan the
	// nonIntlists of each block at each level to look for bigger blocks and  
        // initialize the prec matrices sparsity pattern with the corresponding nodes
	 
    for (unsigned int level = 1; level < comp_dom.num_octree_levels + 1;  level++) // loop over levels

	{           
	std::vector<unsigned int>
	dofs_filled_blocks = comp_dom.dofs_filled_blocks[level];
	unsigned int startBlockLevel = comp_dom.startLevel[level];
	            
		    // we loop over blocks of each level
	
	for (unsigned int jj = 0; jj < dofs_filled_blocks.size();  jj++) 
		{
                
		OctreeBlock<dim> *block1 = comp_dom.blocks[dofs_filled_blocks[jj]];
		const std::vector <unsigned int> &nodesBlk1Ids = block1->GetBlockNodeList();
                
		    // again, no need to perform next operations if block has no nodes
		
		if  (nodesBlk1Ids.size() > 0)
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
		                OctreeBlock<dim>* block2 = comp_dom.blocks[*pos];
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
	                        for (std::set<unsigned int>::iterator pos = directNodes.begin(); pos != directNodes.end(); pos++)
	                            init_prec_sparsity_pattern.add(nodesBlk1Ids[i],*pos);	
			
			    } // end loop over sublevels
		    } // end if: is there any node in the block? 	
		}// end loop over block of a level
	}//end loop over octree levels
	
  
         // sparsity pattern is ready and can be compressed; the direct matrices
	 // and the preconditioner one are the initialized with the sparsity
	 // pattern just computed
	 
  init_prec_sparsity_pattern.compress();
  double filling_percentage = double(init_prec_sparsity_pattern.n_nonzero_elements())/double(comp_dom.dh.n_dofs()*comp_dom.dh.n_dofs())*100.;     
  std::cout<<init_prec_sparsity_pattern.n_nonzero_elements()<<" Nonzeros out of "<<comp_dom.dh.n_dofs()*comp_dom.dh.n_dofs()<<":  "<<filling_percentage<<"%"<<std::endl;
  
  prec_neumann_matrix.reinit(init_prec_sparsity_pattern);
  prec_dirichlet_matrix.reinit(init_prec_sparsity_pattern);
  init_preconditioner.reinit(init_prec_sparsity_pattern);
  
  Point<dim> D;
  double s;
  
         // here we finally start computing the direct integrals: we
	 // first loop among the childless blocks 
  
  for (unsigned int kk = 0; kk < comp_dom.childlessList.size(); kk++)

	{
        //std::cout<<"processing block "<<kk <<"  of  "<<cMesh->GetNumChildlessBlocks()<<std::endl;
        //std::cout<<"block "<<cMesh->GetChildlessBlockId(kk) <<"  of  "<<cMesh->GetNumBlocks()<<"  in block list"<<std::endl;
	 
	 // this is the Id of the block	
	 unsigned int blockId = comp_dom.childlessList[kk];
	 // and this is the block pointer
	 OctreeBlock<dim>* block1 = comp_dom.blocks[blockId];
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
		 OctreeBlock<dim>* block2 = comp_dom.blocks[*pos];
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
		           comp_dom.integralCheck[block1Nodes[kk]][(*it).first] += 1;
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
            
	             for (unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j) 
	                if (comp_dom.double_nodes_set[nodeIndex].count(local_dof_indices[j]) > 0)
	                   {
		           singular_index = j;
		           is_singular = true;
		           break;
	                   }
		     // first case: the current node does not belong to the current cell:
		     // we use regular quadrature
		     if (is_singular == false)
	                { 
	                //std::cout<<"Node "<<i<<"  Elem "<<cell<<" (Direct) Nodes: ";
	                //for(unsigned int j=0; j<fe.dofs_per_cell; ++j) std::cout<<" "<<local_dof_indices[j];
	                //std::cout<<std::endl;
	                
			// we start looping on the quad points of the cell: *pos will be the
			// index of the quad point
	                for (std::set<unsigned int>::iterator pos=cellQuadPoints.begin(); pos!=cellQuadPoints.end(); pos++)
		            {
			    // here we compute the distance R between the node and the quad point
			    const Tensor<1, dim> R = comp_dom.quadPoints[cell][*pos] - support_points[nodeIndex];
                            LaplaceKernel::kernels(R, D, s);
			    
			    // and here are the integrals for each of the degrees of freedom of the cell: note
			    // how the quadrature values (position, normals, jacobianXweight, shape functions) 
			    // are taken from the precomputed ones in ComputationalDomain class
		            for(unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j)
			       {
		               local_neumann_matrix_row_i(j) += ( ( D * 
							            comp_dom.quadNormals[cell][*pos] ) *
						                    comp_dom.quadShapeFunValues[cell][*pos][j] *
						                    comp_dom.quadJxW[cell][*pos] );
		               local_dirichlet_matrix_row_i(j) += ( s * 
							            comp_dom.quadShapeFunValues[cell][*pos][j] *
							            comp_dom.quadJxW[cell][*pos] );
				//std::cout<<D<<" "<<comp_dom.quadNormals[cell][*pos]<<" ";
				//std::cout<<comp_dom.quadShapeFunValues[cell][*pos][j]<<" ";
				//std::cout<<comp_dom.quadJxW[cell][*pos]<<std::endl;			 
		               }
		           }
	                } // end if
		     else
		        {
			// after some checks, we have to create the singular quadrature:
			// here the quadrature points of the cell will be IGNORED,
			// and the singular quadrature points are instead used.
			// the 3d and 2d quadrature rules are different
   		        Assert(singular_index != numbers::invalid_unsigned_int,
		        ExcInternalError());

	                const Quadrature<dim-1> *
	                      singular_quadrature
	                      = (dim == 2
		                 ?
		                dynamic_cast<Quadrature<dim-1>*>(
		                new QGaussLogR<1>(comp_dom.singular_quadrature_order,
				Point<1>((double)singular_index),
				1./cell->measure(), true))
		              :
		                (dim == 3
		                 ?
		                 dynamic_cast<Quadrature<dim-1>*>(
		                 &sing_quadratures_3d[singular_index])
		              :
		              0));
	                Assert(singular_quadrature, ExcInternalError());
			
			// once the singular quadrature has been created, we employ it
			// to create the corresponding fe_values 
                        
	                FEValues<dim-1,dim> fe_v_singular (*comp_dom.mapping, comp_dom.fe, *singular_quadrature, 
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
		            for (unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j)
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
					   
    	                for (unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j)
	                    {
	                    prec_neumann_matrix.add(nodeIndex,local_dof_indices[j],local_neumann_matrix_row_i(j));
	                    prec_dirichlet_matrix.add(nodeIndex,local_dof_indices[j],local_dirichlet_matrix_row_i(j));
	                    if (cell->material_id() == comp_dom.free_sur_ID1 ||
	                       cell->material_id() == comp_dom.free_sur_ID2 ||
		               cell->material_id() == comp_dom.free_sur_ID3    )
	      	               init_preconditioner.add(nodeIndex,local_dof_indices[j],-local_dirichlet_matrix_row_i(j));
	                    else
	                       init_preconditioner.add(nodeIndex,local_dof_indices[j], local_neumann_matrix_row_i(j));
	                    }

		     } // end loop on cells of the intList		 		 
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
	 
    for (unsigned int level = 1; level < comp_dom.num_octree_levels + 1;  level++) // loop over levels

	{
	unsigned int startBlockLevel = comp_dom.startLevel[level];
	for (unsigned int jj = 0; jj < comp_dom.dofs_filled_blocks[level].size();  jj++) // loop over blocks of each level
		{
             
		OctreeBlock<dim> *block1 = comp_dom.blocks[comp_dom.dofs_filled_blocks[level][jj]];
		const std::vector <unsigned int> &nodesBlk1Ids = block1->GetBlockNodeList();
                
		for (unsigned int i = 0; i < nodesBlk1Ids.size(); i++)
		    { 
		    // for each block containing nodes, loop over all sublevels in his NN list (this is because if a 
		    // block remains childless BEFORE the last level, at this point we need to compute
		    // all its contributions up to the bottom level)
		    unsigned int nodeIndex = nodesBlk1Ids[i];
		    std::map <cell_it,std::set<unsigned int> > directQuadPoints;
		    for (unsigned int subLevel = 0; subLevel < block1->NumNearNeighLevels();  subLevel++) 
			{
			const std::set <unsigned int> &nonIntList = block1->GetNonIntList(subLevel);			
			
			// loop over well separated blocks of higher size (level)-----> in this case
			//we must use direct evaluation (luckily being childless they only contain 1 element)
			for (std::set<unsigned int>::iterator pos = nonIntList.begin(); pos !=nonIntList.lower_bound(startBlockLevel); pos++)
			    {
		            OctreeBlock<dim>* block2 = comp_dom.blocks[*pos];
	                    std::map <cell_it, std::vector<unsigned int> >
	                    blockQuadPointsList = block2->GetBlockQuadPointsList();
	                    typename std::map <cell_it, std::vector<unsigned int> >::iterator it;    
                            for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
	                        {
		                for (unsigned int ii=0; ii<(*it).second.size(); ii++)
		                    {
				    directQuadPoints[(*it).first].insert((*it).second[ii]);
                                    
				    /*////////this is for a check/////////////////////
		                    comp_dom.integralCheck[nodesBlk1Ids[i]][(*it).first] += 1;
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
		          
	                //std::cout<<"Node "<<i<<"  Elem "<<cell<<" (Direct) Nodes: ";
	                //for(unsigned int j=0; j<fe.dofs_per_cell; ++j) std::cout<<" "<<local_dof_indices[j];
	                //std::cout<<std::endl;
	                
			// we start looping on the quad points of the cell: *pos will be the
			// index of the quad point
	                for (std::set<unsigned int>::iterator pos=cellQuadPoints.begin(); pos!=cellQuadPoints.end(); pos++)
		            {
			    // here we compute the distance R between the node and the quad point
			      const Tensor<1,dim> R = comp_dom.quadPoints[cell][*pos] - support_points[nodeIndex];
                            LaplaceKernel::kernels(R, D, s);
		            
			    // and here are the integrals for each of the degrees of freedom of the cell: note
			    // how the quadrature values (position, normals, jacobianXweight, shape functions) 
			    // are taken from the precomputed ones in ComputationalDomain class
			    
		            for (unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j)
			        {
		                local_neumann_matrix_row_i(j) += ( ( D * 
							             comp_dom.quadNormals[cell][*pos] ) *
						                     comp_dom.quadShapeFunValues[cell][*pos][j] *
						                     comp_dom.quadJxW[cell][*pos] );
		                local_dirichlet_matrix_row_i(j) += ( s * 
							             comp_dom.quadShapeFunValues[cell][*pos][j] *
							             comp_dom.quadJxW[cell][*pos] );
							 
		                } // end loop over the dofs in the cell
		            } // end loop over the quad points in a cell
			
			// Finally, we need to add
			// the contributions of the
			// current cell to the
			// global matrix.
					   
    	                for (unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j)
	                    {
	                    prec_neumann_matrix.add(nodeIndex,local_dof_indices[j],local_neumann_matrix_row_i(j));
	                    prec_dirichlet_matrix.add(nodeIndex,local_dof_indices[j],local_dirichlet_matrix_row_i(j));
	                    if (cell->material_id() == comp_dom.free_sur_ID1 ||
	                       cell->material_id() == comp_dom.free_sur_ID2 ||
		               cell->material_id() == comp_dom.free_sur_ID3    )
	      	               init_preconditioner.add(nodeIndex,local_dof_indices[j],-local_dirichlet_matrix_row_i(j));
	                    else
	                       init_preconditioner.add(nodeIndex,local_dof_indices[j], local_neumann_matrix_row_i(j));
	                    }
    
			    
                        } // end loop over quad points in the direct quad points list
		     } // end loop over nodes in a block 	
		}// end loop over block of a level
	}//end loop over octree levels


 
std::cout<<"...done computing direct integrals"<<std::endl;
}



template <int dim>
void BEMFMA<dim>::multipole_integrals()
{

std::cout<<"Computing multipole integrals..."<<std::endl;

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

  const unsigned int dofs_per_cell = comp_dom.fe.dofs_per_cell;
        
   
  AssertThrow(dofs_per_cell == GeometryInfo<dim-1>::vertices_per_cell,
	      ExcMessage("The code in this function can only be used for "
			 "the usual Q1 elements."));

  // now we start looping on the childless blocks to perform the integrals
  for (unsigned int kk = 0; kk < comp_dom.childlessList.size(); kk++)

	{
        //std::cout<<"processing block "<<kk <<"  of  "<<cMesh->GetNumChildlessBlocks()<<std::endl;
        //std::cout<<"block "<<cMesh->GetChildlessBlockId(kk) <<"  of  "<<cMesh->GetNumBlocks()<<"  in block list"<<std::endl;
	
	// we get the current block and its Id, and then we
	// compute its center, which is needed to construct the
	// multipole expansion in which we store the integrals
		
	unsigned int blockId = comp_dom.childlessList[kk];
	OctreeBlock<dim>* block = comp_dom.blocks[blockId];
	double delta = block->GetDelta();
	Point<dim> deltaHalf;
	for (unsigned int i=0; i<dim;i++)
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
	    for (unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j)
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
		for (unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j)
		    {
                    elemMultipoleExpansionsKer1[blockId][cell][j].AddNormDer(comp_dom.quadShapeFunValues[cell][q][j]*comp_dom.quadJxW[cell][q]/4/numbers::PI,comp_dom.quadPoints[cell][q],comp_dom.quadNormals[cell][q]);
		    elemMultipoleExpansionsKer2[blockId][cell][j].Add(comp_dom.quadShapeFunValues[cell][q][j]*comp_dom.quadJxW[cell][q]/4/numbers::PI,comp_dom.quadPoints[cell][q]);
		    }
	        } // end loop on cell quadrature points in the block
	     
	    } // end loop over cells in the block

	} // fine loop sui blocchi childless



std::cout<<"...done computing multipole integrals"<<std::endl;
}



template <int dim>
void BEMFMA<dim>::generate_multipole_expansions(const TrilinosWrappers::MPI::Vector &phi_values, const TrilinosWrappers::MPI::Vector &dphi_dn_values) const
{
std::cout<<"Generating multipole expansions..."<<std::endl;


// also here we clear the structures storing the multipole
// expansions for Dirichlet and Neumann matrices  


blockMultipoleExpansionsKer1.clear(); 
blockMultipoleExpansionsKer2.clear();

// we reisze them: there's going to be an expansion per block

blockMultipoleExpansionsKer1.resize(comp_dom.num_blocks); 
blockMultipoleExpansionsKer2.resize(comp_dom.num_blocks);

// these two variables will be handy in the following

std::vector<unsigned int> local_dof_indices(comp_dom.fe.dofs_per_cell);
double delta;

// we loop on blocks and for each of them we create an empty multipole expansion
// centered in the block center

for (unsigned int ii = 0; ii < comp_dom.num_blocks ; ii++)

    {
    delta = comp_dom.blocks[ii]->GetDelta();
    Point<dim> deltaHalf;
    for (unsigned int i=0; i<dim;i++)
	deltaHalf(i) = delta/2.;
	
    Point<dim> blockCenter = comp_dom.blocks[ii]->GetPMin()+deltaHalf;

    blockMultipoleExpansionsKer1[ii] = MultipoleExpansion(trunc_order, blockCenter, &assLegFunction);
    blockMultipoleExpansionsKer2[ii] = MultipoleExpansion(trunc_order, blockCenter, &assLegFunction);		
    }

// we now begin the rising phase of the algorithm: starting from the lowest block levels (childless blocks)
// we get all the values of the multipole integrals and aggregate them in the multipole expansion for
// each blocks 

for (unsigned int kk = 0; kk < comp_dom.childlessList.size(); kk++)

	{
	 
	 // for each block we get the center and the quad points
	 	
	 unsigned int blockId = comp_dom.childlessList[kk];
	 OctreeBlock<dim>* block = comp_dom.blocks[blockId];
	 
         delta = comp_dom.blocks[blockId]->GetDelta();
         Point<dim> deltaHalf;
         for (unsigned int i=0; i<dim;i++)
	     deltaHalf(i) = delta/2.;
         Point<dim> blockCenter = comp_dom.blocks[blockId]->GetPMin()+deltaHalf;
	 
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
	    
 	    for (unsigned int jj=0; jj < comp_dom.fe.dofs_per_cell; ++jj)
		{
		blockMultipoleExpansionsKer2.at(blockId).Add(&elemMultipoleExpansionsKer2[blockId][cell][jj],dphi_dn_values(local_dof_indices[jj]));
		blockMultipoleExpansionsKer1.at(blockId).Add(&elemMultipoleExpansionsKer1[blockId][cell][jj],phi_values(local_dof_indices[jj]));
		}	 	
	    } //end loop ond block elements	
	}	// end loop on childless blocks


// now all the lower level blocks have a multipole expansion containing the contribution to the integrals
// of all the quad points in them: we now begin summing the child multipole expansion to the the parents
// expansions: to do that we need to translate che children expansions to the parent block center: there
// is a MultipoleExpansion clas for this too  

// we loop the levels starting from the bottom one

for (unsigned int level = comp_dom.num_octree_levels; level > 0; level--)

	{
	
	// for each block we add the (translated) multipole expansion to the the parent expansion
	
	std::cout<<"processing level "<<level <<"  of  "<<comp_dom.num_octree_levels<<std::endl;

	for (unsigned int kk = comp_dom.startLevel[level]; kk < comp_dom.endLevel[level]+1; kk++)

		{
		unsigned int parentId = comp_dom.blocks[kk]->GetParentId();
		
		blockMultipoleExpansionsKer1.at(parentId).Add(&blockMultipoleExpansionsKer1.at(kk));
	 	blockMultipoleExpansionsKer2.at(parentId).Add(&blockMultipoleExpansionsKer2.at(kk));
	 	
		} // end loop over blocks of a level
			
	} // end loop over levels


std::cout<<"...done generating multipole expansions"<<std::endl;

}



template <int dim>
void BEMFMA<dim>::multipole_matr_vect_products(const TrilinosWrappers::MPI::Vector &phi_values, const TrilinosWrappers::MPI::Vector &dphi_dn_values,
		                                         TrilinosWrappers::MPI::Vector &matrVectProdN,    TrilinosWrappers::MPI::Vector &matrVectProdD) const
{
std::cout<<"Computing multipole matrix-vector products... "<<std::endl;

    // we start re-initializing matrix-vector-product vectors
    matrVectProdN = 0;
    matrVectProdD = 0;

    //and here we compute the direct integral contributions (stored in two sparse matrices)
    prec_neumann_matrix.vmult(matrVectProdN, phi_values);
    prec_dirichlet_matrix.vmult(matrVectProdD, dphi_dn_values);

    //from here on, we compute the multipole expansions contributions
    
    // we start cleaning past sessions
    blockLocalExpansionsKer1.clear(); 
    blockLocalExpansionsKer2.clear();

    blockLocalExpansionsKer1.resize(comp_dom.num_blocks); 
    blockLocalExpansionsKer2.resize(comp_dom.num_blocks);

    // we declare some familiar variables that will be useful in the method
    std::vector<Point<dim> > support_points(comp_dom.dh.n_dofs());
    DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping,
                                                  comp_dom.dh, support_points);
    std::vector<unsigned int> local_dof_indices(comp_dom.fe.dofs_per_cell);
    double delta;


    // here we loop on all the blocks and build an empty local expansion for
    // each of them
    for (unsigned int ii = 0; ii < comp_dom.num_blocks ; ii++)

        {
        delta = comp_dom.blocks[ii]->GetDelta();
        Point<dim> deltaHalf;
        for (unsigned int i=0; i<dim;i++)
	    deltaHalf(i) = delta/2.;
        Point<dim> blockCenter = comp_dom.blocks[ii]->GetPMin()+deltaHalf;
	  		
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
    for (unsigned int level = 1; level < comp_dom.num_octree_levels + 1;  level++)

	{
	std::cout<<"processing level "<<level <<"  of  "<<comp_dom.num_octree_levels<<std::endl;
	
	// we get the ids of the first and last block of the level 
	unsigned int startBlockLevel = comp_dom.startLevel[level];
	unsigned int endBlockLevel = comp_dom.endLevel[level];
	
	// to reduce computational cost, we decide to loop on the list of blocks which
	// contain at least one node (the local and multipole expansion will be finally evaluated
	// at the nodes positions)
		
	for (unsigned int k = 0; k < comp_dom.dofs_filled_blocks[level].size();  k++) // loop over blocks of each level
		{
		//std::cout<<"Block "<<jj<<std::endl;
		unsigned int jj = comp_dom.dofs_filled_blocks[level][k];
		OctreeBlock<dim> *block1 = comp_dom.blocks[jj];
		unsigned int block1Parent = block1->GetParentId();
		std::vector <unsigned int> nodesBlk1Ids = block1->GetBlockNodeList();
		

		// here we get parent contribution to local expansion and transfer it in current
		// block: this operation requires a local expansion translation, implemented
		// in a specific LocalExpansion class member
		
		blockLocalExpansionsKer1[jj].Add(&blockLocalExpansionsKer1[block1Parent]);
		blockLocalExpansionsKer2[jj].Add(&blockLocalExpansionsKer2[block1Parent]);
                
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
				//std::cout<<"NonIntListPart2 Blocks: "<<*pos1<<" ";
				unsigned int block2Id = *pos1;
				blockLocalExpansionsKer1[jj].Add(&blockMultipoleExpansionsKer1[block2Id]);
				blockLocalExpansionsKer2[jj].Add(&blockMultipoleExpansionsKer2[block2Id]);	
			        
				/*////////this is for a check///////////////////////
				OctreeBlock<dim>* block2 = comp_dom.blocks[block2Id];
				std::vector<unsigned int> nodesBlk1Ids = block1->GetBlockNodeList();
				std::map <cell_it, std::vector<unsigned int> >
	                        blockQuadPointsList = block2->GetBlockQuadPointsList();
	                        typename std::map <cell_it, std::vector<unsigned int> >::iterator it;    
                                for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
	                            {
		                    for (unsigned int i=0; i<(*it).second.size(); i++)
		                        {
		                        for (unsigned int kk=0; kk<nodesBlk1Ids.size(); kk++)
		                            comp_dom.integralCheck[nodesBlk1Ids[kk]][(*it).first] += 1;
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
				//std::cout<<"NonIntListPart3 Blocks: "<<*pos1<<" ";
				unsigned int block2Id = *pos1;
				//std::vector <cell_it> elemBlk2Ids = block2.GetBlockElementsList();
				for (unsigned int ii = 0; ii < nodesBlk1Ids.size(); ii++) //loop over each node of block1
					{
					Point<dim> &nodeBlk1 = support_points[nodesBlk1Ids.at(ii)];
					matrVectProdD(nodesBlk1Ids[ii]) += blockMultipoleExpansionsKer2[block2Id].Evaluate(nodeBlk1);
					matrVectProdN(nodesBlk1Ids[ii]) += blockMultipoleExpansionsKer1[block2Id].Evaluate(nodeBlk1);			        					
					
					/*//////this is for a check/////////////////////////
				        OctreeBlock<dim>* block2 = comp_dom.blocks[block2Id];
				        std::map <cell_it, std::vector<unsigned int> >
	                                blockQuadPointsList = block2->GetBlockQuadPointsList();
	                                typename std::map <cell_it, std::vector<unsigned int> >::iterator it;    
                                        for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
	                                    {
		                            for (unsigned int i=0; i<(*it).second.size(); i++)
		                                {
		                                comp_dom.integralCheck[nodesBlk1Ids[ii]][(*it).first] += 1;
				                }   
                                            }
                                        ////////////////////////////*/
					}												
				} // end loop over well separated blocks of smaller size (level)

				
			} // end loop over all sublevels in  nonIntlist
			
		} // end loop over blocks of each level
	
	} // end loop over levels
	
	
	// finally, when the loop over levels is done, we need to evaluate local expansions of all
	// childless blocks, at each block node(s)
	for (unsigned int kk = 0; kk < comp_dom.childlessList.size(); kk++)

		{
		unsigned int block1Id = comp_dom.childlessList[kk];	
		OctreeBlock<dim> *block1 = comp_dom.blocks[block1Id];
		std::vector <unsigned int> nodesBlk1Ids = block1->GetBlockNodeList();
		
		// loop over nodes of block
		for (unsigned int ii = 0; ii < nodesBlk1Ids.size(); ii++) //loop over each node of block1
			{
			Point<dim> &nodeBlk1 = support_points[nodesBlk1Ids.at(ii)];
			matrVectProdD(nodesBlk1Ids[ii]) += (blockLocalExpansionsKer2[block1Id]).Evaluate(nodeBlk1);
			matrVectProdN(nodesBlk1Ids[ii]) += (blockLocalExpansionsKer1[block1Id]).Evaluate(nodeBlk1);
			}	// end loop over nodes
			
		} // end loop over childless blocks
		
/*////////this is for a check//////////////////////		
for (unsigned int i = 0; i < comp_dom.dh.n_dofs(); i++)
    {
    for (cell_it cell = comp_dom.dh.begin_active(); cell != comp_dom.dh.end(); ++cell)
        {
	std::cout<<i<<" "<<cell<<" "<<comp_dom.integralCheck[i][cell]<<std::endl;
	comp_dom.integralCheck[i][cell] = 0;
	}
    }						
//////////////////////////////*/	
	



std::cout<<"...done computing multipole matrix-vector products"<<std::endl;

}

         // this method computes the preconditioner needed for the GMRES:
	 // to do it, it needs to receive the alpha vector from the bem_problem
	 // class, along with the constraint matrix of the bem problem

template <int dim>
TrilinosWrappers::PreconditionILU &BEMFMA<dim>::FMA_preconditioner(const TrilinosWrappers::MPI::Vector &alpha, ConstraintMatrix &c)
{

  // the final preconditioner (with constraints) has a slightly different sparsity pattern with respect
  // to the non constrained one. we must here initialize such sparsity pattern
  final_prec_sparsity_pattern.reinit(alpha.vector_partitioner(),125*comp_dom.fe.dofs_per_cell);
  //final_prec_sparsity_pattern.reinit(comp_dom.dh.n_dofs(),comp_dom.dh.n_dofs(),125*comp_dom.fe.dofs_per_cell);
  
  for (unsigned int i=0; i < comp_dom.dh.n_dofs(); i++)
      {
      if (c.is_constrained(i))
         {//cout<<i<<"  (c):"<<endl;
         // constrained nodes entries are taken from the bem problem constraint matrix
         final_prec_sparsity_pattern.add(i,i);
         const std::vector< std::pair < unsigned int, double > >
	         * entries = c.get_constraint_entries (i);
	 for (unsigned int j=0; j< entries->size(); ++j)
             final_prec_sparsity_pattern.add(i,(*entries)[j].first); 
         }
      else
         {//cout<<i<<"  (nc): ";
         // other nodes entries are taken from the unconstrained preconditioner matrix
         for (unsigned int j=0; j<comp_dom.dh.n_dofs(); ++j)
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
  
  final_prec_sparsity_pattern.compress();
  final_preconditioner.reinit(final_prec_sparsity_pattern);


  // now we assemble the final preconditioner matrix: the loop works
  // exactly like the previous one
  for (unsigned int i=0; i < comp_dom.dh.n_dofs(); i++)
      if (c.is_constrained(i))
         {
         final_preconditioner.set(i,i,1);
         // constrainednodes entries are taken from the bem problem constraint matrix
         const std::vector< std::pair < unsigned int, double > >
	         * entries = c.get_constraint_entries (i);
	 for (unsigned int j=0; j< entries->size(); ++j)
             {
             final_preconditioner.set(i,(*entries)[j].first,(*entries)[j].second);
             //cout<<i<<" "<<(*entries)[j].first<<" "<<(*entries)[j].second<<endl;
             } 
         }
      else
         {
         // other nodes entries are taken from the unconstrained preconditioner matrix
         for (unsigned int j=0; j<comp_dom.dh.n_dofs(); ++j)
             {
             if (init_prec_sparsity_pattern.exists(i,j))
                {
                final_preconditioner.set(i,j,init_preconditioner(i,j));
                //cout<<i<<" "<<j<<" "<<init_preconditioner(i,j)<<endl;
                }
             }
         }
  // finally, we have to add the alpha values on the diagonal, whenever dealing with a 
  // neumann (in such nodes the potential phi is an unknown) and non constrained node  

  for (unsigned int i=0; i < comp_dom.dh.n_dofs(); i++)
      if (comp_dom.surface_nodes(i) == 0 && !c.is_constrained(i))
         final_preconditioner.add(i,i,alpha(i));
 
     
  //preconditioner.print_formatted(std::cout,4,true,0," 0 ",1.);
  preconditioner.initialize(final_preconditioner);
  
  return preconditioner;  
}

template class BEMFMA<3>;
