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

				 // This class contains all the methods
				 // of the Fast Multipole Algorithm
				 // applied to the Boundary Element 
				 // Method




#ifndef bem_fma_h
#define bem_fma_h
				 
#include "../include/octree_block.h"
#include "../include/local_expansion.h"
#include "../include/multipole_expansion.h"
#include "../include/ass_leg_function.h"
#include "../include/computational_domain.h"
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>

#include <mpi.h>

using namespace dealii;

template <int dim>
class BEMFMA 
{
  public:
				 // Just renaming the cell iterator type
				 				   
    typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;

				 // Contructor: needs a reference to a
				 // the ComputationalDomain class,
				 // containnig geometry of the problem,
				 // quadrature methods and octree
				 // partitioning methods
    
    BEMFMA(ComputationalDomain<dim> &comp_dom);
    
				 // Parameters declaration    

    void declare_parameters(ParameterHandler &prm);
    
				 // Parametersparsing from input file  

    void parse_parameters(ParameterHandler &prm);
    
				 // Method computing the parts of the
				 // BEM system matrices in which the
				 // integrals have to be performed
				 // directly   

    void direct_integrals();

				 // Method computing the multipole
				 // expansion containing the integrals
				 // values for each bottom level block.
				 //  It is called once for each
				 // GMRES solved.

    void multipole_integrals();
    
                                 // Ascending phase of the FMA method.
				 // Multipole expansions are genarated
				 //  at the bottom level blocks, and then
				 // translated to their parent blocks up
				 // to the highest level. It is 
				 // called once per GMRES iteration.
    
    void generate_multipole_expansions(const TrilinosWrappers::MPI::Vector &phi_values, const TrilinosWrappers::MPI::Vector &dphi_dn_values) const;
    
                                 // Descending phase of the FMA method. Local
				 // Expansions are obtained for each block
				 // starting from the top level and down
				 // to the bottom ones, where they are used
				 // to approximete the values of the
				 // integrals, i.e. the BEM matrix-vector
				 // product values
    
    void multipole_matr_vect_products(const TrilinosWrappers::MPI::Vector &phi_values, const TrilinosWrappers::MPI::Vector &dphi_dn_values,
		                            TrilinosWrappers::MPI::Vector &matrVectProdN,    TrilinosWrappers::MPI::Vector &matrVectProdD) const;
    
                                 // Method for the assembling of the
				 // sparse preconitioning matrix for FMA
    
    TrilinosWrappers::PreconditionILU &FMA_preconditioner(const TrilinosWrappers::MPI::Vector &alpha, ConstraintMatrix &c);
    
    private:
    
                                 // Reference to the ComputationalDomain
				 // class
    
    ComputationalDomain<dim> &comp_dom;

                                 // Truncation order for the multipole
				 // and local expansion series: it is
				 // read from the parameters input file.

    unsigned int trunc_order;
    
                                 // Sparsity pattern for the
				 // initial preconditioning matrix
                                 // (no constraints applied)
    
    TrilinosWrappers::SparsityPattern init_prec_sparsity_pattern;
    
                                 // Sparsity pattern for the
				 // final preconditioning matrix
                                 // (constraints applied)
    
    TrilinosWrappers::SparsityPattern final_prec_sparsity_pattern;

                                 // Sparse Neumann matrix containing
				 // only direct integrals contributions
    TrilinosWrappers::SparseMatrix prec_neumann_matrix;    
    //SparseMatrix<double> prec_neumann_matrix;
    
                                 // Sparse Dirichlet matrix containing
				 // only direct integrals contributions

    TrilinosWrappers::SparseMatrix prec_dirichlet_matrix;    
    //SparseMatrix<double> prec_dirichlet_matrix;
    
                                 // Initial sparse preconditioning matrix (without constraints)
    
    TrilinosWrappers::SparseMatrix init_preconditioner;

                                 // Final sparse preconditioning matrix (with constraints)
    
    TrilinosWrappers::SparseMatrix final_preconditioner;
    
                                 // Structures where the Dirichlet
				 // matrix multipole
				 // integrals are stored: for each cell
				 // there are as many multipole
				 // expansions as the lowest level
				 // blocks in which element's quad
				 // points lie.  
        
    mutable std::map <unsigned int, std::map <cell_it, std::vector <MultipoleExpansion > > > elemMultipoleExpansionsKer1; 

                                 // Structures where the Neumann
				 // matrix multipole
				 // integrals are stored: for each cell
				 // there are as many multipole
				 // expansions as the lowest level
				 // blocks in which element's quad
				 // points lie.  
 

    mutable std::map <unsigned int, std::map <cell_it, std::vector <MultipoleExpansion > > > elemMultipoleExpansionsKer2;

                                 // Vector storing the Dirichlet
				 // integrals multipole expansions 
				 // for each block

    mutable std::vector <MultipoleExpansion > blockMultipoleExpansionsKer1;

                                 // Vector storing the Neumann
				 // integrals multipole expansions 
				 // for each block

    mutable std::vector <MultipoleExpansion > blockMultipoleExpansionsKer2;

                                 // Vector storing the Dirichlet
				 // integrals local expansions 
				 // for each block
    
    mutable std::vector <LocalExpansion > blockLocalExpansionsKer1;
    
                                 // Vector storing the Neumann
				 // integrals local expansions 
				 // for each block
    
    mutable std::vector <LocalExpansion > blockLocalExpansionsKer2;   
 
                                 // Associated Legendre functions class
    
    AssLegFunction assLegFunction;
    
                                 // the preconditioner to be passed to bem_problem

    TrilinosWrappers::PreconditionILU preconditioner;

                                 // mpi related variables

    MPI_Comm mpi_communicator;

    unsigned int n_mpi_processes;

    unsigned int this_mpi_process;


    
};

#endif
