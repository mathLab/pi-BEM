#include "bem_problem.h"

using namespace dealii;
using namespace deal2lkit;

template <int dim>
class BEMProblemAccess
{

   BEMProblemAccess(){};

   void initialize(const BEMProblem<dim> *in){
     bem_problem = in;
   };
   const DoFHandler<dim-1, dim> * get_bem_dh(){
     return bem_problem->dh;
   };

   const FiniteElement<dim-1, dim> * get_bem_fe(){
     return bem_problem->fe;
   }
  //  const Triangulation<dim-1, dim> & get_tria();
   const Mapping<dim-1,dim> * get_bem_mapping(){
     return bem_problem->mapping;
   };
   const Quadrature<dim-1> * get_bem_quadrature(){
     return bem_problem->quadrature;
   };
   unsigned int get_bem_singular_quadrature_order(){
     return bem_problem->singular_quadrature_order;
   };
   ConstraintMatrix & get_bem_constraint_matrix(){
     return bem_problem->constraints;
   }
   std::vector <std::set<types::global_dof_index> > & get_bem_double_nodes_set()
    {
      return bem_problem->double_nodes_set;
    }
   TrilinosWrappers::MPI::Vector & get_bem_dirichlet_nodes(){
     return bem_problem->dirichlet_nodes;
   }

  //  get_pinco();
  //  get_palla();

private:
   BEMProblem<dim> * bem_problem;

};
