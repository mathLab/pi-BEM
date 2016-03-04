#include <bem_problem.h>

using namespace dealii;
// using namespace deal2lkit;

// template <int dim> class BEMProblem;

template <int dim>
class BEMProblemAccess
{
public:
   BEMProblemAccess(){};

   void initialize(BEMProblem<dim> *in){
     bem_problem = in;
   };
  DoFHandler<dim-1, dim> * get_bem_dh() const {
     return bem_problem->dh;
   };

   FiniteElement<dim-1, dim> * get_bem_fe() const{
     return bem_problem->fe;
   }
  //  const Triangulation<dim-1, dim> & get_tria();
   Mapping<dim-1,dim> * get_bem_mapping() const{
     return bem_problem->mapping;
   };
   Quadrature<dim-1> * get_bem_quadrature() const{
     return bem_problem->quadrature;
   };
   unsigned int get_bem_singular_quadrature_order(){
     return bem_problem->singular_quadrature_order;
   };
   ConstraintMatrix & get_bem_constraint_matrix() const{
     return bem_problem->constraints;
   }
   std::vector <std::set<types::global_dof_index> > & get_bem_double_nodes_set()
    {
      return bem_problem->double_nodes_set;
    }
   TrilinosWrappers::MPI::Vector & get_bem_dirichlet_nodes() const{
     return bem_problem->dirichlet_nodes;
   }

  //  get_pinco();
  //  get_palla();

private:
   BEMProblem<dim> * bem_problem;

};
