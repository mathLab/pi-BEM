
#include "step_fma.h"

template <int dim> // REMOVE EVERY DOUBLE THING. BEM_FMA ALREADY KNOWS A LOT!
Operator::MinBEMOperator<dim>::MinBEMOperator(const BEMFMA<dim> &fma_in, const MPI_Comm comm_in, const IndexSet &cpu_set_in, const unsigned int mpi_process_in)
  :
  op_fma(fma_in),
  mpi_communicator(comm_in),
  this_mpi_process(mpi_process_in),
  this_cpu_set(cpu_set_in)
{}


template <int dim>
void Operator::MinBEMOperator<dim>::set_alpha()
{
  static TrilinosWrappers::MPI::Vector ones, zeros, dummy;
  ones.reinit(this_cpu_set,mpi_communicator);
  ones.add(-1.);
  zeros.reinit(this_cpu_set,mpi_communicator);
  dummy.reinit(this_cpu_set,mpi_communicator);
  alpha.reinit(this_cpu_set,mpi_communicator);

  op_fma.generate_multipole_expansions(ones,zeros);
  op_fma.multipole_matr_vect_products(ones,zeros,alpha,dummy);

  std::ofstream ofs;
  ofs.open ("fmm_alpha_"+Utilities::int_to_string(op_fma.trunc_order)+".txt", std::ofstream::out | std::ofstream::app);
  ofs << alpha.linfty_norm() <<" "<< alpha.l2_norm() << std::endl;;
  ofs.close();


}

template <int dim>
void Operator::MinBEMOperator<dim>::vmult(TrilinosWrappers::MPI::Vector &dst, const TrilinosWrappers::MPI::Vector &src) const
{
  TrilinosWrappers::MPI::Vector serv_phi(src);
  TrilinosWrappers::MPI::Vector serv_dphi_dn(src);



  TrilinosWrappers::MPI::Vector matrVectProdN;
  TrilinosWrappers::MPI::Vector matrVectProdD;

  //const unsigned int n_dofs =  op_fma.fma_dh->n_dofs();

  matrVectProdN.reinit(this_cpu_set,mpi_communicator);
  matrVectProdD.reinit(this_cpu_set,mpi_communicator);

  dst = 0;

  // TrilinosWrappers::MPI::Vector dirichlet_nodes(this_cpu_set,mpi_communicator);
  TrilinosWrappers::MPI::Vector other_nodes(this_cpu_set,mpi_communicator);

  // dirichlet_nodes = *op_fma.dirichlet_nodes;
  for ( unsigned int i=0; i<op_fma.dirichlet_nodes->size(); ++i)
    {
      other_nodes(i) = (double)(((int)(*op_fma.dirichlet_nodes)(i)+1)%2);
    }

  serv_phi.scale(other_nodes);
  serv_dphi_dn.scale(*op_fma.dirichlet_nodes);

  op_fma.generate_multipole_expansions(serv_phi,serv_dphi_dn);
  op_fma.multipole_matr_vect_products(serv_phi,serv_dphi_dn,matrVectProdN,matrVectProdD);
  serv_phi.scale(alpha);
  dst += matrVectProdD;
  dst *= -1;
  dst += matrVectProdN;
  dst += serv_phi;


}

template <int dim>
void Operator::MinBEMOperator<dim>::compute_rhs(TrilinosWrappers::MPI::Vector &dst, const TrilinosWrappers::MPI::Vector &src_dir, const TrilinosWrappers::MPI::Vector &src_neum)
{
  TrilinosWrappers::MPI::Vector serv_phi(src_dir);
  TrilinosWrappers::MPI::Vector serv_dphi_dn(src_neum);



  TrilinosWrappers::MPI::Vector matrVectProdN;
  TrilinosWrappers::MPI::Vector matrVectProdD;

  //const unsigned int n_dofs =  op_fma.fma_dh->n_dofs();

  matrVectProdN.reinit(this_cpu_set,mpi_communicator);
  matrVectProdD.reinit(this_cpu_set,mpi_communicator);

  dst = 0;

  // COSTRUIRE SURFACE E OTHER NODES!!!
  // TrilinosWrappers::MPI::Vector dirichlet_nodes(*op_fma.dirichlet_nodes);
  TrilinosWrappers::MPI::Vector other_nodes(this_cpu_set,mpi_communicator);

  // std::cout<<dirichlet_nodes.size()<< " " <<op_fma.dirichlet_nodes->size()<<std::endl;
  for ( unsigned int i=0; i<op_fma.dirichlet_nodes->size(); ++i)
    {
      other_nodes(i) = (double)(((int)(*op_fma.dirichlet_nodes)(i)+1)%2);
    }

  serv_phi.scale(other_nodes);
  serv_dphi_dn.scale(*op_fma.dirichlet_nodes);

  op_fma.generate_multipole_expansions(serv_phi,serv_dphi_dn);
  op_fma.multipole_matr_vect_products(serv_phi,serv_dphi_dn,matrVectProdN,matrVectProdD);
  serv_phi.scale(alpha);
  dst += matrVectProdD;
  dst *= -1;
  dst += matrVectProdN;
  dst += serv_phi;



}

template<int dim>
void Operator::MinBEMOperator<dim>::increase_fma_order()
{
  op_fma.trunc_order +=1;
}

template class Operator::MinBEMOperator<3>;
