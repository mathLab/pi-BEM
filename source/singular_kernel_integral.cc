#include "../include/singular_kernel_integral.h"

template<>
SingularKernelIntegral<3>::
SingularKernelIntegral(DoFHandler<2, 3>::active_cell_iterator in_cell,
                       FiniteElement<2, 3> &in_fe,
                       Mapping<2, 3> &in_mapping,
                       Point<2> &in_eta)
 : cell(in_cell)
 , fe(in_fe)
 , mapping(in_mapping)
 , eta(in_eta)
{
};


template<>
SingularKernelIntegral<2>::
SingularKernelIntegral(DoFHandler<1, 2>::active_cell_iterator in_cell,
                       FiniteElement<1, 2> &in_fe,
                       Mapping<1, 2> &in_mapping,
                       Point<1> &in_eta)
 : cell(in_cell)
 , fe(in_fe)
 , mapping(in_mapping)
 , eta(in_eta)
{
};


template<>
Tensor<1,2> SingularKernelIntegral<2>::evaluate_Vk_integrals()
{
ExcNotImplemented();
Tensor<1,2> dummy;
return dummy;
};

template<>
Tensor<1,3> SingularKernelIntegral<3>::evaluate_Vk_integrals()
{


  Tensor<1,3> I_0;
  Tensor<1,3> I_1;
  Tensor<1,3> I_2;

  // we now need to compute all the mapping first and second order derivatives
  // in correspondence with the singularity
  // we have to pass through a FEValues object, but to build one we will 
  // need std::vectors with location of quadrature points and quadrature nodes

  // we have a single quadrature point (eta) located at P
  std::vector<Point<2> > eta_q_points(1);
  eta_q_points[0] = eta;
  // and a single quadrature weight set to one
  std::vector<double> eta_q_weights(1);
  eta_q_weights[0] = 1;

  // here's the quadrature rule obtained with the point and weight generated
  Quadrature<2> eta_quadrature(eta_q_points, eta_q_weights);
  // and here's the FEValues class resulting by it
  FEValues<2,3> eta_fe_values(mapping,
                              fe,
                              eta_quadrature,
                              update_values | update_gradients | update_quadrature_points |
                              update_JxW_values | update_normal_vectors | update_jacobians |
                              update_jacobian_grads);
  // we initialize FEValues object on the current cell
  eta_fe_values.reinit(cell);
  // the single quadrature point is the point in the three dimensional domain corresponding to eta/P
  const std::vector<Point<3>> &eta_q_points_spacedim    = eta_fe_values.get_quadrature_points();
  // we also get the acobian and jacobian gradient at such a location
  auto eta_jacobian = eta_fe_values.jacobian(0);
  auto eta_jacobian_grad = eta_fe_values.jacobian_grad(0);
  // we now have all the first and second order mapping derivaives
  // the product of the jacobian by the normal vector Jxn is the surface normal vector
  // let us obtain it in terms of the mapping derivatives
  Tensor<1,3> eta_Jxn;
  eta_Jxn[0] = eta_jacobian[1][0]*eta_jacobian[2][1]-eta_jacobian[1][1]*eta_jacobian[2][0];
  eta_Jxn[1] = eta_jacobian[0][1]*eta_jacobian[2][0]-eta_jacobian[0][0]*eta_jacobian[2][1];
  eta_Jxn[2] = eta_jacobian[0][0]*eta_jacobian[1][1]-eta_jacobian[0][1]*eta_jacobian[1][0];
  // we also want the derivative of the previous vector with respect to the first parametric plane coordinate u
  Tensor<1,3> d_eta_Jxn_du;
  d_eta_Jxn_du[0] = eta_jacobian_grad[1][0][0]*eta_jacobian[2][1]+eta_jacobian[1][0]*eta_jacobian_grad[2][1][0]-
                    eta_jacobian_grad[1][1][0]*eta_jacobian[2][0]-eta_jacobian[1][1]*eta_jacobian_grad[2][0][0];
  d_eta_Jxn_du[1] = eta_jacobian_grad[0][1][0]*eta_jacobian[2][0]+eta_jacobian[0][1]*eta_jacobian_grad[2][0][0]-
                    eta_jacobian_grad[0][0][0]*eta_jacobian[2][1]-eta_jacobian[0][0]*eta_jacobian_grad[2][1][0];
  d_eta_Jxn_du[2] = eta_jacobian_grad[0][0][0]*eta_jacobian[1][1]+eta_jacobian[0][0]*eta_jacobian_grad[1][1][0]-
                    eta_jacobian_grad[0][1][0]*eta_jacobian[1][0]-eta_jacobian[0][1]*eta_jacobian_grad[1][0][0];
  // we also want the derivative of the previous vector with respect to the second parametric plane coordinate v
  Tensor<1,3> d_eta_Jxn_dv;
  d_eta_Jxn_dv[0] = eta_jacobian_grad[1][0][1]*eta_jacobian[2][1]+eta_jacobian[1][0]*eta_jacobian_grad[2][1][1]-
                    eta_jacobian_grad[1][1][1]*eta_jacobian[2][0]-eta_jacobian[1][1]*eta_jacobian_grad[2][0][1];
  d_eta_Jxn_dv[1] = eta_jacobian_grad[0][1][1]*eta_jacobian[2][0]+eta_jacobian[0][1]*eta_jacobian_grad[2][0][1]-
                    eta_jacobian_grad[0][0][1]*eta_jacobian[2][1]-eta_jacobian[0][0]*eta_jacobian_grad[2][1][1];
  d_eta_Jxn_dv[2] = eta_jacobian_grad[0][0][1]*eta_jacobian[1][1]+eta_jacobian[0][0]*eta_jacobian_grad[1][1][1]-
                    eta_jacobian_grad[0][1][1]*eta_jacobian[1][0]-eta_jacobian[0][1]*eta_jacobian_grad[1][0][1];
  //std::cout<<"dim singularity location: "<<eta<<std::endl;
  //std::cout<<"Spacedim singularity location: "<<eta_q_points_spacedim[0]<<std::endl;
  //std::cout<<"d_eta_Jxn_du "<<d_eta_Jxn_du<<std::endl;
  //std::cout<<"d_eta_Jxn_dv "<<d_eta_Jxn_dv<<std::endl;
//      std::cout<<"H1r: "<<eta_jacobian_grad[0][0][0]<<" "<<eta_jacobian_grad[0][0][1]<<std::endl;
//      std::cout<<"     "<<eta_jacobian_grad[0][1][0]<<" "<<eta_jacobian_grad[0][1][1]<<std::endl;
//      std::cout<<"H2r: "<<eta_jacobian_grad[1][0][0]<<" "<<eta_jacobian_grad[1][0][1]<<std::endl;
//      std::cout<<"     "<<eta_jacobian_grad[1][1][0]<<" "<<eta_jacobian_grad[1][1][1]<<std::endl;
//      std::cout<<"H3r: "<<eta_jacobian_grad[2][0][0]<<" "<<eta_jacobian_grad[2][0][1]<<std::endl;
//      std::cout<<"     "<<eta_jacobian_grad[2][1][0]<<" "<<eta_jacobian_grad[2][1][1]<<std::endl;

  
  // we now loop over the cell faces to carry out the integrals
  for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
      {
      // (to the best of my knowledge)
      // deal.ii does not offer an easy way to get the values of the
      // face vertices in the cell parametric plane
      // thus, we generate them here with a switch, given
      // the deal.ii numbering of faces and vertices
      // in generating them, we already center them in the
      // singularity point, subtracting P
      Point<2> csi_0;
      Point<2> csi_1;
      double coeff;
      switch (f)
          {
          case 0:
               csi_0(0) =  0-eta(0); csi_1(0) =  0-eta(0);
               csi_0(1) =  1-eta(1); csi_1(1) =  0-eta(1);
               coeff = csi_0(0);
          break;
          case 1:
               csi_0(0) =  1-eta(0); csi_1(0) =  1-eta(0);
               csi_0(1) =  0-eta(1); csi_1(1) =  1-eta(1);
               coeff = csi_0(0);
          break;
          case 2:
               csi_0(0) =  0-eta(0); csi_1(0) =  1-eta(0);
               csi_0(1) =  0-eta(1); csi_1(1) =  0-eta(1);
               coeff = csi_1(1);
          break;
          case 3:
               csi_0(0) =  1-eta(0); csi_1(0) =  0-eta(0);
               csi_0(1) =  1-eta(1); csi_1(1) =  1-eta(1);
               coeff = csi_1(1);
          break;
          default:
               csi_0(0) =  0-eta(0); csi_1(0) =  0-eta(0);
               csi_0(1) =  0-eta(1); csi_1(1) =  0-eta(1);
          break;
          }
          
      // if the singularity is located on the face, do not consider this face in
      // the integration
      if (!(((csi_0(0) == 0) && (csi_0(1) == 0)) || ((csi_1(0) == 0) && (csi_1(1) == 0))))
          { 
          
          // once the coordinates in the parametric plane of the
          // face vertices are known, we can compute the polar coordinate
          // theta angle, taking care of having always growing angles
          // on each face f
          double theta_0 = atan2((csi_0)[1],(csi_0)[0]);
          if (theta_0 < 0)
             theta_0 += 2*dealii::numbers::PI;
          double theta_1 = atan2((csi_1)[1],(csi_1)[0]);
          if (theta_1 < 0)
             theta_1 += 2*dealii::numbers::PI;
          //std::cout<<"PRE:   "<<f<<" ->    theta_0: "<<theta_0*180/dealii::numbers::PI<<"   theta_1: "<<theta_1*180/dealii::numbers::PI<<std::endl;   
          if (theta_1<theta_0)
             theta_0 -= 2*dealii::numbers::PI;
          //std::cout<<f<<" ->    theta_0: "<<theta_0*180/dealii::numbers::PI<<"   theta_1: "<<theta_1*180/dealii::numbers::PI<<std::endl;
          

          // we now create a 1D triangulation for the integration
          // in the theta direction of the polar coordinates
          Triangulation<1> theta_triangulation;

          std::vector<Point<1> > theta_vertices;
          std::vector<CellData<1> > theta_cells;
          SubCellData theta_subcelldata;
          
          // we only have one cell that goes from the minimum to maximum
          // theta on the face f 
          theta_vertices.resize(2);
          theta_vertices[0](0)=theta_0;
          theta_vertices[1](0)=theta_1;

          theta_cells.resize(1);

          theta_cells[0].vertices[0]=0;
          theta_cells[0].vertices[1]=1;
          
          GridTools::delete_unused_vertices (theta_vertices, theta_cells, theta_subcelldata);
          GridTools::consistently_order_cells(theta_cells);
          // here the triangulation is initialized and ready
          theta_triangulation.create_triangulation(theta_vertices, theta_cells, theta_subcelldata );

          // we of course also create a dof handler finite element and mapping
          DoFHandler<1> theta_dof_handler(theta_triangulation);
          const FE_Q<1> theta_finite_element(fe.degree);
          theta_dof_handler.distribute_dofs(theta_finite_element);
          MappingQ<1> theta_mapping(theta_finite_element.degree);

          // we create a 1D quadature with the user selected number of quadrature nodes
          QGauss<1> theta_quadrature(theta_quadrature_order);

          // we are now ready to create a FEValues objects that will allow us
          // to obtain all we need for the integration over theta on the face f
          FEValues<1> theta_fe_v(theta_mapping,
                                theta_finite_element,
                                theta_quadrature,
                                update_values |
                                update_quadrature_points | update_JxW_values);

          // we initialize theta_fe_v with the only face available
          auto theta_cell = theta_dof_handler.begin_active(); //this is the only cell
          theta_fe_v.reinit(theta_cell);


          // because we also will need to carry out the surface integral in polar coordinates
          // for I_0, we will also need a 1D triangulation for rho
          Triangulation<1> rho_triangulation;

          std::vector<Point<1> >  rho_vertices;
          std::vector<CellData<1> >  rho_cells;
          SubCellData  rho_subcelldata;
          
          // in this case we will just consider a cell that goes from 0 to 1 and
          // take care of the transformation and jacobian by ourself
          rho_vertices.resize(2);
          rho_vertices[0](0)=0;
          rho_vertices[1](0)=1;

          rho_cells.resize(1);

          rho_cells[0].vertices[0]=0;
          rho_cells[0].vertices[1]=1;
          
          GridTools::delete_unused_vertices (rho_vertices, rho_cells, rho_subcelldata);
          GridTools::consistently_order_cells(rho_cells);
          // here the triangulation is initialized and ready
          rho_triangulation.create_triangulation(rho_vertices, rho_cells, rho_subcelldata );
          // we of course also create a dof handler finite element and mapping
          DoFHandler<1> rho_dof_handler(rho_triangulation);
          const FE_Q<1> rho_finite_element(fe.degree);
          rho_dof_handler.distribute_dofs(rho_finite_element);
          MappingQ<1> rho_mapping(rho_finite_element.degree);
          // we create a 1D quadature with the user selected number of quadrature nodes
          QGauss<1> rho_quadrature(rho_quadrature_order);
          // we are now ready to create a FEValues object that will allow us
          // to obtain all we need for the integration over rho on the face f
          FEValues<1> rho_fe_v(rho_mapping,
                               rho_finite_element,
                               rho_quadrature,
                               update_values |
                               update_quadrature_points | update_JxW_values);
          // we initialize theta_fe_v with the only face available
          auto rho_cell = rho_dof_handler.begin_active(); //this is the only cell
          rho_fe_v.reinit(rho_cell);
          // we get the quadrature points for rho (remember that they go from 0 to 1)
          const std::vector<Point<1>> &rho_q_points    = rho_fe_v.get_quadrature_points();

          // we get the quadrature points for theta (remember that they go from
          // theta_min to theta_max in the cell)
          const std::vector<Point<1>> &theta_q_points    = theta_fe_v.get_quadrature_points();
          // we want to create a vector with the \hat(rho) values corresponding to the face boundary
          // at each theta value
          std::vector<double> q_rho(theta_q_points.size());
          // using deal.ii numbering, these lines do the trick
          for (unsigned int q=0; q<theta_q_points.size(); ++q)
              {
              double a = (f<2)?0.0:1.0;
              double b = (f<2)?1.0:0.0;
              q_rho[q] = coeff/(a*sin(theta_q_points[q](0))+b*(cos(theta_q_points[q](0))));
              }
          // we now create a set of quadrature points and weights to be passed to a
          // FEValues<2,3> to take care of the mapping to the 3D space
          const std::vector<double> &uv_q_weights = theta_fe_v.get_JxW_values();
          std::vector<Point<2> > uv_q_points(theta_q_points.size());
          for (unsigned int q=0; q<theta_q_points.size(); ++q)
              {
              // the quadrature points are the points of the face f in the parametric space
              uv_q_points[q] = Point<2>(q_rho[q]*cos(theta_q_points[q](0))+eta(0),
                                        q_rho[q]*sin(theta_q_points[q](0))+eta(1));
              }
          Quadrature<2> uv_quadrature(uv_q_points, uv_q_weights);
          //here we create the FEValues and initialize it with the <2,3> dof handler cell
          FEValues<2,3> face_fe_values(fe,
                                       uv_quadrature,
                                       update_values | update_gradients | update_quadrature_points |
                                       update_JxW_values | update_normal_vectors);
          face_fe_values.reinit(cell);

          // we now create a set of vectors that will contain useful theta dependent coefficients
          // for each value of theta on the face f
          std::vector<Tensor<1,3> > A(theta_q_points.size());
          std::vector<Tensor<1,3> > B(theta_q_points.size());
          // values of A.norm()
          std::vector<double> q_A(theta_q_points.size());
          // values of A*B
          std::vector<double> q_C(theta_q_points.size());
          // values of jacobian 0-th order Taylor series coefficient
          std::vector<Tensor<1,3> > Jk0(theta_q_points.size());
          // values of jacobian 1-st order Taylor series coefficient
          std::vector<Tensor<1,3> > Jk1(theta_q_points.size());
          for (unsigned int q=0; q<theta_q_points.size(); ++q)
              {
              for (unsigned int i=0; i<3; ++i)
                  {
                  A[q][i] = eta_jacobian[i][0]*cos(theta_q_points[q](0)) +
                         eta_jacobian[i][1]*sin(theta_q_points[q](0));
                  B[q][i] = eta_jacobian_grad[i][0][0]*pow(cos(theta_q_points[q](0)),2)/2 +
                            eta_jacobian_grad[i][0][1]*cos(theta_q_points[q](0))*sin(theta_q_points[q](0)) +
                            eta_jacobian_grad[i][1][1]*pow(sin(theta_q_points[q](0)),2)/2;
                  }
              q_A[q] = A[q].norm();
              q_C[q] = A[q]*B[q];
              Jk0[q] = eta_Jxn;
              Jk1[q] = d_eta_Jxn_du*cos(theta_q_points[q](0))+d_eta_Jxn_dv*sin(theta_q_points[q](0));
              }

          // so let's start computing the integrals on the face f
          // for I_1 and I_2
          for (unsigned int q=0; q<theta_q_points.size(); ++q)
              {
              // the taylor expansion coefficients F_1 and F_2 are finally computed
              double beta = 1/q_A[q];
              double gamma = -q_C[q]/pow(q_A[q],4);
              Tensor<1,3> F_2 = Jk0[q]/pow(q_A[q],3);
              Tensor<1,3> F_1 = -3*q_C[q]*Jk0[q]/pow(q_A[q],5)
                                -3*A[q]/pow(q_A[q],5)*(Jk0[q]*B[q]+Jk1[q]*A[q])
                                +Jk1[q]/pow(q_A[q],3);
              //std::cout<<"f: "<<f<<"  q: "<<q<<" theta: "<<theta_q_points[q](0)*180/dealii::numbers::PI<<"  F_1: "<<F_1<<"  F_2: "<<F_2<<"  rho: "<<q_rho[q]<<"  beta: "<<beta<<"  gamma: "<<gamma<<"  w: "<<uv_q_weights[q]<<"  Jk0: "<<Jk0[q]<<"  Jk1: "<<Jk1[q]<<"  A: "<<A[q]<<"  B: "<<B[q]<<std::endl;
              // and the integral contribution of the present quadrature node is here added
              I_2 += -F_2*(gamma/pow(beta,2)+1/q_rho[q])*uv_q_weights[q];
              I_1 += F_1*log(fabs(q_rho[q]/beta))*uv_q_weights[q];

              //I_1 += F_1*log(fabs(q_rho[q]))*uv_q_weights[q];
              
              
              // we still miss the surface integral I_0
              // for this, we will need to also consider the
              // quadrature nodes of the rho_fe_values
              // remembering they go from 0 to 1, we must rescale them from
              // 0 to q_rho[q]. Here is the obvious scaling factor
              double rho_jac_fact = q_rho[q];

              // we must again create a set of quadrature points and weights to be passed to a
              // FEValues<2,3> to take care of the mapping to the 3D space
              
              // we only want the mapping jacobian from such a FEValues, so we set all the
              // weights to 1
              std::vector<double> inner_uv_q_weights(rho_q_points.size(),1.0);
              
              std::vector<Point<2> > inner_uv_q_points(rho_q_points.size());
              for (unsigned int p=0; p<rho_q_points.size(); ++p)
                  {
                  // based on rho and theta, the quadrature points are generated in polar coordinates
                  // centered at the singularity eta_q_points[0]
                  double rho = rho_jac_fact*rho_q_points[p](0);
                  inner_uv_q_points[p] = Point<2>(rho*cos(theta_q_points[q](0))+eta_q_points[0](0),
                                                  rho*sin(theta_q_points[q](0))+eta_q_points[0](1));
                  }
              
              // we can then create the quadrature and the FEValues
              Quadrature<2> inner_uv_quadrature(inner_uv_q_points, inner_uv_q_weights);
              FEValues<2,3> inner_face_fe_values(fe,
                                                 inner_uv_quadrature,
                                                 update_values | update_gradients | update_quadrature_points |
                                                 update_JxW_values | update_normal_vectors);
              inner_face_fe_values.reinit(cell);
              // let's obtain the quadrature points in the 3D space from the inner_face_fe_values
              const std::vector<Point<3>> &inner_uv_q_points_spacedim    = inner_face_fe_values.get_quadrature_points();
              const std::vector<Tensor<1, 3 >> &inner_uv_q_normals = inner_face_fe_values.get_normal_vectors();
              for (unsigned int p=0; p<rho_q_points.size(); ++p)
                  {
                  // this is the rho obtained scaling the quadrature point in the [0,1] interval
                  double rho = rho_jac_fact*rho_q_points[p](0);
                  // this is the distance in the 3D space between quadrature point and singularity
                  Tensor<1,3> R = inner_uv_q_points_spacedim[p]-eta_q_points_spacedim[0];
                  double r = R.norm();
                  // finally, this is the integral argument
                  // the first term also needs the jacobian of the 3D mapping, which was instead
                  // already included in the F_2 and F_1 factors of the latter terms, and is therefore
                  // not needed there
                  Tensor<1,3> fun = -1/pow(r,3)*(3*(R/r)*(R*inner_uv_q_normals[p]/r)-inner_uv_q_normals[p]);
                  Tensor<1,3> arg = fun*rho*inner_face_fe_values.JxW(p)-F_2/pow(rho,2)-F_1/rho;

    //                  std::cout<<"rho: "<<rho<<"  r: "<<r<<"  arg1: "<<1/pow(r,3)*rho*inner_face_fe_values.JxW(p)<<
    //                  "   arg2: "<<-F_2/pow(rho,2)-F_1/rho<<"  diff: "<<1/pow(r,3)*rho*inner_face_fe_values.JxW(p)-F_2/pow(rho,2)-F_1/rho<<std::endl;

                  //std::cout<<rho<<"  "<<r<<"  "<<1/pow(r,3)*rho*inner_face_fe_values.JxW(p)<<
                  //"  "<<-F_2/pow(rho,2)-F_1/rho<<"  "<<1/pow(r,3)*rho*inner_face_fe_values.JxW(p)-F_2/pow(rho,2)-F_1/rho<<std::endl;

    //                  std::cout<<rho<<"  "<<r<<"  "<<Jk0[q][2]+rho*Jk1[q][2]<<"  "<<inner_uv_q_normals[p][2]*inner_face_fe_values.JxW(p)<<"   "<<inner_uv_q_normals[p][2]*inner_face_fe_values.JxW(p)-Jk0[q][2]-rho*Jk1[q][2]<<std::endl;

    //                  std::cout<<rho<<"  "<<R<<"  "<<rho*A[q]+rho*rho*B[q]<<"  "<<(R-rho*A[q]-rho*rho*B[q]).norm()<<std::endl;
                  
                  // we add the conribution of the present node to the surface integral
                  // the inner weights are given by rho_fe_v.JxW(p) but the jacobian is given by the scaling
                  // factor rho_jac_fact
                  I_0 += arg*rho_jac_fact*rho_fe_v.JxW(p)*theta_fe_v.JxW(q);
                  }
              }
          /*
          std::cout<<"temp I_0 = "<<I_0<<std::endl;
          std::cout<<"temp I_1 = "<<I_1<<std::endl;
          std::cout<<"temp I_2 = "<<I_2<<std::endl;
          */
          }
      }



 return (I_0 + I_1 + I_2)/numbers::PI/4.0;
};

template<>
std::vector<Tensor<1,2> > SingularKernelIntegral<2>::evaluate_VkNj_integrals()
{
ExcNotImplemented();
std::vector<Tensor<1,2> > dummy;
return dummy;
}

template<>
std::vector<Tensor<1,3> > SingularKernelIntegral<3>::evaluate_VkNj_integrals()
{
  Tensor<1,3> zero_tensor;
  std::vector<Tensor<1,3> > I_0(fe.dofs_per_cell,zero_tensor);
  std::vector<Tensor<1,3> > I_1(fe.dofs_per_cell,zero_tensor);
  std::vector<Tensor<1,3> > I_2(fe.dofs_per_cell,zero_tensor);

  // We must obtain the shape function derivatives at the singularity point eta in
  // the parametric plane.
  // In the deal.II framework the (to my knowledge) best way is to
  // create a reference Triangulation<dim,spacedim> and FE<dim,spacedim> and use it
  // to get the shape function values and derivatives
  // building reference cell
  Triangulation<2,3> ref_triangulation;

  std::vector<Point<3> > ref_vertices;
  std::vector<CellData<2> > ref_cells;
  SubCellData ref_subcelldata;

  ref_vertices.resize(4);
  ref_vertices[0](0)= 0.0;
  ref_vertices[0](1)= 0.0;
  ref_vertices[0](2)=0.0;
  ref_vertices[1](0)= 1.0;
  ref_vertices[1](1)= 0.0;
  ref_vertices[1](2)=0.0;
  ref_vertices[2](0)= 0.0;
  ref_vertices[2](1)= 1.0;
  ref_vertices[2](2)=0.0;
  ref_vertices[3](0)= 1.0;
  ref_vertices[3](1)= 1.0;
  ref_vertices[3](2)=0.0;

  ref_cells.resize(1);

  ref_cells[0].vertices[0]=0;
  ref_cells[0].vertices[1]=1;
  ref_cells[0].vertices[2]=2;
  ref_cells[0].vertices[3]=3;
  ref_cells[0].material_id = 1;

  GridTools::delete_unused_vertices (ref_vertices, ref_cells, ref_subcelldata);
  GridTools::consistently_order_cells(ref_cells);

  // here's the triangulation set up
  ref_triangulation.create_triangulation(ref_vertices, ref_cells, ref_subcelldata );

  // we will need a codimension one finite element and
  // the corresponding dof handler
  DoFHandler<2,3> ref_dof_handler(ref_triangulation);
  ref_dof_handler.distribute_dofs(fe);


  // we now need to compute all the mapping first and second order derivatives
  // in correspondence with the singularity
  // we have to pass through a FEValues object, but to build one we will 
  // need std::vectors with location of quadrature points and quadrature nodes

  // we have a single quadrature point (eta) located at P
  std::vector<Point<2> > eta_q_points(1);
  eta_q_points[0] = eta;
  // and a single quadrature weight set to one
  std::vector<double> eta_q_weights(1);
  eta_q_weights[0] = 1;

  // here's the quadrature rule obtained with the point and weight generated
  Quadrature<2> eta_quadrature(eta_q_points, eta_q_weights);
  // and here are the reference and spacedim FEValues class resulting by it
  // (both FEValues are initialized with the correct cell)
  
  FEValues<2,3> ref_eta_fe_values(fe,
                                  eta_quadrature,
                                  update_values | update_gradients | update_quadrature_points |
                                  update_JxW_values | update_normal_vectors | update_jacobians |
                                  update_jacobian_grads);
  // we initialize rference FEValues object on the reference cell
  ref_eta_fe_values.reinit(ref_dof_handler.begin_active());
  
  FEValues<2,3> eta_fe_values(mapping,
                              fe,
                              eta_quadrature,
                              update_values | update_gradients | update_quadrature_points |
                              update_JxW_values | update_normal_vectors | update_jacobians |
                              update_jacobian_grads);
  // we initialize FEValues object on the current cell
  eta_fe_values.reinit(cell);
  
  // from the reference FEValues we get the n_dofs shape function values and gradients at the
  // singularity point
  std::vector<double> eta_shape_values(fe.dofs_per_cell);
  std::vector<Tensor<1,3> > eta_shape_grads(fe.dofs_per_cell);
  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      {
      eta_shape_values[i] = ref_eta_fe_values.shape_value(i,0);
      eta_shape_grads[i] = ref_eta_fe_values.shape_grad(i,0);
      }
  
  // the single quadrature point is the point in the three dimensional domain corresponding to eta/P
  const std::vector<Point<3>> &eta_q_points_spacedim    = eta_fe_values.get_quadrature_points();
  // we also get the acobian and jacobian gradient at such a location
  auto eta_jacobian = eta_fe_values.jacobian(0);
  auto eta_jacobian_grad = eta_fe_values.jacobian_grad(0);
  // we now have all the first and second order mapping derivaives
  // the product of the jacobian by the normal vector Jxn is the surface normal vector
  // let us obtain it in terms of the mapping derivatives
  Tensor<1,3> eta_Jxn;
  eta_Jxn[0] = eta_jacobian[1][0]*eta_jacobian[2][1]-eta_jacobian[1][1]*eta_jacobian[2][0];
  eta_Jxn[1] = eta_jacobian[0][1]*eta_jacobian[2][0]-eta_jacobian[0][0]*eta_jacobian[2][1];
  eta_Jxn[2] = eta_jacobian[0][0]*eta_jacobian[1][1]-eta_jacobian[0][1]*eta_jacobian[1][0];
  // we also want the derivative of the previous vector with respect to the first parametric plane coordinate u
  Tensor<1,3> d_eta_Jxn_du;
  d_eta_Jxn_du[0] = eta_jacobian_grad[1][0][0]*eta_jacobian[2][1]+eta_jacobian[1][0]*eta_jacobian_grad[2][1][0]-
                    eta_jacobian_grad[1][1][0]*eta_jacobian[2][0]-eta_jacobian[1][1]*eta_jacobian_grad[2][0][0];
  d_eta_Jxn_du[1] = eta_jacobian_grad[0][1][0]*eta_jacobian[2][0]+eta_jacobian[0][1]*eta_jacobian_grad[2][0][0]-
                    eta_jacobian_grad[0][0][0]*eta_jacobian[2][1]-eta_jacobian[0][0]*eta_jacobian_grad[2][1][0];
  d_eta_Jxn_du[2] = eta_jacobian_grad[0][0][0]*eta_jacobian[1][1]+eta_jacobian[0][0]*eta_jacobian_grad[1][1][0]-
                    eta_jacobian_grad[0][1][0]*eta_jacobian[1][0]-eta_jacobian[0][1]*eta_jacobian_grad[1][0][0];
  // we also want the derivative of the previous vector with respect to the second parametric plane coordinate v
  Tensor<1,3> d_eta_Jxn_dv;
  d_eta_Jxn_dv[0] = eta_jacobian_grad[1][0][1]*eta_jacobian[2][1]+eta_jacobian[1][0]*eta_jacobian_grad[2][1][1]-
                    eta_jacobian_grad[1][1][1]*eta_jacobian[2][0]-eta_jacobian[1][1]*eta_jacobian_grad[2][0][1];
  d_eta_Jxn_dv[1] = eta_jacobian_grad[0][1][1]*eta_jacobian[2][0]+eta_jacobian[0][1]*eta_jacobian_grad[2][0][1]-
                    eta_jacobian_grad[0][0][1]*eta_jacobian[2][1]-eta_jacobian[0][0]*eta_jacobian_grad[2][1][1];
  d_eta_Jxn_dv[2] = eta_jacobian_grad[0][0][1]*eta_jacobian[1][1]+eta_jacobian[0][0]*eta_jacobian_grad[1][1][1]-
                    eta_jacobian_grad[0][1][1]*eta_jacobian[1][0]-eta_jacobian[0][1]*eta_jacobian_grad[1][0][1];
  //std::cout<<"dim singularity location: "<<eta<<std::endl;
  //std::cout<<"Spacedim singularity location: "<<eta_q_points_spacedim[0]<<std::endl;
  //std::cout<<"d_eta_Jxn_du "<<d_eta_Jxn_du<<std::endl;
  //std::cout<<"d_eta_Jxn_dv "<<d_eta_Jxn_dv<<std::endl;
//      std::cout<<"H1r: "<<eta_jacobian_grad[0][0][0]<<" "<<eta_jacobian_grad[0][0][1]<<std::endl;
//      std::cout<<"     "<<eta_jacobian_grad[0][1][0]<<" "<<eta_jacobian_grad[0][1][1]<<std::endl;
//      std::cout<<"H2r: "<<eta_jacobian_grad[1][0][0]<<" "<<eta_jacobian_grad[1][0][1]<<std::endl;
//      std::cout<<"     "<<eta_jacobian_grad[1][1][0]<<" "<<eta_jacobian_grad[1][1][1]<<std::endl;
//      std::cout<<"H3r: "<<eta_jacobian_grad[2][0][0]<<" "<<eta_jacobian_grad[2][0][1]<<std::endl;
//      std::cout<<"     "<<eta_jacobian_grad[2][1][0]<<" "<<eta_jacobian_grad[2][1][1]<<std::endl;

  
  // we now loop over the cell faces to carry out the integrals
  for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
      {
      // (to the best of my knowledge)
      // deal.ii does not offer an easy way to get the values of the
      // face vertices in the cell parametric plane
      // thus, we generate them here with a switch, given
      // the deal.ii numbering of faces and vertices
      // in generating them, we already center them in the
      // singularity point, subtracting P
      Point<2> csi_0;
      Point<2> csi_1;
      double coeff;
      switch (f)
          {
          case 0:
               csi_0(0) =  0-eta(0); csi_1(0) =  0-eta(0);
               csi_0(1) =  1-eta(1); csi_1(1) =  0-eta(1);
               coeff = csi_0(0);
          break;
          case 1:
               csi_0(0) =  1-eta(0); csi_1(0) =  1-eta(0);
               csi_0(1) =  0-eta(1); csi_1(1) =  1-eta(1);
               coeff = csi_0(0);
          break;
          case 2:
               csi_0(0) =  0-eta(0); csi_1(0) =  1-eta(0);
               csi_0(1) =  0-eta(1); csi_1(1) =  0-eta(1);
               coeff = csi_1(1);
          break;
          case 3:
               csi_0(0) =  1-eta(0); csi_1(0) =  0-eta(0);
               csi_0(1) =  1-eta(1); csi_1(1) =  1-eta(1);
               coeff = csi_1(1);
          break;
          default:
               csi_0(0) =  0-eta(0); csi_1(0) =  0-eta(0);
               csi_0(1) =  0-eta(1); csi_1(1) =  0-eta(1);
          break;
          }

      // if the singularity is located on the face, do not consider this face in
      // the integration
      if (!(((csi_0(0) == 0) && (csi_0(1) == 0)) || ((csi_1(0) == 0) && (csi_1(1) == 0))))
          { 
          // once the coordinates in the parametric plane of the
          // face vertices are known, we can compute the polar coordinate
          // theta angle, taking care of having always growing angles
          // on each face f
          double theta_0 = atan2((csi_0)[1],(csi_0)[0]);
          if (theta_0 < 0)
             theta_0 += 2*dealii::numbers::PI;
          double theta_1 = atan2((csi_1)[1],(csi_1)[0]);
          if (theta_1 < 0)
             theta_1 += 2*dealii::numbers::PI;
          //std::cout<<"PRE:   "<<f<<" ->    theta_0: "<<theta_0*180/dealii::numbers::PI<<"   theta_1: "<<theta_1*180/dealii::numbers::PI<<std::endl;   
          if (theta_1<theta_0)
             theta_0 -= 2*dealii::numbers::PI;
          //std::cout<<f<<" ->    theta_0: "<<theta_0*180/dealii::numbers::PI<<"   theta_1: "<<theta_1*180/dealii::numbers::PI<<std::endl;
          

          // we now create a 1D triangulation for the integration
          // in the theta direction of the polar coordinates
          Triangulation<1> theta_triangulation;

          std::vector<Point<1> > theta_vertices;
          std::vector<CellData<1> > theta_cells;
          SubCellData theta_subcelldata;
          
          // we only have one cell that goes from the minimum to maximum
          // theta on the face f 
          theta_vertices.resize(2);
          theta_vertices[0](0)=theta_0;
          theta_vertices[1](0)=theta_1;

          theta_cells.resize(1);

          theta_cells[0].vertices[0]=0;
          theta_cells[0].vertices[1]=1;
          
          GridTools::delete_unused_vertices (theta_vertices, theta_cells, theta_subcelldata);
          GridTools::consistently_order_cells(theta_cells);
          // here the triangulation is initialized and ready
          theta_triangulation.create_triangulation(theta_vertices, theta_cells, theta_subcelldata );

          // we of course also create a dof handler finite element and mapping
          DoFHandler<1> theta_dof_handler(theta_triangulation);
          const FE_Q<1> theta_finite_element(fe.degree);
          theta_dof_handler.distribute_dofs(theta_finite_element);
          MappingQ<1> theta_mapping(theta_finite_element.degree);

          // we create a 1D quadature with the user selected number of quadrature nodes
          QGauss<1> theta_quadrature(theta_quadrature_order);

          // we are now ready to create a FEValues objects that will allow us
          // to obtain all we need for the integration over theta on the face f
          FEValues<1> theta_fe_v(theta_mapping,
                                theta_finite_element,
                                theta_quadrature,
                                update_values |
                                update_quadrature_points | update_JxW_values);

          // we initialize theta_fe_v with the only face available
          auto theta_cell = theta_dof_handler.begin_active(); //this is the only cell
          theta_fe_v.reinit(theta_cell);


          // because we also will need to carry out the surface integral in polar coordinates
          // for I_0, we will also need a 1D triangulation for rho
          Triangulation<1> rho_triangulation;

          std::vector<Point<1> >  rho_vertices;
          std::vector<CellData<1> >  rho_cells;
          SubCellData  rho_subcelldata;
          
          // in this case we will just consider a cell that goes from 0 to 1 and
          // take care of the transformation and jacobian by ourself
          rho_vertices.resize(2);
          rho_vertices[0](0)=0;
          rho_vertices[1](0)=1;

          rho_cells.resize(1);

          rho_cells[0].vertices[0]=0;
          rho_cells[0].vertices[1]=1;
          
          GridTools::delete_unused_vertices (rho_vertices, rho_cells, rho_subcelldata);
          GridTools::consistently_order_cells(rho_cells);
          // here the triangulation is initialized and ready
          rho_triangulation.create_triangulation(rho_vertices, rho_cells, rho_subcelldata );
          // we of course also create a dof handler finite element and mapping
          DoFHandler<1> rho_dof_handler(rho_triangulation);
          const FE_Q<1> rho_finite_element(fe.degree);
          rho_dof_handler.distribute_dofs(rho_finite_element);
          MappingQ<1> rho_mapping(rho_finite_element.degree);
          // we create a 1D quadature with the user selected number of quadrature nodes
          QGauss<1> rho_quadrature(rho_quadrature_order);
          // we are now ready to create a FEValues object that will allow us
          // to obtain all we need for the integration over rho on the face f
          FEValues<1> rho_fe_v(rho_mapping,
                               rho_finite_element,
                               rho_quadrature,
                               update_values |
                               update_quadrature_points | update_JxW_values);
          // we initialize theta_fe_v with the only face available
          auto rho_cell = rho_dof_handler.begin_active(); //this is the only cell
          rho_fe_v.reinit(rho_cell);
          // we get the quadrature points for rho (remember that they go from 0 to 1)
          const std::vector<Point<1>> &rho_q_points    = rho_fe_v.get_quadrature_points();

          // we get the quadrature points for theta (remember that they go from
          // theta_min to theta_max in the cell)
          const std::vector<Point<1>> &theta_q_points    = theta_fe_v.get_quadrature_points();
          // we want to create a vector with the \hat(rho) values corresponding to the face boundary
          // at each theta value
          std::vector<double> q_rho(theta_q_points.size());
          // using deal.ii numbering, these lines do the trick
          for (unsigned int q=0; q<theta_q_points.size(); ++q)
              {
              double a = (f<2)?0.0:1.0;
              double b = (f<2)?1.0:0.0;
              q_rho[q] = coeff/(a*sin(theta_q_points[q](0))+b*(cos(theta_q_points[q](0))));
              }
          // we now create a set of quadrature points and weights to be passed to a
          // FEValues<2,3> to take care of the mapping to the 3D space
          const std::vector<double> &uv_q_weights = theta_fe_v.get_JxW_values();
          std::vector<Point<2> > uv_q_points(theta_q_points.size());
          for (unsigned int q=0; q<theta_q_points.size(); ++q)
              {
              // the quadrature points are the points of the face f in the parametric space
              uv_q_points[q] = Point<2>(q_rho[q]*cos(theta_q_points[q](0))+eta(0),
                                        q_rho[q]*sin(theta_q_points[q](0))+eta(1));
              }
          Quadrature<2> uv_quadrature(uv_q_points, uv_q_weights);
          //here we create the FEValues and initialize it with the <2,3> dof handler cell
          FEValues<2,3> face_fe_values(fe,
                                       uv_quadrature,
                                       update_values | update_gradients | update_quadrature_points |
                                       update_JxW_values | update_normal_vectors);
          face_fe_values.reinit(cell);

          // we now create a set of vectors that will contain useful theta dependent coefficients
          // for each value of theta on the face f
          std::vector<Tensor<1,3> > A(theta_q_points.size());
          std::vector<Tensor<1,3> > B(theta_q_points.size());
          // values of A.norm()
          std::vector<double> q_A(theta_q_points.size());
          // values of A*B
          std::vector<double> q_C(theta_q_points.size());
          // values of jacobian 0-th order Taylor series coefficient
          std::vector<Tensor<1,3> > Jk0(theta_q_points.size());
          // values of jacobian 1-st order Taylor series coefficient
          std::vector<Tensor<1,3> > Jk1(theta_q_points.size());
          // values of shape function 0-th order Taylor series coefficient
          std::vector<double> serv(fe.dofs_per_cell,0.0);
          std::vector<std::vector<double> > N0(theta_q_points.size(),serv);
          // values of shape function 1-st order Taylor series coefficient
          std::vector<std::vector<double> > N1(theta_q_points.size(),serv);
          for (unsigned int q=0; q<theta_q_points.size(); ++q)
              {
              for (unsigned int i=0; i<3; ++i)
                  {
                  A[q][i] = eta_jacobian[i][0]*cos(theta_q_points[q](0)) +
                         eta_jacobian[i][1]*sin(theta_q_points[q](0));
                  B[q][i] = eta_jacobian_grad[i][0][0]*pow(cos(theta_q_points[q](0)),2)/2 +
                            eta_jacobian_grad[i][0][1]*cos(theta_q_points[q](0))*sin(theta_q_points[q](0)) +
                            eta_jacobian_grad[i][1][1]*pow(sin(theta_q_points[q](0)),2)/2;
                  }
              q_A[q] = A[q].norm();
              q_C[q] = A[q]*B[q];
              Jk0[q] = eta_Jxn;
              Jk1[q] = d_eta_Jxn_du*cos(theta_q_points[q](0))+d_eta_Jxn_dv*sin(theta_q_points[q](0));
              for (unsigned int ii=0; ii<fe.dofs_per_cell; ++ii)
                  {
                  N0[q][ii] = eta_shape_values[ii];
                  N1[q][ii] = eta_shape_grads[ii][0]*cos(theta_q_points[q](0)) +
                             eta_shape_grads[ii][1]*sin(theta_q_points[q](0));
                  }
              }

          // so let's start computing the integrals on the face f
          // for I_1 and I_2
          for (unsigned int q=0; q<theta_q_points.size(); ++q)
              {
              // the taylor expansion coefficients F_1 and F_2 are finally computed
              double beta = 1/q_A[q];
              double gamma = -q_C[q]/pow(q_A[q],4);
              std::vector<Tensor<1,3> > F_2(fe.dofs_per_cell,zero_tensor);
              std::vector<Tensor<1,3> > F_1(fe.dofs_per_cell,zero_tensor);
              for (unsigned int ii=0; ii<fe.dofs_per_cell; ++ii)
                  {
                  F_2[ii] = Jk0[q]/pow(q_A[q],3)*N0[q][ii];
                  F_1[ii] = (-3*q_C[q]*Jk0[q]/pow(q_A[q],5)
                               -3*A[q]/pow(q_A[q],5)*(Jk0[q]*B[q]+Jk1[q]*A[q])
                               +Jk1[q]/pow(q_A[q],3))*N0[q][ii]+Jk0[q]/pow(q_A[q],3)*N1[q][ii];
                  //std::cout<<"f: "<<f<<"  q: "<<q<<" theta: "<<theta_q_points[q](0)*180/dealii::numbers::PI<<"  F_1: "<<F_1<<"  F_2: "<<F_2<<"  rho: "<<q_rho[q]<<"  beta: "<<beta<<"  gamma: "<<gamma<<"  w: "<<uv_q_weights[q]<<"  Jk0: "<<Jk0[q]<<"  Jk1: "<<Jk1[q]<<"  A: "<<A[q]<<"  B: "<<B[q]<<std::endl;
                  // and the integral contribution of the present quadrature node is here added
                  I_2[ii] += -F_2[ii]*(gamma/pow(beta,2)+1/q_rho[q])*uv_q_weights[q];
                  I_1[ii] += F_1[ii]*log(fabs(q_rho[q]/beta))*uv_q_weights[q];
                  //I_1 += F_1*log(fabs(q_rho[q]))*uv_q_weights[q];
                  }
              
              // we still miss the surface integral I_0
              // for this, we will need to also consider the
              // quadrature nodes of the rho_fe_values
              // remembering they go from 0 to 1, we must rescale them from
              // 0 to q_rho[q]. Here is the obvious scaling factor
              double rho_jac_fact = q_rho[q];

              // we must again create a set of quadrature points and weights to be passed to a
              // FEValues<2,3> to take care of the mapping to the 3D space
              
              // we only want the mapping jacobian from such a FEValues, so we set all the
              // weights to 1
              std::vector<double> inner_uv_q_weights(rho_q_points.size(),1.0);
              
              std::vector<Point<2> > inner_uv_q_points(rho_q_points.size());
              for (unsigned int p=0; p<rho_q_points.size(); ++p)
                  {
                  // based on rho and theta, the quadrature points are generated in polar coordinates
                  // centered at the singularity eta_q_points[0]
                  double rho = rho_jac_fact*rho_q_points[p](0);
                  inner_uv_q_points[p] = Point<2>(rho*cos(theta_q_points[q](0))+eta_q_points[0](0),
                                                  rho*sin(theta_q_points[q](0))+eta_q_points[0](1));
                  }
              
              // we can then create the quadrature and the FEValues
              Quadrature<2> inner_uv_quadrature(inner_uv_q_points, inner_uv_q_weights);
              FEValues<2,3> inner_face_fe_values(fe,
                                                 inner_uv_quadrature,
                                                 update_values | update_gradients | update_quadrature_points |
                                                 update_JxW_values | update_normal_vectors);
              inner_face_fe_values.reinit(cell);
              // let's obtain the quadrature points in the 3D space from the inner_face_fe_values
              const std::vector<Point<3>> &inner_uv_q_points_spacedim    = inner_face_fe_values.get_quadrature_points();
              const std::vector<Tensor<1, 3 >> &inner_uv_q_normals = inner_face_fe_values.get_normal_vectors();
              for (unsigned int p=0; p<rho_q_points.size(); ++p)
                  {
                  // this is the rho obtained scaling the quadrature point in the [0,1] interval
                  double rho = rho_jac_fact*rho_q_points[p](0);
                  // this is the distance in the 3D space between quadrature point and singularity
                  Tensor<1,3> R = inner_uv_q_points_spacedim[p]-eta_q_points_spacedim[0];
                  double r = R.norm();
                  // finally, this is the integral argument
                  // the first term also needs the jacobian of the 3D mapping, which was instead
                  // already included in the F_2 and F_1 factors of the latter terms, and is therefore
                  // not needed there
                  for (unsigned int ii=0; ii<fe.dofs_per_cell; ++ii)
                      {
                      Tensor<1,3> fun = -1/pow(r,3)*(3*(R/r)*(R*inner_uv_q_normals[p]/r)-inner_uv_q_normals[p])*
                                        inner_face_fe_values.shape_value(ii,p);
                      Tensor<1,3> arg = fun*rho*inner_face_fe_values.JxW(p)-F_2[ii]/pow(rho,2)-F_1[ii]/rho;

//                          std::cout<<"rho: "<<rho<<"  r: "<<r<<"  arg1: "<<fun*rho*inner_face_fe_values.JxW(p)<<
//                      "   arg2: "<<-F_2[0]/pow(rho,2)-F_1[0]/rho<<"  diff: "<<fun*rho*inner_face_fe_values.JxW(p)-F_2[0]/pow(rho,2)-F_1[0]/rho<<std::endl;

                      //std::cout<<rho<<"  "<<r<<"  "<<1/pow(r,3)*rho*inner_face_fe_values.JxW(p)<<
                      //"  "<<-F_2/pow(rho,2)-F_1/rho<<"  "<<1/pow(r,3)*rho*inner_face_fe_values.JxW(p)-F_2/pow(rho,2)-F_1/rho<<std::endl;

    //                      std::cout<<rho<<"  "<<r<<"  "<<Jk0[q][2]+rho*Jk1[q][2]<<"  "<<inner_uv_q_normals[p][2]*inner_face_fe_values.JxW(p)<<"   "<<inner_uv_q_normals[p][2]*inner_face_fe_values.JxW(p)-Jk0[q][2]-rho*Jk1[q][2]<<std::endl;

    //                  std::cout<<rho<<"  "<<R<<"  "<<rho*A[q]+rho*rho*B[q]<<"  "<<(R-rho*A[q]-rho*rho*B[q]).norm()<<std::endl;
                  
                      // we add the conribution of the present node to the surface integral
                      // the inner weights are given by rho_fe_v.JxW(p) but the jacobian is given by the scaling
                      // factor rho_jac_fact
                      I_0[ii] += arg*rho_jac_fact*rho_fe_v.JxW(p)*theta_fe_v.JxW(q);
                      }
                  }
              }
          /*
          for (unsigned int ii=0; ii<fe.dofs_per_cell; ++ii)
              {
              std::cout<<ii<<"  temp I_0 = "<<I_0[ii]<<std::endl;
              std::cout<<ii<<"  temp I_1 = "<<I_1[ii]<<std::endl;
              std::cout<<ii<<"  temp I_2 = "<<I_2[ii]<<std::endl;
              }
          */
          }
      }


std::vector<Tensor<1,3> > II(fe.dofs_per_cell,zero_tensor);
for (unsigned int ii=0; ii<fe.dofs_per_cell; ++ii)
    {
    II[ii] = (I_0[ii] + I_1[ii] + I_2[ii])/numbers::PI/4.0;
    }

 return II;
};

template<>
std::vector<Tensor<1,2> > SingularKernelIntegral<2>::evaluate_WkNj_integrals()
{
ExcNotImplemented();
std::vector<Tensor<1,2> > dummy;
return dummy;
};

template<>
std::vector<Tensor<1,3> > SingularKernelIntegral<3>::evaluate_WkNj_integrals()
{
  Tensor<1,3> zero_tensor;
  std::vector<Tensor<1,3> > I_0(fe.dofs_per_cell,zero_tensor);
  std::vector<Tensor<1,3> > I_1(fe.dofs_per_cell,zero_tensor);
  std::vector<Tensor<1,3> > I_2(fe.dofs_per_cell,zero_tensor);

  // We must obtain the shape function derivatives at the singularity point eta in
  // the parametric plane.
  // In the deal.II framework the (to my knowledge) best way is to
  // create a reference Triangulation<dim,spacedim> and FE<dim,spacedim> and use it
  // to get the shape function values and derivatives
  // building reference cell
  Triangulation<2,3> ref_triangulation;

  std::vector<Point<3> > ref_vertices;
  std::vector<CellData<2> > ref_cells;
  SubCellData ref_subcelldata;

  ref_vertices.resize(4);
  ref_vertices[0](0)= 0.0;
  ref_vertices[0](1)= 0.0;
  ref_vertices[0](2)=0.0;
  ref_vertices[1](0)= 1.0;
  ref_vertices[1](1)= 0.0;
  ref_vertices[1](2)=0.0;
  ref_vertices[2](0)= 0.0;
  ref_vertices[2](1)= 1.0;
  ref_vertices[2](2)=0.0;
  ref_vertices[3](0)= 1.0;
  ref_vertices[3](1)= 1.0;
  ref_vertices[3](2)=0.0;

  ref_cells.resize(1);

  ref_cells[0].vertices[0]=0;
  ref_cells[0].vertices[1]=1;
  ref_cells[0].vertices[2]=2;
  ref_cells[0].vertices[3]=3;
  ref_cells[0].material_id = 1;

  GridTools::delete_unused_vertices (ref_vertices, ref_cells, ref_subcelldata);
  GridTools::consistently_order_cells(ref_cells);

  // here's the triangulation set up
  ref_triangulation.create_triangulation(ref_vertices, ref_cells, ref_subcelldata );

  // we will need a codimension one finite element and
  // the corresponding dof handler
  DoFHandler<2,3> ref_dof_handler(ref_triangulation);
  ref_dof_handler.distribute_dofs(fe);


  // we now need to compute all the mapping first and second order derivatives
  // in correspondence with the singularity
  // we have to pass through a FEValues object, but to build one we will 
  // need std::vectors with location of quadrature points and quadrature nodes

  // we have a single quadrature point (eta) located at P
  std::vector<Point<2> > eta_q_points(1);
  eta_q_points[0] = eta;
  // and a single quadrature weight set to one
  std::vector<double> eta_q_weights(1);
  eta_q_weights[0] = 1;

  // here's the quadrature rule obtained with the point and weight generated
  Quadrature<2> eta_quadrature(eta_q_points, eta_q_weights);
  // and here are the reference and spacedim FEValues class resulting by it
  // (both FEValues are initialized with the correct cell)
  
  FEValues<2,3> ref_eta_fe_values(fe,
                                  eta_quadrature,
                                  update_values | update_gradients | update_quadrature_points |
                                  update_JxW_values | update_normal_vectors | update_jacobians |
                                  update_jacobian_grads);
  // we initialize rference FEValues object on the reference cell
  ref_eta_fe_values.reinit(ref_dof_handler.begin_active());
  
  FEValues<2,3> eta_fe_values(mapping,
                              fe,
                              eta_quadrature,
                              update_values | update_gradients | update_quadrature_points |
                              update_JxW_values | update_normal_vectors | update_jacobians |
                              update_jacobian_grads);
  // we initialize FEValues object on the current cell
  eta_fe_values.reinit(cell);
  
  // from the reference FEValues we get the n_dofs shape function values and gradients at the
  // singularity point
  std::vector<double> eta_shape_values(fe.dofs_per_cell);
  std::vector<Tensor<1,3> > eta_shape_grads(fe.dofs_per_cell);
  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      {
      eta_shape_values[i] = ref_eta_fe_values.shape_value(i,0);
      eta_shape_grads[i] = ref_eta_fe_values.shape_grad(i,0);
      }
  
  // the single quadrature point is the point in the three dimensional domain corresponding to eta/P
  const std::vector<Point<3>> &eta_q_points_spacedim    = eta_fe_values.get_quadrature_points();
  // we also get the acobian and jacobian gradient at such a location
  auto eta_jacobian = eta_fe_values.jacobian(0);
  auto eta_jacobian_grad = eta_fe_values.jacobian_grad(0);
  // we now have all the first and second order mapping derivaives
  // the product of the jacobian by the normal vector Jxn is the surface normal vector
  // let us obtain it in terms of the mapping derivatives
  Tensor<1,3> eta_Jxn;
  eta_Jxn[0] = eta_jacobian[1][0]*eta_jacobian[2][1]-eta_jacobian[1][1]*eta_jacobian[2][0];
  eta_Jxn[1] = eta_jacobian[0][1]*eta_jacobian[2][0]-eta_jacobian[0][0]*eta_jacobian[2][1];
  eta_Jxn[2] = eta_jacobian[0][0]*eta_jacobian[1][1]-eta_jacobian[0][1]*eta_jacobian[1][0];
  // we also want the derivative of the previous vector with respect to the first parametric plane coordinate u
  Tensor<1,3> d_eta_Jxn_du;
  d_eta_Jxn_du[0] = eta_jacobian_grad[1][0][0]*eta_jacobian[2][1]+eta_jacobian[1][0]*eta_jacobian_grad[2][1][0]-
                    eta_jacobian_grad[1][1][0]*eta_jacobian[2][0]-eta_jacobian[1][1]*eta_jacobian_grad[2][0][0];
  d_eta_Jxn_du[1] = eta_jacobian_grad[0][1][0]*eta_jacobian[2][0]+eta_jacobian[0][1]*eta_jacobian_grad[2][0][0]-
                    eta_jacobian_grad[0][0][0]*eta_jacobian[2][1]-eta_jacobian[0][0]*eta_jacobian_grad[2][1][0];
  d_eta_Jxn_du[2] = eta_jacobian_grad[0][0][0]*eta_jacobian[1][1]+eta_jacobian[0][0]*eta_jacobian_grad[1][1][0]-
                    eta_jacobian_grad[0][1][0]*eta_jacobian[1][0]-eta_jacobian[0][1]*eta_jacobian_grad[1][0][0];
  // we also want the derivative of the previous vector with respect to the second parametric plane coordinate v
  Tensor<1,3> d_eta_Jxn_dv;
  d_eta_Jxn_dv[0] = eta_jacobian_grad[1][0][1]*eta_jacobian[2][1]+eta_jacobian[1][0]*eta_jacobian_grad[2][1][1]-
                    eta_jacobian_grad[1][1][1]*eta_jacobian[2][0]-eta_jacobian[1][1]*eta_jacobian_grad[2][0][1];
  d_eta_Jxn_dv[1] = eta_jacobian_grad[0][1][1]*eta_jacobian[2][0]+eta_jacobian[0][1]*eta_jacobian_grad[2][0][1]-
                    eta_jacobian_grad[0][0][1]*eta_jacobian[2][1]-eta_jacobian[0][0]*eta_jacobian_grad[2][1][1];
  d_eta_Jxn_dv[2] = eta_jacobian_grad[0][0][1]*eta_jacobian[1][1]+eta_jacobian[0][0]*eta_jacobian_grad[1][1][1]-
                    eta_jacobian_grad[0][1][1]*eta_jacobian[1][0]-eta_jacobian[0][1]*eta_jacobian_grad[1][0][1];
  //std::cout<<"dim singularity location: "<<eta<<std::endl;
  //std::cout<<"Spacedim singularity location: "<<eta_q_points_spacedim[0]<<std::endl;
  //std::cout<<"d_eta_Jxn_du "<<d_eta_Jxn_du<<std::endl;
  //std::cout<<"d_eta_Jxn_dv "<<d_eta_Jxn_dv<<std::endl;
//      std::cout<<"H1r: "<<eta_jacobian_grad[0][0][0]<<" "<<eta_jacobian_grad[0][0][1]<<std::endl;
//      std::cout<<"     "<<eta_jacobian_grad[0][1][0]<<" "<<eta_jacobian_grad[0][1][1]<<std::endl;
//      std::cout<<"H2r: "<<eta_jacobian_grad[1][0][0]<<" "<<eta_jacobian_grad[1][0][1]<<std::endl;
//      std::cout<<"     "<<eta_jacobian_grad[1][1][0]<<" "<<eta_jacobian_grad[1][1][1]<<std::endl;
//      std::cout<<"H3r: "<<eta_jacobian_grad[2][0][0]<<" "<<eta_jacobian_grad[2][0][1]<<std::endl;
//      std::cout<<"     "<<eta_jacobian_grad[2][1][0]<<" "<<eta_jacobian_grad[2][1][1]<<std::endl;

  double J0 = eta_Jxn.norm();
  double J0_du = 1/eta_Jxn.norm()*eta_Jxn*d_eta_Jxn_du;
  double J0_dv = 1/eta_Jxn.norm()*eta_Jxn*d_eta_Jxn_dv;
  // we now loop over the cell faces to carry out the integrals
  for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
      {
      // (to the best of my knowledge)
      // deal.ii does not offer an easy way to get the values of the
      // face vertices in the cell parametric plane
      // thus, we generate them here with a switch, given
      // the deal.ii numbering of faces and vertices
      // in generating them, we already center them in the
      // singularity point, subtracting P
      Point<2> csi_0;
      Point<2> csi_1;
      double coeff;
      switch (f)
          {
          case 0:
               csi_0(0) =  0-eta(0); csi_1(0) =  0-eta(0);
               csi_0(1) =  1-eta(1); csi_1(1) =  0-eta(1);
               coeff = csi_0(0);
          break;
          case 1:
               csi_0(0) =  1-eta(0); csi_1(0) =  1-eta(0);
               csi_0(1) =  0-eta(1); csi_1(1) =  1-eta(1);
               coeff = csi_0(0);
          break;
          case 2:
               csi_0(0) =  0-eta(0); csi_1(0) =  1-eta(0);
               csi_0(1) =  0-eta(1); csi_1(1) =  0-eta(1);
               coeff = csi_1(1);
          break;
          case 3:
               csi_0(0) =  1-eta(0); csi_1(0) =  0-eta(0);
               csi_0(1) =  1-eta(1); csi_1(1) =  1-eta(1);
               coeff = csi_1(1);
          break;
          default:
               csi_0(0) =  0-eta(0); csi_1(0) =  0-eta(0);
               csi_0(1) =  0-eta(1); csi_1(1) =  0-eta(1);
          break;
          }

      // if the singularity is located on the face, do not consider this face in
      // the integration
      if (!(((csi_0(0) == 0) && (csi_0(1) == 0)) || ((csi_1(0) == 0) && (csi_1(1) == 0))))
          { 

          // once the coordinates in the parametric plane of the
          // face vertices are known, we can compute the polar coordinate
          // theta angle, taking care of having always growing angles
          // on each face f
          double theta_0 = atan2((csi_0)[1],(csi_0)[0]);
          if (theta_0 < 0)
             theta_0 += 2*dealii::numbers::PI;
          double theta_1 = atan2((csi_1)[1],(csi_1)[0]);
          if (theta_1 < 0)
             theta_1 += 2*dealii::numbers::PI;
          //std::cout<<"PRE:   "<<f<<" ->    theta_0: "<<theta_0*180/dealii::numbers::PI<<"   theta_1: "<<theta_1*180/dealii::numbers::PI<<std::endl;   
          if (theta_1<theta_0)
             theta_0 -= 2*dealii::numbers::PI;
          //std::cout<<f<<" ->    theta_0: "<<theta_0*180/dealii::numbers::PI<<"   theta_1: "<<theta_1*180/dealii::numbers::PI<<std::endl;
          

          // we now create a 1D triangulation for the integration
          // in the theta direction of the polar coordinates
          Triangulation<1> theta_triangulation;

          std::vector<Point<1> > theta_vertices;
          std::vector<CellData<1> > theta_cells;
          SubCellData theta_subcelldata;
          
          // we only have one cell that goes from the minimum to maximum
          // theta on the face f 
          theta_vertices.resize(2);
          theta_vertices[0](0)=theta_0;
          theta_vertices[1](0)=theta_1;

          theta_cells.resize(1);

          theta_cells[0].vertices[0]=0;
          theta_cells[0].vertices[1]=1;
          
          GridTools::delete_unused_vertices (theta_vertices, theta_cells, theta_subcelldata);
          GridTools::consistently_order_cells(theta_cells);
          // here the triangulation is initialized and ready
          theta_triangulation.create_triangulation(theta_vertices, theta_cells, theta_subcelldata );

          // we of course also create a dof handler finite element and mapping
          DoFHandler<1> theta_dof_handler(theta_triangulation);
          const FE_Q<1> theta_finite_element(fe.degree);
          theta_dof_handler.distribute_dofs(theta_finite_element);
          MappingQ<1> theta_mapping(theta_finite_element.degree);

          // we create a 1D quadature with the user selected number of quadrature nodes
          QGauss<1> theta_quadrature(theta_quadrature_order);

          // we are now ready to create a FEValues objects that will allow us
          // to obtain all we need for the integration over theta on the face f
          FEValues<1> theta_fe_v(theta_mapping,
                                theta_finite_element,
                                theta_quadrature,
                                update_values |
                                update_quadrature_points | update_JxW_values);

          // we initialize theta_fe_v with the only face available
          auto theta_cell = theta_dof_handler.begin_active(); //this is the only cell
          theta_fe_v.reinit(theta_cell);


          // because we also will need to carry out the surface integral in polar coordinates
          // for I_0, we will also need a 1D triangulation for rho
          Triangulation<1> rho_triangulation;

          std::vector<Point<1> >  rho_vertices;
          std::vector<CellData<1> >  rho_cells;
          SubCellData  rho_subcelldata;
          
          // in this case we will just consider a cell that goes from 0 to 1 and
          // take care of the transformation and jacobian by ourself
          rho_vertices.resize(2);
          rho_vertices[0](0)=0;
          rho_vertices[1](0)=1;

          rho_cells.resize(1);

          rho_cells[0].vertices[0]=0;
          rho_cells[0].vertices[1]=1;
          
          GridTools::delete_unused_vertices (rho_vertices, rho_cells, rho_subcelldata);
          GridTools::consistently_order_cells(rho_cells);
          // here the triangulation is initialized and ready
          rho_triangulation.create_triangulation(rho_vertices, rho_cells, rho_subcelldata );
          // we of course also create a dof handler finite element and mapping
          DoFHandler<1> rho_dof_handler(rho_triangulation);
          const FE_Q<1> rho_finite_element(fe.degree);
          rho_dof_handler.distribute_dofs(rho_finite_element);
          MappingQ<1> rho_mapping(rho_finite_element.degree);
          // we create a 1D quadature with the user selected number of quadrature nodes
          QGauss<1> rho_quadrature(rho_quadrature_order);
          // we are now ready to create a FEValues object that will allow us
          // to obtain all we need for the integration over rho on the face f
          FEValues<1> rho_fe_v(rho_mapping,
                               rho_finite_element,
                               rho_quadrature,
                               update_values |
                               update_quadrature_points | update_JxW_values);
          // we initialize theta_fe_v with the only face available
          auto rho_cell = rho_dof_handler.begin_active(); //this is the only cell
          rho_fe_v.reinit(rho_cell);
          // we get the quadrature points for rho (remember that they go from 0 to 1)
          const std::vector<Point<1>> &rho_q_points    = rho_fe_v.get_quadrature_points();

          // we get the quadrature points for theta (remember that they go from
          // theta_min to theta_max in the cell)
          const std::vector<Point<1>> &theta_q_points    = theta_fe_v.get_quadrature_points();
          // we want to create a vector with the \hat(rho) values corresponding to the face boundary
          // at each theta value
          std::vector<double> q_rho(theta_q_points.size());
          // using deal.ii numbering, these lines do the trick
          for (unsigned int q=0; q<theta_q_points.size(); ++q)
              {
              double a = (f<2)?0.0:1.0;
              double b = (f<2)?1.0:0.0;
              q_rho[q] = coeff/(a*sin(theta_q_points[q](0))+b*(cos(theta_q_points[q](0))));
              }
          // we now create a set of quadrature points and weights to be passed to a
          // FEValues<2,3> to take care of the mapping to the 3D space
          const std::vector<double> &uv_q_weights = theta_fe_v.get_JxW_values();
          std::vector<Point<2> > uv_q_points(theta_q_points.size());
          for (unsigned int q=0; q<theta_q_points.size(); ++q)
              {
              // the quadrature points are the points of the face f in the parametric space
              uv_q_points[q] = Point<2>(q_rho[q]*cos(theta_q_points[q](0))+eta(0),
                                        q_rho[q]*sin(theta_q_points[q](0))+eta(1));
              }
          Quadrature<2> uv_quadrature(uv_q_points, uv_q_weights);
          //here we create the FEValues and initialize it with the <2,3> dof handler cell
          FEValues<2,3> face_fe_values(fe,
                                       uv_quadrature,
                                       update_values | update_gradients | update_quadrature_points |
                                       update_JxW_values | update_normal_vectors);
          face_fe_values.reinit(cell);

          // we now create a set of vectors that will contain useful theta dependent coefficients
          // for each value of theta on the face f
          std::vector<Tensor<1,3> > A(theta_q_points.size());
          std::vector<Tensor<1,3> > B(theta_q_points.size());
          // values of A.norm()
          std::vector<double> q_A(theta_q_points.size());
          // values of A*B
          std::vector<double> q_C(theta_q_points.size());
          // values of jacobian-k 0-th order Taylor series coefficient
          std::vector<Tensor<1,3> > Jk0(theta_q_points.size());
          // values of jacobian-k 1-st order Taylor series coefficient
          std::vector<Tensor<1,3> > Jk1(theta_q_points.size());
          // values of jacobian 1-st order Taylor series coefficient
          std::vector<double > J1(theta_q_points.size());
          // values of shape function 0-th order Taylor series coefficient
          std::vector<double> serv(fe.dofs_per_cell,0.0);
          std::vector<std::vector<double> > N0(theta_q_points.size(),serv);
          // values of shape function 1-st order Taylor series coefficient
          std::vector<std::vector<double> > N1(theta_q_points.size(),serv);
          for (unsigned int q=0; q<theta_q_points.size(); ++q)
              {
              for (unsigned int i=0; i<3; ++i)
                  {
                  A[q][i] = eta_jacobian[i][0]*cos(theta_q_points[q](0)) +
                         eta_jacobian[i][1]*sin(theta_q_points[q](0));
                  B[q][i] = eta_jacobian_grad[i][0][0]*pow(cos(theta_q_points[q](0)),2)/2 +
                            eta_jacobian_grad[i][0][1]*cos(theta_q_points[q](0))*sin(theta_q_points[q](0)) +
                            eta_jacobian_grad[i][1][1]*pow(sin(theta_q_points[q](0)),2)/2;
                  }
              q_A[q] = A[q].norm();
              q_C[q] = A[q]*B[q];
              Jk1[q] = d_eta_Jxn_du*cos(theta_q_points[q](0))+d_eta_Jxn_dv*sin(theta_q_points[q](0));
              J1[q] = J0_du*cos(theta_q_points[q](0))+J0_dv*sin(theta_q_points[q](0));
              for (unsigned int ii=0; ii<fe.dofs_per_cell; ++ii)
                  {
                  N0[q][ii] = eta_shape_values[ii];
                  N1[q][ii] = eta_shape_grads[ii][0]*cos(theta_q_points[q](0)) +
                             eta_shape_grads[ii][1]*sin(theta_q_points[q](0));
                  }
              }

          // so let's start computing the integrals on the face f
          // for I_1 and I_2
          for (unsigned int q=0; q<theta_q_points.size(); ++q)
              {
              // the taylor expansion coefficients F_1 and F_2 are finally computed
              double beta = 1/q_A[q];
              //double gamma = -q_C[q]/pow(q_A[q],4);
              std::vector<Tensor<1,3> > F_2(fe.dofs_per_cell,zero_tensor);
              std::vector<Tensor<1,3> > F_1(fe.dofs_per_cell,zero_tensor);
              for (unsigned int ii=0; ii<fe.dofs_per_cell; ++ii)
                  {
                  F_2[ii] = 0;
                  F_1[ii] = A[q]*J0/pow(q_A[q],3)*N0[q][ii];
                  //std::cout<<"f: "<<f<<"  q: "<<q<<" theta: "<<theta_q_points[q](0)*180/dealii::numbers::PI<<"  F_1: "<<F_1[ii]<<"  F_2: "<<F_2[ii]<<"  rho: "<<q_rho[q]<<"  beta: "<<beta<<"  gamma: "<<gamma<<"  w: "<<uv_q_weights[q]<<"  Jk0: "<<Jk0[q]<<"  Jk1: "<<Jk1[q]<<"  A: "<<A[q]<<"  B: "<<B[q]<<std::endl;
                  // and the integral contribution of the present quadrature node is here added
                  I_1[ii] += F_1[ii]*log(fabs(q_rho[q]/beta))*uv_q_weights[q];
                  //I_1 += F_1*log(fabs(q_rho[q]))*uv_q_weights[q];
                  }
              
              // we still miss the surface integral I_0
              // for this, we will need to also consider the
              // quadrature nodes of the rho_fe_values
              // remembering they go from 0 to 1, we must rescale them from
              // 0 to q_rho[q]. Here is the obvious scaling factor
              double rho_jac_fact = q_rho[q];

              // we must again create a set of quadrature points and weights to be passed to a
              // FEValues<2,3> to take care of the mapping to the 3D space
              
              // we only want the mapping jacobian from such a FEValues, so we set all the
              // weights to 1
              std::vector<double> inner_uv_q_weights(rho_q_points.size(),1.0);
              
              std::vector<Point<2> > inner_uv_q_points(rho_q_points.size());
              for (unsigned int p=0; p<rho_q_points.size(); ++p)
                  {
                  // based on rho and theta, the quadrature points are generated in polar coordinates
                  // centered at the singularity eta_q_points[0]
                  double rho = rho_jac_fact*rho_q_points[p](0);
                  inner_uv_q_points[p] = Point<2>(rho*cos(theta_q_points[q](0))+eta_q_points[0](0),
                                                  rho*sin(theta_q_points[q](0))+eta_q_points[0](1));
                  }
              
              // we can then create the quadrature and the FEValues
              Quadrature<2> inner_uv_quadrature(inner_uv_q_points, inner_uv_q_weights);
              FEValues<2,3> inner_face_fe_values(fe,
                                                 inner_uv_quadrature,
                                                 update_values | update_gradients | update_quadrature_points |
                                                 update_JxW_values | update_normal_vectors);
              inner_face_fe_values.reinit(cell);
              // let's obtain the quadrature points in the 3D space from the inner_face_fe_values
              const std::vector<Point<3>> &inner_uv_q_points_spacedim    = inner_face_fe_values.get_quadrature_points();
              //const std::vector<Tensor<1, 3 >> &inner_uv_q_normals = inner_face_fe_values.get_normal_vectors();
              for (unsigned int p=0; p<rho_q_points.size(); ++p)
                  {
                  // this is the rho obtained scaling the quadrature point in the [0,1] interval
                  double rho = rho_jac_fact*rho_q_points[p](0);
                  // this is the distance in the 3D space between quadrature point and singularity
                  Tensor<1,3> R = inner_uv_q_points_spacedim[p]-eta_q_points_spacedim[0];
                  double r = R.norm();
                  // finally, this is the integral argument
                  // the first term also needs the jacobian of the 3D mapping, which was instead
                  // already included in the F_2 and F_1 factors of the latter terms, and is therefore
                  // not needed there
                  for (unsigned int ii=0; ii<fe.dofs_per_cell; ++ii)
                      {
                      Tensor<1,3> fun = 1/pow(r,3)*(R)*
                                        inner_face_fe_values.shape_value(ii,p);
                      Tensor<1,3> arg = fun*rho*inner_face_fe_values.JxW(p)-F_2[ii]/pow(rho,2)-F_1[ii]/rho;

//                          std::cout<<"rho: "<<rho<<"  r: "<<r<<"  arg1: "<<1/pow(r,3)*R*rho*inner_face_fe_values.JxW(p)<<
//                      "   arg2: "<<-F_2[ii]/pow(rho,2)-F_1[ii]/rho<<"  diff: "<<1/pow(r,2)*R*rho*inner_face_fe_values.JxW(p)-F_2[ii]/pow(rho,2)-F_1[ii]/rho<<std::endl;

                      //std::cout<<rho<<"  "<<r<<"  "<<1/pow(r,3)*rho*inner_face_fe_values.JxW(p)<<
                      //"  "<<-F_2/pow(rho,2)-F_1/rho<<"  "<<1/pow(r,3)*rho*inner_face_fe_values.JxW(p)-F_2/pow(rho,2)-F_1/rho<<std::endl;

    //                      std::cout<<rho<<"  "<<r<<"  "<<Jk0[q][2]+rho*Jk1[q][2]<<"  "<<inner_uv_q_normals[p][2]*inner_face_fe_values.JxW(p)<<"   "<<inner_uv_q_normals[p][2]*inner_face_fe_values.JxW(p)-Jk0[q][2]-rho*Jk1[q][2]<<std::endl;

    //                  std::cout<<rho<<"  "<<R<<"  "<<rho*A[q]+rho*rho*B[q]<<"  "<<(R-rho*A[q]-rho*rho*B[q]).norm()<<std::endl;
                  
                      // we add the conribution of the present node to the surface integral
                      // the inner weights are given by rho_fe_v.JxW(p) but the jacobian is given by the scaling
                      // factor rho_jac_fact
                      I_0[ii] += arg*rho_jac_fact*rho_fe_v.JxW(p)*theta_fe_v.JxW(q);
                      }
                  }
              }
          /*for (unsigned int ii=0; ii<fe.dofs_per_cell; ++ii)
              {
              std::cout<<ii<<"  temp I_0 = "<<I_0[ii]<<std::endl;
              std::cout<<ii<<"  temp I_1 = "<<I_1[ii]<<std::endl;
              std::cout<<ii<<"  temp I_2 = "<<I_2[ii]<<std::endl;
              }
          */
          }
      }


std::vector<Tensor<1,3> > II(fe.dofs_per_cell,zero_tensor);
for (unsigned int ii=0; ii<fe.dofs_per_cell; ++ii)
    {
    II[ii] = (I_0[ii] + I_1[ii] + I_2[ii])/numbers::PI/4.0;
    }

 return II;
};

template<>
double SingularKernelIntegral<2>::evaluate_integral()
{
ExcNotImplemented();
return 0.0;
};

template<>
double SingularKernelIntegral<3>::evaluate_integral()
{
  double I_0 = 0.0;
  double I_1 = 0.0;
  double I_2 = 0.0;

  // We must obtain the shape function derivatives at the singularity point eta in
  // the parametric plane.
  // In the deal.II framework the (to my knowledge) best way is to
  // create a reference Triangulation<dim,spacedim> and FE<dim,spacedim> and use it
  // to get the shape function values and derivatives
  // building reference cell
  Triangulation<2,3> ref_triangulation;

  std::vector<Point<3> > ref_vertices;
  std::vector<CellData<2> > ref_cells;
  SubCellData ref_subcelldata;

  ref_vertices.resize(4);
  ref_vertices[0](0)= 0.0;
  ref_vertices[0](1)= 0.0;
  ref_vertices[0](2)=0.0;
  ref_vertices[1](0)= 1.0;
  ref_vertices[1](1)= 0.0;
  ref_vertices[1](2)=0.0;
  ref_vertices[2](0)= 0.0;
  ref_vertices[2](1)= 1.0;
  ref_vertices[2](2)=0.0;
  ref_vertices[3](0)= 1.0;
  ref_vertices[3](1)= 1.0;
  ref_vertices[3](2)=0.0;

  ref_cells.resize(1);

  ref_cells[0].vertices[0]=0;
  ref_cells[0].vertices[1]=1;
  ref_cells[0].vertices[2]=3;
  ref_cells[0].vertices[3]=2;
  ref_cells[0].material_id = 1;

  GridTools::delete_unused_vertices (ref_vertices, ref_cells, ref_subcelldata);
  GridTools::consistently_order_cells(ref_cells);

  // here's the triangulation set up
  ref_triangulation.create_triangulation(ref_vertices, ref_cells, ref_subcelldata );

  // we will need a codimension one finite element and
  // the corresponding dof handler
  DoFHandler<2,3> ref_dof_handler(ref_triangulation);
  ref_dof_handler.distribute_dofs(fe);


  // we now need to compute all the mapping first and second order derivatives
  // in correspondence with the singularity
  // we have to pass through a FEValues object, but to build one we will 
  // need std::vectors with location of quadrature points and quadrature nodes

  // we have a single quadrature point (eta) located at P
  std::vector<Point<2> > eta_q_points(1);
  eta_q_points[0] = eta;
  // and a single quadrature weight set to one
  std::vector<double> eta_q_weights(1);
  eta_q_weights[0] = 1;

  // here's the quadrature rule obtained with the point and weight generated
  Quadrature<2> eta_quadrature(eta_q_points, eta_q_weights);
  // and here are the reference and spacedim FEValues class resulting by it
  // (both FEValues are initialized with the correct cell)
  
  FEValues<2,3> ref_eta_fe_values(fe,
                                  eta_quadrature,
                                  update_values | update_gradients | update_quadrature_points |
                                  update_JxW_values | update_normal_vectors | update_jacobians |
                                  update_jacobian_grads);
  // we initialize rference FEValues object on the reference cell
  ref_eta_fe_values.reinit(ref_dof_handler.begin_active());
  
  FEValues<2,3> eta_fe_values(mapping,
                              fe,
                              eta_quadrature,
                              update_values | update_gradients | update_quadrature_points |
                              update_JxW_values | update_normal_vectors | update_jacobians |
                              update_jacobian_grads);
  // we initialize FEValues object on the current cell
  eta_fe_values.reinit(cell);
  
  // from the reference FEValues we get the n_dofs shape function values and gradients at the
  // singularity point
  std::vector<double> eta_shape_values(fe.dofs_per_cell);
  std::vector<Tensor<1,3> > eta_shape_grads(fe.dofs_per_cell);
  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      {
      eta_shape_values[i] = ref_eta_fe_values.shape_value(i,0);
      eta_shape_grads[i] = ref_eta_fe_values.shape_grad(i,0);
      }
  
  // the single quadrature point is the point in the three dimensional domain corresponding to eta/P
  const std::vector<Point<3>> &eta_q_points_spacedim    = eta_fe_values.get_quadrature_points();
  // we also get the acobian and jacobian gradient at such a location
  auto eta_jacobian = eta_fe_values.jacobian(0);
  auto eta_jacobian_grad = eta_fe_values.jacobian_grad(0);
  // we now have all the first and second order mapping derivaives
  // the product of the jacobian by the normal vector Jxn is the surface normal vector
  // let us obtain it in terms of the mapping derivatives
  Tensor<1,3> eta_Jxn;
  eta_Jxn[0] = eta_jacobian[1][0]*eta_jacobian[2][1]-eta_jacobian[1][1]*eta_jacobian[2][0];
  eta_Jxn[1] = eta_jacobian[0][1]*eta_jacobian[2][0]-eta_jacobian[0][0]*eta_jacobian[2][1];
  eta_Jxn[2] = eta_jacobian[0][0]*eta_jacobian[1][1]-eta_jacobian[0][1]*eta_jacobian[1][0];
  // we also want the derivative of the previous vector with respect to the first parametric plane coordinate u
  Tensor<1,3> d_eta_Jxn_du;
  d_eta_Jxn_du[0] = eta_jacobian_grad[1][0][0]*eta_jacobian[2][1]+eta_jacobian[1][0]*eta_jacobian_grad[2][1][0]-
                    eta_jacobian_grad[1][1][0]*eta_jacobian[2][0]-eta_jacobian[1][1]*eta_jacobian_grad[2][0][0];
  d_eta_Jxn_du[1] = eta_jacobian_grad[0][1][0]*eta_jacobian[2][0]+eta_jacobian[0][1]*eta_jacobian_grad[2][0][0]-
                    eta_jacobian_grad[0][0][0]*eta_jacobian[2][1]-eta_jacobian[0][0]*eta_jacobian_grad[2][1][0];
  d_eta_Jxn_du[2] = eta_jacobian_grad[0][0][0]*eta_jacobian[1][1]+eta_jacobian[0][0]*eta_jacobian_grad[1][1][0]-
                    eta_jacobian_grad[0][1][0]*eta_jacobian[1][0]-eta_jacobian[0][1]*eta_jacobian_grad[1][0][0];
  // we also want the derivative of the previous vector with respect to the second parametric plane coordinate v
  Tensor<1,3> d_eta_Jxn_dv;
  d_eta_Jxn_dv[0] = eta_jacobian_grad[1][0][1]*eta_jacobian[2][1]+eta_jacobian[1][0]*eta_jacobian_grad[2][1][1]-
                    eta_jacobian_grad[1][1][1]*eta_jacobian[2][0]-eta_jacobian[1][1]*eta_jacobian_grad[2][0][1];
  d_eta_Jxn_dv[1] = eta_jacobian_grad[0][1][1]*eta_jacobian[2][0]+eta_jacobian[0][1]*eta_jacobian_grad[2][0][1]-
                    eta_jacobian_grad[0][0][1]*eta_jacobian[2][1]-eta_jacobian[0][0]*eta_jacobian_grad[2][1][1];
  d_eta_Jxn_dv[2] = eta_jacobian_grad[0][0][1]*eta_jacobian[1][1]+eta_jacobian[0][0]*eta_jacobian_grad[1][1][1]-
                    eta_jacobian_grad[0][1][1]*eta_jacobian[1][0]-eta_jacobian[0][1]*eta_jacobian_grad[1][0][1];
  //std::cout<<"dim singularity location: "<<eta<<std::endl;
  //std::cout<<"Spacedim singularity location: "<<eta_q_points_spacedim[0]<<std::endl;
  //std::cout<<"d_eta_Jxn_du "<<d_eta_Jxn_du<<std::endl;
  //std::cout<<"d_eta_Jxn_dv "<<d_eta_Jxn_dv<<std::endl;
//      std::cout<<"H1r: "<<eta_jacobian_grad[0][0][0]<<" "<<eta_jacobian_grad[0][0][1]<<std::endl;
//      std::cout<<"     "<<eta_jacobian_grad[0][1][0]<<" "<<eta_jacobian_grad[0][1][1]<<std::endl;
//      std::cout<<"H2r: "<<eta_jacobian_grad[1][0][0]<<" "<<eta_jacobian_grad[1][0][1]<<std::endl;
//      std::cout<<"     "<<eta_jacobian_grad[1][1][0]<<" "<<eta_jacobian_grad[1][1][1]<<std::endl;
//      std::cout<<"H3r: "<<eta_jacobian_grad[2][0][0]<<" "<<eta_jacobian_grad[2][0][1]<<std::endl;
//      std::cout<<"     "<<eta_jacobian_grad[2][1][0]<<" "<<eta_jacobian_grad[2][1][1]<<std::endl;

  
  // we now loop over the cell faces to carry out the integrals
  for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
      {
      // (to the best of my knowledge)
      // deal.ii does not offer an easy way to get the values of the
      // face vertices in the cell parametric plane
      // thus, we generate them here with a switch, given
      // the deal.ii numbering of faces and vertices
      // in generating them, we already center them in the
      // singularity point, subtracting P
      Point<2> csi_0;
      Point<2> csi_1;
      double coeff;
      switch (f)
          {
          case 0:
               csi_0(0) =  0-eta(0); csi_1(0) =  0-eta(0);
               csi_0(1) =  1-eta(1); csi_1(1) =  0-eta(1);
               coeff = csi_0(0);
          break;
          case 1:
               csi_0(0) =  1-eta(0); csi_1(0) =  1-eta(0);
               csi_0(1) =  0-eta(1); csi_1(1) =  1-eta(1);
               coeff = csi_0(0);
          break;
          case 2:
               csi_0(0) =  0-eta(0); csi_1(0) =  1-eta(0);
               csi_0(1) =  0-eta(1); csi_1(1) =  0-eta(1);
               coeff = csi_1(1);
          break;
          case 3:
               csi_0(0) =  1-eta(0); csi_1(0) =  0-eta(0);
               csi_0(1) =  1-eta(1); csi_1(1) =  1-eta(1);
               coeff = csi_1(1);
          break;
          default:
               csi_0(0) =  0-eta(0); csi_1(0) =  0-eta(0);
               csi_0(1) =  0-eta(1); csi_1(1) =  0-eta(1);
          break;
          }
      // once the coordinates in the parametric plane of the
      // face vertices are known, we can compute the polar coordinate
      // theta angle, taking care of having always growing angles
      // on each face f
      double theta_0 = atan2((csi_0)[1],(csi_0)[0]);
      if (theta_0 < 0)
         theta_0 += 2*dealii::numbers::PI;
      double theta_1 = atan2((csi_1)[1],(csi_1)[0]);
      if (theta_1 < 0)
         theta_1 += 2*dealii::numbers::PI;
      //std::cout<<"PRE:   "<<f<<" ->    theta_0: "<<theta_0*180/dealii::numbers::PI<<"   theta_1: "<<theta_1*180/dealii::numbers::PI<<std::endl;   
      if (theta_1<theta_0)
         theta_0 -= 2*dealii::numbers::PI;
      //std::cout<<f<<" ->    theta_0: "<<theta_0*180/dealii::numbers::PI<<"   theta_1: "<<theta_1*180/dealii::numbers::PI<<std::endl;
      

      // we now create a 1D triangulation for the integration
      // in the theta direction of the polar coordinates
      Triangulation<1> theta_triangulation;

      std::vector<Point<1> > theta_vertices;
      std::vector<CellData<1> > theta_cells;
      SubCellData theta_subcelldata;
      
      // we only have one cell that goes from the minimum to maximum
      // theta on the face f 
      theta_vertices.resize(2);
      theta_vertices[0](0)=theta_0;
      theta_vertices[1](0)=theta_1;

      theta_cells.resize(1);

      theta_cells[0].vertices[0]=0;
      theta_cells[0].vertices[1]=1;
      
      GridTools::delete_unused_vertices (theta_vertices, theta_cells, theta_subcelldata);
      GridTools::consistently_order_cells(theta_cells);
      // here the triangulation is initialized and ready
      theta_triangulation.create_triangulation(theta_vertices, theta_cells, theta_subcelldata );

      // we of course also create a dof handler finite element and mapping
      DoFHandler<1> theta_dof_handler(theta_triangulation);
      const FE_Q<1> theta_finite_element(fe.degree);
      theta_dof_handler.distribute_dofs(theta_finite_element);
      MappingQ<1> theta_mapping(theta_finite_element.degree);

      // we create a 1D quadature with the user selected number of quadrature nodes
      QGauss<1> theta_quadrature(theta_quadrature_order);

      // we are now ready to create a FEValues objects that will allow us
      // to obtain all we need for the integration over theta on the face f
      FEValues<1> theta_fe_v(theta_mapping,
                            theta_finite_element,
                            theta_quadrature,
                            update_values |
                            update_quadrature_points | update_JxW_values);

      // we initialize theta_fe_v with the only face available
      auto theta_cell = theta_dof_handler.begin_active(); //this is the only cell
      theta_fe_v.reinit(theta_cell);


      // because we also will need to carry out the surface integral in polar coordinates
      // for I_0, we will also need a 1D triangulation for rho
      Triangulation<1> rho_triangulation;

      std::vector<Point<1> >  rho_vertices;
      std::vector<CellData<1> >  rho_cells;
      SubCellData  rho_subcelldata;
      
      // in this case we will just consider a cell that goes from 0 to 1 and
      // take care of the transformation and jacobian by ourself
      rho_vertices.resize(2);
      rho_vertices[0](0)=0;
      rho_vertices[1](0)=1;

      rho_cells.resize(1);

      rho_cells[0].vertices[0]=0;
      rho_cells[0].vertices[1]=1;
      
      GridTools::delete_unused_vertices (rho_vertices, rho_cells, rho_subcelldata);
      GridTools::consistently_order_cells(rho_cells);
      // here the triangulation is initialized and ready
      rho_triangulation.create_triangulation(rho_vertices, rho_cells, rho_subcelldata );
      // we of course also create a dof handler finite element and mapping
      DoFHandler<1> rho_dof_handler(rho_triangulation);
      const FE_Q<1> rho_finite_element(fe.degree);
      rho_dof_handler.distribute_dofs(rho_finite_element);
      MappingQ<1> rho_mapping(rho_finite_element.degree);
      // we create a 1D quadature with the user selected number of quadrature nodes
      QGauss<1> rho_quadrature(rho_quadrature_order);
      // we are now ready to create a FEValues object that will allow us
      // to obtain all we need for the integration over rho on the face f
      FEValues<1> rho_fe_v(rho_mapping,
                           rho_finite_element,
                           rho_quadrature,
                           update_values |
                           update_quadrature_points | update_JxW_values);
      // we initialize theta_fe_v with the only face available
      auto rho_cell = rho_dof_handler.begin_active(); //this is the only cell
      rho_fe_v.reinit(rho_cell);
      // we get the quadrature points for rho (remember that they go from 0 to 1)
      const std::vector<Point<1>> &rho_q_points    = rho_fe_v.get_quadrature_points();

      // we get the quadrature points for theta (remember that they go from
      // theta_min to theta_max in the cell)
      const std::vector<Point<1>> &theta_q_points    = theta_fe_v.get_quadrature_points();
      // we want to create a vector with the \hat(rho) values corresponding to the face boundary
      // at each theta value
      std::vector<double> q_rho(theta_q_points.size());
      // using deal.ii numbering, these lines do the trick
      for (unsigned int q=0; q<theta_q_points.size(); ++q)
          {
          double a = (f<2)?0.0:1.0;
          double b = (f<2)?1.0:0.0;
          q_rho[q] = coeff/(a*sin(theta_q_points[q](0))+b*(cos(theta_q_points[q](0))));
          }
      // we now create a set of quadrature points and weights to be passed to a
      // FEValues<2,3> to take care of the mapping to the 3D space
      const std::vector<double> &uv_q_weights = theta_fe_v.get_JxW_values();
      std::vector<Point<2> > uv_q_points(theta_q_points.size());
      for (unsigned int q=0; q<theta_q_points.size(); ++q)
          {
          // the quadrature points are the points of the face f in the parametric space
          uv_q_points[q] = Point<2>(q_rho[q]*cos(theta_q_points[q](0))+eta(0),
                                    q_rho[q]*sin(theta_q_points[q](0))+eta(1));
          }
      Quadrature<2> uv_quadrature(uv_q_points, uv_q_weights);
      //here we create the FEValues and initialize it with the <2,3> dof handler cell
      FEValues<2,3> face_fe_values(fe,
                                   uv_quadrature,
                                   update_values | update_gradients | update_quadrature_points |
                                   update_JxW_values | update_normal_vectors);
      face_fe_values.reinit(cell);

      // we now create a set of vectors that will contain useful theta dependent coefficients
      // for each value of theta on the face f
      std::vector<Tensor<1,3> > A(theta_q_points.size());
      std::vector<Tensor<1,3> > B(theta_q_points.size());
      // values of A.norm()
      std::vector<double> q_A(theta_q_points.size());
      // values of A*B
      std::vector<double> q_C(theta_q_points.size());
      // values of jacobian 0-th order Taylor series coefficient
      std::vector<Tensor<1,3> > Jk0(theta_q_points.size());
      // values of jacobian 1-st order Taylor series coefficient
      std::vector<Tensor<1,3> > Jk1(theta_q_points.size());
      // values of shape function 0-th order Taylor series coefficient
      std::vector<double> serv(fe.dofs_per_cell,0.0);
      std::vector<std::vector<double> > N0(theta_q_points.size(),serv);
      // values of shape function 1-st order Taylor series coefficient
      std::vector<std::vector<double> > N1(theta_q_points.size(),serv);
      for (unsigned int q=0; q<theta_q_points.size(); ++q)
          {
          for (unsigned int i=0; i<3; ++i)
              {
              A[q][i] = eta_jacobian[i][0]*cos(theta_q_points[q](0)) +
                     eta_jacobian[i][1]*sin(theta_q_points[q](0));
              B[q][i] = eta_jacobian_grad[i][0][0]*pow(cos(theta_q_points[q](0)),2)/2 +
                        eta_jacobian_grad[i][0][1]*cos(theta_q_points[q](0))*sin(theta_q_points[q](0)) +
                        eta_jacobian_grad[i][1][1]*pow(sin(theta_q_points[q](0)),2)/2;
              }
          q_A[q] = A[q].norm();
          q_C[q] = A[q]*B[q];
          Jk0[q] = eta_Jxn;
          Jk1[q] = d_eta_Jxn_du*cos(theta_q_points[q](0))+d_eta_Jxn_dv*sin(theta_q_points[q](0));
          for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
              {
              N0[q][i] = eta_shape_values[i];
              N1[q][i] = eta_shape_grads[i][0]*cos(theta_q_points[q](0)) +
                         eta_shape_grads[i][1]*sin(theta_q_points[q](0));
              }
          }

      // so let's start computing the integrals on the face f
      // for I_1 and I_2
      for (unsigned int q=0; q<theta_q_points.size(); ++q)
          {
          // the taylor expansion coefficients F_1 and F_2 are finally computed
          double beta = 1/q_A[q];
          double gamma = -q_C[q]/pow(q_A[q],4);
          double F_2 = Jk0[q][2]/pow(q_A[q],3);
          double F_1 = -3*q_C[q]*Jk0[q][2]/pow(q_A[q],5)
                       -3*A[q][2]/pow(q_A[q],5)*(Jk0[q]*B[q]+Jk1[q]*A[q])
                       +Jk1[q][2]/pow(q_A[q],3);
          //std::cout<<"f: "<<f<<"  q: "<<q<<" theta: "<<theta_q_points[q](0)*180/dealii::numbers::PI<<"  F_1: "<<F_1<<"  F_2: "<<F_2<<"  rho: "<<q_rho[q]<<"  beta: "<<beta<<"  gamma: "<<gamma<<"  w: "<<uv_q_weights[q]<<"  Jk0: "<<Jk0[q]<<"  Jk1: "<<Jk1[q]<<"  A: "<<A[q]<<"  B: "<<B[q]<<std::endl;
          // and the integral contribution of the present quadrature node is here added
          I_2 += -F_2*(gamma/pow(beta,2)+1/q_rho[q])*uv_q_weights[q];
          I_1 += F_1*log(fabs(q_rho[q]/beta))*uv_q_weights[q];

          //I_1 += F_1*log(fabs(q_rho[q]))*uv_q_weights[q];
          
          
          // we still miss the surface integral I_0
          // for this, we will need to also consider the
          // quadrature nodes of the rho_fe_values
          // remembering they go from 0 to 1, we must rescale them from
          // 0 to q_rho[q]. Here is the obvious scaling factor
          double rho_jac_fact = q_rho[q];

          // we must again create a set of quadrature points and weights to be passed to a
          // FEValues<2,3> to take care of the mapping to the 3D space
          
          // we only want the mapping jacobian from such a FEValues, so we set all the
          // weights to 1
          std::vector<double> inner_uv_q_weights(rho_q_points.size(),1.0);
          
          std::vector<Point<2> > inner_uv_q_points(rho_q_points.size());
          for (unsigned int p=0; p<rho_q_points.size(); ++p)
              {
              // based on rho and theta, the quadrature points are generated in polar coordinates
              // centered at the singularity eta_q_points[0]
              double rho = rho_jac_fact*rho_q_points[p](0);
              inner_uv_q_points[p] = Point<2>(rho*cos(theta_q_points[q](0))+eta_q_points[0](0),
                                              rho*sin(theta_q_points[q](0))+eta_q_points[0](1));
              }
          
          // we can then create the quadrature and the FEValues
          Quadrature<2> inner_uv_quadrature(inner_uv_q_points, inner_uv_q_weights);
          FEValues<2,3> inner_face_fe_values(fe,
                                             inner_uv_quadrature,
                                             update_values | update_gradients | update_quadrature_points |
                                             update_JxW_values | update_normal_vectors);
          inner_face_fe_values.reinit(cell);
          // let's obtain the quadrature points in the 3D space from the inner_face_fe_values
          const std::vector<Point<3>> &inner_uv_q_points_spacedim    = inner_face_fe_values.get_quadrature_points();
          const std::vector<Tensor<1, 3 >> &inner_uv_q_normals = inner_face_fe_values.get_normal_vectors();
          for (unsigned int p=0; p<rho_q_points.size(); ++p)
              {
              // this is the rho obtained scaling the quadrature point in the [0,1] interval
              double rho = rho_jac_fact*rho_q_points[p](0);
              // this is the distance in the 3D space between quadrature point and singularity
              Tensor<1,3> R = inner_uv_q_points_spacedim[p]-eta_q_points_spacedim[0];
              double r = R.norm();
              // finally, this is the integral argument
              // the first term also needs the jacobian of the 3D mapping, which was instead
              // already included in the F_2 and F_1 factors of the latter terms, and is therefore
              // not needed there
              double fun = -1/pow(r,3)*(3*(R[2]/r)*(R*inner_uv_q_normals[p]/r)-inner_uv_q_normals[p][2]);
              double arg = fun*rho*inner_face_fe_values.JxW(p)-F_2/pow(rho,2)-F_1/rho;

//                  std::cout<<"rho: "<<rho<<"  r: "<<r<<"  arg1: "<<1/pow(r,3)*rho*inner_face_fe_values.JxW(p)<<
//                  "   arg2: "<<-F_2/pow(rho,2)-F_1/rho<<"  diff: "<<1/pow(r,3)*rho*inner_face_fe_values.JxW(p)-F_2/pow(rho,2)-F_1/rho<<std::endl;

              //std::cout<<rho<<"  "<<r<<"  "<<1/pow(r,3)*rho*inner_face_fe_values.JxW(p)<<
              //"  "<<-F_2/pow(rho,2)-F_1/rho<<"  "<<1/pow(r,3)*rho*inner_face_fe_values.JxW(p)-F_2/pow(rho,2)-F_1/rho<<std::endl;

//                  std::cout<<rho<<"  "<<r<<"  "<<Jk0[q][2]+rho*Jk1[q][2]<<"  "<<inner_uv_q_normals[p][2]*inner_face_fe_values.JxW(p)<<"   "<<inner_uv_q_normals[p][2]*inner_face_fe_values.JxW(p)-Jk0[q][2]-rho*Jk1[q][2]<<std::endl;

//                  std::cout<<rho<<"  "<<R<<"  "<<rho*A[q]+rho*rho*B[q]<<"  "<<(R-rho*A[q]-rho*rho*B[q]).norm()<<std::endl;
              
              // we add the conribution of the present node to the surface integral
              // the inner weights are given by rho_fe_v.JxW(p) but the jacobian is given by the scaling
              // factor rho_jac_fact
              I_0 += arg*rho_jac_fact*rho_fe_v.JxW(p)*theta_fe_v.JxW(q);
              }
          }
      /*
      std::cout<<"temp I_0 = "<<I_0<<std::endl;
      std::cout<<"temp I_1 = "<<I_1<<std::endl;
      std::cout<<"temp I_2 = "<<I_2<<std::endl;
      */

      }



 return (I_0 + I_1 + I_2)/numbers::PI/4.0;
};


template class SingularKernelIntegral<2>;
template class SingularKernelIntegral<3>;
