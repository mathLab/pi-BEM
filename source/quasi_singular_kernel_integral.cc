#include "../include/quasi_singular_kernel_integral.h"


template <>
QuasiSingularKernelIntegral<3>::QuasiSingularKernelIntegral(
  const DoFHandler<2, 3>::active_cell_iterator &in_cell,
  const FiniteElement<2, 3>                    &in_fe,
  const Mapping<2, 3>                          &in_mapping,
  const Point<3>                               &in_external_point)
  : cell(in_cell)
  , fe(in_fe)
  , mapping(in_mapping)
{
  eta = find_closest_reference_cell_point(in_external_point);
};


template <>
QuasiSingularKernelIntegral<2>::QuasiSingularKernelIntegral(
  const DoFHandler<1, 2>::active_cell_iterator &in_cell,
  const FiniteElement<1, 2>                    &in_fe,
  const Mapping<1, 2>                          &in_mapping,
  const Point<2> & /*in_external_point*/)
  : cell(in_cell)
  , fe(in_fe)
  , mapping(in_mapping){};



template <int dim>
Point<dim - 1>
QuasiSingularKernelIntegral<dim>::find_closest_reference_cell_point(
  const Point<dim> &external_point)
{
  Point<dim - 1> eta_sol;

  if (dim == 2)
    {
      ExcNotImplemented();
    }
  else
    {
      Vector<double> solution(2);
      solution[0] = 0.5;
      solution[1] = 0.5;
      Vector<double> residual(2);
      residual                  = 1.0;
      unsigned int iter_counter = 0;
      double       min_dist;
      while (residual.l2_norm() > 1e-7 && iter_counter < 30)
        {
          // here's the quadrature rule obtained with the solution point only
          // point
          Point<2> solution_point;
          solution_point[0] = solution[0];
          solution_point[1] = solution[1];

          std::vector<Point<2>> sol_q_points;
          sol_q_points.push_back(solution_point);
          Quadrature<2> sol_quadrature(sol_q_points);

          // the quadrature is needed to declare an FEValues object
          // to give us all we need in terms of point position
          // on the real cell, and its derivatives w/r to the
          // reference cell coordinates
          FEValues<2, 3> sol_fe_values(mapping,
                                       fe,
                                       sol_quadrature,
                                       //update_values | update_gradients |
                                       update_quadrature_points |
                                         update_normal_vectors |
                                         update_jacobians |
                                         update_jacobian_grads);


          sol_fe_values.reinit(cell);

          // the single quadrature point is the point in the three dimensional
          // domain corresponding to eta/P
          const std::vector<Point<3>> &sol_q_points_spacedim =
            sol_fe_values.get_quadrature_points();
          // we also get the jacobian and jacobian gradient at such a location
          auto sol_jacobian      = sol_fe_values.jacobian(0);
          auto sol_jacobian_grad = sol_fe_values.jacobian_grad(0);

          // this is the distance vector between the external point and the
          // current point on the cell
          Vector<double> dist_vect(3);
          dist_vect[0] = external_point[0] - sol_q_points_spacedim[0][0];
          dist_vect[1] = external_point[1] - sol_q_points_spacedim[0][1];
          dist_vect[2] = external_point[2] - sol_q_points_spacedim[0][2];
          min_dist     = dist_vect.l2_norm();
          // std::cout<<"Min Dist: "<<min_dist<<std::endl;


          // this is the jacobian of the mapping
          // at the current solution point
          FullMatrix<double> dq_du(3, 2);
          dq_du[0][0] = sol_jacobian[0][0];
          dq_du[1][0] = sol_jacobian[1][0];
          dq_du[2][0] = sol_jacobian[2][0];
          dq_du[0][1] = sol_jacobian[0][1];
          dq_du[1][1] = sol_jacobian[1][1];
          dq_du[2][1] = sol_jacobian[2][1];

          // the residual of the nonlinear function of
          // which we are searching the root
          // it is the gradient of the distance vector
          residual = 0.0;
          dq_du.Tvmult(residual, dist_vect);
          residual *= 2.0;

          // we need the jacobian of the residual just assembled
          // so we collect the second derivatives
          // of the mapping in two 3x2 matrices (its basically
          // a triple tensor)
          FullMatrix<double> ddq_du_1(3, 2);
          ddq_du_1[0][0] = sol_jacobian_grad[0][0][0];
          ddq_du_1[1][0] = sol_jacobian_grad[1][0][0];
          ddq_du_1[2][0] = sol_jacobian_grad[2][0][0];
          ddq_du_1[0][1] = sol_jacobian_grad[0][1][0];
          ddq_du_1[1][1] = sol_jacobian_grad[1][1][0];
          ddq_du_1[2][1] = sol_jacobian_grad[2][1][0];

          FullMatrix<double> ddq_du_2(3, 2);
          ddq_du_2[0][0] = sol_jacobian_grad[0][0][1];
          ddq_du_2[1][0] = sol_jacobian_grad[1][0][1];
          ddq_du_2[2][0] = sol_jacobian_grad[2][0][1];
          ddq_du_2[0][1] = sol_jacobian_grad[0][1][1];
          ddq_du_2[1][1] = sol_jacobian_grad[1][1][1];
          ddq_du_2[2][1] = sol_jacobian_grad[2][1][1];

          // we now have all the ingredients to assemble
          // the jacobian. the first component is given by
          // the point position gradient multiplied by its
          // transpose
          FullMatrix<double> jacobian(2, 2);
          dq_du.Tmmult(jacobian, dq_du);
          jacobian *= -2.0;

          // then we need to multiply the distance
          // vector by the triple tensor (here split into
          // two matrices)
          Vector<double> col_1(2);
          ddq_du_1.Tvmult(col_1, dist_vect);
          col_1 *= 2.0;

          Vector<double> col_2(2);
          ddq_du_2.Tvmult(col_2, dist_vect);
          col_2 *= 2.0;

          jacobian[0][0] += col_1[0];
          jacobian[0][1] += col_2[0];
          jacobian[1][0] += col_1[1];
          jacobian[1][1] += col_2[1];

          // let's finally carry out newton's step,
          // inverting the jabobian and
          // multiplying it by -residual
          Vector<double>     delta_solution(2);
          FullMatrix<double> jacobian_inverse(2, 2);
          jacobian_inverse.invert(jacobian);
          residual *= -1.0;
          jacobian_inverse.vmult(delta_solution, residual);

          solution += delta_solution;
          // std::cout<<"Iteration: "<<iter_counter<<"  Residual:
          // "<<residual.l2_norm()<<std::endl;
          iter_counter++;
        }

      if (iter_counter == 30)
        {
          // std::cout<<"30 Iteratons"<<std::endl;
          eta_sol[0] = std::numeric_limits<double>::infinity();
          eta_sol[1] = std::numeric_limits<double>::infinity();
        }
      else
        {
          // std::cout<<"Solved"<<std::endl;
          // std::cout<<"Min Dist: "<<min_dist<<std::endl;
          eta_sol[0] = solution[0];
          eta_sol[1] = solution[1];
          // std::cout<<"Cell minimum distance point: "<<eta_sol<<std::endl;
        }

      if (eta_sol[0] < 0 || eta_sol[0] > 1 || eta_sol[1] < 0 || eta_sol[1] > 1)
        {
          double minimum_distance_among_faces =
            std::numeric_limits<double>::infinity();
          unsigned int face_with_minimum_distance =
            GeometryInfo<2>::faces_per_cell;
          // in this case we must look for the closest point, but on the
          // boundary we'll use bisection in the interval [0,1]
          for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
            {
              unsigned int iter_counter_line = 0;
              double       min_dist_line;
              double       residual_line = 1.0;
              // initial guess
              double solution_line = 0.5;
              while (fabs(residual_line) > 1e-7 && iter_counter_line < 30)
                {
                  // here's the quadrature rule obtained with the solution point
                  Point<1> solution_point_line;
                  solution_point_line[0] = solution_line;

                  std::vector<Point<1>> sol_q_points_line;
                  sol_q_points_line.push_back(solution_point_line);
                  Quadrature<1> sol_quadrature_line(sol_q_points_line);

                  FEFaceValues<2, 3> sol_fe_values_line(
                    mapping,
                    fe,
                    sol_quadrature_line,
                    update_quadrature_points | update_jacobians |
                      update_jacobian_grads);
                  sol_fe_values_line.reinit(cell, f);

                  // the single quadrature point is the point in the three
                  // dimensional domain corresponding to eta/P
                  const std::vector<Point<3>> &sol_q_points_spacedim_line =
                    sol_fe_values_line.get_quadrature_points();
                  // we also get the jacobian and jacobian gradient at such a
                  // location
                  auto sol_jacobian_line = sol_fe_values_line.jacobian(0);
                  auto sol_jacobian_grad_line =
                    sol_fe_values_line.jacobian_grad(0);

                  // this is the distance vector between the external point and
                  // the current point on the cell
                  Tensor<1, 3> dist_vect_line;
                  dist_vect_line[0] =
                    external_point[0] - sol_q_points_spacedim_line[0][0];
                  dist_vect_line[1] =
                    external_point[1] - sol_q_points_spacedim_line[0][1];
                  dist_vect_line[2] =
                    external_point[2] - sol_q_points_spacedim_line[0][2];
                  min_dist_line = dist_vect_line.norm();
                  // std::cout<<"Min Dist Line "<<f<<":
                  // "<<min_dist_line<<std::endl;


                  // this is the jacobian of the mapping
                  // at the current solution point
                  Tensor<1, 3> dq_du_line;
                  if (f < 2)
                    {
                      dq_du_line[0] = sol_jacobian_line[0][1];
                      dq_du_line[1] = sol_jacobian_line[1][1];
                      dq_du_line[2] = sol_jacobian_line[2][1];
                    }
                  else
                    {
                      dq_du_line[0] = sol_jacobian_line[0][0];
                      dq_du_line[1] = sol_jacobian_line[1][0];
                      dq_du_line[2] = sol_jacobian_line[2][0];
                    }

                  // the residual of the nonlinear function of
                  // which we are searching the root
                  // it is the gradient of the distance vector
                  residual_line = -2.0 * dq_du_line * dist_vect_line;

                  // we need the jacobian of the residual just assembled
                  // so we collect the second derivatives
                  Tensor<1, 3> ddq_duu;
                  if (f < 2)
                    {
                      ddq_duu[0] = sol_jacobian_grad_line[0][1][1];
                      ddq_duu[1] = sol_jacobian_grad_line[1][1][1];
                      ddq_duu[2] = sol_jacobian_grad_line[2][1][1];
                    }
                  else
                    {
                      ddq_duu[0] = sol_jacobian_grad_line[0][0][0];
                      ddq_duu[1] = sol_jacobian_grad_line[1][0][0];
                      ddq_duu[2] = sol_jacobian_grad_line[2][0][0];
                    }
                  // the jecobian is then
                  double jacobian_line = 2.0 * dq_du_line * dq_du_line -
                                         2.0 * ddq_duu * dist_vect_line;

                  // we have all the ingredients now
                  // to implement the newton step
                  double delta_solution_line = -residual_line / jacobian_line;
                  solution_line += delta_solution_line;

                  // std::cout<<"Face "<<f<<"  Iteration:
                  // "<<iter_counter_line<<"  Residual: "<<residual_line<<"
                  // Distance: "<<min_dist_line<<"  Solution:
                  // "<<solution_line<<" 3D Point:
                  // "<<sol_q_points_spacedim_line[0]<<std::endl;
                  // std::cout<<"JAC: "<<dq_du_line<<std::endl;
                  // std::cout<<"Deriv form: "<<sol_jacobian_line<<std::endl;


                  iter_counter_line++;
                }
              // std::cout<<"Face "<<f<<std::endl;
              if (iter_counter_line == 30)
                {
                  // std::cout<<"30 Iteratons"<<std::endl;
                }
              else
                {
                  // std::cout<<"Solved"<<std::endl;
                  // std::cout<<"Min Dist: "<<min_dist_line<<"  Solution line:
                  // "<<solution_line<<std::endl;
                  if (solution_line > 0 && solution_line < 1)
                    if (min_dist_line < minimum_distance_among_faces)
                      {
                        minimum_distance_among_faces = min_dist_line;
                        face_with_minimum_distance   = f;
                        switch (f)
                          {
                            case 0:
                              eta_sol[0] = 0.0;
                              eta_sol[1] = solution_line;
                              break;
                            case 1:
                              eta_sol[0] = 1.0;
                              eta_sol[1] = solution_line;
                              break;
                            case 2:
                              eta_sol[0] = solution_line;
                              eta_sol[1] = 0.0;
                              break;
                            case 3:
                              eta_sol[0] = solution_line;
                              eta_sol[1] = 1.0;
                              break;
                            default:
                              AssertThrow(
                                true,
                                ExcMessage(
                                  "Cell has four faces and we are not in face 0,1,2,3"));
                              break;
                          }
                      }
                }
            }
          // finally, if no one has touched minimum_distance_among_faces,
          // it means the minimum distance is not on one of the edges,
          // and must be found among the vertices
          if (minimum_distance_among_faces ==
              std::numeric_limits<double>::infinity())
            {
              double minimum_distance_among_vertices =
                minimum_distance_among_faces;
              // we loop on the 4 vertices and get the closest one
              for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell;
                   ++v)
                {
                  double v_distance = external_point.distance(cell->vertex(v));
                  if (v_distance < minimum_distance_among_vertices)
                    {
                      minimum_distance_among_vertices = v_distance;
                      eta_sol = GeometryInfo<2>::unit_cell_vertex(v);
                    }
                }
              // std::cout<<"Minimum distance found among vertices:
              // "<<minimum_distance_among_vertices<<std::endl;
            }
        }
    }

  // std::cout<<"External point coordinates: "<<external_point<<std::endl;
  // std::cout<<"Ref cell min distance point coordinates: "<<eta_sol<<std::endl;
  Point<3> real_min_distance_point =
    mapping.transform_unit_to_real_cell(cell, eta_sol);
  this->min_distance = external_point.distance(real_min_distance_point);
  // std::cout<<"Real cell min distance point coordinates:
  // "<<real_min_distance_point<<std::endl; std::cout<<"Min distance:
  // "<<min_distance<<std::endl; std::cout<<"Cell size:
  // "<<cell->diameter()<<std::endl;

  return eta_sol;
}
