#include <cached_nurbs_patch_manifold.h>
#include <GCPnts_AbscissaPoint.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_CompCurve.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <ShapeAnalysis_Curve.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <TopoDS.hxx>
#include <BRepExtrema_DistShapeShape.hxx>

#include <Standard_Version.hxx>
#include <boost/container/small_vector.hpp>
DEAL_II_NAMESPACE_OPEN

template<int spacedim>
std::string point_to_string(const Point<spacedim> &point)
{
  std::string result = std::to_string(point[0])+std::to_string(point[1])+std::to_string(point[2]);
  return result;
}

using namespace OpenCASCADE;
//  Point<3>
//  CachedNURBSPatchManifold::get_new_point(const ArrayView<const Point<3>> &surrounding_points,
//                  const ArrayView<const double>          &weights) const
//  {
//   // std::cout<<"USING"<<std::endl;
// const std::size_t n_points = surrounding_points.size();
// boost::container::small_vector<Point<2>, 200> chart_points(n_points);
//    for (unsigned int i=0; i<n_points; ++i)
//    {
//     auto it= projections_cache.find(point_to_string(surrounding_points[i]));
//     if (it != projections_cache.end())
//     {
//       chart_points[i] = it->second;
//     }
//     else
//     {
//       chart_points[i] = pull_back(surrounding_points[i]);
//       projections_cache.insert({point_to_string(surrounding_points[i]), chart_points[i]});
//     }
//    }
//     const Point<2> p_chart = sub_manifold.get_new_point
//                    (make_array_view(chart_points.begin(), chart_points.end()),
//                                     weights);

//     return push_forward(p_chart);


//  }


// void
//  CachedNURBSPatchManifold::get_new_points(const ArrayView<const Point<3>> &surrounding_points,
//                  const Table<2,double>                  &weights,
//                   ArrayView<Point<3>>              new_points) const
//  {
//   // std::cout<<"USING"<<std::endl;
// const std::size_t n_points = surrounding_points.size();
// boost::container::small_vector<Point<2>, 200> chart_points(n_points);
//    for (unsigned int i=0; i<n_points; ++i)
//    {
//     auto it= projections_cache.find(point_to_string(surrounding_points[i]));
//     if (it != projections_cache.end())
//     {
//       chart_points[i] = it->second;
//     }
//     else
//     {
//       chart_points[i] = pull_back(surrounding_points[i]);
//       projections_cache.insert({point_to_string(surrounding_points[i]), chart_points[i]});
//     }
//    }

//     boost::container::small_vector<Point<2>, 200> new_points_on_chart(weights.size(0));
//     sub_manifold.get_new_points
//                    (make_array_view(chart_points.begin(), chart_points.end()),
//                                     weights,
//                                      make_array_view(new_points_on_chart.begin(),
//                                     new_points_on_chart.end()));

//    for (std::size_t row=0; row<weights.size(0); ++row)
//       new_points[row] = push_forward(new_points_on_chart[row]);

//  }



// Tensor<1,3>
//   CachedNURBSPatchManifold::get_tangent_vector (const Point<3> &x1,
//                       const Point<3> &x2) const
//   {
//    Point<2> p_x1, p_x2;

//    // auto it_1= projections_cache.find(point_to_string(x1));
//    //  auto it_2= projections_cache.find(point_to_string(x2));

//    //  if(it_1 != projections_cache.end())
//    //    p_x1 = it_1->second;
//    //  else
//      p_x1 = pull_back(x1);

//    //  if(it_2 != projections_cache.end())
//    //    p_x2 = it_2->second;
//    //  else
//      p_x2 = pull_back(x2);


//     const DerivativeForm<1,2,3> F_prime = push_forward_gradient((p_x1));
//      Assert (std::pow(std::abs(F_prime.determinant()), 1./2) >= 1e-12 * F_prime.norm(),
//             ExcMessage("The derivative of a chart function must not be singular."));




//     const Tensor<1,2>    delta   = sub_manifold.get_tangent_vector((p_x1),
//                                                         (p_x2));
//     Tensor<1,3> result;
//     for (unsigned int i=0; i<3; ++i)
//       result[i] += F_prime[i] * delta;

//     return result;
//   }
template<int dim, int spacedim>
Point<dim>
CachedNURBSPatchManifold<dim, spacedim>::pull_back(const Point<spacedim> &space_point) const
{
  auto it_1= projections_cache.find(std::hash<std::string>()(point_to_string(space_point)));
  Point<dim> result;
  if (it_1 != projections_cache.end())
    {
      result = it_1->second;
      cache_hits += 1;
    }
  else
    {
      Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);

      ShapeAnalysis_Surface projector(SurfToProj);
      gp_Pnt2d proj_params = projector.ValueOfUV(point(space_point), tolerance);

      result[0] = proj_params.X();
      result[1] = proj_params.Y();
      projections_cache.insert({std::hash<std::string>()(point_to_string(space_point)), result});
      new_pull_backs += 1;
    }

  return result;
}


template class CachedNURBSPatchManifold<2,3>;

DEAL_II_NAMESPACE_CLOSE
