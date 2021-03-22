// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------






// #  include <BRepAdaptor_CompCurve.hxx>
// #  include <BRepAdaptor_Curve.hxx>
// #  include <BRepAdaptor_HCompCurve.hxx>
// #  include <BRepAdaptor_HCurve.hxx>
// #  include <BRepTools.hxx>
// #  include <BRep_Tool.hxx>
// #  include <GCPnts_AbscissaPoint.hxx>
// #  include <ShapeAnalysis_Curve.hxx>
// #  include <ShapeAnalysis_Surface.hxx>
// #  include <Standard_Version.hxx>
// #  include <TopoDS.hxx>


#include <deal.II/opencascade/utilities.h>

#  include <TopExp_Explorer.hxx>
#  include <TopoDS.hxx>
#  include <TopoDS_Edge.hxx>
#  include <TopoDS_Face.hxx>
#  include <TopoDS_Shape.hxx>
#  include <ShapeAnalysis_Surface.hxx>
#  include <BRepAdaptor_Curve.hxx>
#  include <BRepAdaptor_HCompCurve.hxx>
#  include <BRepAdaptor_HCurve.hxx>
#  include <BRepAdaptor_Surface.hxx>
#  include <BRepAlgo_Section.hxx>
#  include <BRepBndLib.hxx>
#  include <BRepBuilderAPI_MakeVertex.hxx>
#  include <BRepBuilderAPI_MakeEdge.hxx>
#  include <BRepBuilderAPI_Sewing.hxx>
#  include <BRepBuilderAPI_Transform.hxx>
#  include <BRepMesh_IncrementalMesh.hxx>
#  include <BRepTools.hxx>
#  include <BRep_Builder.hxx>
#  include <GeomAPI_ProjectPointOnCurve.hxx>
#  include <BRepExtrema_DistShapeShape.hxx>
#  include <Extrema_ExtAlgo.hxx>
#  include <Extrema_ExtFlag.hxx>
#  include <BRepBndLib.hxx>
#  include <BRepMesh_IncrementalMesh.hxx>
#  include <IntCurvesFace_ShapeIntersector.hxx>
#  include <Geom_Line.hxx>



#include "../include/my_manifold_lib.h"
#include "../include/my_utilities.h"

DEAL_II_NAMESPACE_OPEN


using namespace OpenCASCADE;








  /*============================== MyNormalToMeshProjectionManifold
   * ==============================*/
  template <int dim, int spacedim>
  MyNormalToMeshProjectionManifold<dim, spacedim>::MyNormalToMeshProjectionManifold(
    const TopoDS_Shape &sh,
    const double        tolerance)
    : sh(sh)
    , tolerance(tolerance)
  {
    Assert(spacedim == 3, ExcNotImplemented());
    Assert(
      std::get<0>(count_elements(sh)) > 0,
      ExcMessage(
        "MyNormalToMeshProjectionManifold needs a shape containing faces to operate."));
    Standard_Real aDeflection = 0.0001, deflection;
    Standard_Real aXmin, aYmin, aZmin, aXmax, aYmax, aZmax;
    Bnd_Box box;
    BRepBndLib::Add(sh, box);
    box.Get(aXmin, aYmin, aZmin, aXmax, aYmax, aZmax);
    deflection= sqrt( pow(aXmax-aXmin,2)+ pow(aYmax-aYmin,2)+pow(aZmax-aZmin,2))*aDeflection;
    BRepMesh_IncrementalMesh Inc(sh, deflection);
  }

  template <int dim, int spacedim>
  std::unique_ptr<Manifold<dim, spacedim>>
  MyNormalToMeshProjectionManifold<dim, spacedim>::clone() const
  {
    return std::unique_ptr<Manifold<dim, spacedim>>(
      new MyNormalToMeshProjectionManifold<dim, spacedim>(sh, tolerance));
  }


#include "Teuchos_TimeMonitor.hpp"

Teuchos::RCP<Teuchos::Time> project_to_manifold_all          = Teuchos::TimeMonitor::getNewTimer("project_to_manifold all");
Teuchos::RCP<Teuchos::Time> project_to_manifold_line_inters  = Teuchos::TimeMonitor::getNewTimer("project_to_manifold line inters");

Teuchos::RCP<Teuchos::Time> project_to_manifold1 = Teuchos::TimeMonitor::getNewTimer("project_to_manifold1");
Teuchos::RCP<Teuchos::Time> project_to_manifold2 = Teuchos::TimeMonitor::getNewTimer("project_to_manifold2");
Teuchos::RCP<Teuchos::Time> project_to_manifold3 = Teuchos::TimeMonitor::getNewTimer("project_to_manifold3");
Teuchos::RCP<Teuchos::Time> project_to_manifold31 = Teuchos::TimeMonitor::getNewTimer("project_to_manifold31");
Teuchos::RCP<Teuchos::Time> project_to_manifold311 = Teuchos::TimeMonitor::getNewTimer("project_to_manifold311");
Teuchos::RCP<Teuchos::Time> project_to_manifold32 = Teuchos::TimeMonitor::getNewTimer("project_to_manifold32");
Teuchos::RCP<Teuchos::Time> project_to_manifold4 = Teuchos::TimeMonitor::getNewTimer("project_to_manifold4");
Teuchos::RCP<Teuchos::Time> project_to_manifold5 = Teuchos::TimeMonitor::getNewTimer("project_to_manifold5");

template<int spacedim>
std::string point_to_string(const Point<spacedim> &point)
{
  std::string result = std::to_string(point[0]);
  for(unsigned int i=1; i<spacedim; ++i)
    result =result +std::to_string(point[i]);
  return result;
}

     template <int spacedim>
     Point<spacedim>
     internal_project_to_manifold(std::unordered_map <std::size_t, Tensor<1, spacedim> > &,const TopoDS_Shape &,
                                  const double,
                                  const ArrayView<const Point<spacedim>> &,
                                  const Point<spacedim> &)
     {
       Assert(false, ExcNotImplemented());
       return {};
     }

  template <>
  Point<3>
  internal_project_to_manifold(std::unordered_map <std::size_t, Tensor<1, 3> > &projections_cache, const TopoDS_Shape &sh,
                                  const double tolerance,
    const ArrayView<const Point<3>> &surrounding_points,
    const Point<3> &                 candidate) 
  {
//std::cout << "MyNormalToMeshProjectionManifold<dim, spacedim>::project_to_manifold\n";
  Teuchos::TimeMonitor localTimer1(*project_to_manifold_all);

    TopoDS_Shape out_shape;
    Tensor<1, 3> average_normal;
    
    double cell_size = 0.0;
    for (unsigned int i = 0; i < surrounding_points.size(); ++i)
      {
      cell_size += (candidate-surrounding_points[i]).norm()/surrounding_points.size();
      }
    cell_size *= 2.0;
    
#  ifdef DEBUG
{
  Teuchos::TimeMonitor localTimer11(*project_to_manifold1);
  for (unsigned int i = 0; i < surrounding_points.size(); ++i)
      {
        Assert(closest_point(sh, surrounding_points[i], tolerance)
                   .distance(surrounding_points[i]) <
                 std::max(tolerance * surrounding_points[i].norm(), tolerance),
               ExcPointNotOnManifold<3>(surrounding_points[i]));
      }
}
#  endif
{
  Teuchos::TimeMonitor localTimer12(*project_to_manifold2);
    switch (surrounding_points.size())
      {
        case 2:
          {

  Teuchos::TimeMonitor localTimer13(*project_to_manifold3);
{
  Teuchos::TimeMonitor localTimer131(*project_to_manifold31);
            for (unsigned int i = 0; i < surrounding_points.size(); ++i)
              {
                Tensor<1, 3> surface_normal;
                auto key = std::hash<std::string>()(point_to_string(surrounding_points[i]));
                auto it_1 = projections_cache.find(key);
                if (it_1 != projections_cache.end())
                {
                    surface_normal = it_1->second;
                }
                else
                {
  Teuchos::TimeMonitor localTimer1311(*project_to_manifold311);
                    std::tuple<Point<3>, Tensor<1, 3>, double, double> tmp = my_closest_point_and_differential_forms(sh,surrounding_points[i],tolerance);
                    surface_normal = std::get<1>(tmp);
                    projections_cache.insert({key,surface_normal});
                }
                average_normal += surface_normal;
              }
}
{
  Teuchos::TimeMonitor localTimer132(*project_to_manifold32);
            average_normal /= 2.0;

            Assert(
              average_normal.norm() > 1e-4,
              ExcMessage(
                "Failed to refine cell: the average of the surface normals at the surrounding edge turns out to be a null vector, making the projection direction undetermined."));

            Tensor<1, 3> T = surrounding_points[0] - surrounding_points[1];
            T /= T.norm();
            average_normal = average_normal - (average_normal * T) * T;
            average_normal /= average_normal.norm();
}
            break;
          }
        case 4:
          {
  Teuchos::TimeMonitor localTimer14(*project_to_manifold4);
            Tensor<1, 3> u = surrounding_points[1] - surrounding_points[0];
            Tensor<1, 3> v = surrounding_points[2] - surrounding_points[0];
            const double n1_coords[3] = {u[1] * v[2] - u[2] * v[1],
                                         u[2] * v[0] - u[0] * v[2],
                                         u[0] * v[1] - u[1] * v[0]};
            Tensor<1, 3> n1(n1_coords);
            n1 = n1 / n1.norm();
            u  = surrounding_points[2] - surrounding_points[3];
            v  = surrounding_points[1] - surrounding_points[3];
            const double n2_coords[3] = {u[1] * v[2] - u[2] * v[1],
                                         u[2] * v[0] - u[0] * v[2],
                                         u[0] * v[1] - u[1] * v[0]};
            Tensor<1, 3> n2(n2_coords);
            n2 = n2 / n2.norm();

            average_normal = (n1 + n2) / 2.0;

            Assert(
              average_normal.norm() > tolerance,
              ExcMessage(
                "Failed to refine cell: the normal estimated via the surrounding points turns out to be a null vector, making the projection direction undetermined."));

            average_normal /= average_normal.norm();
            break;
          }
        case 8:
          {
  Teuchos::TimeMonitor localTimer15(*project_to_manifold5);
            Tensor<1, 3> u = surrounding_points[1] - surrounding_points[0];
            Tensor<1, 3> v = surrounding_points[2] - surrounding_points[0];
            const double n1_coords[3] = {u[1] * v[2] - u[2] * v[1],
                                         u[2] * v[0] - u[0] * v[2],
                                         u[0] * v[1] - u[1] * v[0]};
            Tensor<1, 3> n1(n1_coords);
            n1 = n1 / n1.norm();
            u  = surrounding_points[2] - surrounding_points[3];
            v  = surrounding_points[1] - surrounding_points[3];
            const double n2_coords[3] = {u[1] * v[2] - u[2] * v[1],
                                         u[2] * v[0] - u[0] * v[2],
                                         u[0] * v[1] - u[1] * v[0]};
            Tensor<1, 3> n2(n2_coords);
            n2 = n2 / n2.norm();
            u  = surrounding_points[4] - surrounding_points[7];
            v  = surrounding_points[6] - surrounding_points[7];
            const double n3_coords[3] = {u[1] * v[2] - u[2] * v[1],
                                         u[2] * v[0] - u[0] * v[2],
                                         u[0] * v[1] - u[1] * v[0]};
            Tensor<1, 3> n3(n3_coords);
            n3 = n3 / n3.norm();
            u  = surrounding_points[6] - surrounding_points[7];
            v  = surrounding_points[5] - surrounding_points[7];
            const double n4_coords[3] = {u[1] * v[2] - u[2] * v[1],
                                         u[2] * v[0] - u[0] * v[2],
                                         u[0] * v[1] - u[1] * v[0]};
            Tensor<1, 3> n4(n4_coords);
            n4 = n4 / n4.norm();

            average_normal = (n1 + n2 + n3 + n4) / 4.0;

            Assert(
              average_normal.norm() > tolerance,
              ExcMessage(
                "Failed to refine cell: the normal estimated via the surrounding points turns out to be a null vector, making the projection direction undetermined."));

            average_normal /= average_normal.norm();
            break;
          }
        default:
          {
            AssertThrow(false, ExcNotImplemented());
            break;
          }
      }
}
  Teuchos::TimeMonitor localTimer2(*project_to_manifold_line_inters);
    return my_line_intersection(sh, candidate, average_normal, tolerance);
  }


   template <int dim, int spacedim>
   Point<spacedim>
   MyNormalToMeshProjectionManifold<dim, spacedim>::project_to_manifold(
     const ArrayView<const Point<spacedim>> &surrounding_points,
     const Point<spacedim> &                 candidate) const
   {
     return internal_project_to_manifold(projections_cache,sh,
                                         tolerance,
                                         surrounding_points,
                                         candidate);
   }
// Explicit instantiations
// template class MyNormalToMeshProjectionManifold<3, 3>;
template class MyNormalToMeshProjectionManifold<1, 2>;
template class MyNormalToMeshProjectionManifold<2, 3>;
//template class MyNormalToMeshProjectionManifold<2, 2>;
//template class MyNormalToMeshProjectionManifold<1, 2>;


DEAL_II_NAMESPACE_CLOSE

