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

DEAL_II_NAMESPACE_OPEN


using namespace OpenCASCADE;



/* 
Fast Ray-Box Intersection
by Andrew Woo
from "Graphics Gems", Academic Press, 1990
*/


#define NUMDIM	3
#define RIGHT	0
#define LEFT	1
#define MIDDLE	2

bool HitBoundingBox(const Point<3> minB, const Point<3> maxB, const Point<3> origin, const Tensor<1,3> dir, Point<3> coord)
//double minB[NUMDIM], maxB[NUMDIM];		/*box */
//double origin[NUMDIM], dir[NUMDIM];		/*ray */
//double coord[NUMDIM];				/* hit point */
{
	bool inside = true;
	char quadrant[NUMDIM];
	register int i;
	int whichPlane;
	double maxT[NUMDIM];
	double candidatePlane[NUMDIM];

	/* Find candidate planes; this loop can be avoided if
   	rays cast all from the eye(assume perpsective view) */
	for (i=0; i<NUMDIM; i++)
		if(origin(i) < minB(i)) {
			quadrant[i] = LEFT;
			candidatePlane[i] = minB(i);
			inside = false;
		}else if (origin(i) > maxB(i)) {
			quadrant[i] = RIGHT;
			candidatePlane[i] = maxB(i);
			inside = false;
		}else	{
			quadrant[i] = MIDDLE;
		}

	/* Ray origin inside bounding box */
	if(inside)	{
		coord = origin;
		return (true);
	}


	/* Calculate T distances to candidate planes */
	for (i = 0; i < NUMDIM; i++)
		if (quadrant[i] != MIDDLE && dir[i] !=0.)
			maxT[i] = (candidatePlane[i]-origin(i)) / dir[i];
		else
			maxT[i] = -1.;

	/* Get largest of the maxT's for final choice of intersection */
	whichPlane = 0;
	for (i = 1; i < NUMDIM; i++)
		if (maxT[whichPlane] < maxT[i])
			whichPlane = i;

	/* Check final candidate actually inside box */
	if (maxT[whichPlane] < 0.) return (false);
	for (i = 0; i < NUMDIM; i++)
		if (whichPlane != i) {
			coord[i] = origin(i) + maxT[whichPlane] *dir[i];
			if (coord[i] < minB(i) || coord[i] > maxB(i))
				return (false);
		} else {
			coord[i] = candidatePlane[i];
		}
	return (true);				/* ray hits box */
}	





  template <int dim>
  std::tuple<Point<dim>, TopoDS_Shape, double, double>
  my_project_point_and_pull_back(const TopoDS_Shape &in_shape,
                              const Point<dim> &  origin,
                              const double        tolerance)
  {
    TopExp_Explorer exp;
    gp_Pnt          Pproj = point(origin);
//    std::tuple<double,double,double> shape_count = count_elements(in_shape);
//    
//    if (std::get<0>(shape_count) > 0)
//        {
//        TopoDS_Vertex Vproj = BRepBuilderAPI_MakeVertex(Pproj).Vertex();
//        //Extrema_ExtAlgo_Tree
//        BRepExtrema_DistShapeShape vertex_projector(in_shape,Vproj, 1e-4, Extrema_ExtFlag_MINMAX, Extrema_ExtAlgo_Tree);
//        vertex_projector.Perform();
//        //cout<<"***** "<<vertex_projector.Value()<<endl;
//        //cout<<"##### "<<vertex_projector.SupportTypeShape1(1)<<endl;
//        double uV,vV;
//        if (vertex_projector.SupportTypeShape1(1) == 2)
//           {
//           Point<dim> proj_point = point<dim>(vertex_projector.PointOnShape1(1));
//           vertex_projector.ParOnFaceS1(1, uV, vV);
//           TopoDS_Shape support_face = vertex_projector.SupportOnShape1(1);
//           //cout<<"@@@@@ "<<uV<<" "<<vV<<endl;
//           return std::tuple<Point<dim>, TopoDS_Shape, double, double>(
//                  proj_point, support_face, uV, vV);
//           }
//        }

    double minDistance = 1e7;
    gp_Pnt tmp_proj(0.0, 0.0, 0.0);

    unsigned int counter      = 0;
    unsigned int face_counter = 0;

    TopoDS_Shape out_shape;
    double       u = 0;
    double       v = 0;

    for (exp.Init(in_shape, TopAbs_FACE); exp.More(); exp.Next())
      {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        Bnd_Box box;
        BRepBndLib::Add(face, box);
//        Standard_Real aXmin, aYmin, aZmin, aXmax, aYmax, aZmax;
//        box.Get(aXmin, aYmin, aZmin, aXmax, aYmax, aZmax);
//        cout<<origin<<endl;
//        cout<<"X ----> ["<<aXmin<<","<<aXmax<<"]"<<endl;
//        cout<<"Y ----> ["<<aYmin<<","<<aYmax<<"]"<<endl;
//        cout<<"Z ----> ["<<aZmin<<","<<aZmax<<"]"<<endl;
        if (!box.IsOut(Pproj))
           {
            // the projection function needs a surface, so we obtain the
            // surface upon which the face is defined
            Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);
    
            ShapeAnalysis_Surface projector(SurfToProj);
            gp_Pnt2d proj_params = projector.ValueOfUV(point(origin), tolerance);

            SurfToProj->D0(proj_params.X(), proj_params.Y(), tmp_proj);

            double distance = point<dim>(tmp_proj).distance(origin);
//            cout<<face_counter<<" distance: "<<distance<<endl;
            if (distance < minDistance)
              {
                minDistance = distance;
                Pproj       = tmp_proj;
                out_shape   = face;
                u           = proj_params.X();
                v           = proj_params.Y();
                ++counter;
              }
          } 
        ++face_counter;
      }
//    cout<<counter<<" vs "<<face_counter<<endl;
    // face counter tells us if the shape contained faces: if it does, there is
    // no need to loop on edges. Even if the closest point lies on the boundary
    // of a parametric surface, we need in fact to retain the face and both u
    // and v, if we want to use this method to retrieve the surface normal
    if (face_counter == 0)
      for (exp.Init(in_shape, TopAbs_EDGE); exp.More(); exp.Next())
        {
          TopoDS_Edge edge = TopoDS::Edge(exp.Current());
          if (!BRep_Tool::Degenerated(edge))
            {
              TopLoc_Location L;
              Standard_Real   First;
              Standard_Real   Last;

              // the projection function needs a Curve, so we obtain the
              // curve upon which the edge is defined
              Handle(Geom_Curve) CurveToProj =
                BRep_Tool::Curve(edge, L, First, Last);

              GeomAPI_ProjectPointOnCurve Proj(point(origin), CurveToProj);
              unsigned int                num_proj_points = Proj.NbPoints();
              if ((num_proj_points > 0) && (Proj.LowerDistance() < minDistance))
                {
                  minDistance = Proj.LowerDistance();
                  Pproj       = Proj.NearestPoint();
                  out_shape   = edge;
                  u           = Proj.LowerDistanceParameter();
                  ++counter;
                }
            }
        }

    Assert(counter > 0, ExcMessage("Could not find projection points."));
    return std::tuple<Point<dim>, TopoDS_Shape, double, double>(
      point<dim>(Pproj), out_shape, u, v);
  }




  std::tuple<Point<3>, Tensor<1, 3>, double, double>
  my_closest_point_and_differential_forms(const TopoDS_Shape &in_shape,
                                       const Point<3> &    origin,
                                       const double        tolerance)

  {
    std::tuple<Point<3>, TopoDS_Shape, double, double> shape_and_params =
      my_project_point_and_pull_back(in_shape, origin, tolerance);

    TopoDS_Shape &out_shape = std::get<1>(shape_and_params);
    double &      u         = std::get<2>(shape_and_params);
    double &      v         = std::get<3>(shape_and_params);

    // just a check here: the number of faces in out_shape must be 1, otherwise
    // something is wrong
    std::tuple<unsigned int, unsigned int, unsigned int> numbers =
      count_elements(out_shape);
    (void)numbers;

    Assert(
      std::get<0>(numbers) > 0,
      ExcMessage(
        "Could not find normal: the shape containing the closest point has 0 faces."));
    Assert(
      std::get<0>(numbers) < 2,
      ExcMessage(
        "Could not find normal: the shape containing the closest point has more than 1 face."));


    TopExp_Explorer exp;
    exp.Init(out_shape, TopAbs_FACE);
    TopoDS_Face face = TopoDS::Face(exp.Current());
    return push_forward_and_differential_forms(face, u, v, tolerance);
  }



  template <int dim>
  Point<dim>
  my_line_intersection(const TopoDS_Shape &  in_shape,
                       const Point<dim> &    origin,
                       const Tensor<1, dim> &direction,
                       const double          tolerance)
  {
    // translating original Point<dim> to gp point
    //cout<<"Line orig: "<<origin<<"  Line dir: "<<direction<<endl;
        
    gp_Pnt P0 = point(origin);
    gp_Ax1 gpaxis(P0,
                  gp_Dir(direction[0],
                         dim > 1 ? direction[1] : 0,
                         dim > 2 ? direction[2] : 0));
    gp_Lin line(gpaxis);
//    TopExp_Explorer exp;
//    for (exp.Init(in_shape, TopAbs_FACE); exp.More(); exp.Next())
//      {
//        TopoDS_Face face = TopoDS::Face(exp.Current());
//        Bnd_Box box;
//        BRepBndLib::Add(face, box);
//        Standard_Real aXmin, aYmin, aZmin, aXmax, aYmax, aZmax;
////        box.Get(aXmin, aYmin, aZmin, aXmax, aYmax, aZmax);
////        cout<<origin<<endl;
////        cout<<"X ----> ["<<aXmin<<","<<aXmax<<"]"<<endl;
////        cout<<"Y ----> ["<<aYmin<<","<<aYmax<<"]"<<endl;
////        cout<<"Z ----> ["<<aZmin<<","<<aZmax<<"]"<<endl;
//        if (!box.IsOut(Pproj))
//           {
//           }

    // destination point
    gp_Pnt Pproj(0.0, 0.0, 0.0);
    Point<dim> result;
    double     minDistance = 1e7;
       
    TopExp_Explorer exp;
    unsigned int cross_count = 0;
    unsigned int face_count = 0;
    for (exp.Init(in_shape, TopAbs_FACE); exp.More(); exp.Next())
       {
        face_count++;
        TopoDS_Face face = TopoDS::Face(exp.Current());
        Bnd_Box box_face;
        BRepBndLib::Add(face, box_face);
        //box_face.Enlarge(0.01*box_face.SquareExtent());
        Point<3> corner_min = point<3>(box_face.CornerMin());
        Point<3> corner_max = point<3>(box_face.CornerMax());
        Point<3> hit_point;
        if (HitBoundingBox(corner_min, corner_max, origin, -direction, hit_point) ||
            HitBoundingBox(corner_min, corner_max, origin,  direction, hit_point))
           {
//           cout<<"Face: "<<face_count<<endl;
//           cout<<"BBox min: "<<corner_min<<"  BBox max: "<<corner_max<<endl;
//           cout<<"Line orig: "<<origin<<"  Line dir: "<<direction<<endl;
//           cout<<"Hit? "<<HitBoundingBox(corner_min, corner_max, origin, -direction, hit_point)<<endl;
           cross_count++;
           IntCurvesFace_ShapeIntersector Inters;
           Inters.Load(face, tolerance);
           Inters.Perform(line, -RealLast(), +RealLast());
           Assert(Inters.IsDone(), ExcMessage("Could not project point."));
           for (int i = 0; i < Inters.NbPnt(); ++i)
               {
               const double distance = point(origin).Distance(Inters.Pnt(i + 1));
               // cout<<"Point "<<i<<": "<<point(Inters.Pnt(i+1))<<"  distance:
               // "<<distance<<endl;
               if (distance < minDistance)
                  {
                  minDistance = distance;
                  result      = point<dim>(Inters.Pnt(i + 1));
                  }
               }
           }
       }
    //cout<<"Result: "<<result<<endl;
    //cout<<"Cross count: "<<cross_count<<" out of  "<<face_count<<" faces"<<endl;
    return result;
    }




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

