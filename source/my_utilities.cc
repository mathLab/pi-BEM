


#  include "../include/my_utilities.h" 

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
#  include <GeomAPI_ProjectPointOnCurve.hxx>
#  include <IntCurvesFace_ShapeIntersector.hxx>




DEAL_II_NAMESPACE_OPEN

using namespace OpenCASCADE;


/* 
Fast Ray-Box Intersection
by Andrew Woo
from "Graphics Gems", Academic Press, 1990
*/


enum class Position
{
  RIGHT,
  LEFT,
  MIDDLE
};
template<int dim>
bool HitBoundingBox(const Point<dim> minB, const Point<dim> maxB, const Point<dim> origin, const Tensor<1,dim> dir, Point<dim> coord)
{
  bool inside = true;
  std::array<Position, dim> quadrant;
  unsigned int whichPlane;
  std::array<double, dim> maxT;
  std::array<double, dim> candidatePlane;

  /* Find candidate planes; this loop can be avoided if
    rays cast all from the eye(assume perpsective view) */
  for (unsigned int i=0; i<dim; i++)
  {
    if(origin(i) < minB(i)) {
      quadrant[i] = Position::LEFT;
      candidatePlane[i] = minB(i);
      inside = false;
    }else if (origin(i) > maxB(i)) {
      quadrant[i] = Position::RIGHT;
      candidatePlane[i] = maxB(i);
      inside = false;
    }else {
      quadrant[i] = Position::MIDDLE;
    }
  }
  /* Ray origin inside bounding box */
  if(inside)  {
    coord = origin;
    return (true);
  }


  /* Calculate T distances to candidate planes */
  for (unsigned int i = 0; i < dim; i++)
    if (quadrant[i] != Position::MIDDLE && dir[i] !=0.)
      maxT[i] = (candidatePlane[i]-origin(i)) / dir[i];
    else
      maxT[i] = -1.;

  /* Get largest of the maxT's for final choice of intersection */
  whichPlane = 0;
  for (unsigned int i = 1; i < dim; i++)
    if (maxT[whichPlane] < maxT[i])
      whichPlane = i;

  /* Check final candidate actually inside box */
  if (maxT[dim] < 0.) return (false);
  for (unsigned int i = 0; i < dim; i++)
    if (whichPlane != i) {
      coord[i] = origin(i) + maxT[whichPlane] *dir[i];
      if (coord[i] < minB(i) || coord[i] > maxB(i))
        return (false);
    } else {
      coord[i] = candidatePlane[i];
    }
  return (true);        /* ray hits box */
} 






  template <int dim>
  std::tuple<Point<dim>, TopoDS_Shape, double, double>
  my_project_point_and_pull_back(const TopoDS_Shape &in_shape,
                              const Point<dim> &  origin,
                              const double        tolerance)
  {
    TopExp_Explorer exp;
    gp_Pnt          Pproj = point(origin);

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
        if (!box.IsOut(Pproj))
           {
            // the projection function needs a surface, so we obtain the
            // surface upon which the face is defined
            Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);
    
            ShapeAnalysis_Surface projector(SurfToProj);
            gp_Pnt2d proj_params = projector.ValueOfUV(point(origin), tolerance);

            SurfToProj->D0(proj_params.X(), proj_params.Y(), tmp_proj);

            double distance = point<dim>(tmp_proj).distance(origin);
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
        Point<dim> corner_min = point<dim>(box_face.CornerMin());
        Point<dim> corner_max = point<dim>(box_face.CornerMax());
        Point<dim> hit_point;
        if (HitBoundingBox(corner_min, corner_max, origin, -direction, hit_point) ||
            HitBoundingBox(corner_min, corner_max, origin,  direction, hit_point))
           {
           cross_count++;
           IntCurvesFace_ShapeIntersector Inters;
           Inters.Load(face, tolerance);
           Inters.Perform(line, -RealLast(), +RealLast());
           Assert(Inters.IsDone(), ExcMessage("Could not project point."));
           for (int i = 0; i < Inters.NbPnt(); ++i)
               {
               const double distance = point(origin).Distance(Inters.Pnt(i + 1));
               if (distance < minDistance)
                  {
                  minDistance = distance;
                  result      = point<dim>(Inters.Pnt(i + 1));
                  }
               }
           }
       }
    return result;
    }



  template std::tuple<Point<2>, TopoDS_Shape, double, double> my_project_point_and_pull_back<2>(const TopoDS_Shape &in_shape, const Point<2> &  origin, const double        tolerance);
  template std::tuple<Point<3>, TopoDS_Shape, double, double> my_project_point_and_pull_back<3>(const TopoDS_Shape &in_shape, const Point<3> &  origin, const double        tolerance);
  template Point<2> my_line_intersection<2>(const TopoDS_Shape &in_shape, const Point<2> &  origin, const Tensor<1,2>& direction, const double        tolerance);;
  template Point<3> my_line_intersection<3>(const TopoDS_Shape &in_shape, const Point<3> &  origin, const Tensor<1,3>& direction, const double        tolerance);;

DEAL_II_NAMESPACE_CLOSE


