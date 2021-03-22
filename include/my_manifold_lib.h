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


#ifndef dealii_occ_my_manifold_lib_h
#define dealii_occ_my_manifold_lib_h

#include <deal.II/base/config.h>

#  include <deal.II/grid/manifold.h>
//#  include <deal.II/opencascade/manifold_lib.h>
#  include <deal.II/opencascade/utilities.h>

#include <unordered_map>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup OpenCASCADE
 * @{
 */

  /**
   * A Manifold object based on OpenCASCADE TopoDS_Shape where new points are
   * first computed by averaging the surrounding points in the same way as
   * FlatManifold does, and then projecting them using OpenCASCADE utilities
   * onto the manifold along a direction which is an estimation of the
   * surrounding points (hence mesh cell) normal.
   *
   * The direction normal to the mesh is particularly useful because it is the
   * direction in which the mesh is missing nodes. For instance, during the
   * refinement of a cell a new node is initially created around the
   * baricenter of the cell. This location somehow ensures a uniform distance
   * from the nodes of the old cell. Projecting such cell baricenter onto the
   * CAD surface in the direction normal to the original cell will then retain
   * uniform distance from the points of the original cell. Of course, at the
   * stage of mesh generation, no dof handler nor finite element are defined,
   * and such direction has to be estimated. For the case in which 8
   * surrounding points are present, 4 different triangles are identified with
   * the points assigned, and the normals of such triangles are averaged to
   * obtain the approximation of the normal to the cell.
   *
   * The case in which 2 surrounding points are present (i.e.:a cell edge is
   * being refined) is of course more tricky. The average of the CAD surface
   * normals at the 2 surrounding points is first computed, and then projected
   * onto the plane normal to the segment linking the surrounding points. This
   * again is an attempt to have the new point with equal distance with
   * respect to the surrounding points
   *
   * This class only operates with CAD faces and makes the assumption that the
   * shape you pass to it contains at least one face. If that is not the case,
   * an Exception is thrown. In debug mode there is a sanity check to make
   * sure that the surrounding points (the ones used in project_to_manifold())
   * actually live on the Manifold, i.e., calling OpenCASCADE::closest_point()
   * on those points leaves them untouched. If this is not the case, an
   * ExcPointNotOnManifold is thrown.
   *
   *
   * Notice that this type of Manifold descriptor may fail to give results if
   * the triangulation to be refined is close to the boundary of the given
   * TopoDS_Shape, or when the normal direction estimated from the surrounding
   * points does not intersect the shape.  An exception is thrown when this
   * happens.
   *
   * @author Luca Heltai, Andrea Mola, 2011--2014.
   */
  template <int dim, int spacedim>
  class MyNormalToMeshProjectionManifold : public FlatManifold<dim, spacedim>
  {
  public:
    /**
     * Construct a Manifold object which will project points on the
     * TopoDS_Shape @p sh, along a direction which is approximately normal to
     * the mesh cell.
     */
    MyNormalToMeshProjectionManifold(const TopoDS_Shape &sh,
                                   const double        tolerance = 1e-7);

    /**
     * Clone the current Manifold.
     */
    virtual std::unique_ptr<Manifold<dim, spacedim>>
    clone() const override;

    /**
     * Perform the actual projection onto the manifold. This function, in
     * debug mode, checks that each of the @p surrounding_points is within
     * tolerance from the given TopoDS_Shape. If this is not the case, an
     * exception is thrown.
     */
    virtual Point<spacedim>
    project_to_manifold(
      const ArrayView<const Point<spacedim>> &surrounding_points,
      const Point<spacedim> &                 candidate) const override;

    
  protected:

    mutable std::unordered_map <std::size_t, Tensor<1, spacedim> > projections_cache;

    /**
     * The topological shape which is used internally to project points. You
     * can construct such a shape by calling the OpenCASCADE::read_IGES()
     * function, which will create a TopoDS_Shape with the geometry contained
     * in the IGES file.
     */
    const TopoDS_Shape sh;
    
     /**
     * .
     */

    /**
     * Relative tolerance used by this class to compute distances.
     */
    const double tolerance;
  };


/*@}*/

DEAL_II_NAMESPACE_CLOSE


#endif // dealii_occ_manifold_lib_h
