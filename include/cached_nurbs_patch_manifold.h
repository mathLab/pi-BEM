#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>

#include <unordered_map>
DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
class CachedNURBSPatchManifold
  : public OpenCASCADE::NURBSPatchManifold<dim, spacedim>
{
public:
  CachedNURBSPatchManifold(const TopoDS_Face &face,
                           const double       tolerance = 1e-7) //,
    // const types::global_dof_index cache_size=100000)
    : OpenCASCADE::NURBSPatchManifold<2, 3>(face, tolerance)
    , tolerance(tolerance)
    , face(face){
        // projections_cache.reserve(cache_size);
      };

  virtual Point<dim>
  pull_back(const Point<spacedim> &space_point) const;

  // virtual Point<3> get_new_point(const ArrayView<const Point<3>>
  // &surrounding_points,
  //                const ArrayView<const double>          &weights) const;

  // virtual void get_new_points(const ArrayView<const Point<3>>
  // &surrounding_points,
  //                const Table<2,double>                  &weights,
  //                 ArrayView<Point<3>>              new_points) const;
  // virtual Tensor<1,3> get_tangent_vector (const Point<3> &x1,
  //                    const Point<3> &x2) const;

  inline std::unordered_map<std::size_t, Point<dim>> &
  get_projections_cache()
  {
    return projections_cache;
  };

  mutable types::global_dof_index new_pull_backs, cache_hits;

private:
  mutable std::unordered_map<std::size_t, Point<dim>> projections_cache;
  double                                              tolerance;
  TopoDS_Face                                         face;
  const FlatManifold<dim, dim>                        sub_manifold;
};



DEAL_II_NAMESPACE_CLOSE