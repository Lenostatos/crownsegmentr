#include "spatial.h"

#include <iterator> // for std::back_inserter
#include <numeric>  // for std::accumulate

namespace spatial
{
    namespace geom = boost::geometry;
    
    std::vector<ptr_to_const_3d_point_t> get_points_intersecting_sphere (
        const r_tree_for_ptrs_to_const_3d_points_t &point_cloud,
        const point_3d_t &center, distance_t radius
    ) {
        // Construct a box around the sphere.
        point_3d_t min_corner_point{ center };
        geom::subtract_value(min_corner_point, radius);
        
        point_3d_t max_corner_point{ center };
        geom::add_value(max_corner_point, radius);
        
        box_t sphere_box{ min_corner_point, max_corner_point };
        
        // Set up a unary predicate object.
        // TODO Measure performance difference between lambda and unary
        // predicate.
        internal::within_distance_predicate within_distance{
            center, radius
        };
        
        // Query the point cloud.
        std::vector<ptr_to_const_3d_point_t> points_in_neighborhood;
        
        point_cloud.query(
            geom::index::intersects(sphere_box)
                && geom::index::satisfies(within_distance),
            std::back_inserter(points_in_neighborhood)
        );
        
        return points_in_neighborhood;
    }
    
    /** \note When modifying this function keep in mind that there might be a
     *  significant impact on performance since the function is called very
     *  often while calculating modes.
     */
    std::vector<point_3d_t> get_points_intersecting_vertical_cylinder (
        const r_tree_for_3d_points_t &point_cloud,
        const point_2d_t &xy_center, distance_t radius,
        distance_t bottom_height, distance_t top_height
    ) {
        
        // Construct a box around the cylinder.
        point_3d_t min_corner_point{
            geom::get<0>(xy_center) - radius,
            geom::get<1>(xy_center) - radius,
            bottom_height
        };
        point_3d_t max_corner_point{
            geom::get<0>(xy_center) + radius,
            geom::get<1>(xy_center) + radius,
            top_height
        };
        
        box_t cylinder_box{ min_corner_point, max_corner_point };
        
        // Set up a unary predicate object.
        internal::within_xy_distance_predicate within_xy_distance{
            xy_center, radius
        };
        
        // Query the point cloud.
        std::vector<point_3d_t> intersecting_points;
        
        point_cloud.query(
            geom::index::intersects(cylinder_box)
                && geom::index::satisfies(within_xy_distance),
            std::back_inserter(intersecting_points)
        );
        
        return intersecting_points;
    }
    
    point_3d_t weighted_mean_of_points(
        const std::vector< point_3d_t > &points,
        const std::vector< double > &weights
    ) {        
        point_3d_t mean_point{ 0, 0 ,0 };
        
        for(std::size_t i{ 0 }; i < points.size(); i++)
        {
            geom::add_point(
                mean_point,
                point_3d_t{
                    geom::get<0>(points[i]) * weights[i],
                    geom::get<1>(points[i]) * weights[i],
                    geom::get<2>(points[i]) * weights[i]
                }
            );
        }
        
        geom::divide_value(
            mean_point,
            std::accumulate(weights.begin(), weights.end(), double{ 0 })
        );
        
        return mean_point;
    }
}
