// This file is part of crownsegmentr, an R package for identifying tree crowns
// within 3D point clouds.
//
// Copyright (C) 2020 Leon Steinmeier, Nikolai Knapp, UFZ
// Contact: Lenostatos@gmx.de
//
// crownsegmentr is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// crownsegmentr is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with crownsegmentr in a file called "COPYING". If not,
// see <http://www.gnu.org/licenses/>.

#include "spatial.h"

#include <iterator> // for std::back_inserter
#include <numeric>  // for std::accumulate

namespace spatial
{
    namespace geom = boost::geometry;

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
