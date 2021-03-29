// This file is part of crownsegmentr, an R package for identifying tree crowns
// within 3D point clouds.
//
// Copyright (C) 2020-2021 Leon Steinmeier, Nikolai Knapp, UFZ
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

#ifndef AMS3D_H
#define AMS3D_H

#include "spatial.h"

#include <cmath>            // std::exp and std::pow

/** \brief Contains an implementation of the mean shift algorithm adapted for
 *  the use case of identifying tree crowns in 3D LiDAR point clouds as
 *  described by Ferraz et. al 2012.
 *
 *  \par
 *  ams3d is short for "3D adaptive mean shift algorithm", a term used in Ferraz
 *  et. al 2016.
 *
 *  \par How it works
 *  The algorithm tries to find the mode, i.e. kind of the location of the tree
 *  crown for each point in the point cloud. It does this by following the steps
 *  below for each point:
 *
 *  1. It constructs a so called kernel, i.e. a space with the shape of a
 *  vertical cylinder with the point at it's center.
 *  2. The lower quarter of the kernel is truncated.
 *  3. All points in the point cloud that intersect with the truncated kernel
 *  are collected.
 *  4. A weighted centroid is calculated for the collected points, weighted
 *  meaning here that points that are closer to the kernel's center get a higher
 *  weight. There are two different weight functions. One for the points'
 *  horizontal distance to the kernel's center and one for the points' vertical
 *  distance to that center.
 *  5. Now a new kernel is constructed around the calculated centroid.
 *  6. Steps 2 to 5 are repeated until consecutive centroids converge, i.e. when
 *  the distance between consecutive centroids gets smaller than a certain
 *  threshold.
 *  7. The centroid calculated last is then treated as the original point's
 *  mode.
 *
 *  It has been observed that the modes of points that belong to the same tree
 *  crown cluster shortly below the crown apex. In order to assign crown IDs to
 *  points, these clusters have to be identified with another algorithm, e.g.
 *  DBSCAN.
 *
 *
 *  \par On Compliance with the Equations Published by Ferraz et. al 2012
 *  The equations in Ferraz et. al assume a symmetric kernel with a point of
 *  the point cloud as it's center. For calculating the kernel's centroid,
 *  equations (13) and (14) are designed in such a way that points in the
 *  kernel's lower quarter are ignored.
 *      This implementation deviates from those equations in that it
 *  uses kernel objects that directly model the upper three quarters of a
 *  symmetric kernel. The respective equations are modified accordingly.
 *      The calculation of point weights on the kernel's vertical profile also
 *  differs from equations (13) and (14) in Ferraz et. al. In equation (13) the
 *  calculation of a point's vertical distance to the kernel's boundary is shown
 *  together with  the normalization of that distance. This relative distance is
 *  then subtracted from one in equation (14), which effectively inverts it to
 *  the relative distance to the kernel's center.
 *      In this implementation the relative vertical distance to the kernel's
 *  center is calculated directly.
 */
namespace ams3d
{
    /** \brief Calculates the mode of \p point within \p indexed_point_cloud. */
    spatial::point_3d_t calculate_a_single_mode (
        const spatial::point_3d_t &point,
        const spatial::index_for_3d_points_t &indexed_point_cloud,
        const double crown_diameter_to_tree_height,
        const double crown_height_to_tree_height,
        const double centroid_convergence_distance,
        const int max_num_centroids_per_mode
    );

    /** \brief Calculates the mode of \p point within \p indexed_point_cloud.
     *
     *  Same as calculate_a_single_mode but also returns the centroids.
     *
     *  \return A pair with the calculated mode at the first position and an
     *      array of the calculated centroids at the second position.
     */
    std::pair <
        spatial::point_3d_t, std::vector< spatial::point_3d_t >
    > calculate_a_single_mode_plus_centroids (
        const spatial::point_3d_t &point,
        const spatial::index_for_3d_points_t &indexed_point_cloud,
        const double crown_diameter_to_tree_height,
        const double crown_height_to_tree_height,
        const spatial::distance_t &centroid_convergence_distance,
        const int max_num_centroids_per_mode
    );


    // ===================
    // Internal Components
    // ===================

    /** \brief Models a kernel with the shape of a three-dimensional
     *      vertical cylinder.
     */
    class _Kernel
    {
    private:
        const spatial::point_2d_t _xy_center;
        const spatial::distance_t _radius;
        const spatial::distance_t _height;
        const spatial::distance_t _top_height;
        const spatial::distance_t _middle_height;
        const spatial::distance_t _bottom_height;

        /** Searches for points in \p point_cloud that intersect with the
         *  kernel.
         */
        std::vector<spatial::point_3d_t> _find_intersecting_points_in (
            const spatial::index_for_3d_points_t &point_cloud
        ) const;

        /** Calculate \p point's distance to the kernel's center on the
         *  x-y-plane, normalized with the kernel's radius.
         *
         *  Analogous to the argument to the function g^s in equation (15)
         *  in Ferraz et al. 2012.
         */
        spatial::distance_t _calculate_relative_horizontal_distance_of_center_to (
            const spatial::point_3d_t &point
        ) const;

        /** Calculate \p point's distance to the kernel center along the
         *  z-axis, normalized with half the kernel's height.
         *
         *  Analogous to parts of equation (13) and (14) in
         *  Ferraz et al. 2012.
         */
        spatial::distance_t _calculate_relative_vertical_distance_of_center_to (
            const spatial::point_3d_t &point
        ) const;

        /** Calculate the weight of a \p point inside the kernel according to
         *  the kernel's horizontal and vertical profile.
         */
        double _calculate_point_weight_of (
            const spatial::point_3d_t &point
        ) const;

    public:
        /** \brief Constructs an asymmetric kernel around \p center.
         *
         *  @param center A three-dimensional point.
         *  @param crown_diameter_to_tree_height The estimated relationship
         *      between crown diameter and tree height.
         *  @param crown_height_to_tree_height The estimated relationship
         *      between crown height and tree height.
         */
        _Kernel (
            const spatial::point_3d_t &center,
            const double crown_diameter_to_tree_height,
            const double crown_height_to_tree_height
        );

        /** Returns the kernel's weighted centroid given a point cloud. */
        spatial::point_3d_t calculate_centroid_in (
            const spatial::index_for_3d_points_t &point_cloud
        ) const;

    };

    namespace _math_functions
    {
        inline constexpr double gaussian_gamma{-5};

        /** The gaussian function f(x) = exp(gaussian_gamma * x^2).
         *
         *  Analogous to equation (11) in Ferraz et al. 2012.
         */
        inline double gauss( const double x ) {
            return std::exp( gaussian_gamma * std::pow(x, 2) );
        }

        /** The epanechnikov distribution function f(x) = 1 - x^2 .
         *
         *  Analogous to parts of equation (14) in Ferraz et al. 2012.
         */
        inline double epanechnikov( const double x ) {
            return 1 - std::pow( x, 2 );
        }
    }
}

#endif  // define AMS3D_H
