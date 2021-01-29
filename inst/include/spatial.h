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

#ifndef SPATIAL_H
#define SPATIAL_H

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

#include <memory>  // for std::shared_ptr
#include <vector>  // for std::vector

/** A facade to some boost::geometry functionality. */
namespace spatial
{
    namespace geom = boost::geometry;

    namespace constants
    {
        inline constexpr int max_num_elements_per_r_tree_node{ 16 };
    }

    using coordinate_t = double; /**< Data type for coordinate values. */
    using distance_t = double;   /**< Data type for distance values. */

    using point_2d_t = geom::model::d2::point_xy< coordinate_t >;
    using point_3d_t = geom::model::point<
        coordinate_t, 3, geom::cs::cartesian
    >;

    namespace internal
    {
        inline point_2d_t get_xy_point(const point_3d_t & point_3d)
        {
            return { geom::get<0>(point_3d), geom::get<1>(point_3d) };
        }
    }

    /** Data type for an R-tree storing points. */
    using r_tree_for_3d_points_t = geom::index::rtree<
        point_3d_t,
        geom::index::rstar<constants::max_num_elements_per_r_tree_node>
    >;

    /** Data type for a std::shared_ptr to a point. */
    using ptr_to_const_3d_point_t = std::shared_ptr< const point_3d_t >;

    inline ptr_to_const_3d_point_t make_ptr_to_const_3d_point(
        const point_3d_t &point
    ) {
        return std::make_shared< const point_3d_t >(point);
    }

    namespace internal
    {
        // The function object (aka functor) below was adapted from an example
        // in the boost/geometry documentation at:
        // https://www.boost.org/doc/libs/1_72_0/libs/geometry/doc/html/geometry/spatial_indexes/rtree_examples/using_indexablegetter_function_object___storing_indexes_of_external_container_s_elements.html

        /** \brief A functor for dereferencing pointers to points. */
        struct get_point_from_ptr_predicate
        {
        public:
            // This type declaration is required by the R-tree implementation.
            using result_type = const point_3d_t&;

            result_type operator()(const ptr_to_const_3d_point_t &pointer) const
            {
                return *pointer;
            }
        };
    }

    /** Data type for an R-tree storing pointers to points. */
    using r_tree_for_ptrs_to_const_3d_points_t = geom::index::rtree<
        ptr_to_const_3d_point_t,
        geom::index::rstar< constants::max_num_elements_per_r_tree_node >,
        internal::get_point_from_ptr_predicate
    >;

    using line_t = geom::model::linestring< point_3d_t >;

    using segment_t = geom::model::segment< point_3d_t >;
    using segment_point_pair_t = std::pair< segment_t, point_3d_t >;

    using r_tree_for_segment_point_pairs_t = geom::index::rtree<
        segment_point_pair_t,
        geom::index::rstar< constants::max_num_elements_per_r_tree_node >
    >;

    using box_t = geom::model::box< point_3d_t >;

    inline box_t bounding_box(
        const r_tree_for_ptrs_to_const_3d_points_t &r_tree
    ) {
        return r_tree.bounds();
    }

    inline std::vector< ptr_to_const_3d_point_t > get_points_intersecting_box (
        const box_t &box,
        const r_tree_for_ptrs_to_const_3d_points_t &point_cloud
    ) {
        std::vector< ptr_to_const_3d_point_t > points_in_box;

        point_cloud.query(
            geom::index::intersects(box), std::back_inserter(points_in_box)
        );

        return points_in_box;
    }

    /** \brief Queries \p point_cloud for all points within a sphere with a
     *  certain \p radius that is centered on \p center.
     *
     *  \param radius Any point with a distance to \p center that is equal to or
     *      lesser than this radius is considered to be in the sphere. This
     *      includes \p center itself.
     */
    std::vector< ptr_to_const_3d_point_t > get_points_intersecting_sphere (
        const r_tree_for_ptrs_to_const_3d_points_t &point_cloud,
        const point_3d_t &center, distance_t radius
    );

    namespace internal
    {
        /** \brief A functor for querying wether the distance between any point
         *      and the point of this functor is smaller than a treshold.
         */
        class within_distance_predicate
        {
        private:
            const point_3d_t &_point;
            const distance_t _comparable_distance;

        public:
            within_distance_predicate(
                const point_3d_t &point, const distance_t distance
            ):
                _point{ point },
                _comparable_distance{ geom::comparable_distance(
                    point_3d_t{ 0, 0, 0 }, point_3d_t{ 0, 0, distance }
                ) }
            {}

            bool operator()(const point_3d_t &point) const
            {
                return geom::comparable_distance(point, _point)
                <= _comparable_distance;
            }

            bool operator()(const ptr_to_const_3d_point_t &point) const
            {
                return geom::comparable_distance(*point, _point)
                <= _comparable_distance;
            }
        };
    }

    std::vector< point_3d_t > get_points_intersecting_vertical_cylinder (
        const r_tree_for_3d_points_t &point_cloud,
        const point_2d_t &xy_center, distance_t radius,
        distance_t bottom_height, distance_t top_height
    );

    namespace internal
    {
        /** \brief A functor for querying wether the distance on the x-y-plane
         *      between any point and the point of this functor is smaller than
         *      a treshold.
         */
        class within_xy_distance_predicate
        {
        private:
            const point_2d_t _xy_point;
            const distance_t _comparable_distance;

        public:
            within_xy_distance_predicate(
                const point_2d_t &xy_point, const distance_t distance
            ):
            _xy_point{ xy_point },
            _comparable_distance{ geom::comparable_distance(
                point_2d_t{ 0, 0 }, point_2d_t{ 0, distance }
            ) }
            {}

            bool operator()(const point_3d_t &point) const
            {
                return geom::comparable_distance(get_xy_point(point), _xy_point)
                <= _comparable_distance;
            }

            bool operator()(const ptr_to_const_3d_point_t &point) const
            {
                return geom::comparable_distance(get_xy_point(*point), _xy_point)
                <= _comparable_distance;
            }
        };
    }

    /** \brief Calculates the weighted arithmetic mean of a set of points. */
    point_3d_t weighted_mean_of_points(
        const std::vector< point_3d_t > &points,
        const std::vector< double > &weights
    );

    // The functions below directly forward some of the boost::geometry
    // functionality.
    template< typename Geometry1, typename Geometry2 >
    inline distance_t distance(const Geometry1 &geom1, const Geometry2 &geom2)
    {
        return geom::distance(geom1, geom2);
    }

    template< typename Geometry >
    inline coordinate_t get_x(const Geometry &geom)
    {
        return geom::get<0>(geom);
    }

    template< typename Geometry >
    inline coordinate_t get_y(const Geometry &geom)
    {
        return geom::get<1>(geom);
    }

    template< typename Geometry >
    inline coordinate_t get_z(const Geometry &geom)
    {
        return geom::get<2>(geom);
    }
}

#endif  // define SPATIAL_H
