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
        inline constexpr int r_tree_max_num_elements_per_node{ 16 };
    }
    
    using coordinate_t = double;
    using distance_t = double;
    
    using point_2d_t = geom::model::d2::point_xy<coordinate_t>;
    using point_3d_t = geom::model::point<coordinate_t, 3, geom::cs::cartesian>;
    
    namespace internal
    {
        inline point_2d_t get_xy_point(const point_3d_t & point_3d)
        {
            return { geom::get<0>(point_3d), geom::get<1>(point_3d) };
        }
    }
    
    using r_tree_for_3d_points_t = geom::index::rtree<
        point_3d_t,
        geom::index::rstar<constants::r_tree_max_num_elements_per_node>
    >;
    
    using ptr_to_const_3d_point_t = std::shared_ptr<const point_3d_t>;
    
    inline ptr_to_const_3d_point_t make_ptr_to_const_3d_point(
        const point_3d_t &point
    ) {
        return std::make_shared<const point_3d_t>(point);
    }
    
    namespace internal
    {
        // The function object (aka functor) below was adapted from an example
        // in the boost/geometry documentation at:
        // https://www.boost.org/doc/libs/1_72_0/libs/geometry/doc/html/geometry/spatial_indexes/rtree_examples/using_indexablegetter_function_object___storing_indexes_of_external_container_s_elements.html
        struct point_getter_for_r_tree
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
    
    using r_tree_for_ptrs_to_const_3d_points_t = geom::index::rtree<
        ptr_to_const_3d_point_t,
        geom::index::rstar<constants::r_tree_max_num_elements_per_node>,
        internal::point_getter_for_r_tree
    >;
    
    using box_t = geom::model::box<point_3d_t>;
    
    inline box_t bounding_box(
        const r_tree_for_ptrs_to_const_3d_points_t &r_tree
    ) {
        return r_tree.bounds();
    }
    
    inline std::vector<ptr_to_const_3d_point_t> get_points_intersecting_box (
        const box_t &box,
        const r_tree_for_ptrs_to_const_3d_points_t &point_cloud
    ) {
        std::vector<ptr_to_const_3d_point_t> points_in_box;
        
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
    std::vector<ptr_to_const_3d_point_t> get_points_intersecting_sphere (
        const r_tree_for_ptrs_to_const_3d_points_t &point_cloud,
        const point_3d_t &center, distance_t radius
    );
    
    namespace internal
    {
        /** \brief A functor for the satisfies query to the r-tree. */
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
        
        /** \brief A functor for the satisfies query to the r-tree. */
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
    
    std::vector<point_3d_t> get_points_intersecting_vertical_cylinder (
        const r_tree_for_3d_points_t &point_cloud,
        const point_2d_t &xy_center, distance_t radius,
        distance_t bottom_height, distance_t top_height
    );
    
    point_3d_t weighted_mean_of_points(
        const std::vector<point_3d_t> &points,
        const std::vector<double> &weights
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
