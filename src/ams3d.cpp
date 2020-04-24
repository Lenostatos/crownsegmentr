#include "ams3d.h"

#include "spatial.h"

// TODO update these includes
#include <map>          // for std::map
#include <vector>       // for std::vector
#include <chrono>       // for std::chrono::*
#include <memory>       // for std::make_shared
#include <ostream>      // for std::basic_ostream
#include <algorithm>    // for std:: accumulate
#include <utility>      // std::move
#include <functional>   // TODO 
#include <cmath>        // std::pow, std::fmin
#include <iterator>     // std::back_inserter
#include <numeric>      // std::accumulate

namespace ams3d
{
    std::vector< spatial::point_3d_t > calculate_modes(
        const std::vector< spatial::point_3d_t > &point_cloud,
        const double crown_diameter_2_tree_height,
        const double crown_height_2_tree_height,
        bool verbose,
        std::basic_ostream< char > &log_output
    ) {
        spatial::r_tree_for_3d_points_t r_tree{
            point_cloud.begin(), point_cloud.end()
        };
        
        // Prepare some time measurements.
        std::vector<
            std::chrono::steady_clock::duration
        > mode_calculation_times;
        if (verbose) { mode_calculation_times.reserve(point_cloud.size()); }
        
        // Calculate modes for the points.
        std::vector< spatial::point_3d_t > modes;
        
        if (verbose) { log_output << "Start calculating modes.\n"; }
        
        for (const auto &point : point_cloud)
        {
            std::chrono::steady_clock::time_point begin;
            if (verbose) { begin = std::chrono::steady_clock::now(); }
            
            modes.push_back(
                internal::calculate_point_mode(
                    point, r_tree,
                    crown_diameter_2_tree_height,
                    crown_height_2_tree_height
                )
            );
            
            if (verbose)
            {                
                auto end{ std::chrono::steady_clock::now() };
                mode_calculation_times.push_back(end - begin);
                
                if (modes.size() % 1000 == 0)
                {
                    log_output << "Calculated " << modes.size() << " of "
                    << point_cloud.size() << " modes...\n";
                }
            }
            
            // TODO Call Rcpp::checkUserInterrupt() in long running loops when
            // using this code in an R package.
        }
        
        if (verbose)
        {
            // Evaluate the time measurements.
            std::chrono::steady_clock::duration mean_duration{
                std::accumulate(
                    mode_calculation_times.begin(),
                                mode_calculation_times.end(),
                                std::chrono::steady_clock::duration::zero()
                )
                / mode_calculation_times.size()
            };
            
            log_output << "Finished calculating modes!" << "\n\n"
                << "Calculation of one mode took on average "
                << std::chrono::duration_cast<std::chrono::microseconds>(
                    mean_duration
                ).count()
                << " microseconds.\n";
            
            double estimated_minutes_for_1_mill_points{
                static_cast<double>(
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        mean_duration
                    ).count())
                * 1e6 / 1e9 / 60
            };
            
            log_output << "This scales to about "
                << estimated_minutes_for_1_mill_points << " minutes "
                << "(" << estimated_minutes_for_1_mill_points / 60 << " hours) "
                << "of processing time for 1 Million points.\n\n";
        }
        
        return modes;
    }
    
    std::unordered_map<
        spatial::ptr_to_const_3d_point_t,
        spatial::ptr_to_const_3d_point_t
    > calculate_modes(
        const std::vector<spatial::ptr_to_const_3d_point_t> &point_cloud,
        const double crown_diameter_2_tree_height,
        const double crown_height_2_tree_height,
        std::basic_ostream<char> &log_output
    ) {
        // Store the points in a spatial index.
        // TODO Create an iterator for point_cloud that can be passed directly
        // when filling the spatial index.
        std::vector<spatial::point_3d_t> simple_points;
        simple_points.reserve(point_cloud.size());
        
        for (const auto &ptr_to_const_point : point_cloud)
        {
            simple_points.push_back(*ptr_to_const_point);
        }
        
        spatial::r_tree_for_3d_points_t spatial_index{
            simple_points.begin(), simple_points.end()
        };
        
        // Free memory used by the simple_points vector.
        simple_points.clear();
        simple_points.shrink_to_fit();
        
        // Do some time measurements.
        std::vector<std::chrono::steady_clock::duration>
            mode_calculation_times;
        mode_calculation_times.reserve(point_cloud.size());
        
        // Calculate modes for the points.
        std::unordered_map<
            spatial::ptr_to_const_3d_point_t,
            spatial::ptr_to_const_3d_point_t
        > res_mode_map;
        
        log_output << "Start calculating modes.\n" ;
        
        for (const auto &ptr_to_const_point : point_cloud)
        {
            auto begin{ std::chrono::steady_clock::now() };
            
            res_mode_map.emplace(
                ptr_to_const_point,
                spatial::make_ptr_to_const_3d_point(
                    internal::calculate_point_mode(
                        *ptr_to_const_point, spatial_index,
                        crown_diameter_2_tree_height,
                        crown_height_2_tree_height
                    )
                )
            );
            
            auto end{ std::chrono::steady_clock::now() };
            
            mode_calculation_times.push_back(end - begin);
            
            if (res_mode_map.size() % 1000 == 0)
            {
                log_output << "Calculated " << res_mode_map.size() << " of "
                    << point_cloud.size() << " modes...\n";
            }
            
            // TODO Call Rcpp::checkUserInterrupt() in long running loops when
            // using this code in an R package.
        }
        
        // Evaluate the time measurements.
        std::chrono::steady_clock::duration mean_duration{
            std::accumulate(
                mode_calculation_times.begin(),
                mode_calculation_times.end(),
                std::chrono::steady_clock::duration::zero()
            )
            / mode_calculation_times.size()
        };
        
        log_output << "Finished calculating modes!" << "\n\n"
            << "Calculation of one mode took on average "
            << std::chrono::duration_cast<std::chrono::microseconds>(
                mean_duration
            ).count()
            << " microseconds.\n";
        
        double estimated_minutes_for_1_mill_points{
            static_cast<double>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    mean_duration
                ).count())
            * 1e6 / 1e9 / 60
        };
        
        log_output << "This scales to about "
            << estimated_minutes_for_1_mill_points << " minutes "
            << "(" << estimated_minutes_for_1_mill_points / 60 << " hours) "
            << "of processing time for 1 Million points.\n\n";
            
        return res_mode_map;
    }
}
    
namespace ams3d::internal
{
    spatial::point_3d_t calculate_point_mode(
        const spatial::point_3d_t &point,
        const spatial::r_tree_for_3d_points_t &point_cloud,
        const double crown_diameter_2_tree_height,
        const double crown_height_2_tree_height
    ) {
        spatial::point_3d_t current_centroid{ point };
        spatial::point_3d_t former_centroid;
        
        int current_iteration{ 0 };
        do
        {
            current_iteration++;
            
            Kernel current_kernel{
                current_centroid,
                crown_diameter_2_tree_height,
                crown_height_2_tree_height
            };
            
            // TODO Do some performance tests in release mode to see whether
            // "former_centroid = std::move(current_centroid);" is faster.
            former_centroid = current_centroid;
            current_centroid = current_kernel.calculate_centroid(point_cloud);
        } 
        while(
            spatial::distance(former_centroid, current_centroid)
                > constants::centroid_convergence_distance
            && current_iteration < constants::max_num_centroids_per_mode
        );
        
        return current_centroid;
    }
    
    Kernel::Kernel(
        const spatial::point_3d_t &center,
        const double crown_diameter_2_tree_height,
        const double crown_height_2_tree_height
    ):
        _xy_center{ spatial::get_x(center), spatial::get_y(center) },
        _radius{ crown_diameter_2_tree_height * spatial::get_z(center) / 2 },
        _height{ crown_height_2_tree_height * spatial::get_z(center) * 3.0/4.0 },
        
        // The kernel is positioned vertically asymmetric around center.
        _top_height   { spatial::get_z(center) + _height * 2.0/3.0 },
        _middle_height{ _top_height - _height * 0.5 },
        _bottom_height{ _top_height - _height }
    {}
    
    std::vector<spatial::point_3d_t> Kernel::_find_intersecting_points(
        const spatial::r_tree_for_3d_points_t &point_cloud
    ) const
    {
        return spatial::get_points_intersecting_vertical_cylinder (
            point_cloud,
            this->_xy_center, this->_radius,
            this->_bottom_height, this->_top_height
        );
    }
    
    spatial::distance_t Kernel::_calculate_relative_horizontal_distance_to_center(
        const spatial::point_3d_t &point
    ) const
    {
        spatial::distance_t absolute_distance{
            spatial::distance(
                this->_xy_center,
                spatial::point_2d_t{
                    spatial::get_x(point), spatial::get_y(point)
                }
            )
        };
        
        return absolute_distance / this->_radius;
    }
    
    spatial::distance_t Kernel::_calculate_relative_vertical_distance_to_center(
        const spatial::point_3d_t &point
    ) const
    {
        spatial::distance_t absolute_distance{
            std::abs(this->_middle_height - spatial::get_z(point))
        };
        
        return absolute_distance / (this->_height / 2);
    }
    
    double Kernel::_calculate_point_weight(const spatial::point_3d_t &point) const
    {
        return
            math_functions::gauss(
                _calculate_relative_horizontal_distance_to_center(point)
            )
            * math_functions::epanechnikov(
                _calculate_relative_vertical_distance_to_center(point)
            );
    }
    
    spatial::point_3d_t Kernel::calculate_centroid(
        const spatial::r_tree_for_3d_points_t &point_cloud
    ) const
    {
        std::vector<spatial::point_3d_t> points_in_kernel{
            _find_intersecting_points(point_cloud)
        };
        
        std::vector<double> point_weights;
        point_weights.reserve(points_in_kernel.size());
        
        std::transform(
            points_in_kernel.begin(), points_in_kernel.end(),
            std::back_inserter(point_weights),
            [&](spatial::point_3d_t const& point) -> double
            {
                return this->_calculate_point_weight(point);
            }
        );
        
        return spatial::weighted_mean_of_points(
            points_in_kernel, point_weights
        );
    }
}


// The comments below contain code of the algorithm's old version that was
// commented, formatted and played around with in order to better understand the
// algorithm.


// /** \brief Checks whether a point [PointX, PointY, PointZ] is within a cylinder
//  *  of a given radius and height centered around the point at
//  *  [CtrX, CtrY, CtrZ].
//  */
// bool InCylinder(
//     double PointX, double PointY, double PointZ,
//     double CtrX, double CtrY, double CtrZ,
//     double Radius,
//     double Height
// ){
//     return pow((PointX - CtrX), 2.0) + pow((PointY - CtrY), 2.0) <= pow(Radius, 2.0) &&
//            // i.e. horizontal distance between Point and Center is smaller than radius
//            (PointZ >= (CtrZ - (0.5 * Height))) &&
//            (PointZ <= (CtrZ + (0.5 * Height)));
//            
//     if ((pow((PointX - CtrX), 2.0) + pow((PointY - CtrY), 2.0) <= pow(Radius, 2.0)) && (PointZ >= (CtrZ - (0.5*Height))) && (PointZ <= (CtrZ + (0.5*Height))) == true) {
//         return true;
//     } else {
//         return false;
//     }
// }
// 
// // Help functions for vertical filter
// double VerticalDistance(double Height, double CtrZ, double PointZ){
//     
//   double BottomDistance{ std::abs((CtrZ - Height/4 - PointZ) / (3 * Height/8)) };
//   double TopDistance{ std::abs((CtrZ + Height/2 - PointZ) / (3 * Height/8)) };
//   
//   double MinDistance{ std::min(BottomDistance, TopDistance) };
//   return MinDistance;
// }
// 
//Equivalent R code
//distx <- function(h, CtrZ, PointZ){
//  bottomdist <- abs((CtrZ-h/4-PointZ)/(3*h/8))
//  topdist <- abs((CtrZ+h/2-PointZ)/(3*h/8))
//  mindist <- pmin(bottomdist, topdist)
//  return(mindist)
//}
// 
// // Selects all points in the upper three quarters of the cylinder.
// double VerticalMask(double Height, double CtrZ, double PointZ){
//     
//     if(CtrZ - Height/4 <= PointZ && PointZ <= CtrZ + Height/2)
//     {
//         return 1;
//     } else {
//         return 0;
//     }
//     
//     if((PointZ >= CtrZ - Height/4) && (PointZ <= CtrZ + Height/2))
//     {
//         return 1;
//     } else {
//         return 0;
//     }
// }
// 
//Equivalent R code
//maskx <- function(h, CtrZ, PointZ){
//  maskvec <- ifelse(PointZ >= CtrZ-h/4 & PointZ <= CtrZ+h/2, 1, 0)
//  return(maskvec)
//}
// 
// // Epanechnikov function for vertical filter
// double EpanechnikovFunction(double Height, double CtrZ, double PointZ){
//     
//     double Result {
//         VerticalMask(Height, CtrZ, PointZ) *
//         (1 - pow(1 - VerticalDistance(Height, CtrZ, PointZ), 2.0))
//     };
//     // i.e. for every point in the upper three quarters of the cylinder do:
//     // 1 - (1 - VerticalDistance(Point))²
//         
//     return Result;
// }
// 
//Equivalent R code
//Epanechnikov <- function(h, CtrZ, PointZ){
//  output <- maskx(h, CtrZ, PointZ)*(1-(1-distx(h, CtrZ, PointZ))^2)
//  return(output)
//}
// 
// Gauss function for horizontal filter
// double GaussFunction(
//     double Width,
//     double CtrX, double CtrY,
//     double PointX, double PointY
// ){
//     double Distance = pow(
//         pow((PointX - CtrX), 2.0) +
//         pow((PointY - CtrY), 2.0),
//         0.5
//     );
//     
//     double NormDistance = Distance / Width;
//     
//     double Result = std::exp(-5.0 * pow(NormDistance, 2.0));
//     return Result;
// }
// 
//Equivalent R code
//gauss <- function(w, CtrX, CtrY, PointX, PointY){
//  distance <- ((PointX-CtrX)^2+(PointY-CtrY)^2)^0.5
//  norm.distance <- distance/w
//  output <- exp(-5*norm.distance^2)
//  return(output)
//}

