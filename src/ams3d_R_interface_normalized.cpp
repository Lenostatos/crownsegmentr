// This file is part of crownsegmentr, an R package for identifying tree crowns
// within 3D point clouds.
//
// Copyright (C) 2021 Leon Steinmeier, Nikolai Knapp, UFZ Leipzig
// Contact: Leon.Steinmeier@posteo.net
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

#include "ams3d_R_interface.h"

#include "spatial.h"
#include "ams3d.h"

//' Calculate modes with the AMS3D algorithm for a lidar point cloud of a forest
//'
//' Employs the 3D adaptive mean shift algorithm (Ferraz et al., 2016) to
//' estimate the mode of each point in a point cloud which is assumed to contain
//' trees. In this context the mode is a theoretical "center of mass" of a tree
//' crown that is usually located shortly below the crown apex. All points of
//' the same tree crown are expected to have the same mode.
//'
//' @param coordinate_table A data.frame. The first three columns are treated
//'     as x-, y-, and z-coordinates of a normalized airborne lidar point cloud
//'     in that order.
//' @param min_point_height_above_ground A single positive number. The minimum point height
//'     above ground at which the function will try to calculate a mode.
//' @param crown_diameter_to_tree_height,crown_height_to_tree_height Single
//'     numbers. Crown diameter and crown height to tree height ratios of the
//'     trees expected to be found in the point cloud.
//' @param show_progress_bar Boolean Scalar. Should the function show a progress
//'     bar during the computation?
//' @param centroid_convergence_distance Numeric Scalar. Distance at which it is
//'     assumed that subsequently calculated centroids have converged to the
//'     nearest mode.
//' @param max_num_centroids_per_mode Integer Scalar. Maximum number of
//'     centroids calculated before the search for the nearest mode is aborted.
//'
//' @return \code{calculate_modes} returns a data.frame with three columns that
//'     hold the x-, y-, and z-coordinates of the calculated modes. The rows are
//'     sorted according to the points in the \code{coordinate_table}.
//'
//' @references Ferraz, A., S. Saatchi, C. Mallet, and V. Meyer (2016)
//'     \emph{Lidar detection of individual tree size in tropical forests}.
//'     Remote Sensing of Environment 183:318–333.
//'     \url{https://doi.org/10.1016/j.rse.2016.05.028}.
// [[Rcpp::export]]
Rcpp::DataFrame calculate_modes_normalized (
    const Rcpp::DataFrame &coordinate_table,
    const spatial::coordinate_t &min_point_height_above_ground,
    const double crown_diameter_to_tree_height,
    const double crown_height_to_tree_height,
    const spatial::distance_t &centroid_convergence_distance,
    const int max_num_centroids_per_mode,
    const bool show_progress_bar
) {
    // Convert the coordinate table to an array of point objects.
    std::vector< spatial::point_3d_t > points {
        ams3d_R_interface_util::create_point_objects_from( coordinate_table )
    };

    // Set up a spatial index for the point objects.
    spatial::index_for_3d_points_t point_cloud_index {
        spatial::create_index_of_finite (
            points,
            ams3d::_Kernel::bottom_height_above_ground_with (
                min_point_height_above_ground,
                crown_height_to_tree_height
            )
        )
    };

    // Set up an array for the to-be-calculated modes.
    std::vector< spatial::point_3d_t > modes;
    modes.reserve( points.size() );

    // Optionally set up a progress bar.
    RProgress::RProgress progress_bar;
    if (show_progress_bar) {
        progress_bar = ams3d_R_interface_util::create_progress_bar (
            points.size()
        );
        progress_bar.tick(0);
    };

    // Calculate modes for all points in the point cloud.
    for (const auto &point : points)
    {
        // Calculate the mode.
        modes.push_back (
            ams3d::calculate_a_single_mode (
                point,
                point_cloud_index,
                min_point_height_above_ground,
                crown_diameter_to_tree_height,
                crown_height_to_tree_height,
                centroid_convergence_distance,
                max_num_centroids_per_mode
            )
        );

        if (modes.size() % ams3d_R_interface_constants::num_modes_per_tick == 0)
        {
            // Check whether the R user wants to abort the computation
            Rcpp::checkUserInterrupt();

            // Advance the progress bar.
            if (show_progress_bar)
                { progress_bar.tick(ams3d_R_interface_constants::num_modes_per_tick); }
        }
    }

    // Finish the progress bar.
    if (show_progress_bar) { progress_bar.tick( points.size() ); }


    // Set up R vectors for storing the mode coordinates.
    Rcpp::NumericVector mode_x_coords( modes.size() );
    Rcpp::NumericVector mode_y_coords( modes.size() );
    Rcpp::NumericVector mode_z_coords( modes.size() );
    // According to this:
    // https://stackoverflow.com/questions/41602024/should-i-prefer-rcppnumericvector-over-stdvector
    // converting std::vector from and to NumericVector entails a deep copy.
    // TODO Test whether filling the mode_*_coords vectors within the
    // mode_calculation loop increases performance.

    // Fill the R vectors.
    for (std::size_t i{0}; i < modes.size(); i++)
    {
        mode_x_coords[i] = spatial::get_x( modes[i] );
        mode_y_coords[i] = spatial::get_y( modes[i] );
        mode_z_coords[i] = spatial::get_z( modes[i] );
    }

    // Return the modes as an R data.frame.
    return Rcpp::DataFrame::create (
        Rcpp::Named("x") = mode_x_coords,
        Rcpp::Named("y") = mode_y_coords,
        Rcpp::Named("z") = mode_z_coords
    );
}

//' @describeIn calculate_modes_normalized Also return centroids calculated
//'     during the mode finding.
//'
//' @return The \code{calculate_modes_plus_centroids_*} functions return a list
//'     with two elements. The first element (named "mode_coordinates") is the
//'     \code{data.frame} which would have been returned if
//'     \code{calculate_modes} had been called. The second element (named
//'     "centroid_coordinates") is another \code{data.frame} holding centroid
//'     coordinates. To enable grouping of these centroids by the point they
//'     belong to, there is one additional column (named "point_index") which
//'     holds row indices of the corresponding points in the
//'     \code{coordinate_table}.
// [[Rcpp::export]]
Rcpp::List calculate_modes_plus_centroids_normalized (
    const Rcpp::DataFrame &coordinate_table,
    const spatial::coordinate_t &min_point_height_above_ground,
    const double crown_diameter_to_tree_height,
    const double crown_height_to_tree_height,
    const spatial::distance_t &centroid_convergence_distance,
    const int max_num_centroids_per_mode,
    const bool show_progress_bar
) {
    // Convert the coordinate table to an array of point objects.
    std::vector< spatial::point_3d_t > points {
        ams3d_R_interface_util::create_point_objects_from( coordinate_table )
    };

    // Set up a spatial index for the point objects.
    spatial::index_for_3d_points_t point_cloud_index {
        spatial::create_index_of_finite (
            points,
            ams3d::_Kernel::bottom_height_above_ground_with (
                min_point_height_above_ground,
                crown_height_to_tree_height
            )
        )
    };

    // Set up an array for the to-be-calculated modes.
    std::vector< spatial::point_3d_t > modes;
    modes.reserve( points.size() );

    // Set up arrays for the to-be-calculated centroids and their point indices.
    std::vector< spatial::point_3d_t > centroids;
    std::vector< int > point_indices;

    // Optionally set up a progress bar.
    RProgress::RProgress progress_bar;
    if (show_progress_bar) {
        progress_bar = ams3d_R_interface_util::create_progress_bar (
            points.size()
        );
        progress_bar.tick(0);
    };

    int point_index{1}; // 1-based point index for use in R

    // Iterate over all points in the input point cloud to calculate modes for them
    for (const auto &point : points)
    {
        std::pair <
            spatial::point_3d_t, std::vector< spatial::point_3d_t >
        > mode_and_centroids {
            ams3d::calculate_a_single_mode_plus_centroids (
                point,
                point_cloud_index,
                min_point_height_above_ground,
                crown_diameter_to_tree_height,
                crown_height_to_tree_height,
                centroid_convergence_distance,
                max_num_centroids_per_mode
            )
        };

        // Store the calculated mode.
        modes.push_back( mode_and_centroids.first );

        // Store the calculated centroids.
        centroids.insert (
            centroids.end(), // append at the end of centroids
            mode_and_centroids.second.begin(),
            mode_and_centroids.second.end()
        );

        // Store the current point index as many times as there are centroids.
        point_indices.insert (
            point_indices.end(), // append at the end of mode_indices
            mode_and_centroids.second.size(), // for the number of centroids...
            point_index // ...repeat this value
        );

        point_index++;

        if (modes.size() % ams3d_R_interface_constants::num_modes_per_tick == 0)
        {
            // Check whether the R user wants to abort the computation
            Rcpp::checkUserInterrupt();

            // Advance the progress bar.
            if (show_progress_bar)
                { progress_bar.tick(ams3d_R_interface_constants::num_modes_per_tick); }
        }
    }

    // Finish the progress bar.
    if (show_progress_bar) { progress_bar.tick( points.size() ); }


    // Return the modes and centroids to R

    // Set up mode coordinate arrays.
    Rcpp::NumericVector mode_x_coords( modes.size() );
    Rcpp::NumericVector mode_y_coords( modes.size() );
    Rcpp::NumericVector mode_z_coords( modes.size() );

    // Fill the mode coordinate arrays.
    for (std::size_t i{ 0 }; i < modes.size(); i++)
    {
        mode_x_coords[i] = spatial::get_x( modes[i] );
        mode_y_coords[i] = spatial::get_y( modes[i] );
        mode_z_coords[i] = spatial::get_z( modes[i] );
    }

    // Create a DataFrame with the mode coordinate arrays.
    Rcpp::DataFrame mode_coords{ Rcpp::DataFrame::create (
        Rcpp::Named("x") = mode_x_coords,
        Rcpp::Named("y") = mode_y_coords,
        Rcpp::Named("z") = mode_z_coords
    ) };

    // Set up centroid coordinate arrays.
    Rcpp::NumericVector centroid_x_coords( centroids.size() );
    Rcpp::NumericVector centroid_y_coords( centroids.size() );
    Rcpp::NumericVector centroid_z_coords( centroids.size() );

    // Fill the centroid coordinate arrays.
    for (std::size_t i{ 0 }; i < centroids.size(); i++)
    {
        centroid_x_coords[i] = spatial::get_x( centroids[i] );
        centroid_y_coords[i] = spatial::get_y( centroids[i] );
        centroid_z_coords[i] = spatial::get_z( centroids[i] );
    }

    // Create a DataFrame with the centroid coordinate arrays.
    Rcpp::DataFrame centroid_coords{ Rcpp::DataFrame::create (
        Rcpp::Named("x") = centroid_x_coords,
        Rcpp::Named("y") = centroid_y_coords,
        Rcpp::Named("z") = centroid_z_coords,
        Rcpp::Named("point_index") = point_indices
    ) };

    // Return the mode and centroid DataFrames in a list.
    return Rcpp::List::create (
        Rcpp::Named("mode_coordinates") = mode_coords,
        Rcpp::Named("centroid_coordinates") = centroid_coords
    );
}
