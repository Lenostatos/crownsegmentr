// This file is part of crownsegmentr, an R package for identifying tree crowns
// within 3D point clouds.
//
// Copyright (C) 2020-2021 Leon Steinmeier, Nikolai Knapp, UFZ Leipzig
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

//' @describeIn calculate_modes_normalized Can take either a single value or
//'     raster data for both the ground height and the
//'     \code{crown_diameter_to_tree_height} and
//'     \code{crown_height_to_tree_height} parameters.
//'
//' @param ground_height_data A list containing either a single ground height
//'     value (named "value") or a set of elements that make up a ground height
//'     raster covering the whole area of the point cloud. Such a set has to
//'     consist of the named elements described in the section "Raster argument
//'     structure" below.
//' @param crown_diameter_to_tree_height_data A list containing either a single
//'     numeric value (named "value") or the data for a raster of values (see
//'     section "Raster argument structure" below for how the raster data has to
//'     be stored in the list). The values indicate the estimated ratio of crown
//'     diameter to tree height for the whole plot or individual raster pixels
//'     respectively.
//' @param crown_height_to_tree_height_data A list containing either a single
//'     numeric value (named "value") or the data for a raster of values (see
//'     section "Raster argument structure" below for how the raster data has to
//'     be stored in the list). The values indicate the estimated ratio of crown
//'     height to tree height for the whole plot or individual raster pixels
//'     respectively.
//'
//' @section Raster argument structure:
//'     Raster data has to be passed as a list comprising the following named
//'     elements:
//'     \itemize{
//'         \item values: Numeric vector holding the values.
//'         \item num_rows: Integer number indicating the number of rows.
//'         \item num_cols: Integer number indicating the number of columns.
//'         \item x_min: Number indicating the lowest x coordinate covered.
//'         \item x_max: Number indicating the largest x coordinate covered.
//'         \item y_min: Number indicating the lowest y coordinate covered.
//'         \item y_max: Number indicating the largest y coordinate covered.
//'     }
// [[Rcpp::export]]
Rcpp::DataFrame calculate_modes_flexible (
    const Rcpp::DataFrame &coordinate_table,
    const spatial::coordinate_t &min_point_height_above_ground,
    const Rcpp::List &ground_height_data,
    const Rcpp::List &crown_diameter_to_tree_height_data,
    const Rcpp::List &crown_height_to_tree_height_data,
    const spatial::distance_t &centroid_convergence_distance,
    const int max_num_centroids_per_mode,
    const bool show_progress_bar
) {
    // Convert the coordinate table to an array of point objects.
    std::vector< spatial::point_3d_t > points {
        ams3d_R_interface_util::create_point_objects_from( coordinate_table )
    };

    // Convert the ground height data into a shared pointer to a raster object.
    std::shared_ptr< spatial::I_Raster< spatial::coordinate_t > >
    ground_height_grid_ptr {
            ams3d_R_interface_util::convert_list_argument_to_raster_double (
            ground_height_data
        )
    };

    // Convert the crown diameter to tree height data into a unique pointer to a
    // raster object.
    std::unique_ptr< spatial::I_Raster< double > >
    crown_diameter_to_tree_height_grid_ptr {
        ams3d_R_interface_util::convert_list_argument_to_raster_double (
            crown_diameter_to_tree_height_data
        )
    };

    // Convert the crown height to tree height data into a unique pointer to a
    // raster object.
    std::unique_ptr< spatial::I_Raster< double > >
    crown_height_to_tree_height_grid_ptr {
        ams3d_R_interface_util::convert_list_argument_to_raster_double (
            crown_height_to_tree_height_data
        )
    };

    // Create a raster of kernel bottom heights above ground.
    std::shared_ptr< spatial::I_Raster< spatial::coordinate_t > >
    kernel_bottom_height_above_ground_grid_ptr {
        ams3d::_Kernel::bottom_height_above_ground_grid_with (
            min_point_height_above_ground,
            *crown_height_to_tree_height_grid_ptr
        )
    };

    // Set up a spatial index for the point objects.
    spatial::index_for_3d_points_t point_cloud_index {
        spatial::create_index_of_above_ground (
            points,
            kernel_bottom_height_above_ground_grid_ptr,
            ground_height_grid_ptr
        )
    };

    // Set up an array for the to-be-calculated modes.
    std::vector< spatial::point_3d_t > modes{};
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
                *ground_height_grid_ptr,
                *crown_diameter_to_tree_height_grid_ptr,
                *crown_height_to_tree_height_grid_ptr,
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

//' @describeIn calculate_modes_normalized Can take either a single value or
//'     raster data for both the ground height and the
//'     \code{crown_diameter_to_tree_height} and
//'     \code{crown_height_to_tree_height} parameters and also returns centroids
//'     calculated during the mode finding.
//'
// [[Rcpp::export]]
Rcpp::List calculate_modes_plus_centroids_flexible (
    const Rcpp::DataFrame &coordinate_table,
    const spatial::coordinate_t &min_point_height_above_ground,
    const Rcpp::List &ground_height_data,
    const Rcpp::List &crown_diameter_to_tree_height_data,
    const Rcpp::List &crown_height_to_tree_height_data,
    const spatial::distance_t &centroid_convergence_distance,
    const int max_num_centroids_per_mode,
    const bool show_progress_bar
) {
    // Convert the coordinate table to an array of point objects.
    std::vector< spatial::point_3d_t > points {
        ams3d_R_interface_util::create_point_objects_from( coordinate_table )
    };

    // Convert the ground height data into a shared pointer to a raster object.
    std::shared_ptr< spatial::I_Raster< spatial::coordinate_t > >
    ground_height_grid_ptr {
        ams3d_R_interface_util::convert_list_argument_to_raster_double (
            ground_height_data
        )
    };

    // Convert the crown diameter to tree height data into a unique pointer to a
    // raster object.
    std::unique_ptr< spatial::I_Raster< double > >
    crown_diameter_to_tree_height_grid_ptr {
        ams3d_R_interface_util::convert_list_argument_to_raster_double (
            crown_diameter_to_tree_height_data
        )
    };

    // Convert the crown height to tree height data into a unique pointer to a
    // raster object.
    std::unique_ptr< spatial::I_Raster< double > >
    crown_height_to_tree_height_grid_ptr {
        ams3d_R_interface_util::convert_list_argument_to_raster_double (
            crown_height_to_tree_height_data
        )
    };

    // Create a raster of kernel bottom heights above ground.
    std::shared_ptr< spatial::I_Raster< spatial::coordinate_t > >
    kernel_bottom_height_above_ground_grid_ptr {
        ams3d::_Kernel::bottom_height_above_ground_grid_with (
            min_point_height_above_ground,
            *crown_height_to_tree_height_grid_ptr
        )
    };

    // Set up a spatial index for the point objects.
    spatial::index_for_3d_points_t point_cloud_index {
        spatial::create_index_of_above_ground (
            points,
            kernel_bottom_height_above_ground_grid_ptr,
            ground_height_grid_ptr
        )
    };

    // Set up an array for the to-be-calculated modes.
    std::vector< spatial::point_3d_t > modes{};
    modes.reserve( points.size() );

    // Set up arrays for the to-be-calculated centroids and their point indices.
    std::vector< spatial::point_3d_t > centroids{};
    std::vector< int > point_indices{};

    // Optionally set up a progress bar.
    RProgress::RProgress progress_bar;
    if (show_progress_bar) {
        progress_bar = ams3d_R_interface_util::create_progress_bar (
            points.size()
        );
        progress_bar.tick( 0 );
    };

    int point_index{ 1 }; // 1-based point index for use in R

    // Iterate over all points in the input point cloud to calculate modes for them
    for (const auto &point : points)
    {
        std::pair< spatial::point_3d_t, std::vector< spatial::point_3d_t > >
        mode_and_centroids {
            ams3d::calculate_a_single_mode_plus_centroids (
                point,
                point_cloud_index,
                min_point_height_above_ground,
                *ground_height_grid_ptr,
                *crown_diameter_to_tree_height_grid_ptr,
                *crown_height_to_tree_height_grid_ptr,
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
                { progress_bar.tick( ams3d_R_interface_constants::num_modes_per_tick ); }
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
