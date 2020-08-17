#include "ams3d_R_interface.h"

// [[Rcpp::interfaces(r)]]

#include "spatial.h"
#include "ams3d.h"

#include <chrono>

//' Calculate the modes of points in a forest point cloud
//'
//' Employs the 3D adaptive mean shift algorithm to estimate the mode of each
//' point in a point cloud that contains trees. In this context the mode is a
//' theoretical "center of mass" of a tree crown that is usually located shortly
//' below the crown apex. All points in the same tree crown are expected to have
//' the same mode.
//'
//' @section TODO:
//' Add reference to Ferraz et. al.
//'
//' @param point_cloud A data.frame. The first three columns are treated as
//'     x-, y-, and z-coordinates in that order.
//' @param crown_diameter_2_tree_height,crown_height_2_tree_height Numeric
//'     Scalars. Estimates of the crown diameter and crown height to tree height
//'     ratios that are common for the trees in \code{point_cloud}.
//' @param verbose Boolean Scalar. Should the function print runtime information
//'     to the console?
//'
//' @return A data.frame with three columns that hold the x-, y-, and
//'     z-coordinates of the calculated modes. The rows are sorted according to
//'     the points in \code{point_cloud} that the modes belong to.
// [[Rcpp::export]]
Rcpp::DataFrame calculate_modes(
    const Rcpp::DataFrame &point_cloud,
    const double &crown_diameter_2_tree_height,
    const double &crown_height_2_tree_height,
    bool verbose
) {
    Rcpp::NumericVector x_coords( point_cloud[0] );
    Rcpp::NumericVector y_coords( point_cloud[1] );
    Rcpp::NumericVector z_coords( point_cloud[2] );

    // Create points from the coordinates
    std::vector< spatial::point_3d_t > points;
    points.reserve(point_cloud.nrow());

    for (int i{ 0 }; i < point_cloud.nrow(); i++)
    {
        points.push_back(
            spatial::point_3d_t{ x_coords[i], y_coords[i], z_coords[i] }
        );
    }

    // Calculate modes for the points
    std::vector< spatial::point_3d_t > modes{
        ams3d::calculate_modes(
            points,
            crown_diameter_2_tree_height,
            crown_height_2_tree_height,
            verbose, Rcpp::Rcout
        )
    };

    // Return the modes as an R data.frame
    std::chrono::steady_clock::time_point begin;
    if (verbose)
    {
        Rcpp::Rcout << "Write data from C++ back to Rcpp objects...\n";
        begin = std::chrono::steady_clock::now();
    }

    Rcpp::NumericVector mode_x_coords(modes.size());
    Rcpp::NumericVector mode_y_coords(modes.size());
    Rcpp::NumericVector mode_z_coords(modes.size());

    for (std::size_t i{ 0 }; i < modes.size(); i++)
    {
        mode_x_coords[i] = spatial::get_x(modes[i]);
        mode_y_coords[i] = spatial::get_y(modes[i]);
        mode_z_coords[i] = spatial::get_z(modes[i]);
    }

    if (verbose) {
        auto writing_R_data_millisecond_duration{
            std::chrono::duration_cast< std::chrono::milliseconds >(
                std::chrono::steady_clock::now() - begin
            ).count()
        };
        Rcpp::Rcout << "Finished writing in "
            << writing_R_data_millisecond_duration
            << " milliseconds. Returning to R...\n\n";
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("mode_x") = mode_x_coords,
        Rcpp::Named("mode_y") = mode_y_coords,
        Rcpp::Named("mode_z") = mode_z_coords
    );
}

//' Calculate the modes of points in a forest point cloud
//'
//' Employs the 3D adaptive mean shift algorithm to estimate the mode of each
//' point in a point cloud that contains trees. In this context the mode is a
//' theoretical "center of mass" of a tree crown that is usually located shortly
//' below the crown apex. All points in the same tree crown are expected to have
//' the same mode.
//'
//' @section TODO:
//' Add reference to Ferraz et. al.
//'
//' @param point_cloud A data.frame. The first three columns are treated as
//'     x-, y-, and z-coordinates in that order.
//' @param crown_diameter_2_tree_height,crown_height_2_tree_height Numeric
//'     Scalars. Estimates of the crown diameter and crown height to tree height
//'     ratios that are common for the trees in \code{point_cloud}.
//' @param verbose Boolean Scalar. Should the function print runtime information
//'     to the console?
//'
//' @return A list with two elements. The first one ("mode_coords") is a
//'     data.frame with three columns holding the x-, y-, and z-coordinates of
//'     the calculated modes. The second one ("centroid_paths") is another list
//'     that contains one data.frame of point coordinates for each mode.
//'     The first point in each data.frame is always one of the points from the
//'     input point cloud while the other coordinates represent the centroids
//'     that were calculated while searching for that point's mode.
//'     There is one additional column named "mode_index" in each data.frame that
//'     holds a one-based index of the mode that the centroids belong to.
// [[Rcpp::export]]
Rcpp::List calculate_modes_and_centroid_paths(
    const Rcpp::DataFrame &point_cloud,
    const double &crown_diameter_2_tree_height,
    const double &crown_height_2_tree_height,
    bool verbose
) {
    Rcpp::NumericVector x_coords( point_cloud[0] );
    Rcpp::NumericVector y_coords( point_cloud[1] );
    Rcpp::NumericVector z_coords( point_cloud[2] );

    // Create points from the coordinates
    std::vector< spatial::point_3d_t > points;
    points.reserve(point_cloud.nrow());

    for (int i{ 0 }; i < point_cloud.nrow(); i++)
    {
        points.push_back(
            spatial::point_3d_t{ x_coords[i], y_coords[i], z_coords[i] }
        );
    }

    // Calculate modes for the points and return centroid paths as well
    auto mode_ptrs_and_centroid_path_map{
        ams3d::calculate_modes_w_centroid_paths(
            points,
            crown_diameter_2_tree_height,
            crown_height_2_tree_height,
            verbose, Rcpp::Rcout
        )
    };

    std::chrono::steady_clock::time_point begin;
    if (verbose)
    {
        Rcpp::Rcout << "Write data from C++ back to Rcpp objects...\n";
        begin = std::chrono::steady_clock::now();
    }

    auto num_modes{ mode_ptrs_and_centroid_path_map.first.size() };

    Rcpp::NumericVector mode_x_coords(num_modes);
    Rcpp::NumericVector mode_y_coords(num_modes);
    Rcpp::NumericVector mode_z_coords(num_modes);

    Rcpp::List centroid_paths(mode_ptrs_and_centroid_path_map.second.size());

    for (std::size_t i{ 0 }; i < num_modes; i++)
    {
        const auto &mode_ptr{ mode_ptrs_and_centroid_path_map.first[i] };

        mode_x_coords[i] = spatial::get_x(*mode_ptr);
        mode_y_coords[i] = spatial::get_y(*mode_ptr);
        mode_z_coords[i] = spatial::get_z(*mode_ptr);

        const std::vector< spatial::point_3d_t > &centroids{
            mode_ptrs_and_centroid_path_map.second[mode_ptr]
        };

        Rcpp::NumericVector centroid_x_coords(centroids.size());
        Rcpp::NumericVector centroid_y_coords(centroids.size());
        Rcpp::NumericVector centroid_z_coords(centroids.size());

        for (std::size_t i{ 0 }; i < centroids.size(); i++)
        {
            centroid_x_coords[i] = spatial::get_x(centroids[i]);
            centroid_y_coords[i] = spatial::get_y(centroids[i]);
            centroid_z_coords[i] = spatial::get_z(centroids[i]);
        }

        // TODO This part is super time intensive. Unfortunately I can't figure
        // out a way to reserve the needed capacity in advance.
        centroid_paths[i] = Rcpp::DataFrame::create(
            Rcpp::Named("x") = centroid_x_coords,
            Rcpp::Named("y") = centroid_y_coords,
            Rcpp::Named("z") = centroid_z_coords,
            Rcpp::Named("mode_index") = i + 1
        );
    }

    Rcpp::DataFrame mode_coords{ Rcpp::DataFrame::create(
        Rcpp::Named("mode_x") = mode_x_coords,
        Rcpp::Named("mode_y") = mode_y_coords,
        Rcpp::Named("mode_z") = mode_z_coords
    ) };

    if (verbose) {
        auto writing_R_data_millisecond_duration{
            std::chrono::duration_cast< std::chrono::milliseconds >(
                std::chrono::steady_clock::now() - begin
            ).count()
        };
        Rcpp::Rcout << "Finished writing in "
            << writing_R_data_millisecond_duration
            << " milliseconds. Returning to R...\n\n";
    }

    return Rcpp::List::create(
        Rcpp::Named("mode_coords") = mode_coords,
        Rcpp::Named("centroid_paths") = centroid_paths
    );
}
