#include "ams3d_R_interface.h"

// [[Rcpp::interfaces(r)]]

#include "spatial.h"
#include "ams3d.h"

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
//' @param verbose Boolean. Should the function print runtime information to the
//'     console?
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

    // Create points from the coordinates.
    std::vector< spatial::point_3d_t > points;
    points.reserve(point_cloud.nrow());

    for (int i{ 0 }; i < point_cloud.nrow(); i++)
    {
        points.push_back(
            spatial::point_3d_t{ x_coords[i], y_coords[i], z_coords[i] }
        );
    }

    // Calculate modes for the points.
    std::vector< spatial::point_3d_t > modes{
        ams3d::calculate_modes(
            points,
            crown_diameter_2_tree_height,
            crown_height_2_tree_height,
            verbose, Rcpp::Rcout
        )
    };

    // Return the modes as an R data.frame.
    Rcpp::NumericVector mode_x_coords;
    Rcpp::NumericVector mode_y_coords;
    Rcpp::NumericVector mode_z_coords;

    for (const auto &mode : modes)
    {
        mode_x_coords.push_back(spatial::get_x(mode));
        mode_y_coords.push_back(spatial::get_y(mode));
        mode_z_coords.push_back(spatial::get_z(mode));
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("mode_x") = mode_x_coords,
        Rcpp::Named("mode_y") = mode_y_coords,
        Rcpp::Named("mode_z") = mode_z_coords
    );
}
