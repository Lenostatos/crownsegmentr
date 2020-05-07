#ifndef AMS3D_R_INTERFACE
#define AMS3D_R_INTERFACE

#include <Rcpp.h>

Rcpp::DataFrame calculate_modes(
    const Rcpp::DataFrame &point_cloud,
    const double &crown_diameter_2_tree_height,
    const double &crown_height_2_tree_height,
    bool verbose = true
);

Rcpp::List calculate_modes_and_centroid_paths(
    const Rcpp::DataFrame &point_cloud,
    const double &crown_diameter_2_tree_height,
    const double &crown_height_2_tree_height,
    bool verbose = true
);

#endif  // define AMS3D_R_INTERFACE
