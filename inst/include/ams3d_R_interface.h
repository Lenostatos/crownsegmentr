#ifndef AMS3D_R_INTERFACE
#define AMS3D_R_INTERFACE

#include <Rcpp.h>

Rcpp::DataFrame calculate_modes(
    const Rcpp::DataFrame &point_cloud,
    const double &crown_diameter_2_tree_height,
    const double &crown_height_2_tree_height,
    bool verbose
);

#endif  // define AMS3D_R_INTERFACE
