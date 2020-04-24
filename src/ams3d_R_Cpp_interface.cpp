#include "ams3d_R_Cpp_interface.h"

#include "spatial.h"
#include "ams3d.h"

#include <Rcpp.h>

std::vector< spatial::point_3d_t > calculate_modes(
    const std::vector< spatial::point_3d_t > &point_cloud,
    const double crown_diameter_2_tree_height,
    const double crown_height_2_tree_height,
    bool verbose
) {
    return ams3d::calculate_modes(
        point_cloud,
        crown_diameter_2_tree_height,
        crown_height_2_tree_height,
        verbose, Rcpp::Rcout
    );
}
