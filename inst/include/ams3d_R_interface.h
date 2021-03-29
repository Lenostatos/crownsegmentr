// This file is part of crownsegmentr, an R package for identifying tree crowns
// within 3D point clouds.
//
// Copyright (C) 2020-2021 Leon Steinmeier, Nikolai Knapp, UFZ
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

#ifndef AMS3D_R_INTERFACE_H
#define AMS3D_R_INTERFACE_H

#include "spatial.h"
#include "ams3d.h"

#include <Rcpp.h>

#include <RProgress.h> // This line needs to go after "#include <Rcpp.h>"!
// (RProgress seems to rely on Rinternal.h which is probably loaded implicitly
// by "#include <Rcpp.h>".)

Rcpp::DataFrame calculate_modes (
    const Rcpp::DataFrame &coordinate_table,
    const double &crown_diameter_to_tree_height,
    const double &crown_height_to_tree_height,
    const spatial::distance_t &centroid_convergence_distance,
    const int max_num_centroids_per_mode,
    const bool show_progress_bar = true
);

Rcpp::List calculate_modes_plus_centroids (
    const Rcpp::DataFrame &coordinate_table,
    const double &crown_diameter_to_tree_height,
    const double &crown_height_to_tree_height,
    const spatial::distance_t &centroid_convergence_distance,
    const int max_num_centroids_per_mode,
    const bool show_progress_bar = true
);

namespace ams3d_R_interface_helper
{
    std::vector< spatial::point_3d_t > create_point_objects_from (
        const Rcpp::DataFrame &coordinate_table
    );

    RProgress::RProgress create_progress_bar( double total );
}

#endif  // define AMS3D_R_INTERFACE_H
