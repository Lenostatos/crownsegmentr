// This file is part of crownsegmentr, an R package for identifying tree crowns
// within 3D point clouds.
//
// Copyright (C) 2020 Leon Steinmeier, Nikolai Knapp, UFZ
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
