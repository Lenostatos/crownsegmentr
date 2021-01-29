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
