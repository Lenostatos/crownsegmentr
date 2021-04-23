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

namespace ams3d_R_interface_util
{
    std::vector< spatial::point_3d_t > create_point_objects_from (
        const Rcpp::DataFrame &coordinate_table
    ) {
        // Get the coordinates as vectors because Rcpp::DataFrames are not
        // directly subsettable down to individual values.
        Rcpp::NumericVector x_coords( coordinate_table[0] );
        Rcpp::NumericVector y_coords( coordinate_table[1] );
        Rcpp::NumericVector z_coords( coordinate_table[2] );

        // Set up an array for point objects.
        std::vector< spatial::point_3d_t > point_cloud;
        point_cloud.reserve( coordinate_table.nrow() );

        // Create point objects from the coordinate table.
        for (int i{ 0 }; i < coordinate_table.nrow(); i++)
        {
            point_cloud.push_back (
                spatial::point_3d_t{
                    x_coords[i],
                    y_coords[i],
                    z_coords[i]
                }
            );
        }

        return point_cloud;
    }

    RProgress::RProgress create_progress_bar( double total )
    {
        // Constructor usage inferred from source at:
        // https://github.com/r-lib/progress/blob/master/inst/include/RProgress.h
        return RProgress::RProgress {
            // format:
            "  Calculating modes [:bar] :percent%  eta::eta  elapsed::elapsed",
            total, // total
            Rf_GetOptionWidth() - 2, // width
            '=', // complete_char
            '-', // incomplete_char
            false // clear
        };
    }
}
