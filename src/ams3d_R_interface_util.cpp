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
        std::vector< spatial::point_3d_t > points;
        points.reserve( coordinate_table.nrow() );

        // Create point objects from the coordinate table.
        for (int i{ 0 }; i < coordinate_table.nrow(); i++)
        {
            points.push_back (
                spatial::point_3d_t{
                    x_coords[i],
                    y_coords[i],
                    z_coords[i]
                }
            );
        }

        return points;
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

    std::unique_ptr< spatial::I_Raster< double > >
    convert_list_argument_to_raster_double (
        const Rcpp::List &list
    ) {
        // Get the list elements' names.
        auto list_element_names{ Rcpp::as< std::vector< std::string > >( list.names() ) };

        // The following condition reads:
        // If any of the list elements' names are equal to "value".
        if (std::any_of (
            list_element_names.begin(),
            list_element_names.end(),
            [](const std::string &name) { return name == "value"; } ))
        {
            return std::make_unique< spatial::Single_value_pseudo_raster< double > > (
                Rcpp::as< double >( list["value"] )
            );
        }
        else
        {
            return std::make_unique< spatial::Raster< double > > (
                Rcpp::as< std::vector< double > >( list["values"] ),
                Rcpp::as< std::size_t >( list["num_rows"] ),
                Rcpp::as< std::size_t >( list["num_cols"] ),
                Rcpp::as< spatial::coordinate_t >( list["x_min"] ),
                Rcpp::as< spatial::coordinate_t >( list["x_max"] ),
                Rcpp::as< spatial::coordinate_t >( list["y_min"] ),
                Rcpp::as< spatial::coordinate_t >( list["y_max"] )
            );
        }
    }

    // TODO The following would be nice-to-have but doesn't work for some reason
//     template< typename T >
//     std::unique_ptr< spatial::I_Raster< T > > convert_list_argument_to_raster (
//         const Rcpp::List &list
//     ) {
//         // Get the list elements' names.
//         auto list_element_names{ Rcpp::as< std::vector< std::string > >( list.names() ) };
//
//         // The following condition reads:
//         // If any of the list elements' names are equal to "value".
//         if (std::any_of (
//             list_element_names.begin(),
//             list_element_names.end(),
//             [](const std::string &name) { return name == "value"; } ))
//         {
//             return std::make_unique< spatial::Single_value_pseudo_raster< T > > (
//                 Rcpp::as< T >( list["value"] )
//             );
//         }
//         else
//         {
//             return std::make_unique< spatial::Raster< T > > (
//                 Rcpp::as< std::vector< T > >( list["values"] ),
//                 Rcpp::as< std::size_t >( list["num_rows"] ),
//                 Rcpp::as< std::size_t >( list["num_cols"] ),
//                 Rcpp::as< spatial::coordinate_t >( list["x_min"] ),
//                 Rcpp::as< spatial::coordinate_t >( list["x_max"] ),
//                 Rcpp::as< spatial::coordinate_t >( list["y_min"] ),
//                 Rcpp::as< spatial::coordinate_t >( list["y_max"] )
//             );
//         }
//     }
}
