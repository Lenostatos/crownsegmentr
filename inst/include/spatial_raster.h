// This file is part of crownsegmentr, an R package for identifying tree crowns
// within 3D point clouds.
//
// Copyright (C) 2021 Leon Steinmeier, Nikolai Knapp, UFZ Leipzig
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

#ifndef SPATIAL_RASTER_H
#define SPATIAL_RASTER_H

#include "spatial_types.h"

namespace spatial
{
    /** Models a rectangular, non-rotated raster. */
    template< typename T_value >
    class Raster
    {
    private:
        /** The raster values go from top left to bottom right, row-wise. */
        const std::vector< T_value > _values;

        const std::size_t _num_rows;
        const std::size_t _num_cols;

        const coordinate_t _x_min;
        const coordinate_t _x_max;
        const coordinate_t _y_min;
        const coordinate_t _y_max;

        const coordinate_t _row_height;
        const coordinate_t _col_width;

    public:
        Raster (
            const std::vector< T_value > &values,
            std::size_t num_rows,
            std::size_t num_cols,
            coordinate_t x_min,
            coordinate_t x_max,
            coordinate_t y_min,
            coordinate_t y_max
        ):
        _values{ values },
        _num_rows{ num_rows },
        _num_cols{ num_cols },
        _x_min{ x_min },
        _x_max{ x_max },
        _y_min{ y_min },
        _y_max{ y_max },
        _row_height{ (_y_max - _y_min) / _num_rows },
        _col_width{ (_x_max - _x_min) / _num_cols }
        {}

        /** Indicates whether the raster has a value at the location x|y. */
        bool has_value_at( const coordinate_t x, const coordinate_t y ) const
        {
            return _x_min <= x && x <= _x_max && _y_min <= y && y <= _y_max;
        }

        /** \brief Returns the raster value at the location x|y.
         *
         *  Throws an exception if x or y are NaN or the location x|y lies
         *      outside of the raster's extent.
         */
        T_value value_at( const coordinate_t x, const coordinate_t y ) const
        {
            if (std::isnan( x ) || std::isnan( y ))
            {
                throw std::runtime_error (
                    "Tried to access raster value with NaN xy-coordinates."
                );
            }

            if (!this->has_value_at(x, y))
            {
                throw std::out_of_range (
                    "Tried to access raster value outside of raster extent."
                );
            }

            return no_throw_value_at(x, y);
        }

        /** Same as value_at but does not throw exceptions. For NaN coordinate
         *      values and locations outside the raster, the behavior of this
         *      method is undefined.
         */
        T_value no_throw_value_at (
            const coordinate_t x, const coordinate_t y
        ) const
        {
            std::size_t row_index {
                static_cast< std::size_t >( (_y_max - y) / _row_height )
            };
            // if y == _y_min, the row index is too big by one
            if (row_index == _num_rows) { --row_index; }

            std::size_t col_index {
                static_cast< std::size_t >( (x - _x_min) / _col_width )
            };
            // if x == _x_max, the column index is too big by one
            if (col_index == _num_cols) { --col_index; }

            return _values[_num_cols * row_index + col_index];
        }
    };
}

#endif  // define SPATIAL_RASTER_H
