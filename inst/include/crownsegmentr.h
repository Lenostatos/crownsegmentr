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

#ifndef RCPP_crownsegmentr
#define RCPP_crownsegmentr

//' \p The rationale for this file:
//' I want the C++ function(s) in the dbscan_Rcpp_interface header to be
//' available to users who import this package. There are no C++ functions that
//' I want Rcpp to generate a wrapper for.
//' The other header files in this include directory become available to the
//' respective source code by using the Makevars files with the
//' "PKG_CPPFLAGS =	-I../inst/include/" setting.

#include "ams3d_R_Cpp_interface.h"

#endif // RCPP_crownsegmentr
