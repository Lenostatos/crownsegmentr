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

#ifndef CROWNSEGMENTR_H
#define CROWNSEGMENTR_H

//' \p Rationale for this file:
//' I wanted to provide a direct C++ interface for other R packages. This could
//' also be done with the [[Rcpp::export]] and [[Rcpp::interfaces(cpp)]]
//' attributes. However, I wanted a cleanly separated R-C++ interface so I
//' created the ams3d_R_Cpp_interface files and included the header here. This
//' should have the same effect as using the Rcpp attributes.

#include "ams3d_R_Cpp_interface.h"

#endif // CROWNSEGMENTR_H
