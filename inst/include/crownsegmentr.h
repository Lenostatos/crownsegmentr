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
