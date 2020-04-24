#ifndef AMS3D_R_CPP_INTERFACE
#define AMS3D_R_CPP_INTERFACE

/** \file
 *  This header exposes C++ functions to users who get this code via an R
 *  package.
 * 
 *  \todo Read up on stuff at
 *      <a href="https://stackoverflow.com/questions/42261822/c-interface-with-rcppinterfaces-not-working-for-a-function-returning-stdpa">
 *          this stackoverflow question
 *      </a>
 *      in order to check whether this is the correct approach. Also just test
 *      it.
 */

#include "spatial.h"

/** \brief Calculates the mode for each point in \p point_cloud.
 * 
 *  \return An array of modes with one for each point in \p point_cloud.
 *      The modes are in the same order as the points they belong to.
 */
std::vector< spatial::point_3d_t > calculate_modes(
    const std::vector< spatial::point_3d_t > &point_cloud,
    const double crown_diameter_2_tree_height,
    const double crown_height_2_tree_height,
    bool verbose = true
);

#endif  // define AMS3D_R_CPP_INTERFACE
