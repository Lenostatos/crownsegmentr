# This file is part of crownsegmentr, an R package for identifying tree crowns
# within 3D point clouds.
#
# Copyright (C) 2020 Leon Steinmeier
# Contact: Lenostatos@gmx.de
#
# crownsegmentr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# crownsegmentr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with crownsegmentr.  If not, see <http://www.gnu.org/licenses/>.

#' Segment Tree Crowns in a 3D Point Cloud
#'
#' Employs a variant of the mean shift algorithm and after that the DBSCAN
#' algorithm in order to find tree crowns in the \code{point_cloud}.
#'
#' @inheritParams calculate_modes
#'
#' @param point_cloud A data.table or data.frame within which the function
#'   searches for coordinate columns by first looking for columns named "x",
#'   "y", or "z" (or one of their upper-case variants) and then for the first
#'   numeric columns if it couldn't find one or more of those.
#' @param min_num_neighbors_per_core Integer Scalar. The minimum number of
#'   neighbors that a point needs to have in order to be considered as a core
#'   point by the DBSCAN clustering algorithm.
#' @param neighborhood_radius Numeric Scalar. The radius of the space around a
#'   point that is treated as the point's neighborhood.
#' @param return_modes Boolean Scalar. Should mode coordinates be returned as
#'   well?
#' @param return_centroids Boolean Scalar. Should centroids be returned as well?
#'   \emph{Warning:} Setting this to \code{TRUE} considerably slows down
#'   processing speed. This option is only intended to be used for parameter
#'   estimation and debugging purposes on small point clouds.
#'
#' @return A data.table with the three coordinate columns of \code{point_cloud}
#'   together with an additional cluster ID column called "crown_id".
#'
#'   If \code{return_modes} was set to TRUE, another three columns named
#'   "mode_x", "mode_y", and "mode_z" that hold the coordinates of the modes are
#'   included in the data.table.
#'
#'   If \code{return_centroids} was set to \code{TRUE} returns a list with two
#'   elements:
#'
#'   The first element (named "clustered_data") is the data.table that would
#'   have been returned if \code{return_centroids} was set to \code{FALSE}. The
#'   second one (named "centroid_data") is a list of data.frames that hold point
#'   coordinates with one data.frame for each mode in the "clustered_data".
#'   The first point in each data.frame is always one of the points from the
#'   input point cloud while the other coordinates represent the centroids that
#'   were calculated while searching for that point's mode.
#'   There is one additional column named "crown_id" in each data.frame that
#'   holds the crown ID of the mode that the centroids belong to.
#'
#' @export
segment_tree_crowns <- function(point_cloud,
                                crown_diameter_2_tree_height,
                                crown_height_2_tree_height,
                                min_num_neighbors_per_core,
                                neighborhood_radius,
                                verbose = TRUE,
                                return_modes = FALSE,
                                return_centroids = FALSE) {
  coordinates <- try_to_extract_coordinate_data(
    point_cloud,
    object_name = "point_cloud"
  )

  if (!return_centroids) {
    modes <- calculate_modes(
      coordinates,
      crown_diameter_2_tree_height, crown_height_2_tree_height,
      verbose
    )
  } else {
    modes_and_centroids <- calculate_modes_and_centroid_paths(
      coordinates,
      crown_diameter_2_tree_height, crown_height_2_tree_height,
      verbose
    )
    modes <- modes_and_centroids$mode_coords
  }

  if (verbose) cat("Start clustering...")

  crown_ids <- dbscan::dbscan(
    modes, neighborhood_radius, min_num_neighbors_per_core + 1
  )$cluster

  if (verbose) cat("Finished clustering.")

  if (return_modes) {
    res_clustered_data <-
      data.table::data.table(coordinates, modes, "crown_id" = crown_ids)
  } else {
    res_clustered_data <-
      data.table::data.table(coordinates, "crown_id" = crown_ids)
  }

  res_value <- res_clustered_data

  if (return_centroids) {
    centroid_paths <- lapply(
      modes_and_centroids$centroid_paths,
      FUN = function(centroid_path_coords) {
        centroid_path_coords <- cbind(
          centroid_path_coords[, 1:3],
          crown_id = res_clustered_data[
            centroid_path_coords$mode_index[1], "crown_id"
          ]
        )
      }
    )

    res_value <- list(
      "clustered_data" = res_clustered_data, "centroid_data" = centroid_paths
    )
  }

  res_value
}
