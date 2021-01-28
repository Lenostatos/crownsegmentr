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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with crownsegmentr. If not, see <http://www.gnu.org/licenses/>.

#' Segment Tree Crowns in a 3D Point Cloud
#'
#' Employs a variant of the mean shift algorithm (Ferraz et. al, 2016) and after
#' that the DBSCAN algorithm in order to find tree crowns in a set of
#' xyz-coordinates.
#'
#' @param point_cloud A set of xyz-coordinates. Can be passed as either a
#'   \code{\link[base]{data.frame}}, a \code{\link[data.table]{data.table}}, or
#'   a lidR- \code{\link[lidR]{LAS}} or \code{\link[lidR]{LAScatalog-class}}
#'   object.
#'
#'   If it's a \code{\link[base]{data.frame}} or a
#'   \code{\link[data.table]{data.table}} the function searches for coordinate
#'   columns by first looking for numeric columns named "X", "Y", or "Z" (or one
#'   of their lower-case variants) and for each instance where it can't find one
#'   of those it selects the numeric column that comes next in the table's
#'   column order.
#'
#' @inheritParams calculate_modes
#' @param min_num_neighbors_per_core Integer Scalar. The minimum number of
#'   neighbors that a point needs to have in order to be considered as a core
#'   point by the DBSCAN clustering algorithm.
#' @param neighborhood_radius Numeric Scalar. The radius of the space around a
#'   point that is treated as the point's neighborhood.
#' @param id_attribute_name Character Scalar. The column name at which IDs for
#'   individual trees should be stored.
#' @param return_modes Boolean Scalar. Should mode coordinates be returned as
#'   well? This is usually used for parameter estimation and/or debugging
#'   purposes. If set to \code{TRUE} the returned point cloud will have three
#'   additional columns/attributes with mode coordinates named "mode_x",
#'   "mode_y", and "mode_z".
#' @param ... This is just a placeholder for the method-specific parameters
#'   listed below.
#'
#' @return A point cloud of the same type as passed to the function but with an
#'   additional column/attribute holding a crown ID for each point. IDs with the
#'   value zero indicate non-tree points.
#'
#'   The method for \code{\link[lidR]{LAScatalog-class}} objects is an
#'   exception. It works just like any other lidR function that accepts such
#'   objects, i.e. it returns either an in-memory \code{\link[lidR]{LAS}} object
#'   or writes the processed chunks to individual files and returns those file's
#'   names. Please refer to the documentation of the
#'   \code{\link[lidR]{LAScatalog-class}} for more details.
#'
#' @section Note: This function is not intended for point clouds with ground
#'   points.
#'
#' @section Warning: When the point cloud is a \code{\link[base]{data.frame}},
#'   it will be converted to a \code{\link[data.table]{data.table}} by reference
#'   at the beginning of this function's execution and converted back to a
#'   \code{\link[base]{data.frame}} again before the function returns.
#'
#' @references Ferraz, A., S. Saatchi, C. Mallet, and V. Meyer (2016)
#'   \emph{Lidar detection of individual tree size in tropical forests}. Remote
#'   Sensing of Environment 183:318–333.
#'   \url{https://doi.org/10.1016/j.rse.2016.05.028}.
#'
#' @export
methods::setGeneric("segment_tree_crowns",
  function(point_cloud,
           crown_diameter_2_tree_height,
           crown_height_2_tree_height,
           min_num_neighbors_per_core = 4,
           neighborhood_radius = 0.3,
           id_attribute_name = "crown_id",
           return_modes = FALSE,
           ...
           ) standardGeneric("segment_tree_crowns"),
  signature = "point_cloud"
)

#' @describeIn segment_tree_crowns Segments coordinates stored as three columns
#'   in a \code{\link[base]{data.frame}}.
#'
#' @param return_centroids Boolean Scalar. Should centroids be returned as well?
#'   \emph{Warning:} Setting this to \code{TRUE} considerably slows down
#'   processing speed. This option is only intended to be used with small point
#'   clouds for debugging purposes.
#'
#'   If \code{return_centroids} is set to \code{TRUE} a list with two
#'   elements is returned:
#'
#'   The first element (named "clustered_data") is the
#'   \code{\link[base]{data.frame}}/\code{\link[data.table]{data.table}} that
#'   would have been returned if \code{return_centroids} was set to
#'   \code{FALSE}. The second one (named "centroid_data") is a list of
#'   \code{\link[base]{data.frame}}s that hold point coordinates with one
#'   \code{\link[base]{data.frame}} for each mode in the "clustered_data". The
#'   first point in each \code{\link[base]{data.frame}} is always one of the
#'   points from the input point cloud while the other coordinates represent the
#'   centroids that were calculated while searching for that point's mode. There
#'   is one additional column named "crown_id" in each
#'   \code{\link[base]{data.frame}} that holds the crown ID of the mode that the
#'   centroids belong to.
methods::setMethod("segment_tree_crowns",
  signature(point_cloud = "data.frame"),
  function(point_cloud,
           crown_diameter_2_tree_height,
           crown_height_2_tree_height,
           min_num_neighbors_per_core = 4,
           neighborhood_radius = 0.3,
           id_attribute_name = "crown_id",
           return_modes = FALSE,
           verbose = TRUE,
           return_centroids = FALSE) {

    # TODO This feels uncool. Maybe there is a way without having side-effects?
    data.table::setDT(point_cloud)
    res <- segment_tree_crowns(point_cloud)
    data.table::setDF(point_cloud)

    return(res)
  }
)

#' @describeIn segment_tree_crowns Segments coordinates stored as three columns
#'   in a \code{\link[data.table]{data.table}}.
#'
#' @import assertthat
methods::setMethod("segment_tree_crowns",
  signature(point_cloud = "data.table"),
  function(point_cloud,
           crown_diameter_2_tree_height,
           crown_height_2_tree_height,
           min_num_neighbors_per_core = 4,
           neighborhood_radius = 0.3,
           id_attribute_name = "crown_id",
           return_modes = FALSE,
           verbose = TRUE,
           return_centroids = FALSE) {

    assert_that(is.number(crown_diameter_2_tree_height))
    assert_that(is.number(crown_height_2_tree_height))

    assert_that(is.count(min_num_neighbors_per_core))
    assert_that(min_num_neighbors_per_core >= 4)

    assert_that(is.number(neighborhood_radius))
    assert_that(neighborhood_radius > 0)

    assert_that(is.string(id_attribute_name))
    assert_that(!(id_attribute_name == ""))

    assert_that(is.flag(return_modes))
    assert_that(is.flag(verbose))
    assert_that(is.flag(return_centroids))

    coordinates <- try_to_extract_coordinate_data(point_cloud)

    # TODO give a warning like "Did you forget to normalize?" if height values
    # are unrealistically high

    if (!return_centroids) {
      modes <- calculate_modes(
        point_cloud = coordinates,
        crown_diameter_2_tree_height,
        crown_height_2_tree_height,
        verbose
      )
    } else {
      modes_and_centroids <- calculate_modes_and_centroid_paths(
        point_cloud = coordinates,
        crown_diameter_2_tree_height,
        crown_height_2_tree_height,
        verbose
      )
      modes <- modes_and_centroids$mode_coords
    }

    is_na_mode_row <-
      is.na(modes$mode_x) | is.na(modes$mode_y) | is.na(modes$mode_z)

    if (verbose) cat("Start clustering...")

    crown_ids <- vector(mode = "integer", length = nrow(modes))
    crown_ids[is_na_mode_row] <- NA_integer_

    crown_ids[!is_na_mode_row] <- dbscan::dbscan(
      modes[!is_na_mode_row, ],
      eps = neighborhood_radius,
      minPts = min_num_neighbors_per_core + 1
    )$cluster

    # I have observed that the dbscan algorithm sometimes identifies clusters
    # with less than minPts points. My guess is that a cluster can loose points
    # to neighboring clusters after having been initially identified with minPts
    # points.
    #
    # In any case, replace any crown_ids with less than minPts points with the
    # ID given to regular noise points
    num_points_per_id <- table(crown_ids)
    ids_with_less_than_minPts_points <- as.integer(names(
      num_points_per_id[
        which(num_points_per_id < min_num_neighbors_per_core + 1)]
    ))

    crown_ids[which(crown_ids %in% ids_with_less_than_minPts_points)] <- 0

    if (verbose) cat("Finished clustering.")

    # TODO I tried to use the user-defined column name for the crown IDs here,
    # but I got weird bugs. I would love to know how to do this with the
    # user-defined name but for now I will just change the column name at the
    # end of this function.
    res_clustered_data <- data.table::data.table(coordinates,
      "crown_id" = crown_ids
    )

    if (return_modes) {
      res_clustered_data <- cbind(res_clustered_data, modes)
    }

    if (!return_centroids) {

      col_names = names(res_clustered_data)
      col_names[which(col_names == "crown_id")] <- id_attribute_name
      names(res_clustered_data) <- col_names

      return(res_clustered_data)

    } else {

      centroid_paths <- lapply(
        modes_and_centroids$centroid_paths,
        FUN = function(centroid_path_coords) {
          centroid_path_coords <- cbind(
            centroid_path_coords[, 1:3],
            crown_id = res_clustered_data[
              centroid_path_coords$mode_index[1], crown_id]
          )
          col_names = names(centroid_path_coords)
          col_names[which(col_names == "crown_id")] <- id_attribute_name
          names(centroid_path_coords) <- col_names
        }
      )

      col_names = names(res_clustered_data)
      col_names[which(col_names == "crown_id")] <- id_attribute_name
      names(res_clustered_data) <- col_names

      return(list(
        "clustered_data" = res_clustered_data, "centroid_data" = centroid_paths
      ))
    }
  }
)

#' @describeIn segment_tree_crowns Segments the point cloud data of a
#'   \code{\link[lidR]{LAS}} object.
#'
#' @param make_id_attribute_file_writable Boolean Scalar. When writing the
#'   returned LAS object to disk, should the crown IDs (and if
#'   \code{return_modes} is set to \code{TRUE} also the mode coordinates) be
#'   written into that file as well?
#' @param id_attribute_file_description Character Scalar. If
#'   \code{make_id_attribute_file_writable} is \code{TRUE} then this will be
#'   used as an additional description of the ID values when the LAS object is
#'   written to disk.
#'
#' @importClassesFrom lidR LAS
methods::setMethod("segment_tree_crowns",
  signature(point_cloud = "LAS"),
  function(point_cloud,
           crown_diameter_2_tree_height,
           crown_height_2_tree_height,
           min_num_neighbors_per_core = 4,
           neighborhood_radius = 0.3,
           id_attribute_name = "crown_id",
           return_modes = FALSE,
           verbose = TRUE,
           make_id_attribute_file_writable = FALSE,
           id_attribute_file_description = id_attribute_name) {

    assert_that(is.flag(make_id_attribute_file_writable))
    assert_that(is.string(id_attribute_file_description))
    assert_that(!(id_attribute_file_description == ""))

    segmented_points <- segment_tree_crowns(
      point_cloud = point_cloud@data,
      crown_diameter_2_tree_height,
      crown_height_2_tree_height,
      min_num_neighbors_per_core,
      neighborhood_radius,
      id_attribute_name,
      return_modes,
      verbose
    )

    res_point_cloud <- NULL

    if (make_id_attribute_file_writable) {

      res_point_cloud <- lidR::add_lasattribute(point_cloud,
        x = segmented_points[[id_attribute_name]],
        name = id_attribute_name,
        desc = id_attribute_file_description
      )

      if (return_modes) {
        res_point_cloud <- lidR::add_lasattribute(res_point_cloud,
          x = segmented_points[["mode_x"]],
          name = "mode_x",
          desc = "mode_x"
        )
        res_point_cloud <- lidR::add_lasattribute(res_point_cloud,
          x = segmented_points[["mode_y"]],
          name = "mode_y",
          desc = "mode_y"
        )
        res_point_cloud <- lidR::add_lasattribute(res_point_cloud,
          x = segmented_points[["mode_z"]],
          name = "mode_z",
          desc = "mode_z"
        )
      }

    } else {

      res_point_cloud <- lidR::add_attribute(point_cloud,
        x = segmented_points[[id_attribute_name]],
        name = id_attribute_name
      )

      if (return_modes) {
        res_point_cloud <- lidR::add_attribute(res_point_cloud,
          x = segmented_points[["mode_x"]],
          name = "mode_x"
        )
        res_point_cloud <- lidR::add_attribute(res_point_cloud,
          x = segmented_points[["mode_y"]],
          name = "mode_y"
        )
        res_point_cloud <- lidR::add_attribute(res_point_cloud,
          x = segmented_points[["mode_z"]],
          name = "mode_z"
        )
      }
    }

    return(res_point_cloud)
  }
)

#' @describeIn segment_tree_crowns Segments the point cloud data of a
#'   \code{\link[lidR]{LAScatalog-class}} instance.
#'
#' @import data.table
#' @importClassesFrom lidR LAScatalog
methods::setMethod("segment_tree_crowns",
  signature(point_cloud = "LAScatalog"),
  function(point_cloud,
           crown_diameter_2_tree_height,
           crown_height_2_tree_height,
           min_num_neighbors_per_core = 4,
           neighborhood_radius = 0.3,
           id_attribute_name = "crown_id",
           return_modes = FALSE,
           # TODO toggle based on catalog options:
           make_id_attribute_file_writable = TRUE,
           id_attribute_file_description = id_attribute_name) {

    # TODO Throw an error if the scale and offset values are not the same for
    # all files

    select_apex_helper <- function(x, y, z) {
      max_z_index <- which.max(z)
      return(list(
        X = x[max_z_index],
        Y = y[max_z_index]
      ))
    }

    # TODO make this function an extra internal function
    catalog_function <- function(las, bbox, ...) {

      # segment the chunk
      segmented_las <- segment_tree_crowns(
        point_cloud = las,
        crown_diameter_2_tree_height,
        crown_height_2_tree_height,
        min_num_neighbors_per_core,
        neighborhood_radius,
        id_attribute_name,
        return_modes,
        verbose = FALSE,
        make_id_attribute_file_writable,
        id_attribute_file_description
      )

      # Get the apices in the core area
      core_apices <- segmented_las@data[
        segmented_las@data[[id_attribute_name]] != 0,
        select_apex_helper(X, Y, Z),
        by = id_attribute_name][
          # filter the apices that are not inside the buffer
          bbox@xmin <= X & X < bbox@xmax &
            bbox@ymin <= Y & Y < bbox@ymax, ]
      # -> I could use the buffer atrribute here instead but I have observed
      # that lidR sometimes produces duplicated trees on chunk borders and
      # that's what I am trying to prevent here by using the <= and <. This is
      # just a precaution. I haven't investigated whether the buffers are
      # actually part of the problem.

      # Get the tree points and any points with ID == 0 in the core area
      core_trees_n_points <- lidR::filter_poi(segmented_las,
        # i.e. either the point's ID belongs to one of the trees in the core
        # area...
        segmented_las@data[[id_attribute_name]] %in%
          core_apices[[id_attribute_name]] | (
            # ...or the ID is zero and the point is inside the core area. The
            # near equality test was taken from dplyr::near.
            abs(segmented_las@data[[id_attribute_name]] - 0) <
            .Machine$double.eps^0.5 &
              as.integer(buffer) == 0
            # TODO Instead of the buffer I could use the same filter here as
            # above but I guess a few duplicated points that don't belong to any
            # tree won't be a big problem anytime soon. One might try to account
            # for such artifacts after processing the catalog.
          )
      )

      # Calculate cross-chunk unique crown IDs
      # The following lines were adapted from the source code of the
      # lidR::find_trees function
      # (https://github.com/Jean-Romain/lidR/blob/14fadf2735b6dca425a9c65da65943eda19a4a25/R/find_trees.R#L105)
      x_offset <- las@header@PHB[["X offset"]]
      y_offset <- las@header@PHB[["Y offset"]]

      x_scale  <- las@header@PHB[["X scale factor"]]
      y_scale  <- las@header@PHB[["Y scale factor"]]

      bit_shift_ids <- core_apices[ , .(
        crown_id = core_apices[[id_attribute_name]],
        bit_shift_id =
          round((X - x_offset) / x_scale) * 2^32 +
          round((Y - y_offset) / y_scale)
      )]
      # Adaption end

      # Add an ID with the value zero
      bit_shift_ids <- rbind(
        bit_shift_ids,
        data.table::data.table(crown_id = 0, bit_shift_id = 0)
      )
      # Keep IDs with the value NA at that value
      bit_shift_ids[is.na(crown_id), bit_shift_id := NA]

      # Join the point cloud with the new IDs
      names(bit_shift_ids)[1] <- id_attribute_name

      joined_ids <- bit_shift_ids[
        core_trees_n_points@data, on = (id_attribute_name)][
          , -id_attribute_name, with = FALSE]$bit_shift_id

      core_trees_n_points@data[, (id_attribute_name) := joined_ids]

      core_trees_n_points@data[, buffer := NULL]

      return(core_trees_n_points)
    }

    return(lidR::catalog_apply(point_cloud,
      catalog_function,
      .options = list(need_buffer = TRUE, automerge = TRUE, autoread = TRUE)
    ))
  }
)
