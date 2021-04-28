# This file is part of crownsegmentr, an R package for identifying tree crowns
# within 3D point clouds.
#
# Copyright (C) 2020-2021 Leon Steinmeier, Nikolai Knapp, UFZ Leipzig
# Contact: Leon.Steinmeier@posteo.net
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
# along with crownsegmentr in a file called "COPYING". If not,
# see <http://www.gnu.org/licenses/>.


# Generic S4 Function -------------------------------------------------

#' Segment Tree Crowns in a 3D Point Cloud
#'
#' Employs a variant of the mean shift algorithm (Ferraz et. al, 2016) and after
#' that the DBSCAN algorithm in order to identify tree crowns in airborne lidar
#' data.
#'
#' @param point_cloud A data set containing xyz-coordinates. Can be passed as
#'   either a \code{\link[base]{data.frame}}, a
#'   \code{\link[data.table]{data.table}}, or a lidR- \code{\link[lidR]{LAS}} or
#'   \code{\link[lidR]{LAScatalog-class}} object.
#'
#'   If it's a \code{\link[base]{data.frame}} or a
#'   \code{\link[data.table]{data.table}} the function searches for coordinate
#'   columns by looking for the first numeric columns named "x"/"X", "y"/"Y", or
#'   "z"/"Z". For each instance where it can't find one of those it selects the
#'   next available numeric column in the table and issues a warning.
#' @param crown_diameter_to_tree_height,crown_height_to_tree_height Single
#'   numbers or \code{\link[raster]{raster}} objects covering the area of the
#'   \code{point_cloud}. The values should be approximate crown diameter and
#'   crown height to tree height ratios of the trees expected to be found in the
#'   point cloud. Points will not be segmented wherever a raster contains
#'   \code{NA} values.
#' @param segment_crowns_only_above A single positive number denoting the
#'   minimum height above ground at which crown IDs will be calculated.
#'
#'   Note that points directly below this threshold will still be considered
#'   during the segmentation if they are within reach of search cylinders
#'   constructed at the \code{segment_crowns_only_above} height. See "How the
#'   algorithm works" to learn about the search cylinders.
#' @param ground_height One of
#'   \itemize{
#'     \item \code{NULL}, indicating that \code{point_cloud} is normalized with
#'     ground height at zero.
#'     \item A \code{\link[raster]{raster}} object providing ground heights for
#'     the area of the (not normalized) \code{point_cloud}.
#'     \item A list of (ideally named) arguments to the
#'     \code{\link[lidR]{grid_terrain}} function, which will be used to generate
#'     a ground height grid from \code{point_cloud}. Currently not supported
#'     with point clouds stored in \code{data.frame}s. The list should not
#'     contain an argument to the "las" parameter of grid_terrain.
#'   }
#'   Points will not be segmented wherever ground heights are NA.
#' @param crown_id_column_name A character string. The column or attribute name
#'   under which IDs for segmented bodies should be stored.
#' @param centroid_convergence_distance A single number. Distance at which it is
#'   assumed that subsequently calculated centroids have converged to the
#'   nearest mode. See "How the algorithm works" to learn about centroids and
#'   modes in the context of the AMS3D algorithm.
#' @param max_num_centroids_per_mode A single integer. Maximum number of
#'   centroids calculated before the search for the nearest mode is aborted. See
#'   "How the algorithm works" to learn about centroids and modes in the context
#'   of the AMS3D algorithm.
#' @param dbscan_neighborhood_radius A single number. Radius for the spherical
#'   DBSCAN neighborhood around a mode. See "How the algorithm works" to learn
#'   about neighborhoods in the context of the DBSCAN algorithm.
#' @param min_num_modes_per_neighborhood A single integer. The minimum number of
#'   modes within a DBSCAN neighborhood at which the mode in the neighborhood's
#'   center will be treated as a core point. See "How the algorithm works" to
#'   learn about neighborhoods and core points in the context of the DBSCAN
#'   algorithm.
#' @param ... Unused.
#'
#' @return The point cloud which was passed to the function but extended with a
#'   column/attribute holding for each point the ID of a segmented body. IDs
#'   with the value \code{NA} indicate that a point was not assigned to any
#'   body.
#'
#'   If \code{also_return_modes} and/or \code{also_return_centroids} were set to
#'   \code{TRUE}, a list with at most three named elements in the following
#'   order:
#'   \describe{
#'     \item{segmented_point_cloud}{
#'       The segmented point cloud which would have been returned directly if
#'       \code{also_return_modes} and \code{also_return_centroids} had been set
#'       to \code{FALSE}.
#'     }
#'     \item{modes}{
#'       If \code{also_return_modes} was set to \code{TRUE}, a point cloud of
#'       the same type as the input point cloud holding the modes calculated
#'       with the AMS3D algorithm and two additional columns/attributes. One of
#'       these columns/attributes holds IDs of the segmented bodies that the
#'       modes belong to and the other (named "point_index") holds indices to
#'       the points in the input point cloud.
#'     }
#'     \item{centroids}{
#'       If \code{also_return_centroids} was set to \code{TRUE}, a point cloud
#'       of the same type as the input point cloud holding the centroids
#'       calculated with the AMS3D algorithm and two additional
#'       columns/attributes. One of these columns/attributes holds IDs of the
#'       segmented bodies that the centroids belong to and the other (named
#'       "point_index") holds indices to the points in the input point cloud.
#'     }
#'   }
#'
#'   The method for \code{\link[lidR]{LAScatalog-class}} objects works just like
#'   any other lidR function that accepts LAScatalogs, i.e. it returns either an
#'   in-memory \code{\link[lidR]{LAS}} object or writes the processed chunks to
#'   individual files and returns those file's names. Please refer to the
#'   documentation of the \code{\link[lidR]{LAScatalog-class}} for more details.
#'
#' @section How the algorithm works:
#'
#'   The basic assumption is that tree crowns form local maxima of point density
#'   and height within lidar point clouds. These local maxima are called
#'   *modes*. The algorithm tries to find the nearest mode for each point. This
#'   is done by looking at the surrounding points and moving into the direction
#'   of the highest point density until the nearest mode is reached.
#'
#'   The surrounding points are found with a three-dimensional search window
#'   which is roughly shaped like a crown itself, i.e. it is a vertical cylinder
#'   whose diameter and height are similar to those of the surrounding tree
#'   crowns. The exact diameter and height are dependent on the cylinder's
#'   height above ground and are controlled with the two main parameters of the
#'   algorithm: \code{crown_diameter_to_tree_height} and
#'   \code{crown_height_to_tree_height}.
#'
#'   The direction of the highest point density is found by calculating the
#'   average position of all points within the cylinder, the cylinder's so
#'   called *centroid*. In order to move further into the direction of the
#'   highest point density, a new cylinder is placed on the centroid and a new
#'   centroid is calculated for that cylinder. This goes on until the cylinders
#'   "stop moving", i.e. until two subsequently calculated centroids are closer
#'   to each other than \code{centroid_convergence_distance}. At this point, the
#'   most recently calculated centroid, is taken as the original point's nearest
#'   mode.
#'
#'   It sometimes happens that centroids converge only after a lot of
#'   iterations. In order to prevent situations where an excessive number of
#'   centroids is calculated for just one point, the parameter
#'   \code{max_num_centroids_per_mode} is used to abort the centroid
#'   calculations after a certain number of them has been performed.
#'   Nonetheless, the last centroid found before the abortion is still taken as
#'   a good enough guess of the nearest mode's position.
#'
#'   After the modes of the individual points have been calculated, it can be
#'   seen that modes of the same tree crown are positioned very close to each
#'   other, shortly below their crown's apex. These dense clusters of modes are
#'   identified with the DBSCAN algorithm which assigns a cluster ID to every
#'   mode. The cluster IDs are then finally connected back to the points of the
#'   point cloud and used as crown IDs.
#'
#'   The DBSCAN clustering is explained nicely in
#'   \href{https://en.wikipedia.org/wiki/DBSCAN#Preliminary}{Wikipedia} but here
#'   is a quick sketch of what it does: The DBSCAN algorithm classifies points
#'   as either core points, border points, or noise and assigns core and border
#'   points to the same cluster if they are close enough to at least one other
#'   core point of the cluster.
#'
#'   In order to be core points, points need to have enough neighbors. The
#'   parameter \code{dbscan_neighborhood_radius} determines the radius of the
#'   neighborhood and the parameter \code{min_num_modes_per_neighborhood}
#'   determines the minimum number of points in the neighborhood (including the
#'   to-be-classified one), which are needed for a core point.
#'
#'   Border points are within the neighborhood of core points but don't have
#'   enough neighbors to be core points themselves. Noise points are not within
#'   the neighborhood of any core point and also don't have enough neighbors to
#'   be core points.
#'
#'   Clusters are identified by iterating over the points and classifying them
#'   one by one. For each point the neighborhood is scanned and the point is
#'   classified accordingly. If the point is a core or border point, the
#'   neighboring points are classified next. As long as it is possible to
#'   directly connect to new core or border points in this way, the same cluster
#'   ID is assigned to each encountered point.
#'
#' @references Ferraz, A., S. Saatchi, C. Mallet, and V. Meyer (2016)
#'   \emph{Lidar detection of individual tree size in tropical forests}. Remote
#'   Sensing of Environment 183:318–333.
#'   \url{https://doi.org/10.1016/j.rse.2016.05.028}.
#'
#' @example R/examples/segment_tree_crowns_examples.R
#'
#' @export
methods::setGeneric("segment_tree_crowns",
  function(point_cloud,
           crown_diameter_to_tree_height,
           crown_height_to_tree_height,
           segment_crowns_only_above = 0,
           ground_height = NULL,
           crown_id_column_name = "crown_id",
           centroid_convergence_distance = 0.01,
           max_num_centroids_per_mode = 500,
           dbscan_neighborhood_radius = 0.3,
           min_num_modes_per_neighborhood = 5,
           ...
           ) standardGeneric("segment_tree_crowns"),
  signature = "point_cloud"
)


# Method for data.frames and data.tables ----------------------------------

#' @describeIn segment_tree_crowns Segments coordinates stored as three columns
#'   in a \code{\link[base]{data.frame}} or
#'   \code{\link[data.table]{data.table}}.
#'
#' @param verbose \code{TRUE} or \code{FALSE}. Should the function show a
#'   progress bar and other runtime information in the console?
#' @param also_return_modes \code{TRUE} or \code{FALSE}. Should mode coordinates
#'   be returned as well?
#' @param also_return_centroids \code{TRUE} or \code{FALSE}. Should centroid
#'   coordinates be returned as well? This slows down processing by a little bit
#'   and will return a data set which requires at least ~10 times more memory
#'   than the input point cloud.
methods::setMethod("segment_tree_crowns",
  signature(point_cloud = c("data.frame")),
  function(point_cloud,
           crown_diameter_to_tree_height,
           crown_height_to_tree_height,
           segment_crowns_only_above,
           ground_height,
           crown_id_column_name,
           centroid_convergence_distance,
           max_num_centroids_per_mode,
           dbscan_neighborhood_radius,
           min_num_modes_per_neighborhood,
           verbose = TRUE,
           also_return_modes = FALSE,
           also_return_centroids = FALSE) {

    validate_coordinate_table(point_cloud)
    validate_crown_diameter_to_tree_height(crown_diameter_to_tree_height,
      point_cloud
    )
    validate_crown_height_to_tree_height(crown_height_to_tree_height,
      point_cloud
    )
    validate_segment_crowns_only_above(segment_crowns_only_above)
    validate_ground_height(ground_height, point_cloud)
    validate_crown_id_column_name(crown_id_column_name, point_cloud)
    validate_centroid_convergence_distance(centroid_convergence_distance)
    validate_max_num_centroids_per_mode(max_num_centroids_per_mode)
    validate_dbscan_neighborhood_radius(dbscan_neighborhood_radius)
    validate_min_num_modes_per_neighborhood(min_num_modes_per_neighborhood)
    validate_verbose(verbose)
    validate_also_return_modes(also_return_modes)
    validate_also_return_centroids(also_return_centroids)

    # Segment the point cloud
    ids_modes_n_centroids <- segment_tree_crowns_core(
      point_cloud,
      segment_crowns_only_above,
      ground_height,
      crown_diameter_to_tree_height,
      crown_height_to_tree_height,
      verbose,
      centroid_convergence_distance,
      max_num_centroids_per_mode,
      dbscan_neighborhood_radius,
      min_num_modes_per_neighborhood,
      also_return_modes,
      also_return_centroids
    )

    # Bind the IDs to the point cloud data
    if (data.table::is.data.table(point_cloud)) {
      res <- data.table::data.table(
        point_cloud,
        ids_modes_n_centroids$crown_ids
      )
    } else {
      res <- cbind(
        point_cloud,
        ids_modes_n_centroids$crown_ids
      )
    }
    # Rename the ID column with the requested name
    names(res)[ncol(res)] <- crown_id_column_name

    # Return the point cloud with crown IDs if no additional data was requested
    if (!also_return_modes && !also_return_centroids) {

      return(res)

    } else { # set up a list with the different result data sets

      res <- list(segmented_point_cloud = res)

      if (also_return_modes) {
        # Add the modes to the list
        if (data.table::is.data.table(point_cloud)) {
          res[["modes"]] <- ids_modes_n_centroids$mode_coordinates
        } else {
          res[["modes"]] <- as.data.frame(
            ids_modes_n_centroids$mode_coordinates
          )
        }
        # Rename the ID column with the requested name
        names(res$modes)[4] <- crown_id_column_name
      }

      if (also_return_centroids) {
        # Add the centroids to the list
        if (data.table::is.data.table(point_cloud)) {
          res[["centroids"]] <- ids_modes_n_centroids$centroid_coordinates
        } else {
          res[["centroids"]] <- as.data.frame(
            ids_modes_n_centroids$centroid_coordinates
          )
        }
        # Rename the ID column with the requested name
        names(res$centroids)[4] <- crown_id_column_name
      }

      return(res)
    }
  }
)

# Method for lidR::LAS Objects ----------------------------------------

#' @describeIn segment_tree_crowns Segments the point cloud data of a
#'   \code{\link[lidR]{LAS}} object.
#'
#' @param write_crown_id_also_to_file \code{TRUE} or \code{FALSE}. When writing
#'   the returned LAS object to disk, should the IDs of segmented bodies be
#'   written into that file as well? See \code{\link[lidR]{add_lasattribute}}
#'   for additional details. Will also be used for all attributes of the LAS
#'   object(s) which are returned if \code{also_return_modes} and/or
#'   \code{also_return_centroids} were set to \code{TRUE}.
#'
#'   For \code{\link[lidR]{LAScatalog-class}} objects, this is only used if the
#'   result is returned as a \code{\link[lidR]{LAS}} object in memory. If the
#'   LAScatalog is set up to write the segmented point clouds into files, the
#'   IDs of segmented bodies will always be written to these files as well.
#'
#' @param crown_id_file_description A character string. If
#'   \code{write_crown_id_also_to_file} is set to \code{TRUE} this will be used
#'   as an additional description of the IDs of segmented bodies when the LAS
#'   object is written to disk. See the "desc" parameter of
#'   \code{\link[lidR]{add_lasattribute}} for additional details.
#'
#' @importClassesFrom lidR LAS
methods::setMethod("segment_tree_crowns",
  signature(point_cloud = "LAS"),
  function(point_cloud,
           crown_diameter_to_tree_height,
           crown_height_to_tree_height,
           segment_crowns_only_above,
           ground_height,
           crown_id_column_name,
           centroid_convergence_distance,
           max_num_centroids_per_mode,
           dbscan_neighborhood_radius,
           min_num_modes_per_neighborhood,
           verbose = TRUE,
           also_return_modes = FALSE,
           also_return_centroids = FALSE,
           write_crown_id_also_to_file = FALSE,
           crown_id_file_description = crown_id_column_name) {

    validate_crown_diameter_to_tree_height(crown_diameter_to_tree_height,
      point_cloud
    )
    validate_crown_height_to_tree_height(crown_height_to_tree_height,
      point_cloud
    )
    validate_segment_crowns_only_above(segment_crowns_only_above)
    validate_ground_height(ground_height, point_cloud)
    validate_crown_id_column_name(crown_id_column_name, point_cloud@data)
    validate_centroid_convergence_distance(centroid_convergence_distance)
    validate_max_num_centroids_per_mode(max_num_centroids_per_mode)
    validate_dbscan_neighborhood_radius(dbscan_neighborhood_radius)
    validate_min_num_modes_per_neighborhood(min_num_modes_per_neighborhood)
    validate_verbose(verbose)
    validate_also_return_modes(also_return_modes)
    validate_also_return_centroids(also_return_centroids)
    validate_write_crown_id_also_to_file(write_crown_id_also_to_file)
    validate_crown_id_file_description(crown_id_file_description)

    # If ground_height is a list of arguments, pass them to lidR::grid_terrain
    # like this:
    if (is.list(ground_height)) {
      ground_height <- do.call(lidR::grid_terrain,
        args = c(las = point_cloud, ground_height)
      )
    }

    # Segment the point cloud
    ids_modes_n_centroids <- segment_tree_crowns_core(
      point_cloud@data,
      segment_crowns_only_above,
      ground_height,
      crown_diameter_to_tree_height,
      crown_height_to_tree_height,
      verbose,
      centroid_convergence_distance,
      max_num_centroids_per_mode,
      dbscan_neighborhood_radius,
      min_num_modes_per_neighborhood,
      also_return_modes,
      also_return_centroids
    )

    # Bind the IDs to the point cloud data using the requested name
    if (write_crown_id_also_to_file) { # IDs will be saved with LAS files
      res <- lidR::add_lasattribute(
        las = point_cloud,
        x = ids_modes_n_centroids$crown_ids,
        name = crown_id_column_name,
        desc = crown_id_file_description
      )
    } else { # IDs will only be available in memory
      res <- lidR::add_attribute(
        las = point_cloud,
        x = ids_modes_n_centroids$crown_ids,
        name = crown_id_column_name
      )
    }

    # Return the point cloud with crown IDs if no additional data was requested
    if (!also_return_modes && !also_return_centroids) {

      return(res)

    } else { # set up a list with the different result data sets

      res <- list(segmented_point_cloud = res)

      if (also_return_modes) {
        # Create a LAS object storing the mode coordinates
        res[["modes"]] <- suppressWarnings(
          lidR::las_update(
            lidR::las_quantize(
              lidR::LAS(
                ids_modes_n_centroids$mode_coordinates[
                  , .(X = x, Y = y, Z = z)],
                header = point_cloud@header
        ))))

        # Add the crown IDs and point indices to the modes
        res$modes <- lidR::add_attribute(
          las = lidR::add_attribute(
            las = res$modes,
            x = ids_modes_n_centroids$mode_coordinates$crown_id,
            name = crown_id_column_name
          ),
          x = ids_modes_n_centroids$mode_coordinates$point_index,
          name = "point_index"
        )

        # If requested for the crown IDs, make the crown ID and point index
        # attributes of the modes LAS "file-writable" as well
        if (write_crown_id_also_to_file) {
          res$modes <- lidR::add_lasattribute(
            las = lidR::add_lasattribute(
              las = res$modes,
              name = crown_id_column_name,
              desc = crown_id_file_description
            ),
            name = "point_index",
            desc = "Index of original points"
          )
        }

      } # end if also_return_modes

      if (also_return_centroids) {
        # Create a LAS object storing the centroid coordinates
        res[["centroids"]] <- suppressWarnings(
          lidR::las_update(
            lidR::las_quantize(
              lidR::LAS(
                ids_modes_n_centroids$centroid_coordinates[
                  , .(X = x, Y = y, Z = z)],
                header = point_cloud@header
        ))))

        # Add the crown IDs and point indices to the centroids
        res$centroids <- lidR::add_attribute(
          las = lidR::add_attribute(
            las = res$centroids,
            x = ids_modes_n_centroids$centroid_coordinates$crown_id,
            name = crown_id_column_name
          ),
          x = ids_modes_n_centroids$centroid_coordinates$point_index,
          name = "point_index"
        )

        # If requested for the crown IDs, make the crown ID and point index
        # attributes of the centroids LAS "file-writable" as well
        if (write_crown_id_also_to_file) {
          res$centroids <- lidR::add_lasattribute(
            las = lidR::add_lasattribute(
              las = res$centroids,
              name = crown_id_column_name,
              desc = crown_id_file_description
            ),
            name = "point_index",
            desc = "Index of original points"
          )
        }

      } # end if also_return_centroids

      return(res)

    } # end returning mode and/or centroid coordinates
  }
)


# Method for lidR::LAScatalog Objects ---------------------------------

#' @describeIn segment_tree_crowns Segments the point cloud data of a
#'   \code{\link[lidR]{LAScatalog-class}} object. This method does not support
#'   additionally returning modes and/or centroids. Instead of the verbose
#'   parameter use the LAScatalog's progress option (see the LAScatalog
#'   documentation -> "Processing options" -> "progress").
#'
#' @import data.table
#' @importClassesFrom lidR LAScatalog
methods::setMethod("segment_tree_crowns",
  signature(point_cloud = "LAScatalog"),
  function(point_cloud,
           crown_diameter_to_tree_height,
           crown_height_to_tree_height,
           segment_crowns_only_above,
           ground_height,
           crown_id_column_name,
           centroid_convergence_distance,
           max_num_centroids_per_mode,
           dbscan_neighborhood_radius,
           min_num_modes_per_neighborhood,
           write_crown_id_also_to_file = TRUE,
           crown_id_file_description = crown_id_column_name) {

    validate_scale_n_offset_are_consistent(point_cloud)

    validate_crown_diameter_to_tree_height(crown_diameter_to_tree_height,
      point_cloud
    )
    validate_crown_height_to_tree_height(crown_height_to_tree_height,
      point_cloud
    )
    # validate_crown_id_column_name(crown_id_column_name, point_cloud@data)
    # TODO find out whether there is a way to validate this here instead of in
    # the individual calls to the LAS methods.
    validate_segment_crowns_only_above(segment_crowns_only_above)
    validate_ground_height(ground_height, point_cloud)
    validate_centroid_convergence_distance(centroid_convergence_distance)
    validate_max_num_centroids_per_mode(max_num_centroids_per_mode)
    validate_dbscan_neighborhood_radius(dbscan_neighborhood_radius)
    validate_min_num_modes_per_neighborhood(min_num_modes_per_neighborhood)
    write_crown_id_also_to_file <-
      validate_write_crown_id_also_to_file_for_LAScatalogs(
        write_crown_id_also_to_file,
        point_cloud
      )
    validate_crown_id_file_description(crown_id_file_description)


    # Used for selecting the x and y coordinates of the highest point of a tree
    # crown point cloud. Called in the catalog_function below.
    select_apex_helper <- function(x, y, z) {
      max_z_index <- which.max(z)
      return(list(
        X = x[max_z_index],
        Y = y[max_z_index]
      ))
    }

    # This function is forwarded to the lidR::LAScatalog framework by which it
    # is used to process the individual chunks of the input LAScatalog.
    catalog_function <- function(las, bbox, ...) {

      # segment the chunk with the LAS method
      segmented_las <- segment_tree_crowns(
        point_cloud = las,
        crown_diameter_to_tree_height,
        crown_height_to_tree_height,
        segment_crowns_only_above,
        ground_height,
        crown_id_column_name,
        centroid_convergence_distance,
        max_num_centroids_per_mode,
        dbscan_neighborhood_radius,
        min_num_modes_per_neighborhood,
        verbose = FALSE,
        also_return_modes = FALSE,
        also_return_centroids = FALSE,
        write_crown_id_also_to_file,
        crown_id_file_description
      )

      # Get the apex xy-coordinates and ID of each tree crown
      # Then select only those crowns whose apex lies in the core area
      core_apices <- segmented_las@data[
        !is.na(segmented_las@data[[crown_id_column_name]]),
        select_apex_helper(X, Y, Z),
        by = crown_id_column_name][
          # filter the apices that are inside the core area
          bbox@xmin <= X & X < bbox@xmax &
            bbox@ymin <= Y & Y < bbox@ymax, ]
      # -> I could use the buffer atrribute here instead but I have observed
      # that lidR sometimes produces duplicated trees on chunk borders and
      # that's what I am trying to prevent here by using the <= and <.
      # TODO Test whether this actually makes a difference.

      # Get the trees and any points with ID == NA in the core area
      # I.e. in the end, trees at the chunk edge will extend into the buffer
      # while unsegmented points are only selected inside the core area.
      core_trees_n_points <- lidR::filter_poi(segmented_las,
        # i.e. either the point's ID belongs to one of the trees in the core
        # area...
        segmented_las@data[[crown_id_column_name]] %in%
          core_apices[[crown_id_column_name]] |
          # ...or the ID is NA and the point is inside the core area.
          (
            is.na(segmented_las@data[[crown_id_column_name]]) &
              as.integer(buffer) == 0
          )
          # TODO Instead of the buffer I could use the same filter here as with
          # the apices above but I guess a few duplicated points that don't
          # belong to any tree won't be a big problem anytime soon. One might
          # try to account for such artifacts after processing the catalog.
      )

      # Calculate cross-chunk unique crown IDs
      # The following lines were adapted from the source code of the
      # lidR::find_trees function
      # (https://github.com/Jean-Romain/lidR/blob/14fadf2735b6dca425a9c65da65943eda19a4a25/R/find_trees.R#L105)
      x_offset <- las@header@PHB[["X offset"]]
      y_offset <- las@header@PHB[["Y offset"]]

      x_scale  <- las@header@PHB[["X scale factor"]]
      y_scale  <- las@header@PHB[["Y scale factor"]]

      # Connect the IDs of the core apices with the their respective bitshift
      # IDs
      bit_shift_ids <- core_apices[, .(
        crown_id = core_apices[[crown_id_column_name]],
        bit_shift_id =
          round((X - x_offset) / x_scale) * 2^32 +
          round((Y - y_offset) / y_scale)
      )]
      # Adaption end

      # The ID of unsegmented points is missing from the core apices so let's
      # add it here
      bit_shift_ids <- rbind(
        bit_shift_ids,
        data.table::data.table(crown_id = NA, bit_shift_id = NA)
      )
      # Keep IDs with the value NA at that value
      # -> I think this is unnecessary
      # bit_shift_ids[is.na(crown_id), bit_shift_id := NA]

      # Give the crown_id column in the bit_shift_ids data.table the same name
      # as the crown IDs in the point cloud data
      names(bit_shift_ids)[1] <- crown_id_column_name

      # Join the bitshift IDs to the individual points based on their crown IDs
      joined_ids <- bit_shift_ids[
        core_trees_n_points@data[, ..crown_id_column_name],
        on = (crown_id_column_name)
      ]$bit_shift_id
      # Code before:
      # core_trees_n_points@data, on = (crown_id_column_name)][
      #   , -crown_id_column_name, with = FALSE]$bit_shift_id

      # Overwrite the "normal" crown IDs with the bitshift ones
      core_trees_n_points@data[, (crown_id_column_name) := joined_ids]

      # Remove the buffer column
      core_trees_n_points@data[, buffer := NULL]

      return(core_trees_n_points)
    }

    return(lidR::catalog_apply(
      ctg = point_cloud,
      FUN = catalog_function,
      .options = list(need_buffer = TRUE, automerge = TRUE, autoread = TRUE)
    ))
  }
)


#' Get the scale and offset values of all files referenced by a LAScatalog
#'
#' @return A data.table with 7 columns. The first three columns hold the x, y,
#'   and z scale factors, the next three columns hold the x, y, and z offsets
#'   and the last column holds the file paths. There is one row for each
#'   referenced file.
#'
#' @keywords internal
#'
#' @export
collect_scale_n_offset_of_LAScatalog_files <- function(LAScatalog) {

  # Get the paths to the files accessed by the LAScatalog
  file_paths <- unique(unlist(lapply(
    X = lidR::catalog_makechunks(LAScatalog, plot = FALSE),
    FUN = function(chunk) {
      chunk@files
    }
  )))

  # Read the scale and offset data from the file headers and convert them into a
  # table
  return(data.table::rbindlist(lapply(
    X = file_paths,
    FUN = function(file_path) {
      c(
        lidR::readLASheader(file_path)@PHB[c(
          "X scale factor",
          "Y scale factor",
          "Z scale factor",
          "X offset",
          "Y offset",
          "Z offset"
        )],
        file_path = file_path
      )
    }
  )))
}
