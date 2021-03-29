# This file is part of crownsegmentr, an R package for identifying tree crowns
# within 3D point clouds.
#
# Copyright (C) 2021 Leon Steinmeier, Nikolai Knapp, UFZ
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
# along with crownsegmentr in a file called "COPYING". If not,
# see <http://www.gnu.org/licenses/>.

# Note: The definition of a valid coordinate table used by this function is
# relied on by the segment_tree_crowns_core function. Be aware of that when
# changing this function.
validate_coordinate_table <- function(coordinate_table) {

  assert_that(is.data.frame(coordinate_table), msg = paste0(
    "The coordinate data needs to be stored in a data.frame or another data ",
    "type which can be treated as one."
  ))

  # Get the names of all numeric columns
  numeric_col_names <- names(which(sapply(coordinate_table, is.numeric)))

  # Assert that the table consists of at least three numeric columns.
  assert_that(length(numeric_col_names) >= 3, msg = paste0(
    "The coordinate table needs to have at least three numeric columns for ",
    "x-, y-, and z-coordinates but there are only ", length(numeric_col_names),
    " numeric columns."
  ))
}

validate_crown_diameter_to_tree_height <-
  function(crown_diameter_to_tree_height) {

    assert_that(
      is.number(crown_diameter_to_tree_height),
      crown_diameter_to_tree_height > 0
    )

    if (crown_diameter_to_tree_height > 2) { warning(paste0(
      "A crown diameter to tree height ratio of ", crown_diameter_to_tree_height,
      " is likely to high."
    )) }
  }

validate_crown_height_to_tree_height <- function(crown_height_to_tree_height) {

  assert_that(
    is.number(crown_height_to_tree_height),
    crown_height_to_tree_height > 0
  )

  if (crown_height_to_tree_height > 2) { warning(paste0(
    "A crown height to tree height ratio of ", crown_height_to_tree_height,
    " is likely to high."
  )) }
}

validate_crown_id_column_name <- function(crown_id_column_name,
                                          coordinate_table) {

  # Check the data type and validity of the crown ID column name
  assert_that(
    is.string(crown_id_column_name),
    crown_id_column_name != ""
  )

  # Assert that crown_id_column_name is not already a column name of the
  # coordinate_table
  assert_that(
    !(make.names(crown_id_column_name) %in% colnames(coordinate_table)),
    msg = paste0(
      "The point cloud data already has a column/attribute with the name '",
      crown_id_column_name, "'. Please either choose a different argument for",
      " the crown_id_column_name parameter or modify the point cloud data."
    )
  )
}

validate_verbose <- function(verbose) {
  assert_that(
    is.flag(verbose),
    !is.na(verbose)
  )
}

validate_also_return_modes <- function(also_return_modes) {
  assert_that(
    is.flag(also_return_modes),
    !is.na(also_return_modes)
  )
}

validate_also_return_centroids <- function(also_return_centroids) {
  assert_that(
    is.flag(also_return_centroids),
    !is.na(also_return_centroids)
  )
}

validate_centroid_convergence_distance <-
  function(centroid_convergence_distance) {
    assert_that(
      is.number(centroid_convergence_distance),
      centroid_convergence_distance > 0
    )
  }

validate_max_num_centroids_per_mode <- function(max_num_centroids_per_mode) {
  assert_that(
    is.count(max_num_centroids_per_mode),
    max_num_centroids_per_mode >= 1
  )
}

validate_dbscan_neighborhood_radius <- function(dbscan_neighborhood_radius) {
  assert_that(
    is.number(dbscan_neighborhood_radius),
    dbscan_neighborhood_radius > 0
  )
}

validate_min_num_modes_per_neighborhood <-
  function(min_num_modes_per_neighborhood) {
    assert_that(
      is.count(min_num_modes_per_neighborhood),
      min_num_modes_per_neighborhood >= 1
    )
  }

validate_write_crown_id_also_to_file <- function(write_crown_id_also_to_file) {
  assert_that(
    is.flag(write_crown_id_also_to_file),
    !is.na(write_crown_id_also_to_file)
  )
}

#' Ensures that crown IDs are written to output files of the LAScatalog
#'
#' Issues a warning if the user wanted to write the output to files but not
#' store IDs of segmented bodies.
#'
#' @param write_crown_id_also_to_file The to-be-validated parameter.
#' @param LAScatalog The \code{\link[lidR]{LAScatalog-class}} object whose
#'   settings are compared to the value of \code{write_crown_id_also_to_file}.
#'
#' @return A possibly corrected value for write_crown_id_also_to_file.
validate_write_crown_id_also_to_file_for_LAScatalogs <-
  function(write_crown_id_also_to_file, LAScatalog) {

    validate_write_crown_id_also_to_file(write_crown_id_also_to_file)

    if (!write_crown_id_also_to_file &&
        lidR::opt_output_files(LAScatalog) != "") {
      warning("IDs of segmented bodies will be written to the output files. ",
              "(This warning was generated because you set the ",
              "write_crown_id_also_to_file parameter to FALSE but requested ",
              "output files instead of an in-memory object via the LAScatalog ",
              "options.)")
      return(TRUE)
    } else {
      return(write_crown_id_also_to_file)
    }
  }

validate_crown_id_file_description <- function(crown_id_file_description) {
  assert_that(
    is.string(crown_id_file_description),
    crown_id_file_description != ""
  )
}

#' Asserts that all files referenced by a LAScatalog have the same scale and
#' offset values.
#'
#' @param LAScatalog The \code{\link[lidR]{LAScatalog-class}} object to be
#'   tested.
validate_scale_n_offset_are_consistent <- function(LAScatalog) {

  scales_n_offsets <- collect_scale_n_offset_of_LAScatalog_files(
    LAScatalog
  )[, 1:6] # select only the values and not the last column with the file paths

  # Iterate over the scale and offset values for the x, y, and z coordinates one
  # after the other.
  for (column_name in colnames(scales_n_offsets)) {

    # Get the column
    values <- scales_n_offsets[[column_name]]

    # the following logic was inspired by this discussion:
    # https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-numeric-vector#

    # Sanity check for NA values
    assert_that(!anyNA(values),
      msg = paste0(
        "This shouldn't happen? Your LAScatalog appears to be referencing at ",
        "least one LAS file with their ", column_name, " set to NA."
    ))

    # Only do the following if there is more than one value
    if (length(values) > 1) {
      # Check that all the values are near-equal.
      assert_that(abs(max(values) - min(values)) < .Machine$double.eps ^ 0.5,
        msg = paste0(
          "The scale and offset values of all used LAS files have to be equal.",
          " Your LAScatalog appears to be referencing LAS files with differing",
          " ", column_name, "s. Call ",
          "crownsegmentr::collect_scale_n_offset_of_LAScatalog_files(<your LAScatalog object>)",
          " to get an overview of these values for all files referenced by ",
          "your LAScatalog."
      ))
    }
  }
}
