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

#' Try to find three numeric columns for coordinate data.
#'
#' TODO
#'
#' @param input_object The object that is checked for being a data.frame with
#'   coordinate data.
#' @param object_name The name of the object to be used in warning and error
#'   messages.
#'
#' @section TODO:
#'   Infer \code{object_name} from the parent environment.
#'
#' @return A data.frame with just three columns that are expected to hold the
#'   x-, y-, and z-coordinates in that order.
try_to_extract_coordinate_data <- function(input_object, object_name) {
  # Assert that the input object is either a data.table or data.frame.
  assertthat::assert_that(
    assertthat::has_attr(input_object, "class"),
    "data.frame" %in% attr(input_object, "class"),
    msg = paste0(object_name, " needs to be a data.table or data.frame.")
  )

  # Assert that the data.frame has at least three columns.
  assertthat::assert_that(
    ncol(input_object) >= 3,
    msg = paste0(
      object_name, " needs to have at least three columns for x-, y-, and ",
      "z-coordinates."
    )
  )

  # Get the indices of the first columns whose names are a upper- or lower-case
  # variant of "x", "y", and "z".
  x_y_and_z <- c("x", "y", "z")
  lower_case_column_names <- tolower(colnames(input_object))
  xyz_column_indices <- match(x_y_and_z, lower_case_column_names)
  names(xyz_column_indices) <- x_y_and_z

  # Remove NA indices.
  xyz_column_indices <- xyz_column_indices[!is.na(xyz_column_indices)]

  # Check any found columns for being numeric.
  for (index in xyz_column_indices) {

    # If they are not, remove them with a warning.
    if (!is.numeric(input_object[, index])) {
      xyz_column_indices <- xyz_column_indices[xyz_column_indices != index]
      warning(paste0(
        "Found column ", names(input_object[, index]),
        " but won't use it because it's not numeric."
      ))
    }
  }

  # If any of x, y, or z are missing...
  if (length(xyz_column_indices) < 3) {

    # ...get all column indices that are not among the already found indices...
    all_indices <- seq_len(ncol(input_object))
    non_xyz_column_indices <- all_indices[
      !(all_indices %in% xyz_column_indices)
    ]

    # ...and get the remaining numeric columns among them.
    remaining_numeric_column_indices <- vector(mode = "integer")
    for (index in non_xyz_column_indices) {
      if (is.numeric(input_object[, index])) {
        remaining_numeric_column_indices <-
          append(remaining_numeric_column_indices, index)
      }
    }

    num_numeric_columns <-
      length(xyz_column_indices) + length(remaining_numeric_column_indices)
    assertthat::assert_that(
      num_numeric_columns >= 3,
      msg = paste0(
        "Need at least three numeric columns for coordinates, but found only ",
        num_numeric_columns, " numeric columns in ", object_name, "."
      )
    )

    # Now fill up any missing x, y, or z column indices with other numeric
    # columns' indices.
    for (one_of_xyz in x_y_and_z) {
      if (is.na(xyz_column_indices[one_of_xyz])) {
        xyz_column_indices[one_of_xyz] <- remaining_numeric_column_indices[1]
        remaining_numeric_column_indices <-
          remaining_numeric_column_indices[-1]
        warning(paste0(
          "Could not find any numeric column named ", one_of_xyz,
          ". Using column ", names(input_object)[xyz_column_indices[one_of_xyz]], " instead."
        ))
      }
    }
  }

  input_object[, xyz_column_indices]
}
