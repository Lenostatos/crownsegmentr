# This file is part of crownsegmentr, an R package for identifying tree crowns
# within 3D point clouds.
#
# Copyright (C) 2020 Leon Steinmeier, Nikolai Knapp, UFZ
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

#' Try to find three numeric columns for coordinate data.
#'
#' TODO
#'
#' @param input_object The object that is checked for being a
#'   \code{\link[base]{data.frame}} or a \code{\link[data.table]{data.table}}
#'   and is searched for xyz-coordinate data.
#'
#' @return A \code{\link[data.table]{data.table}} with just three columns that
#'   are expected to hold the x-, y-, and z-coordinates in that order.
try_to_extract_coordinate_data <- function(input_object) {

  # Assert that the input object is either a data.table or data.frame.
  assert_that(
    input_object %has_attr% "class",
    "data.frame" %in% attr(input_object, "class"),
    msg =
      "The coordinate data needs to be stored in a data.table or data.frame."
  )

  # Get the names of all numeric columns
  numeric_col_names <- names(which(sapply(input_object, is.numeric)))

  # Assert that the table consists of at least three numeric columns.
  assert_that(
    length(numeric_col_names) >= 3,
    msg = paste0(
      "The coordinates table needs to have at least three numeric columns for ",
      "x-, y-, and z-coordinates but there are only ",
      length(numeric_col_names), "."
    )
  )

  # Get the names of the first numeric columns whose names are a upper- or
  # lower-case variant of either "x", "y", or "z".
  xyz_chars <- c("x", "y", "z")
  input_col_names <- colnames(input_object)

  xyz_col_names <- input_col_names[match(xyz_chars, tolower(input_col_names))]
  names(xyz_col_names) <- xyz_chars

  numeric_non_xyz_col_names <- numeric_col_names[
    !(tolower(numeric_col_names) %in% xyz_chars)]

  res_col_names <- vector("character")

  for (xyz_col_name in xyz_col_names) {
    if (!is.na(xyz_col_name) && xyz_col_name %in% numeric_col_names) {
      # If there is a numeric column called either x, y, or z, take that
      res_col_names <- append(res_col_names, xyz_col_name)
    } else {
      # Take the first numeric column that is neither called x, y, or z and is
      # also not part of the already selected columns
      next_meaningful_col_name <- numeric_non_xyz_col_names[
        !(numeric_non_xyz_col_names %in% res_col_names)][1]
      warning(paste0(
        "Couldn't find a numeric column with name ", names(xyz_col_name),
        " or ", toupper(xyz_col_name), ". Using column ",
        next_meaningful_col_name, " for ", names(xyz_col_name),
        "-coordinates instead."
      ))
    }
  }

  if (is.data.table(input_object)) {
    return(input_object[, res_col_names, with = FALSE])
  } else {
    return(input_object[res_col_names])
  }
}
