#' Segment Tree Crowns in a 3D Point Cloud
#'
#' Employs a variant of the mean shift algorithm and after that the DBSCAN
#' algorithm in order to find tree crowns in the \code{point_cloud}.
#'
#' @inheritParams calculate_modes
#'
#' @param min_num_neighbors_per_core Integer Scalar. The minimum number of
#'   neighbors that a point needs to have in order to be considered as a core
#'   point by the DBSCAN clustering algorithm.
#' @param neighborhood_radius Numeric Scalar. The radius of the space around a
#'   point that is treated as the point's neighborhood.
#' @param return_modes Boolean Scalar. Should mode coordinates be returned as
#'   well?
#'
#' @return A data.table with the three coordinate columns of \code{point_cloud}
#'   together with an additional cluster ID column called "cluster_id". If
#'   \code{return_modes} is TRUE, another three columns with the coordinates of
#'   modes are included.
#'
#' @export
segment_tree_crowns <- function(
  point_cloud,
  crown_diameter_2_tree_height, crown_height_2_tree_height,
  min_num_neighbors_per_core, neighborhood_radius,
  verbose = TRUE, return_modes = FALSE
) {
  coordinates <- try_to_extract_coordinate_data(
    point_cloud, object_name = "point_cloud"
  )

  modes <- calculate_modes(
    coordinates,
    crown_diameter_2_tree_height, crown_height_2_tree_height,
    verbose
  )

  cluster_ids <- dbscan::dbscan(
    modes, neighborhood_radius, min_num_neighbors_per_core + 1
  )

  if (return_modes) {
    return(data.table::data.table(coordinates, modes, "cluster_id" = cluster_ids))
  } else {
    return(data.table::data.table(coordinates, "cluster_id" = cluster_ids))
  }
}
