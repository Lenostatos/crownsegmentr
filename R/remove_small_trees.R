# Remove small trees
#' remove small trees
#'

#'
#'
#'

# Generic S4 Function -------------------------------------------------

#' Remove small clusters from segmented point cloud
#'
#' The function takes a point cloud in which trees were segmented, and removes
#' tree clusters that are smaller than a certain radius or a certain height
#'
#' @param point_cloud a point cloud, either as data.frame/data.table, or as
#' lidR::LAS object.
#' @param min_radius (num) the threshold for crown radius, below which trees
#' will be removed
#' @param min_height (num) the threshold for crown height, below which trees
#' will be removed. Works only if las is normalized.
#' @param crown_id_column_name the name of the column in which the id of the
#' crown is saved
#' @return lidR LAS
#' @section Details
#' returns the same las object that was given as input, but with
#' altered crown id's. Trees that are considered too small have their crown id
#' set to NA, and all other crown id's are re-assigned so that they are
#' without gaps
#' @export
methods::setGeneric("remove_small_trees",
                    function(point_cloud,
                             min_radius = 1,
                             min_height = -Inf,
                             crown_id_column_name = "crown_id") {
                      standardGeneric("remove_small_trees")
                    },
                    signature = "point_cloud"
)



# Method for data.frames and data.tables ----------------------------------
#' @describeIn remove_small_trees removes small tree clusters in a segmented
#'   [LAS object][lidR::LAS-class].
#'
#'
#' @importClassesFrom lidR LAS
methods::setMethod(
  "remove_small_trees",
  signature(point_cloud = "data.frame"),
  function(point_cloud,
           min_radius,
           min_height,
           crown_id_column_name){

    # TODO: test if it contains a data column of the required name

    # create crowns with crown_metrics
    metrics <- ~list(height = max(Z),
                     npoints = length(Z))
    crowns <- lidR::crown_metrics(lidR::LAS(point_cloud),
                                  metrics,
                                  attribute = crown_id_column_name,
                                  geom = "convex")
    # calculate radius
    crowns$area <- sf::st_area(crowns)
    crowns$radius <- sqrt(crowns$area) / 2


    # identify ids of crowns that are not too small
    right_size_ids <- crowns[(crowns$radius >= min_radius) &
                               (crowns$height >= min_height),
                             crown_id_column_name][[1]]

    # identify for which points the crown is not too small
    crown_is_not_small <- point_cloud[[crown_id_column_name]] %in% right_size_ids

    # set crown_id of small trees to NA
    point_cloud[[crown_id_column_name]][!crown_is_not_small] <- NA

    # create new ids in ascending order without gaps (ugly method, but works)
    point_cloud[[crown_id_column_name]] <-
      as.integer(as.factor(point_cloud[[crown_id_column_name]]))

    return(point_cloud)

  })



# Method for lidR::LAS Objects ----------------------------------------
#' @describeIn remove_small_trees removes small tree clusters in a segmented
#'   [LAS object][lidR::LAS-class].
#'
#'
#' @importClassesFrom lidR LAS
methods::setMethod(
  "remove_small_trees",
  signature(point_cloud = "LAS"),
  function(point_cloud,
           min_radius,
           min_height,
           crown_id_column_name){

  # TODO: add crown_id_column_name as variable
  # TODO: test if it contains a data column of the required name

  # create crowns with crown_metrics
  metrics <- ~list(height = max(Z),
                   npoints = length(Z))
  crowns <- lidR::crown_metrics(point_cloud,
                                metrics,
                                attribute = crown_id_column_name,
                                geom = "convex")
  # calculate radius
  crowns$area <- as.numeric(sf::st_area(crowns))
  crowns$radius <- sqrt(crowns$area) / 2

  # identify ids of crowns that are not too small
  right_size_ids <- crowns[(crowns$radius >= min_radius) &
                             (crowns$height >= min_height),
                           crown_id_column_name][[1]]

  # identify for which points the crown is not too small
  crown_is_not_small <- point_cloud@data[[crown_id_column_name]] %in% right_size_ids

  # set crown_id of small trees to NA
  point_cloud@data[[crown_id_column_name]][!crown_is_not_small] <- NA

  # create new ids in ascending order without gaps
  point_cloud@data[[crown_id_column_name]] <-
    as.integer(as.factor(point_cloud@data[[crown_id_column_name]]))

  return(point_cloud)

})


# TODO: Add case for LAS Catalog




