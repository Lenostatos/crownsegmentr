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
#' @param point_cloud a point cloud, either as data.frame/data.table, as
#' lidR::LAS object, or as LasCatalog. Needs to have a data column with the
#' name
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

    # copy point cloud to a reduced data frame that can be used for crown_polygons
    point_cloud_reduced_df <- point_cloud[,c("X", "Y", "Z")]
    point_cloud_reduced_df$crown_id <- point_cloud[[crown_id_column_name]]

    # create crown polgons
    crown.polygons <- crown_polygons_crown_id(point_cloud_reduced_df)

    # identify ids of small crowns
    small.ids <- crown.polygons[(crown.polygons$radius < min_radius)|
                                  (crown.polygons$height < min_height),]$crown_id

    # identify for which points the crown is too small
    crown_is_small <- point_cloud[[crown_id_column_name]] %in% small.ids

    # set crown_id of small trees to NA
    point_cloud[[crown_id_column_name]][crown_is_small] <- NA

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

  # TODO: test if it contains a data column of the required name

  # copy point cloud to a reduced data frame that can be used for crown_polygons
  point_cloud_reduced_df <- point_cloud@data[,c("X", "Y", "Z")]
  point_cloud_reduced_df$crown_id <- point_cloud@data[[crown_id_column_name]]

  # create crown polgons
  crown.polygons <- crown_polygons_crown_id(point_cloud_reduced_df)

  # identify ids of small crowns
  small.ids <- crown.polygons[(crown.polygons$radius < min_radius)|
                                (crown.polygons$height < min_height),]$crown_id

  # identify for which points the crown is too small
  crown_is_small <- point_cloud@data[[crown_id_column_name]] %in% small.ids

  # set crown_id of small trees to NA
  point_cloud@data[[crown_id_column_name]][crown_is_small] <- NA

  # create new ids in ascending order without gaps
  point_cloud@data[[crown_id_column_name]] <-
    as.integer(as.factor(point_cloud@data[[crown_id_column_name]]))

  return(point_cloud)

})


# TODO: Add case for LAS Catalog



# Helper function

crown_polygons_crown_id <- function(point_cloud){
  # takes a normalized and segmented point cloud with a column "crown_id"
  # and returns a SpatVector object with the convex hull of each tree as polygon

  # remove NA crown_id
  point_cloud <- point_cloud[!is.na(point_cloud$crown_id),]
  # convert into terra vector object
  points.vect <- terra::vect(point_cloud, geom=c("X", "Y"))

  # calculate convex hull for every tree
  crowns <- terra::hull(points.vect, by = "crown_id")

  # calculate area
  crowns$area <- suppressWarnings( terra::expanse(crowns) )
  crowns$radius <- sqrt(crowns$area) / 2
  # add Z as the highest Z value in the cluster
  crowns <- merge(crowns, aggregate(Z ~ crown_id, point_cloud, max), by = "crown_id")
  # rename column
  names(crowns)[names(crowns)=="Z"] <- "height"

  return(crowns)
}



