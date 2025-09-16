# Remove small trees
#' remove small trees
#'
#' The function takes a point cloud in which trees were segmented, and removes
#' tree clusters that are smaller than a certain radius or a certain height
#'
#' @param las a point cloud as lidR LAS object, with a data column "crown_id"
#' @param min_radius (num) the threshold for crown radius, below which trees
#' will be removed
#' @param min_height (num) the threshold for crown height, below which trees
#' will be removed. Works only if las is normalized.
#' @return lidR LAS
#' @section Details
#' returns the same las object that was given as input, but with 
#' altered crown id's. Trees that are considered too small have their crown id
#' set to NA, and all other crown id's are re-assigned so that they are 
#' without gaps
#' @export 

remove_small_trees <- function(las, min_radius = 1, min_height = -Inf){
  # assumes that "las" is a las object with "crown_id" as a data column

  # create crown polgons
  crown.polygons <- crown_polygons_crown_id(las)

  # identify ids of small crowns
  small.ids <- crown.polygons[(crown.polygons$radius < min_radius)|
                                (crown.polygons$height < min_height),]$crown_id

  # set crown_id of small trees to NA
  las@data$crown_id[las@data$crown_id %in% small.ids] <- NA

  # create new ids in ascending order without gaps
  las@data$crown_id <- as.integer(as.factor(las@data$crown_id))

  return(las)

}




# Helper function

crown_polygons_crown_id <- function(las){
  # takes a normalized and segmented point cloud with a column "crown_id"
  # and returns a SpatVector object with the convex hull of each tree as polygon

  # create data frame with points
  points.df <- las@data[,c("X","Y", "Z", "crown_id")]
  # remove NA crown_id
  points.df <- points.df[!is.na(points.df$crown_id)]
  # convert into terra vector object
  points.vect <- vect(points.df, geom=c("X", "Y"))

  # calculate convex hull for every tree
  crowns <- terra::hull(points.vect, by = "crown_id")

  # calculate area
  crowns$area <- suppressWarnings( terra::expanse(crowns) )
  crowns$radius <- sqrt(crowns$area) / 2
  # add Z as the highest Z value in the cluster
  crowns <- merge(crowns, aggregate(Z ~ crown_id, points.df, max), by = "crown_id")
  # rename column
  names(crowns)[names(crowns)=="Z"] <- "height"

  return(crowns)
}

