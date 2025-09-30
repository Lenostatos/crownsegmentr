# This file is part of crownsegmentr, an R package for identifying tree crowns
# within 3D point clouds.
#
# Copyright (C) 2021 Leon Steinmeier, Nikolai Knapp, UFZ Leipzig
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
#
# ----------------------------------------------------------
# watershed kds raster
#' calculate a raster of kernel diameter slope for AMS3D
#'
#' The function calculates a raster with values for the kernel diameter slope
#' as input for the AMS3D algorithm. It segments the tree crowns with
#' the watershed algorithm, calculates a ratio of crown diameter to tree height
#' for each tree, and converts this into a raster
#'
#' @param las a point cloud as lidR LAS object
#' @param bandwidth_intercept a fixed value for kernel_diameter_intercept, which
#' reduces the crown diameters by the given value before calculating the ratio
#' of crown diameter to tree height
#' @param limits a numeric vector with minimum and maximum allowed values for
#' the ratio
#' @param ground_height (optional) either
#' *NULL, indicating that the point cloud is normalized, or
#' *a [SpatRaster][terra::SpatRaster] digital terrain model, or
#' *a list of arguments to the
#' [lidR rasterize_terrain()][lidR::rasterize_terrain()] function to normalize
#' the point cloud.
#' @return a terra SpatRaster
#'
#' @section Details
#'
#' The output raster can serve as input for the parameter
#' "kernel_diameter_slope" (hence "kds") for the function segment_tree_crowns.
#' It averages the ratio of crown diameter to tree height for a 5 m radius,
#' for trees that were detected with watershed segmentation.
#'
#' @export

watershed_kds_raster <- function(las,
                                 bandwidth_intercept = 0,
                                 limits = c(0,1),
                                 ground_height = NULL){
  # TODO: Tests, e.g.: are the input values plausible
  #       - maybe progress bar?
  #       - require libraries inside function

  validate_ground_height(ground_height, las)

  # If ground_height is a list of arguments, pass them to
  # lidR::rasterize_terrain
  if (is.list(ground_height)) {
    ground_height <- do.call(lidR::rasterize_terrain,
                             args = c(las = las, ground_height)
    )
  }

  # if the limits vector is longer, only the min and max values will be taken
  # into account
  my_limits <- range(limits)

  # normalize point cloud if applicable
  if(!is.null(ground_height)){
    err.msg <- "Ground height raster does not cover the area of the point cloud."
    assert_that_raster_covers_las_point_cloud(ground_height, las, err.msg)

    norm.las <- lidR::normalize_height(las = las,
                                       algorithm = kriging(),
                                       dtm = ground_height)
  } else {
   norm.las <- las
  }

  # define the resolution for the chm: if the point density is higher than 16
  # in more than half of the relevant area, use 0.25m resolution. If point
  # density is higher than 5, use 0.5 m resolution, otherwise 1 m.
  dens <- lidR::rasterize_density(norm.las, res = 1)
  if(terra::global(dens, function(x) sum(x >= 16)) >=
     0.5 * terra::global(dens, function(x) sum(x > 0))){
    chm.res <- 0.25
  } else if(terra::global(dens, function(x) sum(x >= 5)) >=
            0.5 * terra::global(dens, function(x) sum(x > 0))){
    chm.res <- 0.5
  } else{
    chm.res <- 1
  }

  # create canopy height model
  chm <- lidR::rasterize_canopy(norm.las,
                                res = chm.res,
                                algorithm = lidR::p2r(subcircle = 0.25))

  # fill in NA values with 0
  chm[is.na(chm)] <- 0

  # watershed delineation
  wsh <- lidR::watershed(chm)()

  # convert to polygon
  poly.wsh <- terra::as.polygons(wsh)

  # tree height as highest chm value
  height.points <- data.table(terra::extract(chm, poly.wsh))
  trees <- height.points[,.(height = max(Z)), keyby=ID]

  # calculate crown area
  poly.wsh$area <- expanse(poly.wsh, transform=F)
  # calculate crown diameter as square root of area
  poly.wsh$diam <- sqrt(poly.wsh$area)

  # calculate diameter to height ratio
  poly.wsh$height <- trees$height
  poly.wsh$diam.height.ratio <- pmin(pmax((poly.wsh$diam - bandwidth_intercept) /
                                            poly.wsh$height,
                                          my_limits[1]),
                                     my_limits[2])

  # build smooth raster of average ratio in 5m Radius
  dhr.rast <- terra::rasterize(poly.wsh, chm, field = "diam.height.ratio")
  ratio.avg <- terra::focal(x=dhr.rast, w = 11, fun = "mean", na.rm = T, pad=T)

  # where there are NA values, fill with 10m average
  ratio.avg[is.na(ratio.avg)] <- terra::focal(x = dhr.rast,
                                              w = 21,
                                              fun = "mean",
                                              na.rm = T,
                                              pad=T)
  # if there are still NA values, arbitrarily fill with 0.5
  ratio.avg[is.na(ratio.avg)] <- 0.5

  return(ratio.avg)
}




# ----------------------------------------------------------
# Li kds raster
#' calculate a raster of kernel diameter slope for AMS3D
#'
#' The function calculates a raster with values for the kernel diameter slope
#' as input for the AMS3D algorithm. It segments the tree crowns with
#' the Li2012 algorithm, calculates a ratio of crown diameter to tree height
#' for each tree, and converts this into a raster
#'
#' @param las a point cloud as lidR LAS object
#' @param bandwidth_intercept a fixed value for kernel_diameter_intercept, which
#' reduces the crown diameters by the given value before calculating the ratio
#' of crown diameter to tree height
#' @param limits a numeric vector with minimum and maximum allowed values for
#' the ratio
#' @param ground_height (optional) either
#' *NULL, indicating that the point cloud is normalized, or
#' *a [SpatRaster][terra::SpatRaster] digital terrain model, or
#' *a list of arguments to the
#' [lidR rasterize_terrain()][lidR::rasterize_terrain()] function to normalize
#' the point cloud.
#' @return terra SpatRaster
#' @section Details
#'
#' The output raster can serve as input for the parameter
#' "kernel_diameter_slope" (hence  "kds") for the function
#' segment_tree_crowns. It averages the ratio of crown diameter to tree height
#' for a 5 m radius, for trees that were detected with the Li2012 tree
#' segmentation algorithm.
#'
#' @export

li_kds_raster <- function(las,
                          bandwidth_intercept = 0,
                          limits = c(0,1),
                          ground_height = NULL# TODO: add parameters for li
                          ){

  # validate input
  validate_ground_height(ground_height, las)





  # If ground_height is a list of arguments, pass them to
  # lidR::rasterize_terrain
  if (is.list(ground_height)) {
    ground_height <- do.call(lidR::rasterize_terrain,
                             args = c(las = point_cloud, ground_height)
    )
  }

  # if the limits vector is longer, only the min and max values will be taken
  # into account
  my_limits <- range(limits)

  # normalize point cloud if applicable
  if(!is.null(ground_height)){
    err.msg <- "Ground height raster does not cover the area of the point cloud."
    assert_that_raster_covers_las_point_cloud(ground_height, las, err.msg)

    las <- lidR::normalize_height(las = las,
                                  algorithm = kriging(),
                                  dtm = ground_height)
  }

  # define the resolution for the chm: if the point density is higher than 16
  # in more than half of the relevant area, use 0.25m resolution. If point
  # density is higher than 5, use 0.5 m resolution, otherwise 1 m.
  dens <- lidR::rasterize_density(las, res = 1)
  if(terra::global(dens, function(x) sum(x >= 16)) >=
     0.5 * terra::global(dens, function(x) sum(x > 0))){
    chm.res <- 0.25
  } else if(terra::global(dens, function(x) sum(x >= 5)) >=
            0.5 * terra::global(dens, function(x) sum(x > 0))){
    chm.res <- 0.5
  } else{
    chm.res <- 1
  }

  # create canopy height model
  chm <- lidR::rasterize_canopy(las,
                                res = chm.res,
                                algorithm = lidR::p2r(subcircle = 0.25))

  # fill in NA values with 0
  chm[is.na(chm)] <- 0


  # segment trees with Li2012 algorithm with default parameters
  segm <- lidR::segment_trees(las, lidR::li2012())



  # create crown polygons from point cloud
  crowns <- crown_polygons_treeID(segm)

  # calculate cdr, and cap it with limits
  crowns$diam.height.ratio <- pmin(pmax((crowns$diameter - bandwidth_intercept) /
                                          crowns$height,
                                        limits[1]),
                                   limits[2])

  # build smooth raster of average ratio in 5m Radius
  dhr.rast <- terra::rasterize(crowns, chm, field = "diam.height.ratio")
  ratio.avg <- terra::focal(x=dhr.rast, w = 11, fun = "mean", na.rm = T, pad=T)
  # where there are NA values, fill with 10m average
  ratio.avg[is.na(ratio.avg)] <- terra::focal(x = dhr.rast,
                                              w = 21,
                                              fun = "mean",
                                              na.rm = T,
                                              pad=T)
  # if there are still NA values, arbitrarily fill with 0.5
  ratio.avg[is.na(ratio.avg)] <- 0.5

  return(ratio.avg)
}



####
# Helper functions

crown_polygons_treeID <- function(las){
  # takes a normalized and segmented point cloud with a column "treeID"
  # and returns a SpatVector object with the convex hull of each tree as polygon

  # create data frame with points
  points.df <- las@data[,c("X","Y", "Z", "treeID")]
  # remove NA treeID
  points.df <- points.df[!is.na(points.df$treeID)]
  # convert into terra vector object
  points.vect <- vect(points.df, geom=c("X", "Y"))

  # calculate convex hull for every tree
  crowns <- terra::hull(points.vect, by = "treeID")

  # calculate area
  crowns$area <- suppressWarnings( terra::expanse(crowns) )
  crowns$diameter <- sqrt(crowns$area)
  # add Z as the highest Z value in the cluster
  crowns <- merge(crowns, aggregate(Z ~ treeID, points.df, max), by = "treeID")
  # rename column
  names(crowns)[names(crowns)=="Z"] <- "height"

  return(crowns)
}


