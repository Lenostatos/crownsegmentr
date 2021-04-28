# Load a point cloud of some trees included in the lidR package
point_cloud <- lidR::readLAS(system.file(
  "extdata/MixedConifer.laz", package = "lidR"
))

# Set up a plotting function for segmented point clouds
plot_segmented_point_cloud <- function(point_cloud,
                                       crown_colors = NULL,
                                       size = 3) {

  # Generate random crown colors
  if (is.null(crown_colors)) {
    crown_colors <- lidR::random.colors(
      n = length(unique(point_cloud@data[["crown_id"]])) - 1
    )
  }

  # Plot the segmented crown bodies
  lidR::plot(
    point_cloud,
    color = "crown_id",
    colorPalette = crown_colors,
    size = size,
    axis = TRUE
  )
}


# Segment Normalized Point Clouds ----------------------------------------

segmented_point_cloud <- crownsegmentr::segment_tree_crowns(
  point_cloud,
  crown_diameter_to_tree_height = 0.25,
  crown_height_to_tree_height = 0.5
)

plot_segmented_point_cloud(segmented_point_cloud)


# Exclude Points Below any Height on the Fly ------------------------------

segmented_point_cloud <- crownsegmentr::segment_tree_crowns(
  point_cloud,
  crown_diameter_to_tree_height = 0.25,
  crown_height_to_tree_height = 0.5,
  segment_crowns_only_above = 15 # exclude points below 15 m
)

plot_segmented_point_cloud(segmented_point_cloud)


# Segment Terraneous Point Clouds -----------------------------------------

terraneous_point_cloud <- lidR::readLAS(system.file(
  "extdata/Topography.laz", package = "lidR"
))

# Either pass arguments for the lidR::grid_terrain() function
segmented_point_cloud <- crownsegmentr::segment_tree_crowns(
  terraneous_point_cloud,
  crown_diameter_to_tree_height = 0.5,
  crown_height_to_tree_height = 1,
  segment_crowns_only_above = 2,
  ground_height = list(algorithm = lidR::tin())
)

# Or pass a ground height raster
ground_height_grid <- lidR::grid_terrain(terraneous_point_cloud,
  algorithm = lidR::tin()
)

segmented_point_cloud <- crownsegmentr::segment_tree_crowns(
  terraneous_point_cloud,
  crown_diameter_to_tree_height = 0.5,
  crown_height_to_tree_height = 1,
  segment_crowns_only_above = 2,
  ground_height = ground_height_grid
)

plot_segmented_point_cloud(segmented_point_cloud)


# Account for Different Crown Shapes with Parameter Rasters ---------------

# Create a raster of crown_diameter_to_tree_height ratios with 0.25 in the
# western half and 0.75 in the eastern half
crown_diameter_to_tree_height_grid <- raster::raster(
  matrix(c(0.25, 0.75), ncol = 2),
  xmn = lidR::extent(point_cloud)@xmin,
  xmx = lidR::extent(point_cloud)@xmax,
  ymn = lidR::extent(point_cloud)@ymin,
  ymx = lidR::extent(point_cloud)@ymax,
  crs = lidR::crs(point_cloud)
)

segmented_point_cloud <- crownsegmentr::segment_tree_crowns(
  point_cloud,
  crown_diameter_to_tree_height = crown_diameter_to_tree_height_grid,
  crown_height_to_tree_height = 0.5
)

# Observe how adjacent crowns are undersegmented in the right half of the plot
# where the crown_diameter_to_tree_height parameter is set too high
plot_segmented_point_cloud(segmented_point_cloud)


# Additionally Return Modes and/or Centroids ------------------------------

# You can also return the modes and centroids which were calculated during the
# segmentation in order to get a more in-depth impression of what happened
# internally.
segmentation_results <- crownsegmentr::segment_tree_crowns(
  point_cloud,
  crown_diameter_to_tree_height = 0.25,
  crown_height_to_tree_height = 0.5,
  also_return_modes = TRUE,
  also_return_centroids = TRUE
)

# Generate crown colors "manually" so that we can use the same colors for
# points, modes, and centroids.
crown_colors <- lidR::random.colors(
  n = length(unique(
    segmentation_results$segmented_point_cloud@data[["crown_id"]]
  )) - 1
)

# Plot the segemented point cloud
plot_segmented_point_cloud(segmentation_results$segmented_point_cloud,
  crown_colors
)
# Plot the modes
plot_segmented_point_cloud(segmentation_results$modes,
  crown_colors
)
# Plot the centroids
plot_segmented_point_cloud(segmentation_results$centroids,
  crown_colors,
  size = 1
)


# Parameter Tweaking ------------------------------------------------------

# Observe how processing gets faster while segmentation accuracy slightly
# decreases when using a larger centroid convergence distance and a larger
# DBSCAN neighborhood.
system.time(
  segmented_point_cloud <- crownsegmentr::segment_tree_crowns(
    point_cloud,
    crown_diameter_to_tree_height = 0.25,
    crown_height_to_tree_height = 0.5,
    centroid_convergence_distance = 0.02,
    dbscan_neighborhood_radius = 0.5
  )
)
lidR::plot(
  segmented_point_cloud,
  color = "crown_id",
  colorPalette = lidR::random.colors(
    n = length(unique(segmented_point_cloud@data$crown_id)) - 1
  ),
  size = 3,
  axis = TRUE
)

system.time(
  segmented_point_cloud <- crownsegmentr::segment_tree_crowns(
    point_cloud,
    crown_diameter_to_tree_height = 0.25,
    crown_height_to_tree_height = 0.5,
    centroid_convergence_distance = 0.01, # default value
    dbscan_neighborhood_radius = 0.3 # default value
  )
)
lidR::plot(
  segmented_point_cloud,
  color = "crown_id",
  colorPalette = lidR::random.colors(
    n = length(unique(segmented_point_cloud@data$crown_id)) - 1
  ),
  size = 3,
  axis = TRUE
)
