# Load a point cloud of some trees included in the lidR package
point_cloud <- lidR::readLAS(system.file(
  "extdata", "MixedConifer.laz", package = "lidR"
))


# Segment the point cloud
segmented_point_cloud <- crownsegmentr::segment_tree_crowns(
  point_cloud,
  crown_diameter_to_tree_height = 0.25,
  crown_height_to_tree_height = 0.5
)

# Generate some crown colors for plotting
crown_colors <- c(
  # the first color is a very light grey used for unsegmented points
  "#DDDDDD",
  # the colors of the actual crowns are randomly generated
  lidR::random.colors(
    n = length(unique(segmented_point_cloud@data$crown_id)) - 1
  )
)

# Plot the segmented crown bodies
lidR::plot(
  segmented_point_cloud,
  color = "crown_id",
  colorPalette = crown_colors,
  size = 3,
  axis = TRUE
)


# Observe how adjacent crowns are undersegmented when the
# crown_diameter_to_tree_height parameter is set too high
segmented_point_cloud <- crownsegmentr::segment_tree_crowns(
  point_cloud,
  crown_diameter_to_tree_height = 0.5,
  crown_height_to_tree_height = 0.5
)

lidR::plot(
  segmented_point_cloud,
  color = "crown_id",
  colorPalette = c(
    "#DDDDDD",
    lidR::random.colors(
      n = length(unique(segmented_point_cloud@data$crown_id)) - 1
    )
  ),
  size = 3,
  axis = TRUE
)


# In order to get a more in-depth impression of how the segmentation works,
# also return and plot modes and centroids.
segmentation_results <- crownsegmentr::segment_tree_crowns(
  point_cloud,
  crown_diameter_to_tree_height = 0.25,
  crown_height_to_tree_height = 0.5,
  also_return_modes = TRUE,
  also_return_centroids = TRUE
)

lidR::plot(
  segmentation_results$segmented_point_cloud,
  color = "crown_id",
  colorPalette = crown_colors,
  size = 3,
  axis = TRUE
)

lidR::plot(
  segmentation_results$modes,
  color = "crown_id",
  colorPalette = crown_colors,
  size = 3,
  axis = TRUE
)

lidR::plot(
  segmentation_results$centroids,
  color = "crown_id",
  colorPalette = crown_colors,
  size = 1,
  axis = TRUE
)


# See how the processing gets faster and segmentation accuracy worse when using
# a larger distance at which centroids are considered to have converged and a
# larger DBSCAN neighborhood.
segmented_point_cloud <- crownsegmentr::segment_tree_crowns(
  point_cloud,
  crown_diameter_to_tree_height = 0.5,
  crown_height_to_tree_height = 0.5,
  centroid_convergence_distance = 0.02,
  dbscan_neighborhood_radius = 0.5
)

lidR::plot(
  segmented_point_cloud,
  color = "crown_id",
  colorPalette = c(
    "#DDDDDD",
    lidR::random.colors(
      n = length(unique(segmented_point_cloud@data$crown_id)) - 1
    )
  ),
  size = 3,
  axis = TRUE
)
