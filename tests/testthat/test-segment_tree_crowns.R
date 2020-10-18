# For the LAScatalog tests
# TODO Check whether enabling parallelization like this is CRAN conform
future::plan(strategy = future::multisession)

test_that("the different methods for the generic work", {

  # Method dispatch for data.frame
  pseudo_point_cloud <- data.table::data.table(X = 1., Y = 1., Z = 1.)
  expect_equal(
    segment_tree_crowns(pseudo_point_cloud,
      crown_diameter_2_tree_height = 0.25,
      crown_height_2_tree_height = 0.5
    ),
    cbind(pseudo_point_cloud, crown_id = 0)
  )

  # Method dispatch for lidR::LAS object
  las_object <- lidR::readLAS(
    system.file("extdata", "MixedConifer.laz", package = "lidR")
  )
  segmented_las <- segment_tree_crowns(las_object,
    crown_diameter_2_tree_height = 0.25,
    crown_height_2_tree_height = 0.5
  )
  expect_s4_class(segmented_las, "LAS")
  expect_named(segmented_las@data,
    expected = c(names(las_object@data), "crown_id")
  )

  # Method dispatch for lidR::LAScatalog object
  las_catalog <- lidR::readLAScatalog(
    system.file("extdata", "MixedConifer.laz", package = "lidR"),
    chunk_size = 50
  )
  segmented_las <- segment_tree_crowns(las_catalog,
    crown_diameter_2_tree_height = 0.25,
    crown_height_2_tree_height = 0.5
  )
  expect_s4_class(segmented_las, "LAS")
  expect_named(segmented_las@data,
    expected = c(names(las_object@data), "crown_id")
  )

  # Testing visually that the output is okay:
  # lidR::plot(
  #   segment_tree_crowns(las_catalog,
  #     crown_diameter_2_tree_height = 0.25,
  #     crown_height_2_tree_height = 0.5
  #   ), color = "crown_id", colorPalette = lidR::random.colors(1000)
  # )
})

test_that("the calculation of unique crown IDs works for LAScatalog", {

  # TODO Stuff seems to work but I haven't found a way to systematically test it
  # yet.

  las_object <- lidR::readLAS(
    system.file("extdata", "MixedConifer.laz", package = "lidR"),
    filter = "-drop_z_below 1"
  )
  segmented_las <- segment_tree_crowns(las_object,
    crown_diameter_2_tree_height = 0.25,
    crown_height_2_tree_height = 0.5
  )

  las_catalog <- lidR::readLAScatalog(
    system.file("extdata", "MixedConifer.laz", package = "lidR"),
    chunk_size = 80,
    filter = "-drop_z_below 1"
  )

  segmented_las_catalog <- segment_tree_crowns(las_catalog,
    crown_diameter_2_tree_height = 0.25,
    crown_height_2_tree_height = 0.5
  )

  expect_equal(
    uniqueN(segmented_las@data$crown_id),
    uniqueN(segmented_las_catalog@data$crown_id)
  )

  # chunks <- lidR::catalog_makechunks(las_catalog)

  # chunk_id_sets <- lapply(chunks, function(chunk) {
  #   chunk_extent <- lidR::extent(chunk)
  #   chunk_points <- lidR::filter_poi(segmentation,
  #     chunk_extent@xmin <= X & X < chunk_extent@xmax &
  #       chunk_extent@ymin <= Y & Y < chunk_extent@ymax
  #   )
  #   unique(chunk_points@data$crown_id)
  # })

  # lidR::las_check(segmentation)
  #
  # lidR::plot(
  #   lidR::filter_poi(segmentation, crown_id > 0),
  #   color = "crown_id",
  #   colorPalette = lidR::random.colors(uniqueN(segmentation@data$crown_id) - 1),
  #   legend = TRUE,
  #   axis = TRUE
  # )
  #
  # lidR::plot(
  #   lidR::filter_poi(segmentation,
  #     crown_id %in% intersect(chunk_id_sets[[4]], chunk_id_sets[[3]])
  #   ),
  #   axis = TRUE,
  #   legend = TRUE,
  #   color = "crown_id"
  # )
})
