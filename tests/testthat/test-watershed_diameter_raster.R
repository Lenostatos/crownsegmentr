test_that("watershed_diameter_raster works",{
  skip_if_not_installed("EBImage")

  test_point_cloud_file_path <- system.file(
    "extdata", "MixedConifer.laz",
    package = "lidR"
  )
  test_point_cloud_header <- lidR::readLASheader(test_point_cloud_file_path)
  test_point_cloud <- lidR::clip_rectangle(
    lidR::readLAS(test_point_cloud_file_path),
    xleft = test_point_cloud_header@PHB$`Min X`,
    ybottom = test_point_cloud_header@PHB$`Min Y`,
    xright = test_point_cloud_header@PHB$`Min X` + 50,
    ytop = test_point_cloud_header@PHB$`Min Y` + 50
  )

  testthat::expect_no_failure({
    watershed_raster <- watershed_diameter_raster(test_point_cloud)

    assert_that_raster_fits_point_cloud(watershed_raster,
                                        test_point_cloud)
  })

  # test with modified parameters
  testthat::expect_no_failure({
    watershed_kds <- watershed_diameter_raster(test_point_cloud,
                                               crown_diameter_constant = 10,
                                               limits = c(0.3,0.4),
                                               smoothing_radius = 0)
  })
})





