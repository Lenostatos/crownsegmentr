test_that("The kds_raster functions work",{

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
  range(test_point_cloud@data$Z)

  testthat::expect_no_failure({
    li_kds <- li_kds_raster(test_point_cloud)

    assert_that_raster_fits_point_cloud(li_kds,
                                        test_point_cloud)

    #ws_kds <- watershed_kds_raster(test_point_cloud)

    #assert_that_raster_fits_point_cloud(ws_kds,
    #                                    test_point_cloud)
  })

})




