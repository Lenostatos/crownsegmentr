# test_that("watershed kds works for normalized raster",{
# # load point cloud data
# test_point_cloud_file_path <- system.file(
#   "extdata", "MixedConifer.laz",
#   package = "lidR"
# )
# test_point_cloud_header <- lidR::readLASheader(test_point_cloud_file_path)
# test_point_cloud <- lidR::clip_rectangle(
#   lidR::readLAS(test_point_cloud_file_path),
#   xleft = test_point_cloud_header@PHB$`Min X`,
#   ybottom = test_point_cloud_header@PHB$`Min Y`,
#   xright = test_point_cloud_header@PHB$`Min X` + 50,
#   ytop = test_point_cloud_header@PHB$`Min Y` + 50
# )
#
# kds_raster <- watershed_kds_raster(test_point_cloud)
#
# expect(terra::ext(kds_raster) == lidR::ext(test_point_cloud))
#
# }
# )
#


test_that("li kds works for normalized raster",{
  # load point cloud data
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

    kds_raster <- li_kds_raster(test_point_cloud)

    assert_that_raster_fits_point_cloud(kds_raster,
                                        test_point_cloud)
  })
}
)
