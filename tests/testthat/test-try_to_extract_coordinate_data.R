testthat::test_that(
  paste0(
    "try_to_extract_coordinate_data() throws error for non-data.frame-like ",
    "objects"
  ),
  {
    numeric_scalar <- 1.0
    numeric_vector <- c(.0, 1.0, 5, 10.0)
    numeric_matrix <- matrix(data = numeric_vector, ncol = 2)

    testthat::expect_error(
      try_to_extract_coordinate_data(numeric_scalar, "numeric_scalar"),
      regexp = " needs to be a data.table or data.frame.$"
    )
    testthat::expect_error(
      try_to_extract_coordinate_data(numeric_vector, "numeric_vector"),
      regexp = " needs to be a data.table or data.frame.$"
    )
    testthat::expect_error(
      try_to_extract_coordinate_data(numeric_matrix, "numeric_matrix"),
      regexp = " needs to be a data.table or data.frame.$"
    )
  }
)

testthat::test_that(
  paste0(
    "try_to_extract_coordinate_data() throws error for data.frames with less ",
    "than three columns"
  ),
  {
    numeric_vector <- c(.0, 1.0, 5, 10.0)
    one_column_data_frame <- data.frame(a = numeric_vector)
    two_column_data_frame <- data.frame(a = numeric_vector, b = numeric_vector)

    testthat::expect_error(
      try_to_extract_coordinate_data(
        one_column_data_frame, "one_column_data_frame"
      ),
      regexp = paste0(
        " needs to have at least three columns for x-, y-, and ",
        "z-coordinates.$"
      )
    )
    testthat::expect_error(
      try_to_extract_coordinate_data(
        two_column_data_frame, "two_column_data_frame"
      ),
      regexp = paste0(
        " needs to have at least three columns for x-, y-, and ",
        "z-coordinates.$"
      )
    )
  }
)

# TODO expand this test into multiple smaller ones.
testthat::test_that(
  paste0(
    "try_to_extract_coordinate_data() returns a data.frame with just the x-, ",
    "y-, and z-coordinates."
  ),
  {
    numeric_vector <- c(.0, 1.0, 5, 10.0)

    valid_data_frame_w_normal_colnames <- data.frame(
      x = numeric_vector, y = numeric_vector, z = numeric_vector
    )
    valid_data_frame_w_normal_colnames_unordered <- data.frame(
      y = numeric_vector, z = numeric_vector, x = numeric_vector
    )

    testthat::expect_identical(
      try_to_extract_coordinate_data(
        valid_data_frame_w_normal_colnames, "valid_data_frame_w_normal_colnames"
      ),
      expected = valid_data_frame_w_normal_colnames
    )

    returned_object <- try_to_extract_coordinate_data(
      valid_data_frame_w_normal_colnames_unordered,
      "valid_data_frame_w_normal_colnames_unordered"
    )
    testthat::expect_equal(names(returned_object)[1], expected = "x")
    testthat::expect_equal(names(returned_object)[2], expected = "y")
    testthat::expect_equal(names(returned_object)[3], expected = "z")
  }
)
