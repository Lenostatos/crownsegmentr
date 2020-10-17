test_that(paste0(
    "try_to_extract_coordinate_data() throws error for non-data.frame-like ",
    "objects"
  ), {
    numeric_scalar <- 1.0
    numeric_vector <- c(.0, 1.0, 5, 10.0)
    numeric_matrix <- matrix(data = numeric_vector, ncol = 2)

    expect_error(
      try_to_extract_coordinate_data(numeric_scalar),
      regexp = paste0(
        "^The coordinate data needs to be stored in a data\\.table or ",
        "data\\.frame\\.$"
      )
    )
    expect_error(
      try_to_extract_coordinate_data(numeric_vector),
      regexp = paste0(
        "^The coordinate data needs to be stored in a data\\.table or ",
        "data\\.frame\\.$"
      )
    )
    expect_error(
      try_to_extract_coordinate_data(numeric_matrix),
      regexp = paste0(
        "^The coordinate data needs to be stored in a data\\.table or ",
        "data\\.frame\\.$"
      )
    )
  }
)

test_that(paste0(
    "try_to_extract_coordinate_data() throws error for data.frames with less ",
    "than three columns"
  ), {
    numeric_vector <- c(.0, 1.0, 5, 10.0)
    one_column_data_frame <- data.frame(a = numeric_vector)
    two_column_data_frame <- data.frame(a = numeric_vector, b = numeric_vector)

    expect_error(
      try_to_extract_coordinate_data(one_column_data_frame),
      regexp = paste0(
        "^The coordinates table needs to have at least three numeric columns ",
        "for x-, y-, and z-coordinates but there are only [0-9]+\\.$"
      )
    )
    expect_error(
      try_to_extract_coordinate_data(two_column_data_frame),
      regexp = paste0(
        "^The coordinates table needs to have at least three numeric columns ",
        "for x-, y-, and z-coordinates but there are only [0-9]+\\.$"
      )
    )
  }
)

# TODO expand this test into multiple smaller ones.
test_that(paste0(
    "try_to_extract_coordinate_data() returns a data.frame or data.table ",
    "(depending on the input type) with just the x-, y-, and z-coordinates."
  ), {
    numeric_vector <- c(.0, 1.0, 5, 10.0)

    valid_data_frame_w_normal_colnames <- data.frame(
      x = numeric_vector, Y = numeric_vector, z = numeric_vector
    )
    valid_data_frame_w_normal_colnames_unordered <- data.frame(
      Y = numeric_vector, z = numeric_vector, x = numeric_vector
    )

    expect_identical(
      try_to_extract_coordinate_data(valid_data_frame_w_normal_colnames),
      expected = valid_data_frame_w_normal_colnames
    )

    expect_equal(
      try_to_extract_coordinate_data(
        valid_data_frame_w_normal_colnames_unordered
      ),
      expected = valid_data_frame_w_normal_colnames
    )

    valid_data_frame_w_normal_colnames <- as.data.table(
      valid_data_frame_w_normal_colnames
    )
    valid_data_frame_w_normal_colnames_unordered <- as.data.table(
      valid_data_frame_w_normal_colnames_unordered
    )

    expect_identical(
      try_to_extract_coordinate_data(valid_data_frame_w_normal_colnames),
      expected = valid_data_frame_w_normal_colnames
    )

    expect_equal(
      try_to_extract_coordinate_data(
        valid_data_frame_w_normal_colnames_unordered
      ),
      expected = valid_data_frame_w_normal_colnames
    )
  }
)
