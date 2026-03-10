bl_dat <- bl_prepare_data(iris, class_col = "Species",
                          target_class = "versicolor")

train_features <- bl_dat$train_data[, bl_dat$var_names, drop = FALSE]

test_that("get_variable_ranges returns list with one entry per column", {
  ranges <- get_variable_ranges(train_features)
  expect_length(ranges, ncol(train_features))
})

test_that("get_filter_logical_vector returns logical vector of nrow(data)", {
  ranges <- get_variable_ranges(train_features)
  mask   <- get_filter_logical_vector(train_features, ranges)
  expect_length(mask, nrow(train_features))
  expect_true(is.logical(mask))
})

test_that("all training rows pass their own range filter", {
  ranges <- get_variable_ranges(train_features)
  mask   <- get_filter_logical_vector(train_features, ranges)
  expect_true(all(mask))
})

test_that("a point outside the range is correctly identified as FALSE", {
  ranges   <- get_variable_ranges(train_features)
  out_row  <- train_features[1L, , drop = FALSE]
  # Force one value far outside range
  out_row[[1L]] <- max(train_features[[1L]]) + 1000
  mask <- get_filter_logical_vector(out_row, ranges)
  expect_false(mask[1L])
})

test_that("list(TRUE) filter returns all TRUE (no filter)", {
  mask <- get_filter_logical_vector(train_features, list(TRUE))
  expect_true(all(mask))
  expect_length(mask, nrow(train_features))
})
