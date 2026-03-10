test_that("bl_prepare_data returns bl_data class", {
  result <- bl_prepare_data(iris, class_col = "Species",
                            target_class = "versicolor")
  expect_s3_class(result, "bl_data")
})

test_that("class column is binary 0/1 after multiclass conversion", {
  result <- bl_prepare_data(iris, class_col = "Species",
                            target_class = "versicolor")
  expect_true(all(result$train_data$class %in% c(0L, 1L)))
  expect_true(all(result$test_data$class %in% c(0L, 1L)))
})

test_that("target_class maps correctly to 1", {
  result <- bl_prepare_data(iris, class_col = "Species",
                            target_class = "versicolor")
  # All versicolor rows should have class = 1
  original_labels <- iris$Species[seq_len(nrow(iris))]
  # At least some 1s should exist
  expect_true(any(result$train_data$class == 1L))
  expect_true(any(result$train_data$class == 0L))
})

test_that("feature_cols defaults to all non-class columns", {
  result <- bl_prepare_data(iris, class_col = "Species",
                            target_class = "versicolor")
  expected_features <- setdiff(names(iris), "Species")
  expect_equal(result$var_names, expected_features)
  expect_equal(result$num_vars, 4L)
})

test_that("custom feature_cols is respected", {
  result <- bl_prepare_data(iris, class_col = "Species",
                            target_class = "versicolor",
                            feature_cols = c("Sepal.Length", "Petal.Length"))
  expect_equal(result$var_names, c("Sepal.Length", "Petal.Length"))
  expect_equal(result$num_vars, 2L)
})

test_that("train_fraction produces correct split sizes", {
  result <- bl_prepare_data(iris, class_col = "Species",
                            target_class = "versicolor",
                            train_fraction = 0.8)
  total <- nrow(result$train_data) + nrow(result$test_data)
  expect_equal(total, nrow(iris))
  expect_equal(nrow(result$train_data), floor(0.8 * nrow(iris)))
})

test_that("same seed produces identical splits", {
  r1 <- bl_prepare_data(iris, class_col = "Species",
                        target_class = "versicolor", seed = 42L)
  r2 <- bl_prepare_data(iris, class_col = "Species",
                        target_class = "versicolor", seed = 42L)
  expect_equal(r1$train_data, r2$train_data)
  expect_equal(r1$test_data, r2$test_data)
})

test_that("error when class_col not in data", {
  expect_error(
    bl_prepare_data(iris, class_col = "NotAColumn",
                    target_class = "versicolor"),
    "not found in data"
  )
})

test_that("error when target_class value not found", {
  expect_error(
    bl_prepare_data(iris, class_col = "Species",
                    target_class = "octopus"),
    "not found in column"
  )
})

test_that("error when train_fraction outside (0,1)", {
  expect_error(
    bl_prepare_data(iris, class_col = "Species",
                    target_class = "versicolor", train_fraction = 1.5),
    "strictly between"
  )
})

test_that("print.bl_data runs without error", {
  result <- bl_prepare_data(iris, class_col = "Species",
                            target_class = "versicolor")
  expect_output(print(result), "<bl_data>")
})
