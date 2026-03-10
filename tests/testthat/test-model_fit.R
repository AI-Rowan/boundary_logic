bl_dat <- bl_prepare_data(iris, class_col = "Species",
                          target_class = "versicolor")

test_that("bl_fit_model returns bl_model class", {
  result <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
                         model_type = "GLM")
  expect_s3_class(result, "bl_model")
})

test_that("model_type stored correctly", {
  result <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
                         model_type = "GLM")
  expect_equal(result$model_type, "GLM")
})

test_that("var_names stored correctly", {
  result <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
                         model_type = "GLM")
  expect_equal(result$var_names, bl_dat$var_names)
})

test_that("accuracy is numeric in [0,1]", {
  result <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
                         model_type = "GLM")
  expect_true(is.numeric(result$accuracy))
  expect_true(result$accuracy >= 0 && result$accuracy <= 1)
})

test_that("gini is numeric in [-1,1]", {
  result <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
                         model_type = "GLM")
  expect_true(is.numeric(result$gini))
  expect_true(result$gini >= -1 && result$gini <= 1)
})

test_that("pred_function returns vector same length as nrow(new_data)", {
  result   <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
                           model_type = "GLM")
  new_data <- bl_dat$train_data[1:10, bl_dat$var_names, drop = FALSE]
  preds    <- boundarylogic:::.pred_function(
    result$model, result$model_type, result$rounding, new_data
  )
  expect_length(preds, 10L)
  expect_true(all(preds >= 0 & preds <= 1))
})

test_that("error for unrecognised model_type", {
  expect_error(
    bl_fit_model(bl_dat$train_data, bl_dat$var_names,
                 model_type = "BANANA"),
    "not supported"
  )
})

test_that("error when class column missing", {
  bad_data <- bl_dat$train_data[, bl_dat$var_names, drop = FALSE]
  expect_error(
    bl_fit_model(bad_data, bl_dat$var_names, model_type = "GLM"),
    "column named 'class'"
  )
})

test_that("floor rounding: predictions have <= rounding decimal places", {
  result <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
                         model_type = "GLM", rounding = 2L)
  preds  <- boundarylogic:::.pred_function(
    result$model, result$model_type, 2L,
    bl_dat$train_data[, bl_dat$var_names, drop = FALSE]
  )
  # floor to 2 d.p.: preds * 100 should have no fractional part
  expect_true(all(abs(preds * 100 - round(preds * 100)) < 1e-9))
})

test_that("print.bl_model runs without error", {
  result <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
                         model_type = "GLM")
  expect_output(print(result), "<bl_model>")
})
