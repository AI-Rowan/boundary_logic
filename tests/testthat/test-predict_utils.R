bl_dat <- bl_prepare_data(iris, class_col = "Species",
                          target_class = "versicolor")
bl_mod <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
                       model_type = "GLM")

test_that(".pred_function returns numeric vector same length as nrow(new_data)", {
  new_data <- bl_dat$train_data[1:5, bl_dat$var_names, drop = FALSE]
  preds    <- boundarylogic:::.pred_function(
    bl_mod$model, bl_mod$model_type, new_data
  )
  expect_length(preds, 5L)
  expect_true(is.numeric(preds))
})

test_that(".pred_function GLM predictions are in [0,1]", {
  new_data <- bl_dat$train_data[, bl_dat$var_names, drop = FALSE]
  preds    <- boundarylogic:::.pred_function(
    bl_mod$model, "GLM", new_data
  )
  expect_true(all(preds >= 0 & preds <= 1))
})

test_that("floor rounding: predictions have no more than 3 decimal places", {
  new_data <- bl_dat$train_data[, bl_dat$var_names, drop = FALSE]
  preds    <- boundarylogic:::.pred_function(
    bl_mod$model, "GLM", new_data
  )
  # floor to 3 d.p. means preds * 1000 is always an integer
  expect_true(all(abs(preds * 1000 - floor(preds * 1000)) < 1e-9))
})

test_that(".pred_function throws error for unknown model_type", {
  # Must use a non-workflow object to reach the legacy switch path
  new_data   <- bl_dat$train_data[1L, bl_dat$var_names, drop = FALSE]
  fake_model <- structure(list(), class = "some_raw_model")
  expect_error(
    boundarylogic:::.pred_function(fake_model, "BANANA", new_data),
    "Unknown model_type"
  )
})
