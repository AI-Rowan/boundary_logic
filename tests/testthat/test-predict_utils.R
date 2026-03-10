bl_dat <- bl_prepare_data(iris, class_col = "Species",
                          target_class = "versicolor")
bl_mod <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
                       model_type = "GLM")

test_that(".pred_function returns numeric vector same length as nrow(new_data)", {
  new_data <- bl_dat$train_data[1:5, bl_dat$var_names, drop = FALSE]
  preds    <- boundarylogic:::.pred_function(
    bl_mod$model, bl_mod$model_type, 2L, new_data
  )
  expect_length(preds, 5L)
  expect_true(is.numeric(preds))
})

test_that(".pred_function GLM predictions are in [0,1]", {
  new_data <- bl_dat$train_data[, bl_dat$var_names, drop = FALSE]
  preds    <- boundarylogic:::.pred_function(
    bl_mod$model, "GLM", 2L, new_data
  )
  expect_true(all(preds >= 0 & preds <= 1))
})

test_that("floor rounding with rounding=2 has no more than 2 decimal places", {
  new_data <- bl_dat$train_data[, bl_dat$var_names, drop = FALSE]
  preds    <- boundarylogic:::.pred_function(
    bl_mod$model, "GLM", 2L, new_data
  )
  # floor to 2 d.p. means preds * 100 is always an integer
  expect_true(all(abs(preds * 100 - floor(preds * 100)) < 1e-9))
})

test_that(".pred_function throws error for unknown model_type", {
  new_data <- bl_dat$train_data[1L, bl_dat$var_names, drop = FALSE]
  expect_error(
    boundarylogic:::.pred_function(bl_mod$model, "BANANA", 2L, new_data),
    "Unknown model_type"
  )
})
