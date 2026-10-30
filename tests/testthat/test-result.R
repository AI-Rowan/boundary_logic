bl_dat  <- bl_prepare_data(iris, class_col = "Species",
                           target_class = "versicolor")
bl_mod  <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
                        model_type = "GLM")
bl_proj <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                               method = "PCA")
bl_grid <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = 30L)

test_that("bl_assemble returns bl_result class", {
  result <- bl_assemble(bl_dat, bl_mod, bl_proj, bl_grid)
  expect_s3_class(result, "bl_result")
})

test_that("all required fields are present", {
  result   <- bl_assemble(bl_dat, bl_mod, bl_proj, bl_grid)
  required <- c("train_data", "test_data", "var_names", "num_vars",
                "model", "model_type", "cutoff", "b_margin",
                "V", "tV", "X_center", "X_sd",
                "method", "standardise", "biplot_obj",
                "polygon", "hull_fraction", "biplot_grid",
                "accuracy", "gini", "call", "created_at")
  expect_true(all(required %in% names(result)))
})

test_that("var_names matches bl_model var_names", {
  result <- bl_assemble(bl_dat, bl_mod, bl_proj, bl_grid)
  expect_equal(result$var_names, bl_mod$var_names)
})

test_that("model_type matches bl_model model_type", {
  result <- bl_assemble(bl_dat, bl_mod, bl_proj, bl_grid)
  expect_equal(result$model_type, bl_mod$model_type)
})

test_that("standardise is FALSE when method is CVA", {
  pred_col    <- as.factor(ifelse(bl_dat$train_data$class == 1, "red", "blue"))
  bl_proj_cva <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                                     method = "CVA",
                                     cva_classes = pred_col)
  bl_grid_cva <- bl_build_grid(bl_dat$train_data, bl_proj_cva, bl_mod,
                               m = 30L)
  result      <- bl_assemble(bl_dat, bl_mod, bl_proj_cva, bl_grid_cva)
  expect_false(result$standardise)
})

test_that("hull_fraction is a numeric value from bl_grid", {
  result <- bl_assemble(bl_dat, bl_mod, bl_proj, bl_grid)
  expect_true(is.numeric(result$hull_fraction))
  expect_true(result$hull_fraction >= 0 && result$hull_fraction <= 1)
})

test_that("print.bl_result runs without error and returns invisibly", {
  result   <- bl_assemble(bl_dat, bl_mod, bl_proj, bl_grid)
  returned <- expect_output(print(result), "<bl_result>")
})

test_that("bl_assemble stops when the grid and model cutoffs differ", {
  mod_04        <- bl_mod
  mod_04$cutoff <- 0.4
  expect_error(bl_assemble(bl_dat, mod_04, bl_proj, bl_grid),
               "cutoff = 0.5 but the model's cutoff is 0.4", fixed = TRUE)
})

test_that("a grid without a cutoff field (older object) still assembles", {
  old_grid        <- bl_grid
  old_grid$cutoff <- NULL
  expect_s3_class(bl_assemble(bl_dat, bl_mod, bl_proj, old_grid), "bl_result")
})

test_that("bl_build_result builds the grid at the model's cutoff", {
  mod_04        <- bl_mod
  mod_04$cutoff <- 0.4
  res <- bl_build_result(bl_dat, mod_04, method = "PCA", m = 30L)
  expect_identical(res$biplot_grid$cutoff, 0.4)
  expect_identical(res$cutoff, 0.4)
})

test_that("summary.bl_result runs without error", {
  result <- bl_assemble(bl_dat, bl_mod, bl_proj, bl_grid)
  expect_output(summary(result), "Confusion matrix")
})
