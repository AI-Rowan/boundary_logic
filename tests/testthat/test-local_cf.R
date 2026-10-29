# Phase 3: bl_select_target() resolves the target's true class (class_obs)
# from the same row that supplies its features, and plot.bl_local_result()
# colours the target from it.

lcf_bl_dat <- bl_prepare_data(iris, class_col = "Species",
                              target_class = "versicolor", verbose = FALSE)
lcf_bl_mod <- bl_fit_model(lcf_bl_dat$train_data, lcf_bl_dat$var_names,
                           model_type = "GLM")
lcf_bl_results <- suppressWarnings(
  bl_build_result(lcf_bl_dat, lcf_bl_mod, method = "PCA",
                  m = 60L, b_margin = 0.01))
lcf_train <- lcf_bl_results$train_data
lcf_test  <- lcf_bl_results$test_data

test_that("fixture: train_data[1] and test_data[1] differ in class", {
  expect_false(identical(lcf_train$class[1L], lcf_test$class[1L]))
})

test_that("an index into train_data gives train_data[1]'s class", {
  tgt <- bl_select_target(lcf_bl_results, target = 1L, data = lcf_train)
  expect_identical(tgt$class_obs, as.integer(lcf_train$class[1L]))
})

test_that("a one-row data-frame target gives the same class_obs", {
  tgt_idx <- bl_select_target(lcf_bl_results, target = 1L, data = lcf_train)
  tgt_df  <- bl_select_target(lcf_bl_results, target = lcf_train[1L, ])
  expect_identical(tgt_df$class_obs, tgt_idx$class_obs)
  expect_true(is.na(tgt_df$row_id))
})

test_that("the default (test_data) is unchanged", {
  tgt <- bl_select_target(lcf_bl_results, target = 1L)
  expect_identical(tgt$class_obs, as.integer(lcf_test$class[1L]))
  expect_identical(tgt$row_id, 1L)
  expect_equal(tgt$x_obs, lcf_test[1L, lcf_bl_results$var_names, drop = FALSE])
})

test_that("a data-frame target without a class column has class_obs NA", {
  x_only <- lcf_train[1L, lcf_bl_results$var_names, drop = FALSE]
  tgt    <- bl_select_target(lcf_bl_results, target = x_only)
  expect_identical(tgt$class_obs, NA_integer_)
})

test_that("plot.bl_local_result draws a train_data target and an old bl_target", {
  on.exit(grDevices::graphics.off(), add = TRUE)
  grDevices::pdf(NULL)
  tgt <- bl_select_target(lcf_bl_results, target = 1L, data = lcf_train)
  loc <- bl_find_local_cf(lcf_bl_results, tgt, verbose = FALSE)
  skip_if_not(isTRUE(loc$solution_found))
  expect_no_error(plot(loc, print_summary = FALSE))

  # A bl_target built before class_obs existed plots by prediction colour.
  loc_old <- loc
  loc_old$bl_target$class_obs <- NULL
  expect_no_error(plot(loc_old, print_summary = FALSE))
})
