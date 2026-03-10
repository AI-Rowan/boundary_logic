bl_dat  <- bl_prepare_data(iris, class_col = "Species",
                           target_class = "versicolor")
bl_mod  <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
                        model_type = "GLM")
bl_proj <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                               method = "PCA")

test_that("bl_build_grid returns bl_grid class", {
  result <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = 30L)
  expect_s3_class(result, "bl_grid")
})

test_that("Zgrid has m^2 rows and 2 columns when calc_hull=FALSE", {
  m      <- 30L
  result <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = m,
                          calc_hull = FALSE)
  expect_equal(nrow(result$Zgrid), m^2L)
  expect_equal(ncol(result$Zgrid), 2L)
})

test_that("Xgrid has same rows as Zgrid and cols = num_vars", {
  m      <- 30L
  result <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = m,
                          calc_hull = FALSE)
  expect_equal(nrow(result$Xgrid), nrow(result$Zgrid))
  expect_equal(ncol(result$Xgrid), bl_dat$num_vars)
})

test_that("grid_prob has length m^2 and all values in [0,1]", {
  m      <- 30L
  result <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = m)
  expect_length(result$grid_prob, nrow(result$Zgrid))
  expect_true(all(result$grid_prob >= 0 & result$grid_prob <= 1))
})

test_that("col_value has same length as grid_prob", {
  result <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = 30L)
  expect_length(result$col_value, length(result$grid_prob))
})

test_that("ct is a list (may be empty)", {
  result <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = 30L)
  expect_true(is.list(result$ct))
})

test_that("xseq and yseq have length m", {
  m      <- 30L
  result <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = m)
  expect_length(result$xseq, m)
  expect_length(result$yseq, m)
})

test_that("supplied polygon is reused unchanged", {
  result1  <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = 30L)
  poly     <- result1$polygon
  result2  <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = 30L,
                            polygon = poly)
  # The polygon object stored should be the one we supplied
  expect_identical(result2$polygon, poly)
})

test_that("print.bl_grid runs without error", {
  result <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = 30L)
  expect_output(print(result), "<bl_grid>")
})
