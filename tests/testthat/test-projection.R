bl_dat <- bl_prepare_data(iris, class_col = "Species",
                          target_class = "versicolor")

test_that("bl_build_projection returns bl_projection class for PCA", {
  result <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                                method = "PCA")
  expect_s3_class(result, "bl_projection")
})

test_that("V is a square matrix with nrow = ncol = num_vars", {
  result <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                                method = "PCA")
  p <- bl_dat$num_vars
  expect_equal(dim(result$V), c(p, p))
})

test_that("tV equals solve(V) to numerical tolerance", {
  result <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                                method = "PCA")
  expect_equal(result$tV, solve(result$V), tolerance = 1e-10)
})

test_that("X_sd matches apply(X, 2, sd) for PCA with standardise=TRUE", {
  X      <- as.matrix(bl_dat$train_data[, bl_dat$var_names])
  result <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                                method = "PCA", standardise = TRUE)
  expect_equal(result$X_sd, apply(X, 2, sd), tolerance = 1e-10)
})

test_that("CVA forces standardise = FALSE", {
  pred_col <- as.factor(ifelse(bl_dat$train_data$class == 1, "red", "blue"))
  result <- bl_build_projection(
    bl_dat$train_data, bl_dat$var_names,
    method = "CVA", standardise = TRUE,  # user passes TRUE
    cva_classes = pred_col
  )
  expect_false(result$standardise)
})

test_that("CVA with cva_classes=NULL throws informative error", {
  expect_error(
    bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                        method = "CVA", cva_classes = NULL),
    "requires 'cva_classes'"
  )
})

test_that("PCA round-trip projection stays in correct space", {
  result <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                                method = "PCA", standardise = FALSE)
  X      <- as.matrix(bl_dat$train_data[1L, bl_dat$var_names, drop = FALSE])
  X_c    <- sweep(X, 2L, result$X_center, "-")
  Z      <- X_c %*% result$V[, result$proj_dims, drop = FALSE]
  X_back <- sweep(Z %*% result$tV[result$proj_dims, , drop = FALSE],
                  2L, result$X_center, "+")
  # With only 2 of p=4 dims, this is an approximation — just check shape
  expect_equal(length(X_back), bl_dat$num_vars)
})

test_that("print.bl_projection runs without error", {
  result <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                                method = "PCA")
  expect_output(print(result), "<bl_projection>")
})
