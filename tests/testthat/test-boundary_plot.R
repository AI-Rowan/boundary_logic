# Phase 2: plot.bl_boundary()'s Mahalanobis check describes the plotted
# per-feature measure, and compares it with the full Mahalanobis distance.

bp_dat <- bl_prepare_data(iris, class_col = "Species",
                          target_class = "versicolor", verbose = FALSE)
bp_mod <- bl_fit_model(bp_dat$train_data, bp_dat$var_names, model_type = "GLM")
bp_res <- suppressWarnings(suppressMessages(
  bl_build_result(bp_dat, bp_mod, method = "CVA", m = 60L, b_margin = 0.01)))
bp_bnd <- suppressWarnings(suppressMessages(bl_find_boundary(bp_res)))

plot_quiet <- function(x, ...) {
  on.exit(grDevices::graphics.off(), add = TRUE)
  grDevices::pdf(NULL)
  out <- NULL
  invisible(capture.output(out <- plot(x, ...)))
  out
}

test_that("the Mahalanobis check's 'plotted' total is the plotted measure, squared", {
  res <- plot_quiet(bp_bnd)
  chk <- res$mahalanobis_check
  expect_type(chk, "list")
  expect_equal(chk$plotted, sum(res$vec_to_boundary_sd^2, na.rm = TRUE))

  W     <- bp_res$metric
  v     <- sweep(res$vec_to_boundary_sd, 2L, sqrt(diag(W)), "*")
  full  <- sum(rowSums((v %*% bp_res$metric_inv) * v), na.rm = TRUE)
  expect_equal(chk$full, full)
  expect_equal(chk$difference, chk$full - chk$plotted)
  expect_equal(chk$pct, 100 * chk$difference / chk$full)
})

test_that("with a diagonal W the plot leaves nothing out", {
  bnd <- bp_bnd
  W   <- bnd$bl_result$metric
  bnd$bl_result$metric     <- diag(diag(W))
  bnd$bl_result$metric_inv <- diag(1 / diag(W))
  chk <- plot_quiet(bnd)$mahalanobis_check
  expect_equal(chk$difference, 0, tolerance = 1e-8 * chk$full)
})

test_that("distance = 'euclidean' gives no Mahalanobis check", {
  expect_null(plot_quiet(bp_bnd, distance = "euclidean")$mahalanobis_check)
})
