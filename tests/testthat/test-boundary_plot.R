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

test_that("robustness is measured on the in-plane gap, as the plot's totals", {
  # The per-variable distances use the in-plane gap (Z_obs - B_z) V_r^-1,
  # not bl_find_boundary()'s B_x - x_obs, which also holds the observation's
  # off-plane part (documentation/algorithms.md, finding F3).
  mh  <- plot_quiet(bp_bnd)
  tVr <- bp_res$tV[bp_res$proj_dims, , drop = FALSE]
  gap <- (bp_bnd$Z_obs - bp_bnd$B_z) %*% tVr      # CVA: not standardised
  expect_equal(unname(mh$vec_to_boundary_sd),
               unname(sweep(gap, 2L, sqrt(diag(bp_res$metric)), "/")))
  # bl_robustness() gives the same totals as the plot.
  eu <- plot_quiet(bp_bnd, distance = "euclidean")
  expect_equal(bl_robustness(bp_bnd)$sum_of_distance, mh$sum_of_distance)
  expect_equal(bl_robustness(bp_bnd, distance = "euclidean")$sum_of_distance,
               eu$sum_of_distance)
  expect_equal(bl_robustness(bp_bnd)$robustness, mh$robustness)
})
