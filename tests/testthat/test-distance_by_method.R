# The distance follows the thesis rule by method (section 3.9.1): Euclidean
# distance in X-space under PCA, Mahalanobis distance under CVA. distance =
# "auto" (the default) resolves it, in bl_find_local_cf()'s pair search and
# integer rounding, and in plot.bl_boundary() / bl_robustness().

dm_dat <- bl_prepare_data(iris, class_col = "Species",
                          target_class = "versicolor", verbose = FALSE)
dm_mod <- bl_fit_model(dm_dat$train_data, dm_dat$var_names, model_type = "GLM")
dm_pca <- suppressWarnings(suppressMessages(
  bl_build_result(dm_dat, dm_mod, method = "PCA", m = 60L, b_margin = 0.01)))
dm_cva <- suppressWarnings(suppressMessages(
  bl_build_result(dm_dat, dm_mod, method = "CVA", m = 60L, b_margin = 0.01)))

# The first test row whose search finds a counterfactual.
local_cf <- function(bl_results, ...) {
  for (i in seq_len(nrow(bl_results$test_data))) {
    tgt <- bl_select_target(bl_results, target = i)
    loc <- suppressMessages(bl_find_local_cf(bl_results, tgt, verbose = FALSE, ...))
    if (isTRUE(loc$solution_found)) return(loc)
  }
  NULL
}

plot_quiet <- function(x, ...) {
  on.exit(grDevices::graphics.off(), add = TRUE)
  grDevices::pdf(NULL)
  out <- NULL
  invisible(capture.output(out <- plot(x, ...)))
  out
}

test_that(".resolve_distance() applies the thesis rule for 'auto'", {
  resolve <- boundarylogic:::.resolve_distance
  expect_identical(resolve("auto", list(method = "PCA")), "euclidean")
  expect_identical(resolve("auto", list(method = "CVA")), "mahalanobis")
  expect_identical(resolve("auto", list()), "mahalanobis")   # older bl_result
  expect_identical(resolve("mahalanobis", list(method = "PCA")), "mahalanobis")
  expect_identical(resolve("euclidean", list(method = "CVA")), "euclidean")
  expect_error(resolve("manhattan", list(method = "PCA")))
})

test_that("'auto' is the default in every distance-measure function", {
  for (f in list(bl_find_local_cf, plot.bl_boundary, bl_robustness))
    expect_identical(eval(formals(f)$distance)[1L], "auto")
})

test_that(".cf_sq_distance() gives Mahalanobis and total-SD Euclidean distances", {
  v    <- matrix(c(1, 2, -1, 0.5), nrow = 2L, byrow = TRUE)
  Winv <- matrix(c(2, 0.5, 0.5, 1), nrow = 2L)
  sdv  <- c(2, 4)
  d2   <- boundarylogic:::.cf_sq_distance
  expect_equal(d2(v, "mahalanobis", Winv, sdv),
               c(drop(v[1, ] %*% Winv %*% v[1, ]), drop(v[2, ] %*% Winv %*% v[2, ])))
  expect_equal(d2(v, "euclidean", Winv, sdv),
               c((1 / 2)^2 + (2 / 4)^2, (-1 / 2)^2 + (0.5 / 4)^2))
})

test_that("under PCA the default selects pairs by Euclidean distance in X-space", {
  loc <- local_cf(dm_pca)
  skip_if(is.null(loc), "no counterfactual in the PCA fixture")
  expect_identical(loc$distance, "euclidean")
  ok <- is.finite(loc$all_distances_x)
  expect_identical(names(which.min(loc$all_distances_x[ok])),
                   sprintf("(%d,%d)", loc$best_pair[1L], loc$best_pair[2L]))
  # dist_x is the total-SD Euclidean distance of the reported counterfactual
  # (no constraints, so it is the search's point) ...
  vn <- dm_pca$var_names
  v  <- as.numeric(loc$B_x[1L, vn]) - as.numeric(loc$bl_target$x_obs[1L, vn])
  expect_equal(loc$dist_x, sqrt(sum((v / dm_pca$X_sd)^2)))
  # ... and, under a standardised PCA, equals the plane distance for every pair.
  expect_equal(unname(loc$all_distances_x[ok]), unname(loc$all_distances[ok]),
               tolerance = 1e-6)
})

test_that("under PCA an explicit 'mahalanobis' selects by the Mahalanobis distance", {
  loc <- local_cf(dm_pca, distance = "mahalanobis")
  skip_if(is.null(loc), "no counterfactual in the PCA fixture")
  expect_identical(loc$distance, "mahalanobis")
  ok <- is.finite(loc$all_distances_mahalanobis)
  expect_identical(names(which.min(loc$all_distances_mahalanobis[ok])),
                   sprintf("(%d,%d)", loc$best_pair[1L], loc$best_pair[2L]))
})

test_that("under CVA the default is the Mahalanobis search, unchanged", {
  tgt  <- bl_select_target(dm_cva, target = 1L)
  auto <- suppressMessages(bl_find_local_cf(dm_cva, tgt, verbose = FALSE))
  mah  <- suppressMessages(bl_find_local_cf(dm_cva, tgt, verbose = FALSE,
                                            distance = "mahalanobis"))
  expect_identical(auto$distance, "mahalanobis")
  expect_identical(auto$solution_found, mah$solution_found)
  expect_identical(auto$best_pair, mah$best_pair)
  expect_identical(auto$B_x, mah$B_x)
  expect_identical(auto$B_pred, mah$B_pred)
})

test_that("Phase 2 standardises by the total SD under PCA and by W under CVA", {
  bnd_pca <- suppressWarnings(suppressMessages(bl_find_boundary(dm_pca)))
  res_pca <- plot_quiet(bnd_pca, type = "boxplot")
  expect_identical(res_pca$distance, "euclidean")
  expect_null(res_pca$mahalanobis_check)
  expect_identical(bl_robustness(bnd_pca)$distance, "euclidean")
  expect_equal(bl_robustness(bnd_pca)$sum_of_distance,
               bl_robustness(bnd_pca, distance = "euclidean")$sum_of_distance)

  bnd_cva <- suppressWarnings(suppressMessages(bl_find_boundary(dm_cva)))
  res_cva <- plot_quiet(bnd_cva, type = "boxplot")
  expect_identical(res_cva$distance, "mahalanobis")
  expect_type(res_cva$mahalanobis_check, "list")
  expect_identical(bl_robustness(bnd_cva)$distance, "mahalanobis")
})
