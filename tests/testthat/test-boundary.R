# The boundary searches take the exact nearest contour vertex, so their results
# do not depend on the session's random-number state. .nearest_idx_block()
# used max.col(-d2), whose default ties.method = "random" treats every vertex
# within 1e-5 of the farthest vertex's squared distance as tied and draws one
# at random. See documentation/nearest-vertex-tie-break.md.

nb <- boundarylogic:::.nearest_idx_block

test_that(".nearest_idx_block returns the true nearest inside the old tie band", {
  Z <- matrix(c(0, 0), nrow = 1L)
  # Squared distances 1 + 1e-4, 1 and 2e4: the first two are within the old
  # band (1e-5 * 2e4 = 0.2), but only the second is the nearest.
  B <- rbind(c(0, sqrt(1 + 1e-4)), c(1, 0), c(100, 100))
  picks <- vapply(1:20, function(s) { set.seed(s); nb(Z, B) }, integer(1))
  expect_true(all(picks == 2L))
})

test_that(".nearest_idx_block equals which.min on every row", {
  set.seed(1)
  Z  <- matrix(rnorm(40), ncol = 2L)
  B  <- matrix(rnorm(600), ncol = 2L)
  d2 <- outer(Z[, 1L], B[, 1L], "-")^2 + outer(Z[, 2L], B[, 2L], "-")^2
  expect_identical(nb(Z, B), apply(d2, 1L, which.min))
})

test_that("an exact tie returns the first nearest row", {
  Z <- matrix(c(0, 0), nrow = 1L)
  B <- rbind(c(5, 5), c(0, 1), c(1, 0))
  expect_identical(nb(Z, B), 2L)
})

test_that(".nearest_idx_block uses no random numbers", {
  set.seed(42)
  before <- .Random.seed
  nb(matrix(c(0, 0), nrow = 1L), rbind(c(1, 0), c(0, 1 + 1e-9), c(50, 50)))
  expect_identical(.Random.seed, before)
})


# ---- The searches on the standardised-iris pipeline -----------------------------
# With the random tie-break, eigenvector pair (1,3) of target 1 returned one of
# two counterfactuals depending on the RNG state.

std_iris <- iris
std_iris$class   <- as.numeric(std_iris$Species == "versicolor")
std_iris$Species <- NULL
bnd_feats <- setdiff(names(std_iris), "class")
std_iris[, bnd_feats] <- scale(std_iris[, bnd_feats])

bnd_dat <- bl_wrap_data(std_iris, std_iris, var_names = bnd_feats)
bnd_mod <- bl_fit_model(bnd_dat$train_data, bnd_dat$var_names, model_type = "GLM")
bnd_res <- suppressWarnings(suppressMessages(
  bl_build_result(bnd_dat, bnd_mod, method = "PCA", m = 60L, b_margin = 0.01)))
bnd_tgt <- bl_select_target(bnd_res, target = 1L)

quiet <- function(expr) {
  out <- NULL
  invisible(capture.output(out <- suppressWarnings(suppressMessages(expr))))
  out
}

test_that("bl_find_local_cf gives the same answer under any RNG state", {
  runs <- lapply(1:10, function(s) {
    set.seed(s)
    quiet(bl_find_local_cf(bnd_res, bnd_tgt, verbose = FALSE))
  })
  for (r in runs[-1L]) {
    expect_identical(r$all_distances_mahalanobis, runs[[1L]]$all_distances_mahalanobis)
    expect_identical(r$B_x, runs[[1L]]$B_x)
  }

  # ... and after other searches have run in the same session.
  invisible(quiet(bl_find_boundary(bnd_res)))
  invisible(quiet(bl_surrogate(bnd_res)))
  later <- quiet(bl_find_local_cf(bnd_res, bnd_tgt, verbose = FALSE))
  expect_identical(later$all_distances_mahalanobis, runs[[1L]]$all_distances_mahalanobis)
  expect_identical(later$B_x, runs[[1L]]$B_x)
})

test_that("bl_find_boundary and bl_surrogate give the same answer under any RNG state", {
  bnd  <- lapply(1:5, function(s) { set.seed(s); quiet(bl_find_boundary(bnd_res)) })
  surr <- lapply(1:5, function(s) { set.seed(s); quiet(bl_surrogate(bnd_res)) })
  for (k in 2:5) {
    expect_identical(bnd[[k]]$B_z, bnd[[1L]]$B_z)
    expect_identical(bnd[[k]]$B_x, bnd[[1L]]$B_x)
    expect_identical(surr[[k]]$surrogate_pred, surr[[1L]]$surrogate_pred)
  }
})


# ---- dist_mahalanobis -------------------------------------------------------------

test_that("bl_find_boundary returns the squared Mahalanobis distance v' W^-1 v", {
  bnd <- quiet(bl_find_boundary(bnd_res))
  v   <- as.matrix(bnd$B_x) - as.matrix(bnd$x_obs)
  ok  <- which(!is.na(bnd$B_x[[1L]]))
  expect_gt(length(ok), 0L)
  for (i in head(ok, 5L)) {
    expect_equal(bnd$dist_mahalanobis[i],
                 as.numeric(t(v[i, ]) %*% bnd_res$metric_inv %*% v[i, ]))
  }
  # A row without a counterfactual stays NA.
  expect_identical(is.na(bnd$dist_mahalanobis), is.na(bnd$B_x[[1L]]))
  # The print shows the mean of its square root.
  txt <- paste(capture.output(print(bnd)), collapse = "\n")
  expect_match(txt, sprintf("Dist M (mean)  : %.4f",
                            mean(sqrt(bnd$dist_mahalanobis), na.rm = TRUE)),
               fixed = TRUE)
})

test_that("bl_find_boundary gives NA dist_mahalanobis, with a warning, without metric_inv", {
  old <- bnd_res
  old$metric_inv <- NULL
  expect_warning(bnd <- suppressMessages(bl_find_boundary(old)), "metric_inv not found")
  expect_true(all(is.na(bnd$dist_mahalanobis)))
  expect_no_match(paste(capture.output(print(bnd)), collapse = "\n"), "Dist M")
})
