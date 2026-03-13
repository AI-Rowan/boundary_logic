############################################################
# bl_filter_outliers(): iterative convex hull polygon filter
# Implements the iterative polygon filter step from Phase 1.
# Refactored from: scripts/3.1 Biplot search - Train data.r (Steps 3.1-3.2)
############################################################

#' Remove outliers from training data using a convex hull polygon filter
#'
#' Projects training data into a lightweight 2-PC space, builds a convex hull
#' enclosing `hull_fraction` of the points, removes observations outside the
#' hull, and reports the number and percentage of rows retained and removed.
#'
#' This function implements the polygon filter step of Phase 1. Call it
#' multiple times with different `hull_fraction` values to explore the
#' trade-off between data coverage and outlier removal. Once satisfied, fit
#' the final model using `bl_fit_model()` or `bl_wrap_model()` and then
#' pass the filtered data to `bl_build_projection()` and `bl_build_grid()`.
#'
#' @section How the polygon filter works:
#' The polygon is constructed in the standardised feature space (all variables
#' centred and scaled), not in a biplot Z-space, so that it can be applied
#' before the projection matrix is computed. The convex hull is computed using
#' `aplpack::plothulls()` on the first two principal components of the raw
#' feature matrix, which is a lightweight proxy for the full biplot projection.
#'
#' @param bl_data       A `"bl_data"` object returned by `bl_prepare_data()`
#'   or `bl_wrap_data()`.
#' @param hull_fraction Numeric in (0, 1]; the `fraction` argument for
#'   `aplpack::plothulls()`. Values below `1` trim the most extreme points.
#'   `1.0` retains all points (no trimming). Default `0.9`.
#' @param verbose       Logical; if `TRUE`, prints a one-line summary showing
#'   rows retained and removed with percentages. Default `TRUE`.
#'
#' @return A list of class `"bl_filter_result"` with components:
#' \describe{
#'   \item{`train_data`}{Filtered training data frame.}
#'   \item{`test_data`}{Test data frame (filtered to training variable ranges).}
#'   \item{`var_names`}{Feature column names.}
#'   \item{`num_vars`}{Number of features.}
#'   \item{`target_class`}{The positive class value (passed through).}
#'   \item{`polygon`}{`SpatialPolygons` convex hull in standardised PCA space.}
#'   \item{`hull_fraction`}{The fraction value used.}
#'   \item{`n_retained`}{Number of training rows retained.}
#'   \item{`n_removed`}{Number of training rows removed.}
#' }
#'
#' @examples
#' bl_dat  <- bl_prepare_data(datasets::iris,
#'                             class_col    = "Species",
#'                             target_class = "versicolor")
#' bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 0.9)
#' print(bl_filt)
#'
#' @export
bl_filter_outliers <- function(bl_data,
                               hull_fraction = 0.9,
                               verbose       = TRUE) {

  # ---- Validation -------------------------------------------------------
  if (!inherits(bl_data, "bl_data"))
    stop("'bl_data' must be a 'bl_data' object from bl_prepare_data().",
         call. = FALSE)
  stop_if_not_in_range(hull_fraction, 0, 1.001, "hull_fraction")  # allow == 1
  if (hull_fraction > 1) hull_fraction <- 1.0

  train_data <- bl_data$train_data
  test_data  <- bl_data$test_data
  var_names  <- bl_data$var_names
  num_vars   <- bl_data$num_vars
  n_original <- nrow(train_data)

  # ---- Standardise features for PCA proxy ------------------------------
  X <- as.matrix(train_data[, var_names, drop = FALSE])
  X_center <- colMeans(X)
  X_sd     <- apply(X, 2L, stats::sd)
  X_sd[X_sd == 0] <- 1  # avoid division by zero for constant columns
  X_st     <- scale(X, center = X_center, scale = X_sd)

  # ---- Project to first 2 PCs (lightweight proxy) ----------------------
  svd_x  <- svd(X_st, nu = 2L, nv = 2L)
  Z_train <- X_st %*% svd_x$v[, 1:2, drop = FALSE]
  colnames(Z_train) <- c("x", "y")

  # ---- Build convex hull polygon in Z-space ----------------------------
  min_val <- min(Z_train)
  max_val <- max(Z_train)

  polygon <- .build_hull_polygon(Z_train, min_val, max_val,
                                 outlie = hull_fraction)

  # ---- Identify points inside polygon ----------------------------------
  inside <- .points_in_polygon(Z_train, polygon)

  train_filtered <- train_data[inside, , drop = FALSE]
  rownames(train_filtered) <- NULL
  n_retained <- nrow(train_filtered)
  n_removed  <- n_original - n_retained

  if (n_retained < 10L)
    warning(sprintf(
      "Only %d training rows remain after hull_fraction = %g. ",
      n_retained, hull_fraction
    ), "Consider increasing hull_fraction.", call. = FALSE)

  # ---- Verbose output --------------------------------------------------
  if (isTRUE(verbose)) {
    pct_removed  <- 100 * n_removed  / n_original
    pct_retained <- 100 * n_retained / n_original
    cat(sprintf(
      "Hull fraction: %.2f | Retained: %d (%.1f%%) | Removed: %d (%.1f%%)\n",
      hull_fraction, n_retained, pct_retained, n_removed, pct_removed
    ))
  }

  # ---- Return ----------------------------------------------------------
  structure(
    list(
      train_data    = train_filtered,
      test_data     = test_data,
      var_names     = var_names,
      num_vars      = num_vars,
      target_class  = bl_data$target_class,
      polygon       = polygon,
      hull_fraction = hull_fraction,
      n_retained    = n_retained,
      n_removed     = n_removed
    ),
    class = "bl_filter_result"
  )
}


# --------------------------------------------------------------------------
# S3 print method
# --------------------------------------------------------------------------

#' @export
print.bl_filter_result <- function(x, ...) {
  n_total      <- x$n_retained + x$n_removed
  pct_retained <- 100 * x$n_retained / n_total
  pct_removed  <- 100 * x$n_removed  / n_total
  cat("<bl_filter_result>\n")
  cat(sprintf("  Hull fraction : %.2f\n", x$hull_fraction))
  cat(sprintf("  Retained      : %d (%.1f%%)\n", x$n_retained, pct_retained))
  cat(sprintf("  Removed       : %d (%.1f%%)\n", x$n_removed,  pct_removed))
  invisible(x)
}
