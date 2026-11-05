############################################################
# bl_assemble(): Phase 1 result object
# Assembles all Phase 1 artifacts into a single reproducible bl_result.
############################################################

#' Assemble Phase 1 artifacts into a single result object
#'
#' Combines the outputs of `bl_prepare_data()` (or `bl_wrap_data()`),
#' `bl_fit_model()` (or `bl_wrap_model()`), `bl_build_projection()`, and
#' `bl_build_grid()` into a single `bl_result` S3 object. This object is the
#' reproducible anchor for the Phase 2 (global: counterfactual search,
#' surrogate, robustness) and Phase 3 (local interpretation) operations.
#'
#' @section Phase 1 workflow:
#' ```r
#' # Step 1: prepare data (split + outlier filter via hull_fraction)
#' bl_dat  <- bl_prepare_data(datasets::iris, class_col = "Species",
#'                             target_class = "versicolor",
#'                             hull_fraction = 0.9)
#'
#' # Step 2: fit the model on the filtered training data
#' bl_mod  <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
#'                          model_type = "GLM")
#'
#' # Step 3: build projection
#' bl_proj <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
#'                                 method = "PCA")
#'
#' # Step 4: build prediction grid
#' bl_grid <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod)
#'
#' # Step 5: assemble
#' bl_results <- bl_assemble(bl_dat, bl_mod, bl_proj, bl_grid)
#' ```
#' [bl_build_result()] runs Steps 3-5 in one call.
#'
#' @param bl_data   A `"bl_filter_result"` from `bl_prepare_data()` (which
#'   applies the outlier filter), or a `"bl_data"` from `bl_wrap_data()`
#'   (optionally passed through `bl_filter_outliers()`, which also returns a
#'   `"bl_filter_result"`). The function accepts either class and extracts the
#'   training and test data.
#' @param bl_model          A `"bl_model"` from `bl_fit_model()`, or `NULL`.
#'   When `NULL`, model-derived fields (`model`, `model_type`, `cutoff`,
#'   `accuracy`, `gini`) will be `NULL` in the returned object.
#' @param bl_projection     A `"bl_projection"` from `bl_build_projection()`.
#' @param bl_grid           A `"bl_grid"` from `bl_build_grid()`. Its
#'   `cutoff` must equal `bl_model$cutoff`; otherwise the function stops.
#'
#' @return A list of class `"bl_result"` containing all Phase 1 artifacts.
#'   See **Value** section for the full field list.
#'
#' @section Value (fields of `bl_result`):
#' \describe{
#'   \item{`train_data`}{Post-filter training data frame.}
#'   \item{`test_data`}{Test data frame, passed through unfiltered. Set to
#'     `train_data` when there is no test data (`train_fraction = 1`).}
#'   \item{`var_names`}{Feature column names.}
#'   \item{`num_vars`}{Number of features.}
#'   \item{`model`}{Fitted model object.}
#'   \item{`model_type`}{Character; model family.}
#'   \item{`cutoff`}{Decision threshold.}
#'   \item{`b_margin`}{Decision boundary contour band half-width.}
#'   \item{`V`}{Loading matrix (p x p).}
#'   \item{`tV`}{Inverse of V (p x p).}
#'   \item{`X_center`}{Column means.}
#'   \item{`X_sd`}{Column standard deviations.}
#'   \item{`method`}{`"PCA"` or `"CVA"`.}
#'   \item{`standardise`}{Effective standardise flag.}
#'   \item{`proj_dims`}{Integer vector of length 2; eigenvector indices used
#'     for the 2D projection plane.}
#'   \item{`biplot_obj`}{The biplotEZ object for Phase 2 plotting.}
#'   \item{`metric`}{Within-class covariance matrix `W` (p x p), or total
#'     covariance `Sigma` when no class info was available. Used as the
#'     Mahalanobis metric in `bl_find_local_cf()` and `plot.bl_boundary()`
#'     (their default under CVA; under PCA they default to Euclidean
#'     distance).}
#'   \item{`metric_inv`}{Inverse of `metric` via Cholesky (p x p).}
#'   \item{`metric_type`}{`"W_cva"`, `"W_binary"`, or `"Sigma"`.}
#'   \item{`polygon`}{Convex hull of the training data in the final biplot
#'     plane, from `bl_grid` (`SpatialPolygons`), or `NULL` without a grid.
#'     This is the grid hull, built with `outlie` (default `1`, every
#'     training point inside); it is not the outlier-filter hull of
#'     `bl_prepare_data()`.}
#'   \item{`hull_fraction`}{The grid hull's `outlie` value, or `NULL`.}
#'   \item{`train_ranges`}{Per-feature `c(min, max)` ranges of the training
#'     data, from `bl_projection`; the feasibility check in the
#'     counterfactual searches.}
#'   \item{`scaling`}{Display-only raw-unit scaling from [bl_set_scaling()]
#'     (`list(center, scale, method)`), or `NULL`. Used by the biplot plot
#'     methods to relabel axis ticks into original units; does not affect any
#'     computation. Distinct from `X_center`/`X_sd` (internal PCA standardisation)
#'     and the `standardise` flag.}
#'   \item{`biplot_grid`}{The full `bl_grid` list.}
#'   \item{`accuracy`}{Training accuracy.}
#'   \item{`gini`}{Training Gini coefficient.}
#'   \item{`call`}{The R call that created this object.}
#'   \item{`created_at`}{`POSIXct` creation timestamp.}
#' }
#'
#' @export
bl_assemble <- function(bl_data,
                        bl_model      = NULL,
                        bl_projection,
                        bl_grid       = NULL) {

  # ---- Validate inputs --------------------------------------------------
  if (!inherits(bl_data, c("bl_data", "bl_filter_result")))
    stop("'bl_data' must be a 'bl_filter_result' object from ",
         "bl_prepare_data() or bl_filter_outliers(), or a 'bl_data' object ",
         "from bl_wrap_data().",
         call. = FALSE)
  if (!is.null(bl_model) && !inherits(bl_model, "bl_model"))
    stop("'bl_model' must be a 'bl_model' object from bl_fit_model().",
         call. = FALSE)
  if (!inherits(bl_projection, "bl_projection"))
    stop("'bl_projection' must be a 'bl_projection' object from bl_build_projection().",
         call. = FALSE)
  if (!is.null(bl_grid) && !inherits(bl_grid, "bl_grid"))
    stop("'bl_grid' must be a 'bl_grid' object from bl_build_grid().",
         call. = FALSE)

  # The grid's contours and every downstream function must use one cutoff.
  # A bl_grid built before the cutoff field existed has none; skip the check.
  if (!is.null(bl_grid) && !is.null(bl_model) && !is.null(bl_grid$cutoff) &&
      !isTRUE(all.equal(bl_grid$cutoff, bl_model$cutoff)))
    stop(sprintf(paste0(
      "The grid was built with cutoff = %g but the model's cutoff is %g. ",
      "Rebuild the grid with bl_build_grid(cutoff = bl_model$cutoff) ",
      "(the default)."), bl_grid$cutoff, bl_model$cutoff),
      call. = FALSE)

  # ---- Resolve filtered vs unfiltered data ------------------------------
  # Accept either bl_data or bl_filter_result — both carry train_data / test_data
  train_data <- bl_data$train_data
  test_data  <- bl_data$test_data

  # When train_fraction = 1 (all data used for training), test_data is a
  # 0-row data frame. Fall back to train_data so all downstream functions
  # (bl_project_points, bl_predict, bl_find_boundary, etc.) work correctly.
  if (!is.null(test_data) && nrow(test_data) == 0L) {
    message("No test data found (train_fraction = 1?). ",
            "bl_result$test_data set to train_data for downstream compatibility.")
    test_data <- train_data
  }

  # train_ranges: always from bl_projection (computed from whatever training
  # data was passed to bl_build_projection, filtered or unfiltered).
  train_ranges <- bl_projection$train_ranges

  # polygon and hull_fraction: always from bl_grid when available, because
  # the grid is computed in the final biplot Z-space (PCA or CVA, possibly
  # rotated). NULL when no model/grid was supplied.
  polygon       <- bl_grid$polygon
  hull_fraction <- bl_grid$hull_fraction

  # var_names: from bl_model when present, else from bl_data
  var_names_use <- if (!is.null(bl_model)) bl_model$var_names else bl_data$var_names

  # CVA always has standardise = FALSE; store the effective value
  standardise_eff <- bl_projection$standardise

  # ---- Assemble ---------------------------------------------------------
  structure(
    list(
      # Data
      train_data    = train_data,
      test_data     = test_data,
      var_names     = var_names_use,
      num_vars      = length(var_names_use),

      # Model (NULL when no model supplied)
      model         = if (!is.null(bl_model)) bl_model$model      else NULL,
      model_type    = if (!is.null(bl_model)) bl_model$model_type else NULL,
      cutoff        = if (!is.null(bl_model)) bl_model$cutoff     else NULL,
      b_margin      = bl_grid$b_margin,

      # Projection
      V             = bl_projection$V,
      tV            = bl_projection$tV,
      X_center      = bl_projection$X_center,
      X_sd          = bl_projection$X_sd,
      method        = bl_projection$method,
      standardise   = standardise_eff,
      proj_dims     = bl_projection$proj_dims,
      biplot_obj    = bl_projection$biplot_obj,

      # Mahalanobis metric (within-class covariance W or Sigma fallback)
      metric        = bl_projection$metric,
      metric_inv    = bl_projection$metric_inv,
      metric_type   = bl_projection$metric_type,

      # Filtering
      polygon       = polygon,
      hull_fraction = hull_fraction,
      train_ranges  = train_ranges,

      # Raw-unit display scaling (NULL unless bl_set_scaling() was used)
      scaling       = bl_data$scaling,

      # Grid (NULL when no model supplied)
      biplot_grid   = bl_grid,

      # Performance (NULL when no model supplied)
      accuracy      = if (!is.null(bl_model)) bl_model$accuracy   else NULL,
      gini          = if (!is.null(bl_model)) bl_model$gini       else NULL,

      # Metadata
      call          = match.call(),
      created_at    = Sys.time()
    ),
    class = c("bl_result", "list")
  )
}


# --------------------------------------------------------------------------
# bl_build_result(): single-call convenience wrapper for Steps 4-6
# --------------------------------------------------------------------------

#' Build and assemble a full boundary logic result in one step
#'
#' Convenience wrapper that combines [bl_build_projection()],
#' [bl_build_grid()], and [bl_assemble()] into a single call. Accepts the
#' output of [bl_prepare_data()], [bl_wrap_data()] or [bl_filter_outliers()]
#' and extracts `train_data` and `var_names` automatically.
#'
#' When `bl_model` is `NULL` the grid step is skipped and a model-free
#' `bl_result` is returned: its `model`, `cutoff` and `biplot_grid` fields are
#' `NULL`, and [plot()] on it renders the biplot coloured by true class,
#' without a prediction surface or decision boundary.
#'
#' **Assumption:** the projection matrix V is always built from the
#' `train_data` slot of the supplied object. Building V from test data or any
#' other source is not supported here -- use [bl_build_projection()] directly
#' if a non-standard data source is needed.
#'
#' The grid's decision boundary contours are extracted around
#' `bl_model$cutoff`, the same threshold every downstream function uses.
#'
#' For advanced control over any individual step use [bl_build_projection()],
#' [bl_build_grid()], and [bl_assemble()] directly.
#'
#' @section Typical usage (with model):
#' ```r
#' # bl_dat from bl_prepare_data() (outlier filter built in) or bl_wrap_data()
#' bl_results <- bl_build_result(bl_dat, bl_mod, method = "CVA")
#' ```
#'
#' @section Exploratory biplot (no model yet):
#' ```r
#' bl_results <- bl_build_result(bl_dat, method = "PCA")
#' plot(bl_results)   # biplot coloured by true class (red/blue), no surface
#' ```
#'
#' @param bl_data       A `"bl_filter_result"` from [bl_prepare_data()] or
#'   [bl_filter_outliers()], or a `"bl_data"` from [bl_wrap_data()]. The
#'   function extracts `train_data` and `var_names` from whichever is supplied.
#' @param bl_model      A `"bl_model"` object from [bl_fit_model()] or
#'   [bl_wrap_model()]. When `NULL` (default), the projection is built and a
#'   model-free `bl_result` (no grid) is returned.
#' @param method        `"PCA"` or `"CVA"`. Default `"PCA"`.
#' @param proj_dims     Integer vector of length 2; eigenvector indices for
#'   the biplot plane. Default `c(1L, 2L)`.
#' @param standardise   Logical; standardise features before PCA. Ignored for
#'   CVA (always `FALSE`). Default `TRUE`.
#' @param cva_classes   Optional factor of class labels (one per training
#'   row) for CVA. When `NULL` and `bl_model` is supplied, four-class
#'   confusion labels (TP/TN/FP/FN) are used. When `NULL` and no model is
#'   supplied, binary `as.factor(train_data$class)` is used automatically.
#' @param title         Character; biplot title. Default `""`.
#' @param m             Integer; grid resolution (m x m). Default `200L`.
#'   Ignored when `bl_model = NULL`.
#' @param outlie        Numeric in (0, 1]; convex hull fraction for the
#'   prediction grid. Default `1`. Ignored when `bl_model = NULL`.
#' @param calc_hull     Logical; trim the prediction surface to the training
#'   hull. Visual only -- does not affect counterfactual search. Default
#'   `TRUE`. Ignored when `bl_model = NULL`.
#' @param b_margin      Numeric scalar; half-width of the band around
#'   `cutoff` used to extract the decision boundary contour. Default
#'   `0.001` (contours at `cutoff +/- 0.001`). Must satisfy
#'   `0 < b_margin < 0.5`. Independent of prediction precision: all
#'   prediction scores are always floor-rounded to 3 decimal places inside
#'   `.pred_function()` regardless of `b_margin`.
#'   Ignored when `bl_model = NULL`.
#'
#' @return A `"bl_result"` object. When `bl_model = NULL` it is model-free:
#'   `model`, `model_type`, `cutoff`, `b_margin`, `biplot_grid`, `polygon`,
#'   `accuracy` and `gini` are `NULL`.
#'
#' @seealso [bl_build_projection()], [bl_build_grid()], [bl_assemble()]
#'
#' @export
bl_build_result <- function(bl_data     = NULL,
                             bl_model    = NULL,
                             method      = "PCA",
                             proj_dims   = c(1L, 2L),
                             standardise = TRUE,
                             cva_classes = NULL,
                             title       = "",
                             m           = 200L,
                             outlie      = 1,
                             calc_hull   = TRUE,
                             b_margin    = 0.001) {

  # ---- Safely resolve bl_model (lazy promise) ---------------------------
  # If the caller wrote bl_model = bl_mod and bl_mod does not yet exist,
  # R would throw "object 'bl_mod' not found" at the first is.null() touch.
  # force() triggers evaluation here so we can catch that error cleanly and
  # fall through to the standard NULL handling below.
  bl_model <- tryCatch(force(bl_model), error = function(e) NULL)

  # ---- Validate data source ---------------------------------------------
  if (is.null(bl_data))
    stop("bl_build_result() requires 'bl_data'. ",
         "Supply a 'bl_filter_result' object from bl_prepare_data() or ",
         "bl_filter_outliers(), or a 'bl_data' object from bl_wrap_data().",
         call. = FALSE)

  if (!inherits(bl_data, c("bl_data", "bl_filter_result")))
    stop("bl_build_result(): 'bl_data' must be a 'bl_filter_result' object ",
         "from bl_prepare_data() or bl_filter_outliers(), or a 'bl_data' ",
         "object from bl_wrap_data().",
         call. = FALSE)

  # The projection matrix V is always built from the training data.
  # bl_build_result() does not support building V from test data — use
  # bl_build_projection() directly if a non-standard data source is needed.
  train_data <- bl_data$train_data
  var_names  <- bl_data$var_names

  # ---- No model: inform user --------------------------------------------
  if (is.null(bl_model)) {
    message("Note: 'bl_model' is NULL. ",
            "Returning a model-free bl_result. ",
            "plot_biplotEZ() will show the biplot coloured by true class only. ",
            "No prediction grid or decision boundary will be available.")
  }

  # ---- CVA without model: fall back to binary class labels --------------
  if (method == "CVA" && is.null(bl_model) && is.null(cva_classes)) {
    if ("class" %in% names(train_data)) {
      cva_classes <- as.factor(train_data$class)
    }
  }

  # ---- Step 4: build projection -----------------------------------------
  bl_proj <- bl_build_projection(
    train_data  = train_data,
    var_names   = var_names,
    method      = method,
    standardise = standardise,
    proj_dims   = proj_dims,
    bl_model    = bl_model,
    cva_classes = cva_classes,
    title       = title
  )

  # Copy any user-set raw-unit scaling onto the projection so the model-free
  # exploratory plot.bl_projection() can relabel its axes too.
  bl_proj$scaling <- bl_data$scaling

  # ---- Step 5: build prediction grid (skipped when no model) -----------
  bl_grid <- if (!is.null(bl_model)) {
    bl_build_grid(
      train_data    = train_data,
      bl_projection = bl_proj,
      bl_model      = bl_model,
      m             = m,
      outlie        = outlie,
      calc_hull     = calc_hull,
      b_margin      = b_margin
    )
  } else {
    NULL
  }

  # ---- Step 6: assemble and return --------------------------------------
  bl_assemble(
    bl_data       = bl_data,
    bl_model      = bl_model,
    bl_projection = bl_proj,
    bl_grid       = bl_grid
  )
}


# --------------------------------------------------------------------------
# S3 methods: print and summary
# --------------------------------------------------------------------------

#' @export
print.bl_result <- function(x, ...) {
  cat("<bl_result>\n")
  cat(sprintf("  Model type   : %s\n",
              if (!is.null(x$model_type)) x$model_type else "(none)"))
  cat(sprintf("  Method       : %s\n", x$method))
  cat(sprintf("  Features     : %d (%s)\n",
              x$num_vars, paste(x$var_names, collapse = ", ")))
  cat(sprintf("  Train rows   : %d  |  Test rows: %d\n",
              nrow(x$train_data), nrow(x$test_data)))
  if (!is.null(x$model)) {
    cat(sprintf("  Cutoff       : %g\n", x$cutoff))
    cat(sprintf("  Accuracy     : %.4f  |  Gini: %.4f\n",
                x$accuracy, x$gini))
    m <- length(x$biplot_grid$xseq)
    cat(sprintf("  Grid size    : %d x %d\n", m, m))
    # hull_fraction here is the grid hull's outlie, not the outlier filter's
    if (!is.null(x$hull_fraction)) {
      cat(sprintf("  Grid hull    : outlie = %.2f\n", x$hull_fraction))
    } else {
      cat("  Grid hull    : none\n")
    }
  } else {
    cat("  Model        : none (exploratory biplot only)\n")
  }
  if (!is.null(x$scaling))
    cat(sprintf("  Axis units   : raw (%s scaling stored)\n", x$scaling$method))
  cat(sprintf("  Created      : %s\n",
              format(x$created_at, "%Y-%m-%d %H:%M:%S")))
  invisible(x)
}


#' @export
summary.bl_result <- function(object, ...) {
  x <- object
  print(x)
  cat("\n--- Model performance (training data) ---\n")

  pred_prob <- .pred_function(
    model_use  = x$model,
    model_type = x$model_type,
    new_data   = x$train_data[, x$var_names, drop = FALSE]
  )
  pred_class <- as.numeric(pred_prob >= x$cutoff)
  actual     <- x$train_data[["class"]]

  cm <- table(Predicted = pred_class, Actual = actual)
  cat("Confusion matrix (train):\n")
  print(cm)

  cat(sprintf("\nGrid probability range: [%.4f, %.4f]\n",
              min(x$biplot_grid$grid_prob, na.rm = TRUE),
              max(x$biplot_grid$grid_prob, na.rm = TRUE)))
  cat(sprintf("Contour lines found  : %d\n", length(x$biplot_grid$ct)))

  # hull_fraction here is the grid hull's outlie, not the outlier filter's
  if (!is.null(x$hull_fraction)) {
    cat(sprintf("Grid hull            : outlie = %.2f (contours clipped to it)\n",
                x$hull_fraction))
  } else {
    cat("Grid hull            : none\n")
  }

  if (x$method == "CVA" && !isTRUE(x$standardise)) {
    cat("Note: CVA forced standardise = FALSE.\n")
  }

  invisible(x)
}
