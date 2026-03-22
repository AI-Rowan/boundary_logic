############################################################
# bl_assemble(): Phase 1 result object
# Assembles all Phase 1 artifacts into a single reproducible bl_result.
############################################################

#' Assemble Phase 1 artifacts into a single result object
#'
#' Combines the outputs of `bl_prepare_data()`, optionally
#' `bl_filter_outliers()`, `bl_fit_model()`, `bl_build_projection()`, and
#' `bl_build_grid()` into a single `bl_result` S3 object. This object is the
#' reproducible anchor for all Phase 2 operations (counterfactual search,
#' local interpretation, global interpretation).
#'
#' @section Phase 1 workflow:
#' ```r
#' # Step 1: prepare data
#' bl_dat  <- bl_prepare_data(datasets::iris, class_col = "Species",
#'                             target_class = "versicolor")
#'
#' # Step 2 (optional): iterative polygon filter
#' bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 0.9)
#'
#' # Step 3: fit final model (on filtered data)
#' bl_mod  <- bl_fit_model(bl_filt$train_data, bl_filt$var_names,
#'                          model_type = "GLM")
#'
#' # Step 4: build projection
#' bl_proj <- bl_build_projection(bl_filt$train_data, bl_filt$var_names,
#'                                 method = "PCA")
#'
#' # Step 5: build prediction grid
#' bl_grid <- bl_build_grid(bl_filt$train_data, bl_proj, bl_mod)
#'
#' # Step 6: assemble — pass bl_filt if filtering was used, bl_dat otherwise
#' bl_results <- bl_assemble(bl_filt, bl_mod, bl_proj, bl_grid)
#' ```
#'
#' @param bl_data   Either a `"bl_data"` object from `bl_prepare_data()` (when
#'   no polygon filtering was applied), or a `"bl_filter_result"` object from
#'   `bl_filter_outliers()` (when filtering was applied). Pass whichever is the
#'   last step in your data-preparation chain — the function detects the type
#'   automatically and extracts the correct training and test data.
#' @param bl_model          A `"bl_model"` from `bl_fit_model()`.
#' @param bl_projection     A `"bl_projection"` from `bl_build_projection()`.
#' @param bl_grid           A `"bl_grid"` from `bl_build_grid()`.
#'
#' @return A list of class `"bl_result"` containing all Phase 1 artifacts.
#'   See **Value** section for the full field list.
#'
#' @section Value (fields of `bl_result`):
#' \describe{
#'   \item{`train_data`}{Post-filter training data frame.}
#'   \item{`test_data`}{Test data frame (filtered to training variable ranges).}
#'   \item{`var_names`}{Feature column names.}
#'   \item{`num_vars`}{Number of features.}
#'   \item{`model`}{Fitted model object.}
#'   \item{`model_type`}{Character; model family.}
#'   \item{`cutoff`}{Decision threshold.}
#'   \item{`rounding`}{Rounding precision.}
#'   \item{`V`}{Loading matrix (p x p).}
#'   \item{`tV`}{Inverse of V (p x p).}
#'   \item{`X_center`}{Column means.}
#'   \item{`X_sd`}{Column standard deviations.}
#'   \item{`method`}{`"PCA"` or `"CVA"`.}
#'   \item{`standardise`}{Effective standardise flag.}
#'   \item{`proj_dims`}{Integer vector of length 2; eigenvector indices used
#'     for the 2D projection plane.}
#'   \item{`biplot_obj`}{The biplotEZ object for Phase 2 plotting.}
#'   \item{`polygon`}{Convex hull polygon (`SpatialPolygons` or `NULL`).}
#'   \item{`hull_fraction`}{Hull fraction used, or `NULL`.}
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
    stop("'bl_data' must be a 'bl_data' object from bl_prepare_data() or ",
         "a 'bl_filter_result' object from bl_filter_outliers().",
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
      model         = bl_model$model,
      model_type    = bl_model$model_type,
      cutoff        = bl_model$cutoff,
      rounding      = bl_grid$rounding,

      # Projection
      V             = bl_projection$V,
      tV            = bl_projection$tV,
      X_center      = bl_projection$X_center,
      X_sd          = bl_projection$X_sd,
      method        = bl_projection$method,
      standardise   = standardise_eff,
      proj_dims     = bl_projection$proj_dims,
      biplot_obj    = bl_projection$biplot_obj,

      # Filtering
      polygon       = polygon,
      hull_fraction = hull_fraction,
      train_ranges  = train_ranges,

      # Grid (NULL when no model supplied)
      biplot_grid   = bl_grid,

      # Performance (NULL when no model supplied)
      accuracy      = bl_model$accuracy,
      gini          = bl_model$gini,

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
#' output of either [bl_prepare_data()] or [bl_filter_outliers()] and
#' extracts `train_data` and `var_names` automatically.
#'
#' When `bl_model` is `NULL` the grid and assembly steps are skipped. A
#' [`bl_projection`][bl_build_projection()] object is returned instead of a
#' `bl_result`, and [plot()] on that object renders the biplot without any
#' prediction surface.
#'
#' **Assumption:** the projection matrix V is always built from the
#' `train_data` slot of the supplied object. Building V from test data or any
#' other source is not supported here — use [bl_build_projection()] directly
#' if a non-standard data source is needed.
#'
#' For advanced control over any individual step use [bl_build_projection()],
#' [bl_build_grid()], and [bl_assemble()] directly.
#'
#' @section Typical usage (with model):
#' ```r
#' # With outlier filtering:
#' bl_results <- bl_build_result(bl_filt, bl_mod, method = "CVA")
#'
#' # Without filtering:
#' bl_results <- bl_build_result(bl_dat, bl_mod, method = "PCA")
#' ```
#'
#' @section Exploratory biplot (no model yet):
#' ```r
#' bl_proj <- bl_build_result(bl_filt, method = "PCA")
#' plot(bl_proj)   # renders biplot coloured by true class (red/blue)
#' ```
#'
#' @param bl_data       A `"bl_data"` object from [bl_prepare_data()] or a
#'   `"bl_filter_result"` from [bl_filter_outliers()]. The function extracts
#'   `train_data` and `var_names` from whichever is supplied.
#' @param bl_model      A `"bl_model"` object from [bl_fit_model()] or
#'   [bl_wrap_model()]. When `NULL` (default), only the projection is built
#'   and a `bl_projection` object is returned.
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
#'   hull. Visual only — does not affect counterfactual search. Default
#'   `TRUE`. Ignored when `bl_model = NULL`.
#' @param rounding      Integer; controls the decision boundary contour band
#'   width only — `b_margin = 1 / (10^rounding)`. Default `3L` gives contours
#'   at `cutoff ± 0.001`. All prediction scores are always floor-rounded to
#'   3 decimal places regardless of this setting.
#'   Ignored when `bl_model = NULL`.
#'
#' @return A `"bl_result"` object (when `bl_model` is supplied) or a
#'   `"bl_projection"` object (when `bl_model = NULL`).
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
                             rounding    = 3L) {

  # ---- Validate data source ---------------------------------------------
  if (is.null(bl_data)) {
    message("bl_build_result() requires 'bl_data'. ",
            "Supply a 'bl_data' object from bl_prepare_data() or ",
            "a 'bl_filter_result' object from bl_filter_outliers().")
    return(invisible(NULL))
  }

  if (!inherits(bl_data, c("bl_data", "bl_filter_result"))) {
    message("bl_build_result(): 'bl_data' must be a 'bl_data' object from ",
            "bl_prepare_data() or a 'bl_filter_result' object from ",
            "bl_filter_outliers().")
    return(invisible(NULL))
  }

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

  # ---- Step 5: build prediction grid (skipped when no model) -----------
  bl_grid <- if (!is.null(bl_model)) {
    bl_build_grid(
      train_data    = train_data,
      bl_projection = bl_proj,
      bl_model      = bl_model,
      m             = m,
      outlie        = outlie,
      calc_hull     = calc_hull,
      rounding      = rounding
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
    if (!is.null(x$hull_fraction)) {
      cat(sprintf("  Hull fraction: %.2f\n", x$hull_fraction))
    } else {
      cat("  Hull fraction: none applied\n")
    }
  } else {
    cat("  Model        : none (exploratory biplot only)\n")
  }
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

  if (!is.null(x$hull_fraction)) {
    cat(sprintf("Polygon filter       : applied (fraction = %.2f)\n",
                x$hull_fraction))
  } else {
    cat("Polygon filter       : not applied\n")
  }

  if (x$method == "CVA" && !isTRUE(x$standardise)) {
    cat("Note: CVA forced standardise = FALSE.\n")
  }

  invisible(x)
}
