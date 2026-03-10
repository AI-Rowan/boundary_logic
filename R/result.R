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
#' bl_dat  <- bl_prepare_data(iris, class_col = "Species",
#'                             target_class = "versicolor")
#'
#' # Step 2 (optional): iterative polygon filter
#' bl_filt <- bl_filter_outliers(bl_dat, model_type = "GLM",
#'                                hull_fraction = 0.9)
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
#' # Step 6: assemble
#' result  <- bl_assemble(bl_dat, bl_filt, bl_mod, bl_proj, bl_grid)
#' ```
#'
#' @param bl_data           A `"bl_data"` object from `bl_prepare_data()`.
#' @param bl_filter_result  A `"bl_filter_result"` from `bl_filter_outliers()`,
#'   or `NULL` if no polygon filtering was applied.
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
                        bl_filter_result = NULL,
                        bl_model,
                        bl_projection,
                        bl_grid) {

  # ---- Validate inputs --------------------------------------------------
  if (!inherits(bl_data, "bl_data"))
    stop("'bl_data' must be a 'bl_data' object from bl_prepare_data().",
         call. = FALSE)
  if (!is.null(bl_filter_result) && !inherits(bl_filter_result, "bl_filter_result"))
    stop("'bl_filter_result' must be a 'bl_filter_result' object or NULL.",
         call. = FALSE)
  if (!inherits(bl_model, "bl_model"))
    stop("'bl_model' must be a 'bl_model' object from bl_fit_model().",
         call. = FALSE)
  if (!inherits(bl_projection, "bl_projection"))
    stop("'bl_projection' must be a 'bl_projection' object from bl_build_projection().",
         call. = FALSE)
  if (!inherits(bl_grid, "bl_grid"))
    stop("'bl_grid' must be a 'bl_grid' object from bl_build_grid().",
         call. = FALSE)

  # ---- Resolve filtered vs unfiltered data ------------------------------
  if (!is.null(bl_filter_result)) {
    train_data    <- bl_filter_result$train_data
    test_data     <- bl_filter_result$test_data
    polygon       <- bl_filter_result$polygon
    hull_fraction <- bl_filter_result$hull_fraction
  } else {
    train_data    <- bl_data$train_data
    test_data     <- bl_data$test_data
    polygon       <- bl_grid$polygon   # polygon from the grid computation
    hull_fraction <- NULL
  }

  # CVA always has standardise = FALSE; store the effective value
  standardise_eff <- bl_projection$standardise

  # ---- Assemble ---------------------------------------------------------
  structure(
    list(
      # Data
      train_data    = train_data,
      test_data     = test_data,
      var_names     = bl_model$var_names,
      num_vars      = length(bl_model$var_names),

      # Model
      model         = bl_model$model,
      model_type    = bl_model$model_type,
      cutoff        = bl_model$cutoff,
      rounding      = bl_model$rounding,

      # Projection
      V             = bl_projection$V,
      tV            = bl_projection$tV,
      X_center      = bl_projection$X_center,
      X_sd          = bl_projection$X_sd,
      method        = bl_projection$method,
      standardise   = standardise_eff,
      biplot_obj    = bl_projection$biplot_obj,

      # Filtering
      polygon       = polygon,
      hull_fraction = hull_fraction,

      # Grid
      biplot_grid   = bl_grid,

      # Performance
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
# S3 methods: print and summary
# --------------------------------------------------------------------------

#' @export
print.bl_result <- function(x, ...) {
  m <- length(x$biplot_grid$xseq)
  cat("<bl_result>\n")
  cat(sprintf("  Model type   : %s\n", x$model_type))
  cat(sprintf("  Method       : %s\n", x$method))
  cat(sprintf("  Features     : %d (%s)\n",
              x$num_vars, paste(x$var_names, collapse = ", ")))
  cat(sprintf("  Train rows   : %d  |  Test rows: %d\n",
              nrow(x$train_data), nrow(x$test_data)))
  cat(sprintf("  Cutoff       : %g\n", x$cutoff))
  cat(sprintf("  Accuracy     : %.4f  |  Gini: %.4f\n",
              x$accuracy, x$gini))
  cat(sprintf("  Grid size    : %d x %d\n", m, m))
  if (!is.null(x$hull_fraction)) {
    cat(sprintf("  Hull fraction: %.2f\n", x$hull_fraction))
  } else {
    cat("  Hull fraction: none applied\n")
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
    rounding   = x$rounding,
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
