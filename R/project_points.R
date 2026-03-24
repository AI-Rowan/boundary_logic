############################################################
# bl_project_points(): project any dataset into the biplot Z-space
############################################################

#' Project observations into the biplot Z-space
#'
#' Projects any data frame — training, test, or new observations — into the
#' 2D biplot plane using the projection matrix V stored in a `"bl_result"`
#' object. Also scores each observation through the fitted model and assigns
#' point colours:
#' \itemize{
#'   \item **True labels known** (`"class"` column present): confusion colours
#'     TP (red), TN (blue), FP (purple), FN (orange).
#'   \item **Labels unknown** (no `"class"` column): red for predicted class 1,
#'     blue for predicted class 0.
#' }
#'
#' The projection formula mirrors the one used during Phase 1:
#' `Z = scale(X, center, scale) %*% V[, proj_dims]`, where `scale` is
#' `X_sd` for PCA with standardisation and `FALSE` for CVA.
#'
#' @param data      Data frame containing the feature columns named in
#'   `bl_result$var_names`. May optionally include a `"class"` column
#'   (numeric 0/1) for confusion colouring.
#' @param bl_result A `"bl_result"` object from `bl_assemble()`.
#' @param filter_to_polygon Logical; if `TRUE`, rows whose projected
#'   coordinates fall outside the training-data convex hull polygon are
#'   removed from the returned object. Default `FALSE`.
#' @param filter_to_train_ranges Logical; if `TRUE`, rows where any feature
#'   value falls outside the min/max range of the training data (as seen by
#'   `bl_build_projection()`) are removed before projection. Default `FALSE`.
#'
#' @return A list of class `"bl_points"` with components:
#' \describe{
#'   \item{`Z`}{Numeric matrix (n x 2); projected coordinates in Z-space.}
#'   \item{`pred_prob`}{Numeric vector (length n); model probabilities in
#'     [0, 1].}
#'   \item{`pred_class`}{Numeric 0/1 vector (length n); predicted class.}
#'   \item{`pred_col`}{Character vector (length n); point colour.}
#'   \item{`class`}{Numeric 0/1 vector if `"class"` was in `data`, else
#'     `NULL`.}
#'   \item{`inside_polygon`}{Logical vector (length n); `TRUE` if within the
#'     convex hull. Always `TRUE` when `bl_result$polygon` is `NULL`.}
#' }
#'
#' @examples
#' bl_dat  <- bl_prepare_data(datasets::iris, class_col = "Species",
#'                             target_class = "versicolor")
#' bl_mod  <- bl_fit_model(bl_dat$train_data, bl_dat$var_names)
#' bl_proj <- bl_build_projection(bl_dat$train_data, bl_dat$var_names)
#' bl_grid <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = 50L)
#' result  <- bl_assemble(bl_dat, bl_model = bl_mod,
#'                         bl_projection = bl_proj, bl_grid = bl_grid)
#'
#' # Project training data
#' train_pts <- bl_project_points(result$train_data, result)
#'
#' # Project test data (class column present — confusion colours)
#' test_pts  <- bl_project_points(result$test_data, result)
#'
#' # Project new observations (no class column — red/blue by prediction)
#' new_obs   <- result$test_data[1:5, result$var_names]
#' new_pts   <- bl_project_points(new_obs, result)
#'
#' @importFrom dplyr case_when
#' @importFrom sp coordinates over
#' @export
bl_project_points <- function(data, bl_result,
                              filter_to_polygon     = FALSE,
                              filter_to_train_ranges = FALSE) {

  if (!inherits(bl_result, "bl_result"))
    stop("'bl_result' must be a 'bl_result' object from bl_assemble().",
         call. = FALSE)
  stop_if_not_data_frame(data, "data")

  var_names   <- bl_result$var_names
  proj_dims   <- bl_result$proj_dims
  V           <- bl_result$V
  X_center    <- bl_result$X_center
  X_sd        <- bl_result$X_sd
  standardise <- bl_result$standardise
  cutoff      <- bl_result$cutoff

  # ---- Optional filter to training variable ranges -----------------------
  if (isTRUE(filter_to_train_ranges)) {
    if (is.null(bl_result$train_ranges)) {
      warning("'filter_to_train_ranges = TRUE' ignored: no train_ranges in bl_result.",
              call. = FALSE)
    } else {
      in_ranges <- get_filter_logical_vector(
        data[, var_names, drop = FALSE], bl_result$train_ranges
      )
      data <- data[in_ranges, , drop = FALSE]
    }
  }

  # ---- Project to Z-space ------------------------------------------------
  sv   <- X_sd
  if (!isTRUE(standardise)) sv <- FALSE
  X_st <- scale(as.matrix(data[, var_names, drop = FALSE]),
                center = X_center, scale = sv)
  Z    <- X_st %*% V[, proj_dims, drop = FALSE]
  colnames(Z) <- c("x", "y")

  # ---- Polygon membership ------------------------------------------------
  inside_polygon <- rep(TRUE, nrow(Z))
  if (!is.null(bl_result$polygon)) {
    Z_df <- as.data.frame(Z)
    sp::coordinates(Z_df) <- ~x + y
    inside_polygon <- !is.na(sp::over(Z_df, bl_result$polygon))
  }

  # ---- Score through model (skipped when no model) -----------------------
  if (!is.null(bl_result$model)) {
    pred_prob  <- .pred_function(
      model_use  = bl_result$model,
      model_type = bl_result$model_type,
      new_data   = data[, var_names, drop = FALSE]
    )
    pred_class <- as.numeric(pred_prob >= cutoff)

    # Confusion colours when true labels available; predicted-class otherwise
    if ("class" %in% names(data)) {
      actual   <- data[["class"]]
      pred_col <- dplyr::case_when(
        actual == pred_class & actual == 1L ~ "red",    # TP
        actual == pred_class & actual == 0L ~ "blue",   # TN
        actual != pred_class & actual == 0L ~ "purple", # FP
        TRUE                               ~ "orange"   # FN
      )
    } else {
      actual   <- NULL
      pred_col <- ifelse(pred_class == 1, "red", "blue")
    }
  } else {
    # No model: colour by true class only (red = 1, blue = 0)
    pred_prob  <- NULL
    pred_class <- NULL
    if ("class" %in% names(data)) {
      actual   <- data[["class"]]
      pred_col <- ifelse(actual == 1L, "red", "blue")
    } else {
      actual   <- NULL
      pred_col <- rep("grey50", nrow(Z))
    }
  }

  # ---- Optional filter to polygon ----------------------------------------
  if (isTRUE(filter_to_polygon) && any(!inside_polygon)) {
    keep       <- inside_polygon
    Z          <- Z[keep, , drop = FALSE]
    pred_prob  <- pred_prob[keep]
    pred_class <- pred_class[keep]
    pred_col   <- pred_col[keep]
    if (!is.null(actual)) actual <- actual[keep]
    inside_polygon <- rep(TRUE, sum(keep))
  }

  structure(
    list(
      Z              = Z,
      pred_prob      = pred_prob,
      pred_class     = pred_class,
      pred_col       = pred_col,
      class          = actual,
      inside_polygon = inside_polygon
    ),
    class = "bl_points"
  )
}


# --------------------------------------------------------------------------
# bl_predict(): tabular prediction summary
# --------------------------------------------------------------------------

#' Tabular prediction summary for a dataset
#'
#' Scores every row of a data frame through the fitted model and returns a
#' tidy data frame suitable for inspection and target selection. This replaces
#' the ad-hoc \code{data.frame()} construction used in example scripts.
#'
#' The returned data frame contains:
#' \itemize{
#'   \item \code{row} — row index (1-based).
#'   \item \code{pred_prob} — model predicted probability (rounded to 4 d.p.).
#'   \item \code{pred_class} — predicted class (0 or 1).
#'   \item \code{true_class} — true class if a \code{"class"} column is present
#'     in \code{data}; otherwise omitted.
#'   \item \code{confusion} — \code{"TP"}, \code{"TN"}, \code{"FP"}, or
#'     \code{"FN"} if true class is known; otherwise omitted.
#'   \item One column per feature variable from \code{bl_result$var_names}.
#' }
#'
#' @param bl_result A \code{"bl_result"} object from \code{\link{bl_assemble}}.
#' @param data Data frame to score. Defaults to \code{bl_result$test_data}.
#'
#' @return A data frame with columns \code{row}, \code{pred_prob},
#'   \code{pred_class}, optionally \code{true_class} and \code{confusion},
#'   followed by all feature columns.
#'
#' @examples
#' bl_dat  <- bl_prepare_data(datasets::iris, class_col = "Species",
#'                             target_class = "versicolor")
#' bl_mod  <- bl_fit_model(bl_dat$train_data, bl_dat$var_names)
#' bl_proj <- bl_build_projection(bl_dat$train_data, bl_dat$var_names)
#' bl_grid <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = 50L)
#' result  <- bl_assemble(bl_dat, bl_model = bl_mod,
#'                         bl_projection = bl_proj, bl_grid = bl_grid)
#'
#' # Inspect test-set predictions
#' bl_predict(result)
#'
#' # Score a custom data frame (no class column)
#' bl_predict(result, data = result$test_data[, result$var_names])
#'
#' @importFrom dplyr case_when
#' @export
bl_predict <- function(bl_result, data = NULL) {
  if (!inherits(bl_result, "bl_result"))
    stop("'bl_result' must be a 'bl_result' object from bl_assemble().",
         call. = FALSE)
  if (is.null(data)) data <- bl_result$test_data
  if (is.null(data))
    stop("No 'data' supplied and 'bl_result$test_data' is NULL.", call. = FALSE)
  stop_if_not_data_frame(data, "data")

  pts       <- bl_project_points(data, bl_result)
  var_names <- bl_result$var_names

  out <- data.frame(
    row        = seq_len(nrow(data)),
    pred_prob  = round(pts$pred_prob, 3L),
    pred_class = as.integer(pts$pred_class),
    stringsAsFactors = FALSE
  )

  if (!is.null(pts$class)) {
    pc <- as.integer(pts$pred_class)
    tc <- as.integer(pts$class)
    out$true_class <- tc
    out$confusion  <- dplyr::case_when(
      pc == 1L & tc == 1L ~ "TP",
      pc == 0L & tc == 0L ~ "TN",
      pc == 1L & tc == 0L ~ "FP",
      pc == 0L & tc == 1L ~ "FN",
      TRUE ~ NA_character_
    )
  }

  cbind(out, data[, var_names, drop = FALSE])
}


# --------------------------------------------------------------------------
# S3 print method
# --------------------------------------------------------------------------

#' @export
print.bl_points <- function(x, ...) {
  cat("<bl_points>\n")
  cat(sprintf("  Points          : %d\n", nrow(x$Z)))
  cat(sprintf("  Labels known    : %s\n", if (!is.null(x$class)) "yes" else "no"))
  cat(sprintf("  Inside polygon  : %d / %d\n",
              sum(x$inside_polygon), length(x$inside_polygon)))
  cat(sprintf("  Pred prob range : [%.4f, %.4f]\n",
              min(x$pred_prob), max(x$pred_prob)))
  invisible(x)
}
