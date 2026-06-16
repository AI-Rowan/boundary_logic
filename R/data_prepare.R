############################################################
# bl_prepare_data(): data ingestion and binary class conversion
# Refactored from: scripts/2.1 Data Setup.R
############################################################

#' Prepare a data frame for boundary logic analysis
#'
#' Converts a multiclass or binary class column to a numeric 0/1 indicator,
#' selects feature columns, produces a reproducible train/test split, and
#' removes outliers from the training set using a convex hull polygon filter.
#' Outlier filtering is integrated: pass `hull_fraction` to control how
#' aggressively the hull trims the outermost training points (default `0.9`
#' removes approximately the outermost 10%).
#'
#' @section Multiclass to binary conversion:
#' If your data has more than two classes (e.g., `datasets::iris` has three Species),
#' supply `target_class` to specify which class maps to `1`. All other classes
#' map to `0`. For example:
#' ```r
#' # datasets::iris: treat "versicolor" as the positive class
#' bl_prepare_data(datasets::iris, class_col = "Species", target_class = "versicolor")
#' ```
#' If `target_class = NULL`, the class column must already be numeric 0/1.
#'
#' @param data           A data frame.
#' @param class_col      Character scalar; name of the column holding class
#'   labels.
#' @param target_class   The value in `class_col` that should be coded as `1`.
#'   All other values become `0`. If `NULL`, the column must already be
#'   numeric 0/1.
#' @param feature_cols   Character vector of column names to use as features.
#'   If `NULL` (default), all columns except `class_col` are used.
#' @param train_fraction Numeric in (0, 1); proportion of rows allocated to
#'   training. Default `0.8`.
#' @param seed           Integer; random seed for the train/test split.
#'   Default `121`.
#' @param hull_fraction  Numeric in (0, 1]; fraction argument for the convex
#'   hull polygon filter applied to training data. Values below `1` trim the
#'   most extreme points; `1` retains all points. Default `0.9`.
#' @param verbose        Logical; if `TRUE`, prints a one-line summary of rows
#'   retained and removed by the hull filter. Default `TRUE`.
#'
#' @return A list of class `"bl_filter_result"` with components:
#' \describe{
#'   \item{`train_data`}{Filtered training data frame with columns
#'     `var_names` + `"class"` (numeric 0/1).}
#'   \item{`test_data`}{Data frame of test rows, same structure (not filtered).}
#'   \item{`var_names`}{Character vector of feature column names.}
#'   \item{`num_vars`}{Integer; number of features.}
#'   \item{`target_class`}{The value that was mapped to `1`, or `NULL` if
#'     no conversion was applied.}
#'   \item{`polygon`}{`SpatialPolygons` convex hull in standardised PCA space.}
#'   \item{`hull_fraction`}{The fraction value used for the hull filter.}
#'   \item{`n_retained`}{Number of training rows retained after filtering.}
#'   \item{`n_removed`}{Number of training rows removed by the filter.}
#' }
#'
#' @examples
#' # Binary class, versicolor vs rest, default hull filter (0.9)
#' bl_dat <- bl_prepare_data(datasets::iris,
#'                            class_col    = "Species",
#'                            target_class = "versicolor",
#'                            hull_fraction = 0.9)
#' cat(bl_dat$n_retained, "training rows retained\n")
#'
#' @export
bl_prepare_data <- function(data,
                            class_col,
                            target_class   = NULL,
                            feature_cols   = NULL,
                            train_fraction = 0.8,
                            seed           = 121L,
                            hull_fraction  = 0.9,
                            verbose        = TRUE) {

  # ---- Input validation ------------------------------------------------
  stop_if_not_data_frame(data, "data")
  stop_if_not_character(class_col, "class_col")
  if (length(class_col) != 1L)
    stop("'class_col' must be a single column name.", call. = FALSE)
  stop_if_col_missing(data, class_col, "class_col")
  stop_if_not_in_range(train_fraction, 0, 1.0001, "train_fraction")

  # ---- Select feature columns ------------------------------------------
  if (is.null(feature_cols)) {
    feature_cols <- setdiff(names(data), class_col)
  } else {
    stop_if_not_character(feature_cols, "feature_cols")
    missing_cols <- setdiff(feature_cols, names(data))
    if (length(missing_cols) > 0L)
      stop(sprintf("feature_cols not found in data: %s.",
                   paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  if (length(feature_cols) == 0L)
    stop("No feature columns found. Check 'feature_cols' and 'class_col'.",
         call. = FALSE)

  # ---- Binary class conversion -----------------------------------------
  class_raw <- data[[class_col]]

  if (!is.null(target_class)) {
    if (!any(class_raw == target_class))
      stop(sprintf(
        "target_class '%s' not found in column '%s'. Available values: %s.",
        target_class, class_col,
        paste(unique(as.character(class_raw)), collapse = ", ")
      ), call. = FALSE)
    class_bin <- as.numeric(class_raw == target_class)
  } else {
    # Must already be 0/1
    class_num <- suppressWarnings(as.numeric(class_raw))
    if (any(is.na(class_num)) || !all(class_num %in% c(0, 1)))
      stop(paste0(
        "When 'target_class' is NULL, '", class_col,
        "' must contain only 0 and 1. ",
        "Supply 'target_class' to perform a multiclass-to-binary conversion."
      ), call. = FALSE)
    class_bin    <- class_num
    target_class <- NULL
  }

  # ---- Assemble clean data frame ---------------------------------------
  clean_data <- data[, feature_cols, drop = FALSE]
  clean_data[["class"]] <- class_bin

  # ---- Train / test split ----------------------------------------------
  set.seed(seed)
  n         <- nrow(clean_data)
  train_idx <- sample(n, size = floor(train_fraction * n), replace = FALSE)

  train_data <- clean_data[ train_idx, , drop = FALSE]
  test_data  <- clean_data[-train_idx, , drop = FALSE]

  rownames(train_data) <- NULL
  rownames(test_data)  <- NULL

  # ---- Return ----------------------------------------------------------
  bl_data_raw <- structure(
    list(
      train_data   = train_data,
      test_data    = test_data,
      var_names    = feature_cols,
      num_vars     = length(feature_cols),
      target_class = target_class
    ),
    class = "bl_data"
  )
  bl_filter_outliers(bl_data_raw, hull_fraction = hull_fraction, verbose = verbose)
}


# --------------------------------------------------------------------------
# bl_wrap_data(): bypass bl_prepare_data() for pre-split, pre-coded data
# --------------------------------------------------------------------------

#' Wrap already-prepared train/test data into a bl_data object
#'
#' Use this function when your data is already split into training and test
#' sets and the class column is already coded as numeric 0/1. This bypasses
#' the splitting and multiclass-conversion steps of `bl_prepare_data()`.
#'
#' @param train_data  Data frame; feature columns plus a column named
#'   `"class"` (numeric 0/1).
#' @param test_data   Data frame; same column structure as `train_data`.
#'   If `NULL`, an empty data frame with matching columns is used.
#' @param var_names   Character vector of feature column names. If `NULL`
#'   (default), all columns except `"class"` are used.
#' @param target_class Optional informational label indicating what the `1`
#'   class represents. Does not affect any computation.
#'
#' @return A list of class `"bl_data"` compatible with all downstream
#'   functions (`bl_filter_outliers()`, `bl_build_projection()`, etc.).
#'
#' @examples
#' \dontrun{
#' df       <- datasets::iris
#' df$class <- as.numeric(df$Species == "versicolor")
#' df$Species <- NULL
#' train_df <- df[1:100, ]
#' test_df  <- df[101:150, ]
#' bl_dat   <- bl_wrap_data(train_df, test_df, target_class = "versicolor")
#' }
#'
#' @export
bl_wrap_data <- function(train_data,
                         test_data    = NULL,
                         var_names    = NULL,
                         target_class = NULL) {

  stop_if_not_data_frame(train_data, "train_data")

  if (!"class" %in% names(train_data))
    stop("'train_data' must contain a column named 'class'.", call. = FALSE)

  class_num <- suppressWarnings(as.numeric(train_data[["class"]]))
  if (any(is.na(class_num)) || !all(class_num %in% c(0, 1)))
    stop("The 'class' column in 'train_data' must contain only 0 and 1.",
         call. = FALSE)

  if (is.null(var_names)) {
    var_names <- setdiff(names(train_data), "class")
  } else {
    stop_if_not_character(var_names, "var_names")
    missing_cols <- setdiff(var_names, names(train_data))
    if (length(missing_cols) > 0L)
      stop(sprintf("var_names not found in train_data: %s.",
                   paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  if (length(var_names) == 0L)
    stop("No feature columns found. Check 'var_names'.", call. = FALSE)

  if (is.null(test_data)) {
    test_data <- train_data[0L, , drop = FALSE]
  } else {
    stop_if_not_data_frame(test_data, "test_data")
    if (!"class" %in% names(test_data))
      stop("'test_data' must contain a column named 'class'.", call. = FALSE)
    missing_test <- setdiff(var_names, names(test_data))
    if (length(missing_test) > 0L)
      stop(sprintf("var_names not found in test_data: %s.",
                   paste(missing_test, collapse = ", ")), call. = FALSE)
  }

  train_data[["class"]] <- as.numeric(train_data[["class"]])
  test_data[["class"]]  <- as.numeric(test_data[["class"]])
  rownames(train_data)  <- NULL
  rownames(test_data)   <- NULL

  structure(
    list(
      train_data   = train_data,
      test_data    = test_data,
      var_names    = var_names,
      num_vars     = length(var_names),
      target_class = target_class
    ),
    class = "bl_data"
  )
}


# --------------------------------------------------------------------------
# bl_set_scaling(): record a pre-standardisation transform for raw-unit axes
# --------------------------------------------------------------------------

#' Record the feature scaling used before model fitting
#'
#' When features are standardised (or otherwise affinely rescaled) *before*
#' being passed to `boundarylogic` -- a common step to improve model fit for
#' scale-sensitive models -- every biplot axis is labelled in the standardised
#' units, which is hard to read. `bl_set_scaling()` records the per-feature
#' affine transform `standardised = (raw - center) / scale` so that the biplot
#' plot methods can relabel their axis ticks back into the original (raw) units.
#'
#' The scaling is **display-only metadata**. It does not change the projection,
#' the prediction grid, or what the model sees: the model (including a custom
#' model wrapped via [bl_wrap_model()]) continues to receive data in the same
#' (standardised) units it was trained on. Only the plot methods -- and, where
#' applicable, the printed counterfactual/target values -- consume the scaling.
#'
#' Attach it as the final data-preparation step, after [bl_prepare_data()] /
#' [bl_wrap_data()] (and after any [bl_filter_outliers()]). The recorded scaling
#' then flows through [bl_build_result()] / [bl_assemble()] into the `bl_result`.
#'
#' @param x      A `"bl_data"` object from [bl_wrap_data()] or a
#'   `"bl_filter_result"` object from [bl_prepare_data()] / [bl_filter_outliers()].
#' @param center Numeric vector of per-feature centring constants. Either named
#'   (names must cover `x$var_names`; reordered to match) or unnamed and of
#'   length `length(x$var_names)` in `var_names` order. The attributes produced
#'   by base [scale()] can be passed directly, e.g.
#'   `center = attr(z, "scaled:center")`.
#' @param scale  Numeric vector of per-feature scale factors, same length and
#'   naming rules as `center`. Must contain no zero or `NA` entries.
#' @param method Character scalar naming the transform family, recorded for
#'   reference only (default `"z-score"`). All features share this single
#'   affine method; per-variable transform families are future work.
#'
#' @return `x` with a `scaling` component added:
#'   `list(center = <named numeric>, scale = <named numeric>, method = <character>)`,
#'   both vectors ordered to match `x$var_names`. The object's class is unchanged.
#'
#' @examples
#' df       <- datasets::iris
#' df$class <- as.numeric(df$Species == "versicolor")
#' df$Species <- NULL
#' feats    <- setdiff(names(df), "class")
#' z        <- scale(df[, feats])
#' df[, feats] <- z
#' bl_dat   <- bl_wrap_data(df[1:100, ], df[101:150, ])
#' bl_dat   <- bl_set_scaling(bl_dat,
#'                            center = attr(z, "scaled:center"),
#'                            scale  = attr(z, "scaled:scale"))
#'
#' @seealso [bl_wrap_data()], [bl_prepare_data()], [bl_build_result()]
#' @export
bl_set_scaling <- function(x, center, scale, method = "z-score") {

  if (!inherits(x, c("bl_data", "bl_filter_result")))
    stop("'x' must be a 'bl_data' object from bl_wrap_data()/bl_prepare_data() ",
         "or a 'bl_filter_result' object from bl_filter_outliers().",
         call. = FALSE)

  var_names <- x$var_names
  if (!is.character(method) || length(method) != 1L)
    stop("'method' must be a single character string.", call. = FALSE)

  center <- .align_scaling_vec(center, var_names, "center")
  scale  <- .align_scaling_vec(scale,  var_names, "scale")

  if (any(scale == 0))
    stop("'scale' must not contain zero (would make the transform non-invertible).",
         call. = FALSE)

  x$scaling <- list(center = center, scale = scale, method = method)
  x
}


#' Validate and order a scaling vector against the feature names
#'
#' Coerces a `center`/`scale` vector to numeric, checks for finiteness, and
#' reorders it to `var_names`. Named vectors are matched by name (must cover
#' every feature); unnamed vectors must already be in `var_names` order.
#'
#' @param v         Numeric vector supplied by the user.
#' @param var_names Character vector of feature names.
#' @param arg_name  Character; the argument name for error messages.
#' @return Numeric vector of length `length(var_names)`, named by `var_names`.
#' @noRd
.align_scaling_vec <- function(v, var_names, arg_name) {
  if (!is.numeric(v))
    stop(sprintf("'%s' must be numeric.", arg_name), call. = FALSE)
  if (any(!is.finite(v)))
    stop(sprintf("'%s' must not contain NA, NaN, or Inf.", arg_name),
         call. = FALSE)

  p <- length(var_names)
  if (!is.null(names(v))) {
    missing_nm <- setdiff(var_names, names(v))
    if (length(missing_nm) > 0L)
      stop(sprintf("'%s' is missing entries for: %s.",
                   arg_name, paste(missing_nm, collapse = ", ")),
           call. = FALSE)
    out <- v[var_names]
  } else {
    if (length(v) != p)
      stop(sprintf(
        "Unnamed '%s' must have length %d (one per feature, in var_names order); got %d.",
        arg_name, p, length(v)), call. = FALSE)
    out <- v
    names(out) <- var_names
  }
  out
}


#' Convert standardised feature values back to original (raw) units
#'
#' Display-only inverse of a recorded `bl_set_scaling()` transform, used by the
#' Shapley / sparse / target print and plot methods to report feature values in
#' the user's original units. For a per-feature transform
#' `standardised = (raw - center) / scale`:
#'   - `kind = "level"` (an observed value, counterfactual, sparse value):
#'     `raw = value * scale + center`.
#'   - `kind = "delta"` (a difference such as `data_to_boundary = end - start`):
#'     `raw = value * scale` (no centre term -- differences are translation-free).
#'
#' @param values    Numeric vector aligned to `var_names` (one per feature).
#' @param scaling   The `bl_result$scaling` list (`center`, `scale`, `method`)
#'   or `NULL`.
#' @param var_names Feature names matching the order/identity of `values`.
#' @param kind      `"level"` or `"delta"`.
#' @return `values` converted to raw units, or unchanged when `scaling` is
#'   `NULL` or does not cover every feature (with a warning in the latter case).
#' @noRd
.scale_to_raw <- function(values, scaling, var_names,
                          kind = c("level", "delta")) {
  if (is.null(scaling)) return(values)
  kind   <- match.arg(kind)
  center <- scaling$center[var_names]
  scale  <- scaling$scale[var_names]
  if (any(is.na(center)) || any(is.na(scale))) {
    warning("bl_result$scaling does not cover all features; ",
            "values left in model (standardised) units.", call. = FALSE)
    return(values)
  }
  if (kind == "delta") unname(values) * unname(scale)
  else                 unname(values) * unname(scale) + unname(center)
}


# --------------------------------------------------------------------------
# S3 print method
# --------------------------------------------------------------------------

#' @export
print.bl_data <- function(x, ...) {
  target_str <- if (is.null(x$target_class)) {
    "pre-coded 0/1"
  } else {
    sprintf("'%s' -> 1, all others -> 0", x$target_class)
  }
  cat("<bl_data>\n")
  cat(sprintf("  Features   : %d (%s)\n",
              x$num_vars, paste(x$var_names, collapse = ", ")))
  cat(sprintf("  Target     : %s\n", target_str))
  cat(sprintf("  Train rows : %d  |  Test rows : %d\n",
              nrow(x$train_data), nrow(x$test_data)))
  if (!is.null(x$scaling))
    cat(sprintf("  Scaling    : %s (raw-unit biplot axes enabled)\n",
                x$scaling$method))
  invisible(x)
}
