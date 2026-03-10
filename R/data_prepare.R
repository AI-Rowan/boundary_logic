############################################################
# bl_prepare_data(): data ingestion and binary class conversion
# Refactored from: scripts/2.1 Data Setup.R
############################################################

#' Prepare a data frame for boundary logic analysis
#'
#' Converts a multiclass or binary class column to a numeric 0/1 indicator,
#' selects feature columns, and produces a reproducible train/test split.
#'
#' @section Multiclass to binary conversion:
#' If your data has more than two classes (e.g., `iris` has three Species),
#' supply `target_class` to specify which class maps to `1`. All other classes
#' map to `0`. For example:
#' ```r
#' # iris: treat "versicolor" as the positive class
#' bl_prepare_data(iris, class_col = "Species", target_class = "versicolor")
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
#'
#' @return A list of class `"bl_data"` with components:
#' \describe{
#'   \item{`train_data`}{Data frame of training rows with columns
#'     `var_names` + `"class"` (numeric 0/1).}
#'   \item{`test_data`}{Data frame of test rows, same structure.}
#'   \item{`var_names`}{Character vector of feature column names.}
#'   \item{`num_vars`}{Integer; number of features.}
#'   \item{`target_class`}{The value that was mapped to `1`, or `NULL` if
#'     no conversion was applied.}
#' }
#'
#' @examples
#' # Binary class, versicolor vs rest
#' bl_dat <- bl_prepare_data(iris,
#'                            class_col    = "Species",
#'                            target_class = "versicolor")
#' table(bl_dat$train_data$class)
#'
#' @export
bl_prepare_data <- function(data,
                            class_col,
                            target_class   = NULL,
                            feature_cols   = NULL,
                            train_fraction = 0.8,
                            seed           = 121L) {

  # ---- Input validation ------------------------------------------------
  stop_if_not_data_frame(data, "data")
  stop_if_not_character(class_col, "class_col")
  if (length(class_col) != 1L)
    stop("'class_col' must be a single column name.", call. = FALSE)
  stop_if_col_missing(data, class_col, "class_col")
  stop_if_not_in_range(train_fraction, 0, 1, "train_fraction")

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
  structure(
    list(
      train_data   = train_data,
      test_data    = test_data,
      var_names    = feature_cols,
      num_vars     = length(feature_cols),
      target_class = target_class
    ),
    class = "bl_data"
  )
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
  invisible(x)
}
