############################################################
# Internal utilities: input validation and Gini calculation
############################################################

# --------------------------------------------------------------------------
# Input validators
# --------------------------------------------------------------------------

#' @keywords internal
stop_if_not_data_frame <- function(x, arg_name) {
  if (!is.data.frame(x)) {
    stop(sprintf("'%s' must be a data.frame, not %s.", arg_name, class(x)[1]),
         call. = FALSE)
  }
}

#' @keywords internal
stop_if_not_character <- function(x, arg_name) {
  if (!is.character(x) || length(x) == 0) {
    stop(sprintf("'%s' must be a non-empty character vector.", arg_name),
         call. = FALSE)
  }
}

#' @keywords internal
stop_if_col_missing <- function(data, col, arg_name) {
  if (!col %in% names(data)) {
    stop(sprintf("Column '%s' (supplied via '%s') not found in data. Available columns: %s.",
                 col, arg_name, paste(names(data), collapse = ", ")),
         call. = FALSE)
  }
}

#' @keywords internal
stop_if_not_in_range <- function(x, lo, hi, arg_name) {
  if (!is.numeric(x) || length(x) != 1 || x <= lo || x >= hi) {
    stop(sprintf("'%s' must be a single numeric value strictly between %g and %g.",
                 arg_name, lo, hi),
         call. = FALSE)
  }
}

#' @keywords internal
stop_if_not_scalar_numeric <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x)) {
    stop(sprintf("'%s' must be a single finite numeric value.", arg_name),
         call. = FALSE)
  }
}

#' @keywords internal
stop_if_not_positive_integer <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1 || x != round(x) || x < 1L) {
    stop(sprintf("'%s' must be a single positive integer.", arg_name),
         call. = FALSE)
  }
}

# --------------------------------------------------------------------------
# Gini coefficient
# Defined as 2 * AUC - 1, computed from predicted probabilities.
# --------------------------------------------------------------------------

#' Calculate the Gini coefficient (2 * AUC - 1)
#'
#' @param actual Numeric 0/1 vector of true class labels.
#' @param predicted_prob Numeric vector of predicted class-1 probabilities.
#'
#' @return Single numeric in [-1, 1].
#' @keywords internal
calc_gini <- function(actual, predicted_prob) {
  n <- length(actual)
  if (n != length(predicted_prob)) {
    stop("'actual' and 'predicted_prob' must have the same length.", call. = FALSE)
  }
  if (length(unique(actual)) < 2L) return(NA_real_)

  # Sort by descending predicted probability
  ord      <- order(predicted_prob, decreasing = TRUE)
  actual_s <- actual[ord]

  # Trapezoidal AUC
  n1   <- sum(actual == 1)
  n0   <- n - n1
  if (n1 == 0L || n0 == 0L) return(NA_real_)

  tp   <- cumsum(actual_s == 1)
  fp   <- cumsum(actual_s == 0)
  tpr  <- tp / n1
  fpr  <- fp / n0

  # Add origin
  tpr <- c(0, tpr)
  fpr <- c(0, fpr)

  auc  <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)])) / 2
  2 * auc - 1
}
