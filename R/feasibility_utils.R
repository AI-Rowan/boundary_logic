############################################################
# Internal feasibility helpers
# Refactored from: scripts/1 Functions for biplot.R (sections 9)
# Determine plausible variable ranges and filter data accordingly.
############################################################

#' Get the min/max range of each variable in a data frame
#'
#' Returns a named list with one element per column. Each element is a
#' length-2 numeric vector `c(min, max)` for that variable. The list can
#' be passed directly to `get_filter_logical_vector()`.
#'
#' @param data A data frame. All columns are included.
#'
#' @return Named list of `c(min, max)` vectors, one per column of `data`.
#' @keywords internal
get_variable_ranges <- function(data) {
  vars <- names(data)
  ranges <- lapply(vars, function(var) range(data[[var]], na.rm = TRUE))
  names(ranges) <- vars
  ranges
}


#' Evaluate a list of variable ranges as a row-wise AND filter
#'
#' @param data        Data frame to evaluate against.
#' @param range_list  Named list of `c(min, max)` vectors as returned by
#'   `get_variable_ranges()`. Variables not present in `data` are silently
#'   skipped.
#'
#' @return Logical vector of length `nrow(data)`. `TRUE` means the row
#'   satisfies all range constraints.
#' @keywords internal
get_filter_logical_vector <- function(data, range_list) {
  if (identical(range_list, list(TRUE)) || identical(range_list, TRUE))
    return(rep(TRUE, nrow(data)))

  keep <- rep(TRUE, nrow(data))
  for (var in names(range_list)) {
    if (!var %in% names(data)) next
    rng  <- range_list[[var]]
    vals <- data[[var]]
    keep <- keep & !is.na(vals) & vals >= rng[1L] & vals <= rng[2L]
  }
  keep
}
