############################################################
# Internal feasibility helpers
# Refactored from: scripts/1 Functions for biplot.R (sections 9)
# Determine plausible variable ranges and filter data accordingly.
############################################################

#' Get the min/max range of each variable in a data frame
#'
#' Returns a list of `rlang` expressions of the form
#' `var >= min & var <= max` for each column. These can be passed
#' directly to `get_filter_logical_vector()`.
#'
#' @param data A data frame. All columns are included.
#'
#' @importFrom rlang expr sym eval_tidy
#' @return A list of expressions, one per column of `data`.
#' @keywords internal
get_variable_ranges <- function(data) {
  vars <- names(data)
  filters <- lapply(vars, function(var) {
    rng <- range(data[[var]], na.rm = TRUE)
    rlang::expr(
      (!!rlang::sym(var)) >= !!rng[1] & (!!rlang::sym(var)) <= !!rng[2]
    )
  })
  filters
}


#' Evaluate a list of range expressions as a row-wise AND filter
#'
#' @param data         Data frame to evaluate against.
#' @param filter_exprs List of expressions as returned by
#'   `get_variable_ranges()`.
#'
#' @return Logical vector of length `nrow(data)`. `TRUE` means the row
#'   passes all filters.
#' @keywords internal
get_filter_logical_vector <- function(data, filter_exprs) {
  # Single TRUE means "no filter": return all rows
  if (identical(filter_exprs, list(TRUE)) ||
      identical(filter_exprs, TRUE)) {
    return(rep(TRUE, nrow(data)))
  }

  # Combine all expressions with AND
  row_filter_expr <- Reduce(
    function(x, y) rlang::expr((!!x) & (!!y)),
    filter_exprs
  )
  rlang::eval_tidy(row_filter_expr, data = data)
}
