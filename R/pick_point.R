############################################################
# bl_pick_point(): interactive biplot point picker
############################################################

#' Interactively identify the nearest data point on an active biplot
#'
#' After rendering a biplot with [plot_biplotEZ()], call this function and
#' click anywhere on the plot. For each click the nearest data point is
#' identified by Z-space Euclidean distance, its row number printed to the
#' console, and the point highlighted with a yellow circle and row label
#' on the active plot device. Press **Escape** to stop clicking.
#'
#' @section Typical usage:
#' ```r
#' plot_biplotEZ(bl_results)
#' picked <- bl_pick_point(bl_results)
#'
#' # Pick from test data
#' plot_biplotEZ(bl_results,
#'               points = bl_project_points(bl_results$test_data, bl_results))
#' picked <- bl_pick_point(bl_results, data = bl_results$test_data)
#' ```
#'
#' @param bl_result  A `"bl_result"` object from [bl_assemble()].
#' @param data       Data frame to search for the nearest point. Defaults to
#'   `bl_result$train_data`. Pass `bl_result$test_data` to search test
#'   points, or any other data frame whose features match
#'   `bl_result$var_names`.
#'
#' @return Invisibly returns a data frame with one row per click containing:
#' \describe{
#'   \item{`row`}{Row index of the nearest point in `data`.}
#'   \item{`dist_to_click`}{Z-space Euclidean distance from the click to the
#'     nearest point.}
#'   \item{`z1`}{Z-space coordinate 1 of the nearest point.}
#'   \item{`z2`}{Z-space coordinate 2 of the nearest point.}
#'   \item{(feature columns)}{All columns named in `bl_result$var_names`.}
#' }
#' Returns `NULL` invisibly if no points were picked.
#'
#' @importFrom graphics locator points text
#' @export
bl_pick_point <- function(bl_result, data = NULL) {

  if (!inherits(bl_result, "bl_result")) {
    message("bl_pick_point(): 'bl_result' must be a 'bl_result' object from bl_assemble().")
    return(invisible(NULL))
  }

  if (is.null(data)) data <- bl_result$train_data

  # Project data to Z-space
  pts <- bl_project_points(data, bl_result)
  Z   <- pts$Z

  cat("Click on the biplot to identify the nearest point.\n")
  cat("Press Escape to finish.\n\n")

  picked <- list()

  repeat {
    click <- graphics::locator(1)
    if (is.null(click)) break

    # Nearest point by Z-space Euclidean distance
    dists   <- sqrt((Z[, 1L] - click$x)^2 + (Z[, 2L] - click$y)^2)
    nearest <- which.min(dists)

    cat(sprintf("Row: %d  |  Z: (%.3f, %.3f)  |  dist to click: %.4f\n",
                nearest, Z[nearest, 1L], Z[nearest, 2L], dists[nearest]))
    print(data[nearest, , drop = FALSE])
    cat("\n")

    # Highlight on the active plot
    graphics::points(Z[nearest, 1L], Z[nearest, 2L],
                     pch = 21, cex = 2, bg = "yellow", col = "black", lwd = 1.5)
    graphics::text(Z[nearest, 1L], Z[nearest, 2L],
                   labels = nearest, cex = 0.7, font = 2, pos = 3)

    picked[[length(picked) + 1L]] <- data.frame(
      row           = nearest,
      dist_to_click = dists[nearest],
      z1            = Z[nearest, 1L],
      z2            = Z[nearest, 2L],
      data[nearest, bl_result$var_names, drop = FALSE],
      row.names     = NULL
    )
  }

  if (length(picked) == 0L) {
    cat("No points picked.\n")
    return(invisible(NULL))
  }

  result <- do.call(rbind, picked)
  cat(sprintf("Total points picked: %d\n", nrow(result)))
  invisible(result)
}
