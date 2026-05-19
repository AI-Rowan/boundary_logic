############################################################
# bl_pick_point(): interactive biplot point picker
############################################################

#' Interactively identify the nearest data point on an active biplot, with
#' optional counterfactual boundary overlay for each picked point
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
#'
#' # Show counterfactual arrow for each picked point
#' bl_bnd <- bl_find_boundary(bl_results)
#' plot_biplotEZ(bl_results)
#' picked <- bl_pick_point(bl_results, bl_boundary = bl_bnd)
#' ```
#'
#' @param bl_result  A `"bl_result"` object from [bl_assemble()].
#' @param data       Data frame to search for the nearest point. Defaults to
#'   `bl_result$train_data`. Pass `bl_result$test_data` to search test
#'   points, or any other data frame whose features match
#'   `bl_result$var_names`.
#' @param bl_boundary  A `"bl_boundary"` object from [bl_find_boundary()], or
#'   `NULL` (default). When supplied, each picked observation's nearest boundary
#'   counterfactual is drawn on the active plot: an `x` cross at `B_z` and (if
#'   `show_arrows = TRUE`) an arrow from the observation to `B_z`. If no
#'   boundary was found for that observation (`B_z = NA`), a message is
#'   printed. The rows of `bl_boundary` must correspond 1-to-1 with the rows
#'   of `data`; a warning is issued if the row counts differ.
#' @param show_arrows  Logical; if `TRUE` (default) and `bl_boundary` is
#'   supplied, draws an arrow from each picked observation to its boundary
#'   counterfactual. Set `FALSE` to draw only the cross at `B_z`.
#' @param arrow_col    Character; colour for the boundary cross and arrow.
#'   Default `"grey30"`.
#'
#' @section Known limitations:
#' If the biplot was rendered with `rotate_deg != 0` in [plot_biplotEZ()],
#' the boundary cross and arrow drawn by `bl_pick_point()` will be misaligned
#' because `bl_boundary$Z_obs` and `B_z` are stored in unrotated Z-space
#' coordinates. Use `rotate_deg = 0` (the default) when combining
#' `bl_pick_point()` with boundary overlay.
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
#' @importFrom graphics locator points text arrows
#' @importFrom grDevices dev.cur
#' @export
bl_pick_point <- function(bl_result,
                           data        = NULL,
                           bl_boundary = NULL,
                           show_arrows = TRUE,
                           arrow_col   = "grey30") {

  if (!inherits(bl_result, "bl_result")) {
    message("bl_pick_point(): 'bl_result' must be a 'bl_result' object from bl_assemble().")
    return(invisible(NULL))
  }

  if (!is.null(bl_boundary) && !inherits(bl_boundary, "bl_boundary"))
    stop("'bl_boundary' must be a 'bl_boundary' object from bl_find_boundary().",
         call. = FALSE)

  if (grDevices::dev.cur() == 1L)
    stop("bl_pick_point(): no active plot found. ",
         "Call plot_biplotEZ() to render the biplot first, then call bl_pick_point().",
         call. = FALSE)

  if (is.null(data)) data <- bl_result$train_data

  if (!is.null(bl_boundary)) {
    n_bnd <- nrow(bl_boundary$B_z)
    n_dat <- nrow(data)
    if (n_dat != n_bnd) {
      warning(sprintf(
        paste0("bl_pick_point(): 'data' has %d rows but 'bl_boundary' has %d rows. ",
               "Row indices will not align. ",
               "Pass the same data used in bl_find_boundary() to avoid this."),
        n_dat, n_bnd
      ), call. = FALSE)
    }
  }

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

    # Boundary counterfactual overlay for the picked point
    if (!is.null(bl_boundary)) {
      if (nearest > nrow(bl_boundary$B_z)) {
        warning(sprintf(
          paste0("bl_pick_point(): nearest index %d exceeds nrow(bl_boundary$B_z) = %d. ",
                 "Data alignment is off -- pass the training data used in bl_find_boundary()."),
          nearest, nrow(bl_boundary$B_z)
        ), call. = FALSE)
      } else {
        bz_row <- bl_boundary$B_z[nearest, ]
        zo_row <- bl_boundary$Z_obs[nearest, ]

        if (all(is.finite(bz_row))) {
          graphics::points(bz_row[1L], bz_row[2L],
                           pch = 4L, col = arrow_col, cex = 0.9, lwd = 1.5)
          if (isTRUE(show_arrows)) {
            graphics::arrows(
              x0     = zo_row[1L],
              y0     = zo_row[2L],
              x1     = bz_row[1L],
              y1     = bz_row[2L],
              length = 0.08,
              col    = arrow_col,
              lwd    = 1.2
            )
          }
        } else {
          message(sprintf("No boundary found for observation %d.", nearest))
        }
      }
    }

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
