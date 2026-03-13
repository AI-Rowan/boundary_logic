############################################################
# plot.bl_boundary() — distance-to-boundary jitter / boxplot
# bl_robustness()    — scalar total boundary distance
#
# Refactored from: scripts/4.1 Global interpretations - Robust and VIP.r
############################################################


#' Distance-to-boundary plot
#'
#' Produces a signed plot showing each observation's standardised distance to
#' the decision boundary, decomposed by variable. Variables are sorted
#' ascending by their total absolute distance (smallest = least important at
#' the bottom).
#'
#' Two chart types are available via the `type` argument:
#' * `"jitter"` (default) — one point per observation per variable, jittered
#'   vertically for legibility.
#' * `"boxplot"` — box-and-whisker summaries per confusion group and variable.
#'
#' @section Colour scheme:
#' Colours match the biplot convention from `bl_project_points()`:
#' * When true labels are present in the source data: four confusion colours —
#'   TP (red), TN (blue), FP (purple), FN (orange).
#' * When no true labels are available: two colours — predicted class 1 (red),
#'   predicted class 0 (blue).
#'
#' @section Distance computation:
#' The distance vector per observation is computed in biplot Z-space
#' (`Z_obs - Z_boundary`), back-projected to X-space via the 2 x p block
#' of `tV`, then standardised by `X_sd` so that variables are comparable
#' regardless of their original scale. If the projection was standardised
#' (PCA with `standardise = TRUE`), the X-space vector is first
#' unstandardised by multiplying by `X_sd` before the final division.
#'
#' The sign indicates which side of the boundary the observation lies on:
#' positive values indicate the observation is above the boundary in that
#' variable's direction. The vertical line at x = 0 marks the boundary.
#'
#' The y-axis label shows `"VarName : <total>"` where `<total>` is the sum
#' of absolute standardised distances across all observations — a
#' variable-level importance proxy. Use `bl_robustness()` for the scalar
#' total across all variables.
#'
#' @param x    A `"bl_boundary"` object from `bl_find_boundary()`.
#' @param type Character; `"jitter"` (default) or `"boxplot"`.
#' @param ...  Unused; for S3 compatibility.
#'
#' @return Invisibly returns a named list with:
#' \describe{
#'   \item{`plot`}{The `ggplot` object.}
#'   \item{`sum_of_distance`}{Named numeric vector; per-variable total
#'     absolute standardised distance.}
#'   \item{`robustness`}{Numeric scalar; `sum(sum_of_distance)`.}
#'   \item{`vec_to_boundary_sd`}{n x p matrix of signed standardised
#'     distances (one row per obs, one column per variable).}
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_jitter geom_boxplot geom_vline labs
#'   scale_color_manual scale_fill_manual theme_light theme element_text
#'   position_dodge
#' @export
plot.bl_boundary <- function(x, type = c("jitter", "boxplot"), ...) {

  type <- match.arg(type)

  bl_result   <- x$bl_result
  var_names   <- bl_result$var_names
  proj_dims   <- bl_result$proj_dims
  tV          <- bl_result$tV
  X_sd        <- bl_result$X_sd
  standardise <- bl_result$standardise

  tVr <- tV[proj_dims, , drop = FALSE]   # 2 x p

  # ---- Compute signed standardised distance per variable -----------
  vec_to_boundary_Z  <- x$Z_obs - x$B_z
  vec_to_boundary    <- vec_to_boundary_Z %*% tVr
  colnames(vec_to_boundary) <- var_names

  if (isTRUE(standardise)) {
    vec_to_boundary <- sweep(vec_to_boundary, 2L, X_sd, "*")
  }
  vec_to_boundary_sd <- sweep(vec_to_boundary, 2L, X_sd, "/")

  # ---- Per-variable total absolute distance ------------------------
  sum_of_distance <- colSums(abs(vec_to_boundary_sd), na.rm = TRUE)

  # ---- Sort variables by total distance (ascending = bottom to top)
  ord        <- order(sum_of_distance)
  var_sorted <- var_names[ord]
  sd_sorted  <- sum_of_distance[ord]
  y_labels   <- paste0(var_sorted, " : ", round(sd_sorted, 0L))

  # ---- Colour group: 4-colour confusion scheme when labels present --
  pred_class_int <- as.integer(x$pred_obs >= bl_result$cutoff)
  has_labels     <- !all(is.na(x$class_obs))

  if (has_labels) {
    actual  <- x$class_obs
    grp     <- dplyr::case_when(
      actual == 1L & pred_class_int == 1L ~ "TP",
      actual == 0L & pred_class_int == 0L ~ "TN",
      actual == 0L & pred_class_int == 1L ~ "FP",
      actual == 1L & pred_class_int == 0L ~ "FN",
      TRUE                               ~ "Unknown"
    )
    col_map     <- c("TP" = "red", "TN" = "blue", "FP" = "purple", "FN" = "orange")
    legend_name <- "Confusion"
  } else {
    grp         <- as.character(pred_class_int)
    col_map     <- c("0" = "blue", "1" = "red")
    legend_name <- "Predicted class"
  }

  # ---- Build long-format data frame for ggplot --------------------
  mat_sorted <- vec_to_boundary_sd[, ord, drop = FALSE]
  colnames(mat_sorted) <- y_labels

  n     <- nrow(mat_sorted)
  n_var <- ncol(mat_sorted)

  dfp <- data.frame(
    values   = round(as.numeric(unlist(mat_sorted, use.names = FALSE)), 2L),
    Variable = factor(rep(y_labels, each = n), levels = y_labels),
    grp      = rep(grp, times = n_var),
    stringsAsFactors = FALSE
  )

  # ---- Build plot --------------------------------------------------
  base_plot <- ggplot2::ggplot(dfp, ggplot2::aes(x = values, y = Variable)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
    ggplot2::labs(
      x     = "Distance to Boundary \u2013 Standardised",
      y     = "Variable : total distance to boundary",
      title = "Variable Distance to Boundary"
    ) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "bottom")

  if (type == "jitter") {
    p <- base_plot +
      ggplot2::geom_jitter(
        ggplot2::aes(colour = grp),
        size = 1.5, height = 0.25
      ) +
      ggplot2::scale_color_manual(name = legend_name, values = col_map)
  } else {
    p <- base_plot +
      ggplot2::geom_boxplot(
        ggplot2::aes(fill = grp, colour = grp),
        position = ggplot2::position_dodge(width = 0.7),
        alpha = 0.4, outlier.size = 0.8
      ) +
      ggplot2::scale_fill_manual(name  = legend_name, values = col_map) +
      ggplot2::scale_color_manual(name = legend_name, values = col_map)
  }

  print(p)

  # ---- Console summary ---------------------------------------------
  robustness <- sum(sum_of_distance, na.rm = TRUE)
  cat(sprintf("\nRobustness (total distance): %.2f\n", robustness))
  cat("\nPer-variable totals (descending):\n")
  print(round(sort(sum_of_distance, decreasing = TRUE), 2L))

  invisible(list(
    plot               = p,
    sum_of_distance    = sum_of_distance,
    robustness         = robustness,
    vec_to_boundary_sd = vec_to_boundary_sd
  ))
}


#' Global robustness score: total distance to boundary
#'
#' Returns per-variable and overall total absolute standardised distances to
#' the decision boundary across all observations. A larger value indicates
#' the population is, on average, further from the boundary.
#'
#' @param bl_boundary A `"bl_boundary"` object from `bl_find_boundary()`.
#'
#' @return Named list with:
#' \describe{
#'   \item{`sum_of_distance`}{Named numeric vector; per-variable totals.}
#'   \item{`robustness`}{Numeric scalar; `sum(sum_of_distance)`.}
#' }
#'
#' @export
bl_robustness <- function(bl_boundary) {

  if (!inherits(bl_boundary, "bl_boundary"))
    stop("'bl_boundary' must be a 'bl_boundary' object from bl_find_boundary().",
         call. = FALSE)

  bl_result   <- bl_boundary$bl_result
  var_names   <- bl_result$var_names
  proj_dims   <- bl_result$proj_dims
  tV          <- bl_result$tV
  X_sd        <- bl_result$X_sd
  standardise <- bl_result$standardise

  tVr <- tV[proj_dims, , drop = FALSE]

  vec_to_boundary_Z  <- bl_boundary$Z_obs - bl_boundary$B_z
  vec_to_boundary    <- vec_to_boundary_Z %*% tVr
  colnames(vec_to_boundary) <- var_names

  if (isTRUE(standardise)) {
    vec_to_boundary <- sweep(vec_to_boundary, 2L, X_sd, "*")
  }
  vec_to_boundary_sd <- sweep(vec_to_boundary, 2L, X_sd, "/")
  sum_of_distance    <- colSums(abs(vec_to_boundary_sd), na.rm = TRUE)

  list(
    sum_of_distance = sum_of_distance,
    robustness      = sum(sum_of_distance, na.rm = TRUE)
  )
}
