############################################################
# plot.bl_boundary() — distance-to-boundary jitter / boxplot
# bl_robustness()    — scalar total boundary distance
#
# Refactored from: scripts/4.1 Global interpretations - Robust and VIP.r
############################################################

# Bare column-name symbols referenced via NSE in plot.bl_boundary()'s
# ggplot2::aes() calls -- declared so R CMD check does not flag them as
# undefined globals.
utils::globalVariables(c("values", "Variable"))


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
#' of `tV`, then standardised so variables are comparable regardless of
#' their original scale. The denominator depends on `distance`. The default,
#' `"auto"`, follows the thesis rule (section 3.9.1): `"euclidean"` for a PCA
#' biplot and `"mahalanobis"` for a CVA biplot.
#' \describe{
#'   \item{`"mahalanobis"` (CVA default)}{`sqrt(diag(W))` -- within-class
#'     standard deviation per feature. This is the diagonal of the
#'     Mahalanobis metric matrix stored on `bl_result$metric`. Features
#'     that are good class separators (large between-class variance, small
#'     within-class variance) are correctly *amplified* in the importance
#'     measure. A Mahalanobis check is printed to the console: it compares
#'     the squared length of the plotted per-feature distances,
#'     `sum_j v_j^2 / W_jj`, with the full squared Mahalanobis distance
#'     `v' W^{-1} v`, so users can judge how much the plot leaves out by
#'     ignoring the correlations in `W`.}
#'   \item{`"euclidean"` (PCA default)}{`X_sd` -- total per-feature standard
#'     deviation, as the thesis standardises the variable importance
#'     (section 4.2.2). Conflates within- and between-class variance.}
#' }
#' If the projection was standardised (PCA with `standardise = TRUE`),
#' the X-space vector is first unstandardised by multiplying by `X_sd`
#' before the final division.
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
#' @section Why not full Mahalanobis with per-feature decomposition:
#' Full Mahalanobis distance `v' W^{-1} v` is a quadratic form with cross-terms
#' that do not decompose cleanly per feature. The theoretically correct
#' per-feature attribution is the Shapley value decomposition of the full
#' Mahalanobis distance, but its computational cost is `O(n * 2^p)`, which is
#' infeasible for the global Phase 2 use case. See
#' `2 implementation_summary.txt` Section 4.6.1 for the full discussion.
#'
#' @param x        A `"bl_boundary"` object from `bl_find_boundary()`.
#' @param type     Character; `"jitter"` (default) or `"boxplot"`.
#' @param distance Character; per-feature standardisation metric. One of
#'   `"auto"` (default: `"euclidean"` for a PCA biplot, `"mahalanobis"` for a
#'   CVA biplot), `"mahalanobis"` (uses `sqrt(diag(W))`) or `"euclidean"`
#'   (uses `X_sd`).
#' @param ...      Unused; for S3 compatibility.
#'
#' @return Invisibly returns a named list with:
#' \describe{
#'   \item{`plot`}{The `ggplot` object.}
#'   \item{`sum_of_distance`}{Named numeric vector; per-variable total
#'     absolute standardised distance.}
#'   \item{`robustness`}{Numeric scalar; `sum(sum_of_distance)`.}
#'   \item{`vec_to_boundary_sd`}{n x p matrix of signed standardised
#'     distances (one row per obs, one column per variable).}
#'   \item{`distance`}{Character; the denominator used, `"mahalanobis"`
#'     (within-class SD) or `"euclidean"` (total SD), after `"auto"` is
#'     resolved.}
#'   \item{`mahalanobis_check`}{For `"mahalanobis"`, a list of squared
#'     distances summed across observations: `plotted` (the per-feature
#'     measure, `sum(vec_to_boundary_sd^2)`), `full` (`v' W^{-1} v`),
#'     `difference` (`full - plotted`) and `pct` (the difference as a
#'     percentage of `full`). `NULL` for `"euclidean"` or when
#'     `bl_result$metric_inv` is missing.}
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_jitter geom_boxplot geom_vline labs
#' @importFrom ggplot2 scale_color_manual scale_fill_manual theme_light theme
#' @importFrom ggplot2 element_text position_dodge
#' @export
plot.bl_boundary <- function(x, type = c("jitter", "boxplot"),
                              distance = c("auto", "mahalanobis", "euclidean"),
                              ...) {

  type      <- match.arg(type)
  bl_result <- x$bl_result
  # "auto": total SD under PCA, within-class SD under CVA (the thesis rule)
  distance  <- .resolve_distance(distance[1L], bl_result)

  var_names   <- bl_result$var_names
  proj_dims   <- bl_result$proj_dims
  tV          <- bl_result$tV
  X_sd        <- bl_result$X_sd
  standardise <- bl_result$standardise
  metric      <- bl_result$metric
  metric_inv  <- bl_result$metric_inv

  # Fall back to Euclidean if metric is missing (older bl_result objects)
  if (distance == "mahalanobis" && is.null(metric)) {
    warning(
      "bl_result$metric not found; falling back to distance = 'euclidean'. ",
      "Rebuild bl_result with the current package version to enable Mahalanobis.",
      call. = FALSE
    )
    distance <- "euclidean"
  }

  tVr <- tV[proj_dims, , drop = FALSE]   # 2 x p

  # ---- Compute signed standardised distance per variable -----------
  vec_to_boundary_Z  <- x$Z_obs - x$B_z
  vec_to_boundary    <- vec_to_boundary_Z %*% tVr
  colnames(vec_to_boundary) <- var_names

  if (isTRUE(standardise)) {
    vec_to_boundary <- sweep(vec_to_boundary, 2L, X_sd, "*")
  }

  denom <- if (distance == "mahalanobis") sqrt(diag(metric)) else X_sd
  vec_to_boundary_sd <- sweep(vec_to_boundary, 2L, denom, "/")

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
    col_map     <- .confusion_cols
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
  cat(sprintf("\nDistance metric: %s\n", distance))
  cat(sprintf("Robustness (total distance): %.2f\n", robustness))
  cat("\nPer-variable totals (descending):\n")
  print(round(sort(sum_of_distance, decreasing = TRUE), 2L))

  # ---- Mahalanobis check (Mahalanobis only) ------------------------
  # Compares the plotted per-feature measure with the full Mahalanobis
  # distance. The plot divides each feature by sqrt(W_jj), so its squared
  # length is sum_j v_j^2 / W_jj = rowSums(vec_to_boundary_sd^2); the full
  # squared distance is v' W^-1 v. The difference is what the per-feature
  # plot leaves out by ignoring the correlations in W.
  mahalanobis_check <- NULL
  if (distance == "mahalanobis" && !is.null(metric_inv)) {
    full_per_row    <- rowSums((vec_to_boundary %*% metric_inv) * vec_to_boundary)
    plotted_per_row <- rowSums(vec_to_boundary_sd^2)
    total_full      <- sum(full_per_row, na.rm = TRUE)
    total_plotted   <- sum(plotted_per_row, na.rm = TRUE)
    difference      <- total_full - total_plotted
    pct             <- if (total_full != 0) 100 * difference / total_full else 0
    mahalanobis_check <- list(plotted = total_plotted, full = total_full,
                              difference = difference, pct = pct)
    cat(sprintf(paste0(
      "\nMahalanobis check (squared distances, summed across obs):\n",
      "  plotted (per feature) = %.2f | full d_M^2 = %.2f | difference = %.2f (%+.1f%% of full)\n"),
      total_plotted, total_full, difference, pct
    ))
    if (abs(pct) >= 25) {
      cat("  Note: |difference| >= 25% -- the per-feature plot leaves out a large\n")
      cat("        part of the full Mahalanobis distance (the correlations in W).\n")
      cat("        See implementation_summary.txt Section 4.6.1 for alternatives.\n")
    }
  }

  invisible(list(
    plot               = p,
    sum_of_distance    = sum_of_distance,
    robustness         = robustness,
    vec_to_boundary_sd = vec_to_boundary_sd,
    distance           = distance,
    mahalanobis_check  = mahalanobis_check
  ))
}


#' Global robustness score: total distance to boundary
#'
#' Returns per-variable and overall total absolute standardised distances to
#' the decision boundary across all observations. A larger value indicates
#' the population is, on average, further from the boundary.
#'
#' @param bl_boundary A `"bl_boundary"` object from `bl_find_boundary()`.
#' @param distance    Character; per-feature standardisation metric. One of
#'   `"auto"` (default: `"euclidean"` for a PCA biplot, `"mahalanobis"` for a
#'   CVA biplot), `"mahalanobis"` (uses `sqrt(diag(W))` from
#'   `bl_result$metric`) or `"euclidean"` (uses `X_sd`). See
#'   [plot.bl_boundary()] for full details.
#'
#' @return Named list with:
#' \describe{
#'   \item{`sum_of_distance`}{Named numeric vector; per-variable totals.}
#'   \item{`robustness`}{Numeric scalar; `sum(sum_of_distance)`.}
#'   \item{`distance`}{Character; the metric used.}
#' }
#'
#' @export
bl_robustness <- function(bl_boundary,
                           distance = c("auto", "mahalanobis", "euclidean")) {

  if (!inherits(bl_boundary, "bl_boundary"))
    stop("'bl_boundary' must be a 'bl_boundary' object from bl_find_boundary().",
         call. = FALSE)

  bl_result   <- bl_boundary$bl_result
  distance    <- .resolve_distance(distance[1L], bl_result)
  var_names   <- bl_result$var_names
  proj_dims   <- bl_result$proj_dims
  tV          <- bl_result$tV
  X_sd        <- bl_result$X_sd
  standardise <- bl_result$standardise
  metric      <- bl_result$metric

  if (distance == "mahalanobis" && is.null(metric)) {
    warning(
      "bl_result$metric not found; falling back to distance = 'euclidean'. ",
      "Rebuild bl_result with the current package version to enable Mahalanobis.",
      call. = FALSE
    )
    distance <- "euclidean"
  }

  tVr <- tV[proj_dims, , drop = FALSE]

  vec_to_boundary_Z  <- bl_boundary$Z_obs - bl_boundary$B_z
  vec_to_boundary    <- vec_to_boundary_Z %*% tVr
  colnames(vec_to_boundary) <- var_names

  if (isTRUE(standardise)) {
    vec_to_boundary <- sweep(vec_to_boundary, 2L, X_sd, "*")
  }
  denom              <- if (distance == "mahalanobis") sqrt(diag(metric)) else X_sd
  vec_to_boundary_sd <- sweep(vec_to_boundary, 2L, denom, "/")
  sum_of_distance    <- colSums(abs(vec_to_boundary_sd), na.rm = TRUE)

  list(
    sum_of_distance = sum_of_distance,
    robustness      = sum(sum_of_distance, na.rm = TRUE),
    distance        = distance
  )
}
