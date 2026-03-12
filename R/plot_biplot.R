############################################################
# plot_biplotEZ(): Phase 1 biplot visualisation
# Based on: scripts/1.6 Plot_biplotEZ.R
# Uses the biplotEZ pipeline for the base plot and axis rendering.
############################################################

#' Plot the Phase 1 biplot with decision boundary
#'
#' Creates a biplot visualisation of the Phase 1 result using the biplotEZ
#' rendering pipeline. The plot shows:
#' \itemize{
#'   \item **Base axes** — biplotEZ variable axes drawn as a canvas.
#'   \item **Coloured prediction grid** — each grid point coloured blue
#'     (class 0) through white to red (class 1) by model probability.
#'   \item **Training data points** — coloured by confusion category.
#'   \item **Variable axes redrawn on top** — so axes are visible over
#'     the grid and points.
#'   \item **Decision boundary contour lines.**
#'   \item Optionally, a **target observation** highlighted in yellow.
#' }
#'
#' @section Point colour legend:
#' \describe{
#'   \item{Red}{True positive (actual = 1, predicted = 1).}
#'   \item{Blue}{True negative (actual = 0, predicted = 0).}
#'   \item{Purple}{False positive (actual = 0, predicted = 1).}
#'   \item{Orange}{False negative (actual = 1, predicted = 0).}
#'   \item{Yellow}{Target observation (if supplied).}
#' }
#'
#' @param bl_result        A `"bl_result"` object from `bl_assemble()`.
#' @param target_point     Named numeric vector of feature values (in original
#'   X-space) for an observation to highlight. Names must match
#'   `bl_result$var_names`. If `NULL` (default), no target is plotted.
#' @param target_label     Character or integer label displayed inside the
#'   yellow target marker. If `NULL` and `target_point` is provided, the
#'   point is shown without a label.
#' @param no_grid          Logical; if `TRUE`, the prediction grid is hidden.
#'   Default `FALSE`.
#' @param points           A `"bl_points"` object from `bl_project_points()`.
#'   Controls which observations are plotted as coloured points. If `NULL`
#'   (default), training data is projected automatically. Pass
#'   `bl_project_points(result$test_data, result)` to show test data instead.
#' @param no_points        Logical; if `TRUE`, data points are hidden.
#'   Default `FALSE`.
#' @param no_contour       Logical; if `TRUE`, decision boundary contour lines
#'   are hidden. Default `FALSE`.
#' @param new_title        Character; overrides the biplot title stored in
#'   `bl_result`. Use `NA` (default) to keep the stored title.
#' @param cex_z            Numeric; size of training data points (`cex`).
#'   Default `0.5`.
#' @param label_dir        Character; biplotEZ axis label direction.
#'   `"Hor"` (horizontal, default) or `"Rad"` (radial).
#' @param tick_label_cex   Numeric; axis tick label size. Default `0.6`.
#' @param ticks_v          Integer; number of ticks per variable axis.
#'   Default `1L`.
#' @param which            Integer vector; indices of variables to draw axes
#'   for. Defaults to all variables.
#' @param X_names          Character vector; custom variable names for axis
#'   labels. Defaults to `bl_result$var_names`.
#' @param label_offset_var Integer or integer vector; index/indices of variables
#'   whose axis labels should be shifted outward. `0` (default) means no
#'   offset. Supply a vector (e.g., `c(1L, 3L)`) to offset multiple variables.
#' @param label_offset_dist Numeric or numeric vector; outward offset distance(s)
#'   for `label_offset_var`. If a single value, the same distance is applied to
#'   all variables listed in `label_offset_var`. If a vector, must be the same
#'   length as `label_offset_var`. Default `0.5`.
#' @param rotate_deg       Numeric; angle in degrees to rotate the entire plot
#'   clockwise. Rotates all plotted elements — grid, training points, contour
#'   lines, variable axes, and the target point — without rerunning the
#'   pipeline. Default `0` (no rotation).
#' @param contour_col      Character; colour for decision boundary contour
#'   lines. Default `"black"`.
#' @param contour_lwd      Numeric; line width for contour lines. Default
#'   `1.5`.
#' @param contour_lty      Integer; line type for contour lines (1 = solid).
#'   Default `1L`.
#'
#' @return Invisibly returns `bl_result`. Called for its side-effect (plot).
#'
#' @importFrom graphics lines points text
#' @importFrom grDevices col2rgb colorRampPalette rgb
#'
#' @examples
#' bl_dat  <- bl_prepare_data(datasets::iris, class_col = "Species",
#'                             target_class = "versicolor")
#' bl_mod  <- bl_fit_model(bl_dat$train_data, bl_dat$var_names)
#' bl_proj <- bl_build_projection(bl_dat$train_data, bl_dat$var_names)
#' bl_grid <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = 100L)
#' result  <- bl_assemble(bl_dat, bl_model = bl_mod,
#'                         bl_projection = bl_proj, bl_grid = bl_grid)
#'
#' plot_biplotEZ(result)
#'
#' # Highlight a target observation
#' target <- unlist(result$test_data[1, result$var_names])
#' plot_biplotEZ(result, target_point = target, target_label = 1)
#'
#' @export
plot_biplotEZ <- function(bl_result,
                           points            = NULL,
                           target_point      = NULL,
                           target_label      = NULL,
                           no_grid           = FALSE,
                           no_points         = FALSE,
                           no_contour        = FALSE,
                           new_title         = NA,
                           cex_z             = 0.5,
                           label_dir         = "Hor",
                           tick_label_cex    = 0.6,
                           ticks_v           = 1L,
                           which             = NULL,
                           X_names           = NULL,
                           label_offset_var  = 0L,
                           label_offset_dist = 0.5,
                           rotate_deg        = 0,
                           contour_col       = "black",
                           contour_lwd       = 1.5,
                           contour_lty       = 1L) {

  if (!inherits(bl_result, "bl_result"))
    stop("'bl_result' must be a 'bl_result' object from bl_assemble().",
         call. = FALSE)

  # ---- Resolve points to plot -------------------------------------------
  if (is.null(points)) {
    points <- bl_project_points(bl_result$train_data, bl_result)
  } else if (!inherits(points, "bl_points")) {
    stop("'points' must be a 'bl_points' object from bl_project_points().",
         call. = FALSE)
  }

  gr          <- bl_result$biplot_grid
  biplot_plot <- bl_result$biplot_obj
  proj_dims   <- bl_result$proj_dims
  num_vars    <- bl_result$num_vars
  var_names   <- bl_result$var_names

  # ---- Defaults for axis arguments -------------------------------------
  if (is.null(which))   which   <- seq_len(num_vars)
  if (is.null(X_names)) X_names <- var_names

  # ---- Optional title override -----------------------------------------
  if (!is.na(new_title)) biplot_plot$Title <- new_title

  # ---- Label offset vector (shift one or more axis labels outward) ----
  label_line_vec <- rep(0, num_vars)
  valid_idx <- label_offset_var[label_offset_var >= 1L &
                                label_offset_var <= num_vars]
  if (length(valid_idx) > 0L) {
    dist_vec <- rep_len(label_offset_dist, length(valid_idx))
    label_line_vec[valid_idx] <- dist_vec
  }

  # ---- Optional rotation (clockwise by rotate_deg degrees) -------------
  if (!is.null(rotate_deg) && rotate_deg != 0) {
    theta <- -rotate_deg * pi / 180   # negative = clockwise
    R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2L)

    gr$Zgrid   <- gr$Zgrid %*% R
    points$Z   <- points$Z %*% R

    gr$ct <- lapply(gr$ct, function(cl) {
      pts      <- cbind(cl$x, cl$y) %*% R
      cl$x     <- pts[, 1L]
      cl$y     <- pts[, 2L]
      cl
    })

    # Rotate the axis arrow directions in the biplotEZ object
    biplot_plot$Lmat[, proj_dims] <- biplot_plot$Lmat[, proj_dims] %*% R
  }

  # ---- Step 1: biplotEZ base plot (axes only, no samples) --------------
  # Establishes the coordinate system and draws the axis skeleton.
  biplot_plot |>
    biplotEZ::samples(opacity = 0, which = NULL) |>
    biplotEZ::axes(col            = "grey",
                   label.dir      = label_dir,
                   which          = which,
                   X.names        = X_names,
                   tick.label.cex = tick_label_cex,
                   ticks          = ticks_v,
                   label.line     = label_line_vec) |>
    plot()

  # ---- Step 2: prediction grid (coloured squares) ----------------------
  if (!isTRUE(no_grid)) {
    graphics::points(gr$Zgrid,
                     type = "p",
                     col  = gr$col_value,
                     pch  = 15L,
                     cex  = 0.5)
  }

  # ---- Step 3: data points ---------------------------------------------
  if (!isTRUE(no_points)) {
    graphics::points(x   = points$Z[, 1L],
                     y   = points$Z[, 2L],
                     col = points$pred_col,
                     pch = 16L,
                     cex = cex_z)
  }

  # ---- Step 4: axes redrawn on top (visible over grid and points) ------
  biplot_plot |>
    biplotEZ::samples(opacity = 0, which = NULL) |>
    biplotEZ::axes(col            = "grey22",
                   label.dir      = label_dir,
                   which          = which,
                   X.names        = X_names,
                   tick.label.cex = tick_label_cex,
                   ticks          = ticks_v,
                   label.line     = label_line_vec) |>
    plot(add = TRUE)

  # ---- Step 5: decision boundary contour lines -------------------------
  if (!isTRUE(no_contour)) {
    for (cl in gr$ct) {
      graphics::lines(cl$x, cl$y,
                      col = contour_col,
                      lwd = contour_lwd,
                      lty = contour_lty)
    }
  }

  # ---- Step 6: target observation (optional) ---------------------------
  if (!is.null(target_point)) {
    target_vec <- as.numeric(target_point[var_names])
    if (any(is.na(target_vec)))
      warning("Some var_names not found in target_point; those entries are NA.",
              call. = FALSE)

    target_st <- target_vec - bl_result$X_center
    if (isTRUE(bl_result$standardise))
      target_st <- target_st / bl_result$X_sd

    target_z <- matrix(target_st, nrow = 1L) %*%
      bl_result$V[, proj_dims, drop = FALSE]

    graphics::points(x   = target_z[1L, 1L],
                     y   = target_z[1L, 2L],
                     pch = 21L,
                     bg  = "yellow",
                     col = "black",
                     cex = 1.8)

    if (!is.null(target_label)) {
      graphics::text(x      = target_z[1L, 1L],
                     y      = target_z[1L, 2L],
                     labels = as.character(target_label),
                     cex    = 0.75,
                     font   = 2L)
    }
  }

  # ---- Step 7: prediction summary -------------------------------------
  n_total  <- length(points$pred_class)
  n_pos    <- sum(points$pred_class == 1L)
  n_neg    <- n_total - n_pos

  cat("\n--- Prediction summary ---\n")
  cat(sprintf("  Points plotted : %d\n", n_total))
  cat(sprintf("  Predicted 1    : %d  (%.1f %%)\n",
              n_pos, 100 * n_pos / n_total))
  cat(sprintf("  Predicted 0    : %d  (%.1f %%)\n",
              n_neg, 100 * n_neg / n_total))

  if (!is.null(points$class)) {
    correct  <- sum(points$pred_class == points$class)
    accuracy <- correct / n_total
    cat(sprintf("  Accuracy       : %d / %d  (%.1f %%)\n",
                correct, n_total, 100 * accuracy))
    tp <- sum(points$pred_class == 1L & points$class == 1L)
    tn <- sum(points$pred_class == 0L & points$class == 0L)
    fp <- sum(points$pred_class == 1L & points$class == 0L)
    fn <- sum(points$pred_class == 0L & points$class == 1L)
    cat(sprintf("  TP / TN / FP / FN : %d / %d / %d / %d\n", tp, tn, fp, fn))
  } else {
    cat("  (True labels unknown — accuracy not available)\n")
  }
  cat("--------------------------\n")

  invisible(bl_result)
}
