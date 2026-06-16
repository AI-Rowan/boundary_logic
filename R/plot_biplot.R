############################################################
# plot_biplotEZ(): Phase 1 biplot visualisation
# Based on: scripts/1.6 Plot_biplotEZ.R
# Uses the biplotEZ pipeline for the base plot and axis rendering.
############################################################

# ---- Private helpers -------------------------------------------------------

.make_label_line_vec <- function(label_offset_var, label_offset_dist,
                                  num_vars, var_names) {
  if (is.character(label_offset_var)) {
    idx <- match(label_offset_var, var_names)
    bad <- label_offset_var[is.na(idx)]
    if (length(bad) > 0L)
      warning(sprintf(
        "label_offset_var: variable(s) not found: %s",
        paste(bad, collapse = ", ")), call. = FALSE)
    label_offset_var <- idx[!is.na(idx)]
  }
  vec       <- rep(0.1, num_vars)
  valid_idx <- label_offset_var[label_offset_var >= 1L &
                                label_offset_var <= num_vars]
  if (length(valid_idx) > 0L)
    vec[valid_idx] <- rep_len(label_offset_dist, length(valid_idx))
  vec
}

.make_ticks_vec <- function(ticks_v, ticks_var, ticks_n, num_vars, var_names) {
  if (is.character(ticks_var)) {
    idx <- match(ticks_var, var_names)
    bad <- ticks_var[is.na(idx)]
    if (length(bad) > 0L)
      warning(sprintf(
        "ticks_var: variable(s) not found: %s",
        paste(bad, collapse = ", ")), call. = FALSE)
    ticks_var <- idx[!is.na(idx)]
  }
  vec       <- rep(ticks_v, num_vars)
  valid_idx <- ticks_var[ticks_var >= 1L & ticks_var <= num_vars]
  if (length(valid_idx) > 0L)
    vec[valid_idx] <- rep_len(ticks_n, length(valid_idx))
  vec
}

.apply_biplot_rotation <- function(biplot_obj, rotate_deg, proj_dims) {
  if (is.null(rotate_deg) || rotate_deg == 0)
    return(list(biplot_obj = biplot_obj, R_mat = NULL))
  theta              <- -rotate_deg * pi / 180
  R_mat              <- matrix(c(cos(theta), sin(theta),
                                 -sin(theta), cos(theta)), nrow = 2L)
  V_rot              <- biplot_obj$Lmat
  V_rot[, proj_dims] <- V_rot[, proj_dims] %*% R_mat
  tV_rot             <- solve(V_rot)
  tVr_rot            <- tV_rot[proj_dims, , drop = FALSE]
  biplot_obj$Lmat[, proj_dims] <- V_rot[, proj_dims]
  biplot_obj$ax.one.unit       <- (1 / diag(t(tVr_rot) %*% tVr_rot)) * t(tVr_rot)
  biplot_obj$Z[, proj_dims]    <- biplot_obj$Z[, proj_dims] %*% R_mat
  list(biplot_obj = biplot_obj, R_mat = R_mat)
}

#' Relabel biplot axis ticks into the original (raw) feature units
#'
#' biplotEZ computes axis tick *labels* as an affine function of the object's
#' `$means` and `$sd` (see `.calibrate.axis()` in biplotEZ); the 2D geometry
#' (`$Z`, `$Lmat`, `$ax.one.unit`) is independent of those two fields. When the
#' data fed to the package was standardised as `std = (raw - center) / scale`,
#' the labels read in standardised units. Transforming
#' `means_new = means * scale + center` and `sd_new = sd * scale`
#' (and forcing `$scaled`/`$center` TRUE so biplotEZ reconstructs the tick range
#' on the raw scale) makes every axis read in raw units while leaving the plotted
#' geometry byte-identical. This holds uniformly for PCA (`scaled` TRUE/FALSE)
#' and CVA (`scaled` FALSE, where `sd = 1` so `sd_new = scale`).
#'
#' @param biplot_obj A biplotEZ S3 object (possibly already rotated).
#' @param scaling    The `bl_result$scaling` list (`center`, `scale`, `method`)
#'   or `NULL`.
#' @param var_names  Feature names in the column order of the biplot's data
#'   matrix; `scaling$center`/`scale` are aligned to this order.
#' @return The biplot object with `$means`/`$sd` relabelled, or unchanged when
#'   `scaling` is `NULL`.
#' @noRd
.bl_rescale_biplot_axes <- function(biplot_obj, scaling, var_names) {
  if (is.null(scaling)) return(biplot_obj)

  center <- scaling$center[var_names]
  scale  <- scaling$scale[var_names]
  if (any(is.na(center)) || any(is.na(scale))) {
    warning("bl_result$scaling does not cover all biplot features; ",
            "axes left in model (standardised) units.", call. = FALSE)
    return(biplot_obj)
  }
  center <- unname(center)
  scale  <- unname(scale)

  biplot_obj$means  <- biplot_obj$means * scale + center
  biplot_obj$sd     <- biplot_obj$sd * scale
  biplot_obj$scaled <- TRUE
  biplot_obj$center <- TRUE
  biplot_obj
}

# ---------------------------------------------------------------------------

#' Plot the Phase 1 biplot
#'
#' Creates a biplot visualisation of the Phase 1 result using the biplotEZ
#' rendering pipeline. The plot shows:
#' \itemize{
#'   \item **Base axes** -- biplotEZ variable axes drawn as a canvas.
#'   \item **Coloured prediction grid** -- each grid point coloured blue
#'     (class 0) through white to red (class 1) by model probability.
#'   \item **Training data points** -- coloured by confusion category.
#'   \item **Variable axes redrawn on top** -- so axes are visible over
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
#' @param plot_points      Logical; if `TRUE` (default), training data points
#'   are drawn. Set to `FALSE` to hide them.
#' @param confusion_cols   Logical; if `TRUE` (default), points are coloured
#'   by confusion category (TP=red, TN=blue, FP=purple, FN=orange) when true
#'   labels are available. Set to `FALSE` to colour points by predicted class
#'   only (class 1 = red, class 0 = blue), ignoring true labels.
#' @param no_contour       Logical; if `TRUE`, decision boundary contour lines
#'   are hidden. Default `FALSE`.
#' @param new_title        Character; overrides the biplot title stored in
#'   `bl_result`. Use `NA` (default) to keep the stored title.
#' @param cex_z            Numeric; size of training data points (`cex`).
#'   Default `0.5`.
#' @param label_dir        Character; biplotEZ axis label direction.
#'   `"Paral"` (default) follows the plot border -- labels on the left/right
#'   borders are vertical, top/bottom are horizontal. `"Hor"` forces all labels
#'   horizontal. `"Orthog"` draws labels perpendicular to the axis direction.
#' @param label_cex        Numeric; variable name (axis label) size. Default `1`.
#' @param tick_label_cex   Numeric; axis tick label size. Default `0.6`.
#' @param ticks_v          Integer; number of ticks per variable axis.
#'   Default `1L`.
#' @param ticks_var        Integer, integer vector, character, or character
#'   vector; variable(s) whose tick-mark count should override `ticks_v`.
#'   Supply variable names (e.g. `"Sepal.Length"`) or integer indices
#'   (e.g. `c(1L, 3L)`). `0` (default) means no override -- every axis uses
#'   `ticks_v`.
#' @param ticks_n          Integer or integer vector; tick-mark count(s) for
#'   the variable(s) named in `ticks_var`. A single value applies to all
#'   listed variables; a vector must match the length of `ticks_var`.
#'   Default `5L`.
#' @param which            Integer vector; indices of variables to draw axes
#'   for. Defaults to all variables.
#' @param X_names          Character vector; custom variable names for axis
#'   labels. Defaults to `bl_result$var_names`.
#' @param label_offset_var Integer, integer vector, character, or character
#'   vector; variable(s) whose axis labels should be shifted outward from the
#'   border. Supply variable names (e.g., `"Sepal.Length"`) or integer indices
#'   (e.g., `c(1L, 3L)`). `0` (default) means no offset.
#' @param label_offset_dist Numeric or numeric vector; outward offset distance(s)
#'   in margin lines for `label_offset_var`. A single value applies to all
#'   listed variables; a vector must match the length of `label_offset_var`.
#'   Useful range 1--3. Default `1.5`.
#' @param rotate_deg       Numeric; angle in degrees to rotate the entire plot
#'   clockwise. Rotates all plotted elements -- grid, training points, contour
#'   lines, variable axes, and the target point -- without rerunning the
#'   pipeline. Default `0` (no rotation).
#' @param grid_col         `NULL` (default) or a length-2 character vector
#'   `c(class0_col, class1_col)`. When supplied, the prediction grid is
#'   rendered as hard binary colours -- `grid_col[1]` for grid points with
#'   predicted probability < 0.5 (class 0) and `grid_col[2]` for >= 0.5
#'   (class 1) -- instead of the default blue-white-red gradient. Example:
#'   `grid_col = c("blue", "red")`.
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
#' @importFrom grDevices col2rgb colorRampPalette rgb dev.cur dev.new
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
                           plot_points       = TRUE,
                           confusion_cols    = TRUE,
                           no_contour        = FALSE,
                           new_title         = NA,
                           cex_z             = 0.5,
                           label_dir         = "Paral",
                           label_cex         = 1,
                           tick_label_cex    = 0.6,
                           ticks_v           = 1L,
                           ticks_var         = 0L,
                           ticks_n           = 5L,
                           which             = NULL,
                           X_names           = NULL,
                           label_offset_var  = 0L,
                           label_offset_dist = 1.5,
                           rotate_deg        = 0,
                           grid_col          = NULL,
                           contour_col       = "black",
                           contour_lwd       = 1.5,
                           contour_lty       = 1L) {

  if (!inherits(bl_result, "bl_result")) {
    message("plot_biplotEZ(): 'bl_result' must be a 'bl_result' object from bl_assemble(). ",
            "Did bl_build_result() return NULL due to a missing 'bl_data' argument?")
    return(invisible(NULL))
  }

  label_dir <- match.arg(label_dir, c("Paral", "Hor", "Orthog"))

  if (grDevices::dev.cur() == 1L) grDevices::dev.new()

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
  label_line_vec <- .make_label_line_vec(label_offset_var, label_offset_dist,
                                         num_vars, var_names)
  ticks_vec <- .make_ticks_vec(ticks_v, ticks_var, ticks_n, num_vars, var_names)

  # ---- Optional rotation (clockwise by rotate_deg degrees) -------------
  rot         <- .apply_biplot_rotation(biplot_plot, rotate_deg, proj_dims)
  biplot_plot <- rot$biplot_obj
  R_mat       <- rot$R_mat

  # ---- Optional raw-unit axis relabel (display only) -------------------
  biplot_plot <- .bl_rescale_biplot_axes(biplot_plot, bl_result$scaling, var_names)
  if (!is.null(R_mat)) {
    gr$Zgrid <- gr$Zgrid %*% R_mat
    points$Z <- points$Z %*% R_mat
    gr$ct    <- lapply(gr$ct, function(cl) {
      pts  <- cbind(cl$x, cl$y) %*% R_mat
      cl$x <- pts[, 1L]; cl$y <- pts[, 2L]; cl
    })
  }

  # ---- Step 1: biplotEZ base plot (axes only, no samples) --------------
  # Establishes the coordinate system and draws the axis skeleton.
  biplot_plot |>
    biplotEZ::samples(opacity = 0, which = NULL) |>
    biplotEZ::axes(col            = "grey",
                   label.dir      = label_dir,
                   label.cex      = label_cex,
                   which          = which,
                   X.names        = X_names,
                   tick.label.cex = tick_label_cex,
                   ticks          = ticks_vec,
                   label.line     = label_line_vec) |>
    plot()

  # ---- Step 2: prediction grid (coloured squares) ----------------------
  if (!isTRUE(no_grid) && !is.null(gr)) {
    if (!is.null(grid_col) && length(grid_col) == 2L) {
      grid_draw_col <- ifelse(gr$grid_prob >= 0.5, grid_col[2L], grid_col[1L])
    } else {
      grid_draw_col <- gr$col_value
    }
    graphics::points(gr$Zgrid,
                     type = "p",
                     col  = grid_draw_col,
                     pch  = 15L,
                     cex  = 0.5)
  }

  # ---- Step 3: data points ---------------------------------------------
  if (isTRUE(plot_points)) {
    pt_col <- if (!isTRUE(confusion_cols)) {
      ifelse(points$pred_class == 1L, "red", "blue")
    } else {
      points$pred_col
    }
    graphics::points(x   = points$Z[, 1L],
                     y   = points$Z[, 2L],
                     col = pt_col,
                     pch = 16L,
                     cex = cex_z)
  }

  # ---- Step 4: axes redrawn on top (visible over grid and points) ------
  bp_overlay <- biplot_plot |>
    biplotEZ::samples(opacity = 0, which = NULL) |>
    biplotEZ::axes(col            = "grey22",
                   label.dir      = label_dir,
                   label.cex      = label_cex,
                   which          = which,
                   X.names        = X_names,
                   tick.label.cex = tick_label_cex,
                   ticks          = ticks_vec,
                   label.line     = label_line_vec)
  graphics::par(new = TRUE)
  plot(bp_overlay)

  # ---- Step 5: decision boundary contour lines -------------------------
  if (!isTRUE(no_contour) && !is.null(gr)) {
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
    if (!is.null(R_mat)) target_z <- target_z %*% R_mat

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

  # ---- Step 7: summary -------------------------------------------------
  n_total <- nrow(points$Z)

  if (!is.null(points$pred_class)) {
    # Model available: prediction summary
    n_pos <- sum(points$pred_class == 1L)
    n_neg <- n_total - n_pos
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
      cat("  (True labels unknown -- accuracy not available)\n")
    }
    cat("--------------------------\n")
  } else if (!is.null(points$class)) {
    # No model: class balance summary from true labels
    n_pos <- sum(points$class == 1L)
    n_neg <- n_total - n_pos
    cat("\n--- Class summary (no model) ---\n")
    cat(sprintf("  Points plotted : %d\n", n_total))
    cat(sprintf("  Class 1        : %d  (%.1f %%)\n",
                n_pos, 100 * n_pos / n_total))
    cat(sprintf("  Class 0        : %d  (%.1f %%)\n",
                n_neg, 100 * n_neg / n_total))
    cat("--------------------------------\n")
  }

  invisible(bl_result)
}


#' Plot method for bl_result objects
#'
#' Thin wrapper around \code{\link{plot_biplotEZ}}. All parameters are
#' passed through via \code{...}.
#'
#' @param x   A \code{"bl_result"} object from \code{bl_assemble()}.
#' @param ... Arguments forwarded to \code{\link{plot_biplotEZ}}.
#'
#' @return Invisibly returns \code{x}. Called for its side-effect (plot).
#' @export
plot.bl_result <- function(x, ...) plot_biplotEZ(x, ...)
