############################################################
# bl_build_grid(): m x m prediction grid in Z-space
# Refactored from: scripts/1.1 biplot_plane_optimized.R
#
# Design notes:
#   - All computation is separated from plotting.
#   - The biplotEZ object is used only to read plot axis bounds
#     (the par("usr") trick). This is isolated in .bl_get_plot_bounds().
#   - Grid scoring uses 50,000-row chunks for memory safety.
############################################################

# --------------------------------------------------------------------------
# Private: extract square plot bounds from a biplotEZ object
# --------------------------------------------------------------------------

#' Read axis limits from a rendered biplotEZ plot
#'
#' Opens a temporary off-screen device, renders the biplot axes only
#' (no points, no grid), reads `par("usr")`, then closes the device.
#' Returns `c(min_val, max_val)` for a square plot.
#'
#' @param biplot_obj A `biplotEZ` S3 object.
#' @param Z_train    n x 2 matrix of projected training coordinates (fallback
#'   if biplotEZ rendering fails).
#'
#' @importFrom biplotEZ samples axes
#' @importFrom grDevices png dev.off
#' @importFrom graphics par
#' @return Numeric vector of length 2: `c(min_val, max_val)`.
#' @keywords internal
.bl_get_plot_bounds <- function(biplot_obj, Z_train) {
  tmp_file <- tempfile(fileext = ".png")
  grDevices::png(tmp_file, width = 500L, height = 500L)
  on.exit({
    grDevices::dev.off()
    unlink(tmp_file)
  }, add = TRUE)

  tryCatch({
    biplot_obj |>
      biplotEZ::samples(opacity = 0, which = NULL) |>
      biplotEZ::axes(col = "grey", which = NULL) |>
      plot()
    usr <- graphics::par("usr")
    c(min(usr[c(1, 3)]), max(usr[c(2, 4)]))
  }, error = function(e) {
    # Fallback to data range if biplotEZ render fails
    warning("Could not render biplot to get axis bounds; using data range.",
            call. = FALSE)
    c(min(Z_train), max(Z_train))
  })
}


# --------------------------------------------------------------------------
# Main exported function
# --------------------------------------------------------------------------

#' Build the biplot prediction grid
#'
#' Creates an m x m grid of points in the 2D biplot Z-space, back-projects
#' each point to the original feature space, scores every grid point through
#' the fitted model, and extracts contour lines that approximate the decision
#' boundary. Also computes the convex hull polygon bounding the training data.
#'
#' The grid and contours are the foundation for all Phase 2 boundary search
#' and counterfactual identification operations.
#'
#' @section Memory efficiency:
#' Grid scoring is performed in chunks of 50,000 rows to avoid large peak
#' memory allocations, matching the behaviour of the original optimised
#' research script.
#'
#' @section Contour lines:
#' Two sets of contour lines are always computed and stored:
#'
#' * **`ct`** — standard boundary contours at `cutoff ± b_margin`. These
#'   are used by `bl_find_boundary()` for counterfactual search.
#'
#' * **`ct_surrogate`** — hull-clipped contours for the surrogate model.
#'   Grid cells outside the training-data convex hull are set to `cutoff`
#'   before extraction, so contour lines are fully enclosed within the hull.
#'   Used exclusively by `bl_surrogate()`. Do not use these for biplot
#'   rendering or boundary search.
#'
#' @param train_data    Data frame with feature columns (post-filter training
#'   data). Used to project points into Z-space and to define the hull
#'   polygon when `polygon` is `NULL`.
#' @param bl_projection A `"bl_projection"` object from
#'   `bl_build_projection()`.
#' @param bl_model      A `"bl_model"` object from `bl_fit_model()`.
#' @param m             Integer; grid resolution (m x m points per side).
#'   Default `200L`. Increase for finer boundary resolution at the cost of
#'   computation time.
#' @param cutoff        Numeric; the decision threshold. Default `0.5`.
#' @param b_margin      Numeric; half-width of the band around `cutoff` used
#'   to extract contour lines. Default `NULL`, which resolves to
#'   `1 / (10^rounding)`.
#' @param rounding      Integer; floor-based rounding precision for grid
#'   scores. Default `2`.
#' @param polygon       An `sp::SpatialPolygons` object to use as the hull
#'   boundary. If `NULL` (default), a new polygon is computed from
#'   `train_data` using `aplpack::plothulls()` with `fraction = outlie`.
#' @param outlie        Numeric in (0, 1]; hull fraction used only when
#'   `polygon` is `NULL`. Default `0.9`.
#' @param calc_hull     Logical; if `TRUE`, grid points outside the polygon
#'   are removed from the returned grid. Default `FALSE`.
#'
#' @return A list of class `"bl_grid"` with components:
#' \describe{
#'   \item{`Zgrid`}{Numeric matrix (m² x 2); grid coordinates in Z-space.}
#'   \item{`Xgrid`}{Data frame (m² x p); grid back-projected to X-space.}
#'   \item{`grid_prob`}{Numeric vector (length m²); model scores per grid
#'     point, in [0, 1].}
#'   \item{`col_value`}{Character vector (length m²); colour per grid point
#'     on a blue-white-red scale.}
#'   \item{`min_val`}{Numeric; lower bound of the square plot region.}
#'   \item{`max_val`}{Numeric; upper bound of the square plot region.}
#'   \item{`polygon`}{`sp::SpatialPolygons` convex hull.}
#'   \item{`ct`}{List of contour line objects (boundary search contours).}
#'   \item{`ct_surrogate`}{List of contour line objects with hull-clipped
#'     boundaries for `bl_surrogate()`. Outside-hull cells are set to
#'     `cutoff` before extraction. Not used in biplot rendering.}
#'   \item{`xseq`}{Numeric vector (length m); x-axis grid sequence.}
#'   \item{`yseq`}{Numeric vector (length m); y-axis grid sequence.}
#' }
#'
#' @importFrom grDevices contourLines colorRampPalette
#' @importFrom dplyr case_when
#' @importFrom sp coordinates over
#'
#' @examples
#' bl_dat  <- bl_prepare_data(datasets::iris,
#'                             class_col    = "Species",
#'                             target_class = "versicolor")
#' bl_mod  <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
#'                          model_type = "GLM")
#' bl_proj <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
#'                                 method = "PCA")
#' bl_grid <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = 100L)
#'
#' @export
bl_build_grid <- function(train_data,
                          bl_projection,
                          bl_model,
                          m        = 200L,
                          cutoff   = 0.5,
                          b_margin = NULL,
                          rounding = 2L,
                          polygon  = NULL,
                          outlie   = 0.9,
                          calc_hull = FALSE) {

  # ---- Validation -------------------------------------------------------
  stop_if_not_data_frame(train_data, "train_data")
  if (!inherits(bl_projection, "bl_projection"))
    stop("'bl_projection' must be a 'bl_projection' object from bl_build_projection().",
         call. = FALSE)
  if (!inherits(bl_model, "bl_model"))
    stop("'bl_model' must be a 'bl_model' object from bl_fit_model().",
         call. = FALSE)
  stop_if_not_positive_integer(m, "m")
  stop_if_not_scalar_numeric(cutoff, "cutoff")

  if (is.null(b_margin)) b_margin <- 1 / (10^rounding)

  # ---- Unpack projection -----------------------------------------------
  V           <- bl_projection$V
  tV          <- bl_projection$tV
  proj_dims   <- bl_projection$proj_dims
  X_center    <- bl_projection$X_center
  X_sd        <- bl_projection$X_sd
  standardise <- bl_projection$standardise
  biplot_obj  <- bl_projection$biplot_obj
  var_names   <- bl_model$var_names

  Vr  <- V[, proj_dims, drop = FALSE]   # p x 2
  tVr <- tV[proj_dims, , drop = FALSE]  # 2 x p

  # ---- Project training data to Z-space --------------------------------
  sv   <- X_sd
  if (!isTRUE(standardise)) sv <- FALSE
  X_st <- scale(as.matrix(train_data[, var_names, drop = FALSE]),
                center = X_center, scale = sv)
  Z_train <- X_st %*% Vr
  colnames(Z_train) <- c("x", "y")

  # ---- Get plot bounds from biplotEZ -----------------------------------
  bounds  <- .bl_get_plot_bounds(biplot_obj, Z_train)
  min_val <- bounds[1L]
  max_val <- bounds[2L]

  # ---- Build m x m grid in Z-space -------------------------------------
  xseq  <- seq(from = min_val, to = max_val, length.out = m)
  yseq  <- xseq
  Zgrid <- as.matrix(expand.grid(xseq, yseq))
  colnames(Zgrid) <- c("x", "y")

  # ---- Back-project grid to X-space ------------------------------------
  Xgrid <- as.data.frame(Zgrid %*% tVr)
  if (isTRUE(standardise)) Xgrid <- sweep(Xgrid, 2L, X_sd, "*")
  Xgrid <- sweep(Xgrid, 2L, X_center, "+")
  colnames(Xgrid) <- var_names

  # ---- Build or reuse convex hull polygon ------------------------------
  if (is.null(polygon)) {
    polygon <- .build_hull_polygon(Z_train, min_val, max_val, outlie = outlie)
  }

  # Determine inside/outside for grid and training points
  G_df <- as.data.frame(Zgrid)
  sp::coordinates(G_df) <- ~x + y
  inside_grid <- !is.na(sp::over(G_df, polygon))

  Z_df <- as.data.frame(Z_train)
  sp::coordinates(Z_df) <- ~x + y
  inside_data <- !is.na(sp::over(Z_df, polygon))

  # ---- Score each grid point (chunked for memory safety) ---------------
  ##A : make use of the .pred_function to determine grid_prob
  chunk     <- 50000L
  grid_prob <- numeric(nrow(Xgrid))
  for (i0 in seq(1L, nrow(Xgrid), by = chunk)) {
    j0 <- min(i0 + chunk - 1L, nrow(Xgrid))
    grid_prob[i0:j0] <- .pred_function(
      model_use  = bl_model$model,
      model_type = bl_model$model_type,
      rounding   = rounding,
      new_data   = Xgrid[i0:j0, , drop = FALSE]
    )
  }

  # ---- Colour palette (blue -> white -> red) ---------------------------
  col_vec   <- grDevices::colorRampPalette(
    c("deepskyblue", "white", "lightsalmon")
  )(101L)
  col_value <- col_vec[floor(grid_prob * 100) + 1L]

  # ---- Contour lines ---------------------------------------------------
  grid_contour <- matrix(grid_prob, ncol = m, byrow = FALSE)

  # Standard boundary contours: used by bl_find_boundary()
  ct <- grDevices::contourLines(
    xseq, yseq, grid_contour,
    levels = c(cutoff - b_margin, cutoff + b_margin)
  )

  # Surrogate contours: outside-hull cells set to cutoff so contour lines
  # are fully enclosed within the training-data hull. Used only by
  # bl_surrogate() — do not use for biplot rendering or boundary search.
  grid_contour_surr <- grid_contour
  grid_contour_surr[!inside_grid] <- cutoff
  ct_surrogate <- grDevices::contourLines(
    xseq, yseq, grid_contour_surr,
    levels = c(cutoff - b_margin, cutoff, cutoff + b_margin)
  )
  ct_surrogate <- ct_surrogate[
    sapply(ct_surrogate, function(cl) cl$level != cutoff)
  ]

  # ---- Optional: remove grid points outside polygon --------------------
  # Data point colouring is handled by bl_project_points(); only the
  # prediction surface (grid) is filtered here.

  if (isTRUE(calc_hull)) {
    grid_prob <- grid_prob[inside_grid]
    col_value <- col_value[inside_grid]
    Zgrid     <- Zgrid[inside_grid, , drop = FALSE]
    Xgrid     <- Xgrid[inside_grid, , drop = FALSE]
  }

  # ---- Return ----------------------------------------------------------
  structure(
    list(
      Zgrid         = Zgrid,
      Xgrid         = Xgrid,
      grid_prob     = grid_prob,
      col_value     = col_value,
      min_val       = min_val,
      max_val       = max_val,
      polygon       = polygon,
      hull_fraction = outlie,
      ct            = ct,
      ct_surrogate  = ct_surrogate,
      xseq          = xseq,
      yseq          = yseq
    ),
    class = "bl_grid"
  )
}


# --------------------------------------------------------------------------
# S3 print method
# --------------------------------------------------------------------------

#' @export
print.bl_grid <- function(x, ...) {
  m <- length(x$xseq)
  cat("<bl_grid>\n")
  cat(sprintf("  Grid size      : %d x %d (%d points)\n",
              m, m, nrow(x$Zgrid)))
  cat(sprintf("  Prob range     : [%.4f, %.4f]\n",
              min(x$grid_prob, na.rm = TRUE),
              max(x$grid_prob, na.rm = TRUE)))
  cat(sprintf("  Contour lines  : %d  (surrogate: %d)\n",
              length(x$ct), length(x$ct_surrogate)))
  invisible(x)
}
