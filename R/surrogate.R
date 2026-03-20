############################################################
# bl_surrogate(): biplot-region surrogate model
# plot.bl_surrogate(): surrogate region biplot
#
# Refactored from: scripts/4.2.1 Create Grouping.R
#                  scripts/4.2.2 Biplot Grouping CT.R
############################################################


#' Biplot-region surrogate model
#'
#' Assigns each observation's predicted class using the nearest contour line
#' in Z-space rather than the fitted classifier. The surrogate provides a
#' purely spatial approximation of the decision boundary confined to the
#' training-data hull.
#'
#' @section How the surrogate works:
#' The surrogate uses `ct_surrogate` contour lines always computed and stored
#' in `bl_result$biplot_grid$ct_surrogate` by `bl_build_grid()`. These
#' contours are generated with cells outside the training-data convex hull
#' set to `cutoff` before extraction, so all contour lines are fully enclosed
#' within the hull.
#'
#' Each observation is first projected to Z-space. Observations that fall
#' **outside** the training-data convex hull polygon are excluded from the
#' surrogate classification — they receive `surrogate_pred = NA` — because
#' the surrogate regions are defined only within the hull. A message reports
#' how many observations were excluded. Accuracy is computed on in-hull
#' observations only.
#'
#' In-hull observations are assigned to the nearest contour segment in
#' Z-space. The contour's probability level determines the surrogate class:
#' level < `cutoff` → class 0; level >= `cutoff` → class 1.
#'
#' Accuracy is reported against (in-hull observations only):
#' * The fitted model's predictions (`accuracy_vs_model`)
#' * True class labels when present in `data` (`accuracy_vs_labels`)
#'
#' @param bl_result A `"bl_result"` object from `bl_assemble()`.
#' @param data      Data frame to classify. Defaults to
#'   `bl_result$train_data`. Must contain all columns in `var_names`.
#'
#' @return A list of class `"bl_surrogate"` with components:
#' \describe{
#'   \item{`surrogate_pred`}{Integer vector; surrogate class per obs (0 or 1).}
#'   \item{`model_pred`}{Integer vector; model-predicted class per obs.}
#'   \item{`class_obs`}{Integer vector; true labels (`NA` if absent in `data`).}
#'   \item{`accuracy_vs_model`}{Numeric; agreement rate with model predictions.}
#'   \item{`accuracy_vs_labels`}{Numeric; accuracy vs true labels
#'     (`NA` if no labels available).}
#'   \item{`boundary_set_nr`}{Integer vector; contour index assigned per obs.}
#'   \item{`z_boundary_type`}{Numeric vector; probability level per retained contour.}
#'   \item{`Z_obs`}{n x 2 matrix; Z-space coordinates of each observation.}
#'   \item{`bl_result`}{Reference to the parent `bl_result`.}
#' }
#'
#' @importFrom sp point.in.polygon
#' @export
bl_surrogate <- function(bl_result, data = NULL) {

  if (!inherits(bl_result, "bl_result"))
    stop("'bl_result' must be a 'bl_result' object from bl_assemble().",
         call. = FALSE)

  ct_surr <- bl_result$biplot_grid$ct_surrogate
  if (length(ct_surr) == 0L)
    stop("No surrogate contours found in bl_result$biplot_grid$ct_surrogate.",
         call. = FALSE)

  # ---- Unpack --------------------------------------------------------
  V           <- bl_result$V
  proj_dims   <- bl_result$proj_dims
  X_center    <- bl_result$X_center
  X_sd        <- bl_result$X_sd
  standardise <- bl_result$standardise
  cutoff      <- bl_result$cutoff
  var_names   <- bl_result$var_names

  polygon     <- bl_result$polygon

  Vr <- V[, proj_dims, drop = FALSE]

  # ---- Data ----------------------------------------------------------
  if (is.null(data)) data <- bl_result$train_data
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)
  n <- nrow(data)

  # ---- Project to Z-space -------------------------------------------
  sv   <- X_sd
  if (!isTRUE(standardise)) sv <- FALSE
  X_st  <- scale(as.matrix(data[, var_names, drop = FALSE]),
                 center = X_center, scale = sv)
  Z_obs <- X_st %*% Vr
  colnames(Z_obs) <- c("x", "y")

  # ---- Model predictions --------------------------------------------
  model_pred_prob <- .pred_function(
    model_use  = bl_result$model,
    model_type = bl_result$model_type,
    new_data   = data[, var_names, drop = FALSE]
  )
  model_pred <- as.integer(model_pred_prob >= cutoff)
  class_obs  <- if ("class" %in% names(data)) as.integer(data[["class"]]) else
    rep(NA_integer_, n)

  # ---- Filter to observations inside the training polygon ------------
  # The surrogate regions are defined only within the training-data convex
  # hull. Observations projected outside the hull have no valid surrogate
  # region and are excluded from predictions (surrogate_pred = NA).
  in_hull <- rep(TRUE, n)
  if (inherits(polygon, "SpatialPolygons")) {
    in_hull <- sp::point.in.polygon(
      Z_obs[, 1L], Z_obs[, 2L],
      polygon@polygons[[1L]]@Polygons[[1L]]@coords[, 1L],
      polygon@polygons[[1L]]@Polygons[[1L]]@coords[, 2L]
    ) > 0L
    n_excluded <- sum(!in_hull)
    if (n_excluded > 0L)
      message(sprintf(
        "%d observation(s) lie outside the training-data polygon and will have surrogate_pred = NA.",
        n_excluded
      ))
  }
  idx_in  <- which(in_hull)
  Z_in    <- Z_obs[idx_in, , drop = FALSE]

  # ---- Build surrogate boundary list from ct_surrogate --------------
  nr_raw            <- length(ct_surr)
  z_boundaries_list <- vector("list", nr_raw)
  z_boundary_type   <- numeric(nr_raw)

  for (i in seq_len(nr_raw)) {
    Mi <- cbind(ct_surr[[i]]$x, ct_surr[[i]]$y)
    colnames(Mi) <- c("x", "y")
    # Clip to training-data polygon
    if (inherits(polygon, "SpatialPolygons")) {
      inside <- sp::point.in.polygon(
        Mi[, 1L], Mi[, 2L],
        polygon@polygons[[1L]]@Polygons[[1L]]@coords[, 1L],
        polygon@polygons[[1L]]@Polygons[[1L]]@coords[, 2L]
      )
      Mi <- Mi[inside > 0L, , drop = FALSE]
    }
    z_boundaries_list[[i]] <- Mi
    z_boundary_type[i]     <- ct_surr[[i]]$level
  }

  # Retain non-empty boundaries only
  non_empty <- vapply(z_boundaries_list,
                      function(M) !is.null(M) && nrow(M) > 0L, logical(1L))
  z_boundaries_list <- z_boundaries_list[non_empty]
  z_boundary_type   <- z_boundary_type[non_empty]
  nr_boundaries     <- length(z_boundaries_list)

  if (nr_boundaries == 0L)
    stop("No surrogate contour segments remain after polygon clipping.",
         call. = FALSE)

  # ---- For each in-hull obs, find nearest contour segment -----------
  # Following 4.2.2 logic: find nearest point on each boundary, then
  # select the boundary whose nearest point is globally closest.
  # Only in-hull observations are processed; out-of-hull get NA.
  n_in <- length(idx_in)
  nearest_pts_list <- vector("list", nr_boundaries)
  for (i in seq_len(nr_boundaries)) {
    Mi  <- z_boundaries_list[[i]]
    idx <- .nearest_idx_block(Z_in, Mi)
    idx[is.na(idx)] <- 1L   # guard against empty boundary
    nearest_pts_list[[i]] <- Mi[idx, , drop = FALSE]
  }

  z_boundary_row      <- matrix(0, nrow = nr_boundaries, ncol = 2L)
  boundary_set_nr_in  <- integer(n_in)

  for (j in seq_len(n_in)) {
    for (i in seq_len(nr_boundaries)) {
      z_boundary_row[i, ] <- nearest_pts_list[[i]][j, ]
    }
    d2 <- rowSums(
      (z_boundary_row -
         matrix(Z_in[j, ], nrow = nr_boundaries, ncol = 2L, byrow = TRUE))^2
    )
    boundary_set_nr_in[j] <- which.min(d2)
  }

  # ---- Surrogate class from contour level ---------------------------
  boundary_class     <- as.integer(z_boundary_type >= cutoff)
  surrogate_pred_in  <- boundary_class[boundary_set_nr_in]

  # Map back to full n-length vector (NA for out-of-hull)
  surrogate_pred  <- rep(NA_integer_, n)
  boundary_set_nr <- rep(NA_integer_, n)
  surrogate_pred[idx_in]  <- surrogate_pred_in
  boundary_set_nr[idx_in] <- boundary_set_nr_in

  # ---- Accuracy (in-hull observations only) --------------------------
  accuracy_vs_model  <- mean(surrogate_pred == model_pred, na.rm = TRUE)
  accuracy_vs_labels <- if (!all(is.na(class_obs))) {
    mean(surrogate_pred == class_obs, na.rm = TRUE)
  } else {
    NA_real_
  }

  # ---- Return --------------------------------------------------------
  structure(
    list(
      surrogate_pred     = surrogate_pred,
      model_pred         = model_pred,
      class_obs          = class_obs,
      accuracy_vs_model  = accuracy_vs_model,
      accuracy_vs_labels = accuracy_vs_labels,
      boundary_set_nr    = boundary_set_nr,
      z_boundary_type    = z_boundary_type,
      Z_obs              = Z_obs,
      bl_result          = bl_result
    ),
    class = "bl_surrogate"
  )
}


# --------------------------------------------------------------------------
# S3 print method
# --------------------------------------------------------------------------

#' @export
print.bl_surrogate <- function(x, ...) {
  cat("<bl_surrogate>\n")
  cat(sprintf("  Observations        : %d\n", nrow(x$Z_obs)))
  cat(sprintf("  Surrogate regions   : %d\n", length(x$z_boundary_type)))
  cat(sprintf("  Accuracy vs model   : %.4f\n", x$accuracy_vs_model))
  if (!is.na(x$accuracy_vs_labels)) {
    cat(sprintf("  Accuracy vs labels  : %.4f\n", x$accuracy_vs_labels))
  } else {
    cat("  Accuracy vs labels  : (no labels in data)\n")
  }
  invisible(x)
}


# --------------------------------------------------------------------------
# S3 plot method
# --------------------------------------------------------------------------

#' Plot surrogate model region assignments on the biplot
#'
#' Renders the prediction surface background (coloured squares, clipped to
#' the training-data convex hull), overlays observation points coloured by
#' their surrogate region assignment, and draws the hull-clipped surrogate
#' contour lines. No background colour is shown outside the hull.
#'
#' The rendering order mirrors `plot_biplotEZ()`:
#' 1. biplotEZ axis skeleton (axes only, no sample points).
#' 2. Hull-clipped prediction grid (blue-white-red squares, inside hull only).
#' 3. Surrogate-coloured observation points (red = class 1, blue = class 0).
#' 4. Axes redrawn on top so labels remain visible over the grid.
#' 5. Hull-clipped surrogate contour lines in black.
#'
#' @param x   A `"bl_surrogate"` object from `bl_surrogate()`.
#' @param ... Unused; for S3 compatibility.
#'
#' @importFrom graphics points lines legend
#' @importFrom sp over coordinates
#' @export
plot.bl_surrogate <- function(x, ...) {

  bl_result <- x$bl_result
  gr        <- bl_result$biplot_grid
  polygon   <- bl_result$polygon

  # ---- Step 1: biplotEZ axis skeleton ----------------------------------
  bl_result$biplot_obj |>
    biplotEZ::samples(opacity = 0, which = NULL) |>
    biplotEZ::axes(col = "grey", which = NULL) |>
    plot()

  # ---- Step 2: prediction grid, clipped to hull -----------------------
  Zgrid     <- gr$Zgrid
  col_value <- gr$col_value

  if (inherits(polygon, "SpatialPolygons")) {
    G_df <- as.data.frame(Zgrid)
    sp::coordinates(G_df) <- ~x + y
    inside_grid <- !is.na(sp::over(G_df, polygon))
    Zgrid     <- Zgrid[inside_grid, , drop = FALSE]
    col_value <- col_value[inside_grid]
  }

  graphics::points(Zgrid, type = "p", col = col_value, pch = 15L, cex = 0.5)

  # ---- Step 3: surrogate-coloured observation points ------------------
  surr_col <- ifelse(x$surrogate_pred == 1L, "red", "blue")
  graphics::points(x$Z_obs[, 1L], x$Z_obs[, 2L],
                   col = surr_col, pch = 16L, cex = 0.7)

  # ---- Step 4: axes redrawn on top ------------------------------------
  bl_result$biplot_obj |>
    biplotEZ::samples(opacity = 0, which = NULL) |>
    biplotEZ::axes(col = "grey22", which = NULL) |>
    plot(add = TRUE)

  # ---- Step 5: hull-clipped surrogate contour lines -------------------
  # ct_surrogate contours are already bounded within the training hull
  for (cl in gr$ct_surrogate) {
    graphics::lines(cl$x, cl$y, col = "black", lwd = 1.5)
  }

  # ---- Legend ----------------------------------------------------------
  graphics::legend(
    "topright",
    legend = c("Surrogate class 1", "Surrogate class 0"),
    col    = c("red", "blue"),
    pch    = 16L, bty = "n", cex = 0.85
  )

  # ---- Accuracy --------------------------------------------------------
  cat(sprintf("Surrogate accuracy vs model : %.4f\n", x$accuracy_vs_model))
  if (!is.na(x$accuracy_vs_labels)) {
    cat(sprintf("Surrogate accuracy vs labels: %.4f\n", x$accuracy_vs_labels))
  }

  invisible(x)
}
