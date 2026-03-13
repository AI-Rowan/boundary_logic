############################################################
# bl_find_boundary(): global counterfactual search
# Refactored from: scripts/1.4 Biplot Boundary Search_optimized.R
############################################################


# --------------------------------------------------------------------------
# Private helpers (ported directly from research code)
# --------------------------------------------------------------------------

#' Block-wise nearest-neighbour indices (no n x m distance matrix)
#'
#' @param Z n x 2 matrix of query points.
#' @param B m x 2 matrix of candidate boundary points.
#' @param block Integer; block size for Z (default 5000).
#'
#' @return Integer vector of length n; row index of the nearest B row.
#' @keywords internal
.nearest_idx_block <- function(Z, B, block = 5000L) {
  Z <- as.matrix(Z); B <- as.matrix(B)
  if (ncol(Z) != 2L || ncol(B) != 2L)
    stop(".nearest_idx_block: Z and B must have 2 columns.", call. = FALSE)
  n <- nrow(Z); m <- nrow(B)
  if (m == 0L || n == 0L) return(rep(NA_integer_, n))
  best_idx <- rep.int(NA_integer_, n)
  best_d2  <- rep.int(Inf, n)

  from <- 1L
  while (from <= n) {
    to  <- min(from + block - 1L, n)
    Zb  <- Z[from:to, , drop = FALSE]

    b_from <- 1L
    stepB  <- max(1000L, as.integer(1e7 / max(1L, to - from + 1L)))
    while (b_from <= m) {
      b_to <- min(b_from + stepB - 1L, m)
      Bb   <- B[b_from:b_to, , drop = FALSE]
      if (nrow(Bb) == 0L) { b_from <- b_to + 1L; next }

      dx <- outer(Zb[, 1L], Bb[, 1L], "-")
      dy <- outer(Zb[, 2L], Bb[, 2L], "-")
      d2 <- dx * dx + dy * dy

      min_loc <- max.col(-d2)
      d2_min  <- d2[cbind(seq_len(nrow(Zb)), min_loc)]

      upd <- d2_min < best_d2[from:to]
      if (any(upd)) {
        best_d2[from:to][upd]  <- d2_min[upd]
        best_idx[from:to][upd] <- (b_from - 1L) + min_loc[upd]
      }
      b_from <- b_to + 1L
    }
    from <- to + 1L
  }
  best_idx
}

#' Filter polyline points to those inside a SpatialPolygons object
#'
#' @param Mxy n x 2 numeric matrix of (x, y) coordinates.
#' @param poly A `SpatialPolygons` object, or any other object (returned as-is).
#'
#' @return n' x 2 matrix with only the rows inside `poly`.
#' @keywords internal
.poly_clip <- function(Mxy, poly) {
  if (is.null(Mxy) || length(Mxy) == 0L)
    return(matrix(numeric(0L), ncol = 2L))
  if (!inherits(poly, "SpatialPolygons"))
    return(as.matrix(Mxy))
  inside <- sp::point.in.polygon(
    Mxy[, 1L], Mxy[, 2L],
    poly@polygons[[1L]]@Polygons[[1L]]@coords[, 1L],
    poly@polygons[[1L]]@Polygons[[1L]]@coords[, 2L]
  )
  Mxy[inside > 0L, , drop = FALSE]
}

#' Detect whether a polyline is closed (first row equals last row)
#'
#' @param Mxy Numeric matrix with at least 2 rows and 2 columns.
#' @return Logical scalar.
#' @keywords internal
.is_closed <- function(Mxy) {
  Mxy <- as.matrix(Mxy)
  if (nrow(Mxy) < 2L) return(FALSE)
  (Mxy[1L, 1L] == Mxy[nrow(Mxy), 1L]) && (Mxy[1L, 2L] == Mxy[nrow(Mxy), 2L])
}


# --------------------------------------------------------------------------
# Main exported function
# --------------------------------------------------------------------------

#' Find the nearest decision boundary point for each observation
#'
#' For each observation in `data`, finds the nearest point on the decision
#' boundary in biplot Z-space, then back-projects it to the original feature
#' space to produce a counterfactual. The search uses the contour lines
#' stored in `bl_result$biplot_grid$ct`.
#'
#' @section Counterfactual logic:
#' A class-1 observation's counterfactual is the nearest point on a class-0
#' contour (and vice versa). Boundary segments are filtered to the training
#' data variable ranges (`bl_result$train_ranges`) to ensure counterfactuals
#' remain within the observed feature space.
#'
#' @param bl_result A `"bl_result"` object from `bl_assemble()`.
#' @param data      Data frame to process. Defaults to
#'   `bl_result$test_data`. Must contain all columns in `var_names`. Can
#'   also be `bl_result$train_data` for a global training-data search.
#' @param tdp       Integer vector of row indices to process. `NULL`
#'   (default) processes all rows of `data`.
#'
#' @return A list of class `"bl_boundary"` with components:
#' \describe{
#'   \item{`B_z`}{n x 2 matrix; Z-space boundary point per observation.}
#'   \item{`B_x`}{Data frame (n x p); X-space counterfactual.}
#'   \item{`B_pred`}{Numeric vector; model probability at counterfactual
#'     (should be close to `cutoff` — validates back-projection quality).}
#'   \item{`dist_z`}{Numeric vector; Euclidean distance obs to boundary
#'     in Z-space.}
#'   \item{`dist_x`}{Numeric vector; Euclidean distance obs to boundary
#'     in standardised X-space (each dimension divided by `X_sd`).}
#'   \item{`Z_obs`}{n x 2 matrix; Z-space coordinates of each observation.}
#'   \item{`x_obs`}{Data frame; original observations (feature columns).}
#'   \item{`pred_obs`}{Numeric vector; model probability at each observation.}
#'   \item{`class_obs`}{Integer vector; actual 0/1 labels (`NA` if absent).}
#'   \item{`z_boundaries_list`}{List of clipped contour segment matrices.}
#'   \item{`z_boundary_type`}{Numeric vector; probability level per contour.}
#'   \item{`nr_boundaries`}{Integer; number of active contour segments.}
#'   \item{`position_store`}{n x 2 integer matrix; contour index used per obs.}
#'   \item{`bl_result`}{Reference to the parent `bl_result` object.}
#' }
#'
#' @importFrom sp point.in.polygon
#' @export
bl_find_boundary <- function(bl_result, data = NULL, tdp = NULL) {

  if (!inherits(bl_result, "bl_result"))
    stop("'bl_result' must be a 'bl_result' object from bl_assemble().",
         call. = FALSE)

  # ---- Unpack --------------------------------------------------------
  V            <- bl_result$V
  tV           <- bl_result$tV
  proj_dims    <- bl_result$proj_dims
  X_center     <- bl_result$X_center
  X_sd         <- bl_result$X_sd
  standardise  <- bl_result$standardise
  cutoff       <- bl_result$cutoff
  rounding     <- bl_result$rounding
  var_names    <- bl_result$var_names
  polygon      <- bl_result$polygon
  train_ranges <- bl_result$train_ranges
  ct           <- bl_result$biplot_grid$ct

  Vr  <- V[, proj_dims, drop = FALSE]   # p x 2
  tVr <- tV[proj_dims, , drop = FALSE]  # 2 x p

  # ---- Resolve data --------------------------------------------------
  if (is.null(data)) data <- bl_result$test_data
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)
  missing_vars <- setdiff(var_names, names(data))
  if (length(missing_vars))
    stop(sprintf("'data' is missing columns: %s.",
                 paste(missing_vars, collapse = ", ")), call. = FALSE)

  if (is.null(tdp)) tdp <- seq_len(nrow(data))
  data_sub <- data[tdp, , drop = FALSE]
  n        <- nrow(data_sub)

  # ---- Project observations to Z-space -------------------------------
  sv   <- X_sd
  if (!isTRUE(standardise)) sv <- FALSE
  X_st  <- scale(as.matrix(data_sub[, var_names, drop = FALSE]),
                 center = X_center, scale = sv)
  Z_obs <- X_st %*% Vr
  colnames(Z_obs) <- c("x", "y")

  # ---- Model predictions per obs ------------------------------------
  pred_obs       <- .pred_function(
    model_use  = bl_result$model,
    model_type = bl_result$model_type,
    rounding   = rounding,
    new_data   = data_sub[, var_names, drop = FALSE]
  )
  pred_class_obs <- as.integer(pred_obs >= cutoff)

  class_obs <- if ("class" %in% names(data_sub)) {
    as.integer(data_sub[["class"]])
  } else {
    rep(NA_integer_, n)
  }

  # ---- Helper: build empty return ----------------------------------
  .empty_return <- function() {
    structure(list(
      B_z               = matrix(0, nrow = n, ncol = 2L),
      B_x               = as.data.frame(matrix(NA_real_, nrow = n,
                            ncol = length(var_names),
                            dimnames = list(NULL, var_names))),
      B_pred            = rep(NA_real_, n),
      dist_z            = rep(NA_real_, n),
      dist_x            = rep(NA_real_, n),
      Z_obs             = Z_obs,
      x_obs             = data_sub[, var_names, drop = FALSE],
      pred_obs          = pred_obs,
      class_obs         = class_obs,
      z_boundaries_list = list(),
      z_boundary_type   = numeric(0L),
      nr_boundaries     = 0L,
      position_store    = matrix(0L, nrow = n, ncol = 2L),
      bl_result         = bl_result
    ), class = "bl_boundary")
  }

  if (length(ct) == 0L) {
    message("No contour lines found in bl_result. Returning empty bl_boundary.")
    return(.empty_return())
  }

  # ---- Clip contours to polygon ------------------------------------
  nr_raw            <- length(ct)
  z_boundaries_list <- vector("list", nr_raw)
  z_boundary_type   <- numeric(nr_raw)

  for (i in seq_len(nr_raw)) {
    Mi <- cbind(ct[[i]]$x, ct[[i]]$y)
    colnames(Mi) <- c("x", "y")
    if (inherits(polygon, "SpatialPolygons")) Mi <- .poly_clip(Mi, polygon)
    z_boundaries_list[[i]] <- Mi
    z_boundary_type[i]     <- ct[[i]]$level
  }

  # ---- Retain boundaries that contain/enclose data -----------------
  boundary_used <- integer(0L)
  for (i in seq_len(nr_raw)) {
    Mi <- z_boundaries_list[[i]]
    if (is.null(Mi) || nrow(Mi) == 0L) next
    keep <- TRUE
    if (.is_closed(Mi) && inherits(polygon, "SpatialPolygons")) {
      inside_ct <- sp::point.in.polygon(
        Z_obs[, 1L], Z_obs[, 2L], Mi[, 1L], Mi[, 2L]
      )
      keep <- any(inside_ct > 0L)
    }
    if (keep) boundary_used <- c(boundary_used, i)
  }
  if (length(boundary_used) == 0L) boundary_used <- seq_len(nr_raw)

  # Rebuild with retained boundaries only
  z_boundary_type   <- z_boundary_type[boundary_used]
  z_boundaries_list <- lapply(boundary_used, function(j) {
    Mi <- cbind(ct[[j]]$x, ct[[j]]$y)
    colnames(Mi) <- c("x", "y")
    if (inherits(polygon, "SpatialPolygons")) Mi <- .poly_clip(Mi, polygon)
    Mi
  })
  nr_boundaries <- length(z_boundaries_list)

  # ---- Consistency pruning: keep segments where model agrees --------
  # Back-project each boundary segment to X-space, apply train_ranges
  # filter, score through the model, keep only segments whose model
  # prediction is consistent with the contour's probability level.
  for (i in seq_len(nr_boundaries)) {
    Mi <- z_boundaries_list[[i]]
    if (is.null(Mi) || nrow(Mi) == 0L) next

    Bx <- as.data.frame(Mi %*% tVr)
    if (isTRUE(standardise)) Bx <- sweep(Bx, 2L, X_sd, "*")
    Bx <- sweep(Bx, 2L, X_center, "+")
    colnames(Bx) <- var_names

    keep_rows <- if (is.list(train_ranges)) {
      get_filter_logical_vector(Bx, train_ranges)
    } else {
      rep(TRUE, nrow(Bx))
    }

    if (!any(keep_rows)) {
      z_boundaries_list[[i]] <- Mi[0L, , drop = FALSE]
      next
    }

    Bx_keep <- Bx[keep_rows, , drop = FALSE]
    Mi_keep <- Mi[keep_rows, , drop = FALSE]

    p_bnd <- .pred_function(
      model_use  = bl_result$model,
      model_type = bl_result$model_type,
      rounding   = rounding,
      new_data   = Bx_keep
    )
    ok <- if (z_boundary_type[i] < cutoff) (p_bnd < cutoff) else (p_bnd >= cutoff)
    z_boundaries_list[[i]] <- if (any(ok)) {
      Mi_keep[ok, , drop = FALSE]
    } else {
      Mi[0L, , drop = FALSE]
    }
  }

  # Drop emptied boundaries
  non_empty <- vapply(z_boundaries_list,
                      function(M) !is.null(M) && nrow(M) > 0L, logical(1L))
  z_boundaries_list <- z_boundaries_list[non_empty]
  z_boundary_type   <- z_boundary_type[non_empty]
  nr_boundaries     <- length(z_boundaries_list)

  if (nr_boundaries == 0L) {
    message("No valid boundary segments remain after pruning. ",
            "Returning empty bl_boundary.")
    return(.empty_return())
  }

  # ---- Nearest boundary point per obs, per boundary ----------------
  B_boundary_z <- vector("list", nr_boundaries)
  for (i in seq_len(nr_boundaries)) {
    Mi  <- z_boundaries_list[[i]]
    idx <- .nearest_idx_block(Z_obs, Mi)
    # Guard against NA indices (empty boundary)
    idx[is.na(idx)] <- 1L
    B_boundary_z[[i]] <- Mi[idx, , drop = FALSE]
  }

  # ---- Choose opposing-side contour per obs ------------------------
  # B0 = nearest class-1 contour point per obs
  # B1 = nearest class-0 contour point per obs
  # class-1 predicted obs → counterfactual is on class-0 side (B1)
  # class-0 predicted obs → counterfactual is on class-1 side (B0)
  boundary_class_1 <- z_boundary_type >= cutoff
  boundary_class_0 <- !boundary_class_1

  z_boundary_row <- matrix(0, nrow = nr_boundaries, ncol = 2L)
  B0 <- matrix(0, nrow = n, ncol = 2L)
  B1 <- matrix(0, nrow = n, ncol = 2L)
  position_store <- matrix(0L, nrow = n, ncol = 2L)

  for (j in seq_len(n)) {
    for (i in seq_len(nr_boundaries)) {
      Bij <- B_boundary_z[[i]][j, ]
      z_boundary_row[i, ] <- if (length(Bij) == 0L || any(!is.finite(Bij))) {
        c(1e9, 1e9)
      } else {
        Bij
      }
    }

    d2_all <- rowSums(
      (z_boundary_row -
         matrix(Z_obs[j, ], nrow = nr_boundaries, ncol = 2L, byrow = TRUE))^2
    )
    d2_1 <- d2_all[boundary_class_1]
    d2_0 <- d2_all[boundary_class_0]

    pos1 <- if (length(d2_1)) which.min(d2_1) else NA_integer_
    pos0 <- if (length(d2_0)) which.min(d2_0) else NA_integer_

    if (is.na(pos1)) {
      position_store[j, 1L] <- 0L; B0[j, ] <- c(NA_real_, NA_real_)
    } else {
      k1 <- which(boundary_class_1)[pos1]
      position_store[j, 1L] <- k1
      B0[j, ] <- B_boundary_z[[k1]][j, ]
    }

    if (is.na(pos0)) {
      position_store[j, 2L] <- 0L; B1[j, ] <- c(NA_real_, NA_real_)
    } else {
      k0 <- which(boundary_class_0)[pos0]
      position_store[j, 2L] <- k0
      B1[j, ] <- B_boundary_z[[k0]][j, ]
    }
  }

  # Combine: class-1 obs gets class-0 boundary (B1); class-0 obs gets B0
  group1 <- pred_class_obs == 1L
  group0 <- !group1
  B0_out <- B0; B1_out <- B1
  B0_out[group1, ] <- 0
  B1_out[group0, ] <- 0
  B_z <- B0_out + B1_out
  # Replace rows that are still exactly zero with NA (no boundary found)
  zero_rows <- rowSums(B_z^2) == 0
  if (any(zero_rows)) B_z[zero_rows, ] <- NA_real_

  # ---- Back-project boundary to X-space ----------------------------
  B_x <- as.data.frame(B_z %*% tVr)
  if (isTRUE(standardise)) B_x <- sweep(B_x, 2L, X_sd, "*")
  B_x <- sweep(B_x, 2L, X_center, "+")
  colnames(B_x) <- var_names

  # Score counterfactuals (should be near cutoff)
  B_pred <- tryCatch(
    .pred_function(
      model_use  = bl_result$model,
      model_type = bl_result$model_type,
      rounding   = rounding,
      new_data   = B_x
    ),
    error = function(e) rep(NA_real_, n)
  )

  # ---- Distances ---------------------------------------------------
  dist_z <- sqrt(rowSums((Z_obs - B_z)^2, na.rm = FALSE))

  x_obs_mat  <- as.matrix(data_sub[, var_names, drop = FALSE])
  B_x_mat    <- as.matrix(B_x)
  diff_x_std <- sweep(x_obs_mat - B_x_mat, 2L, X_sd, "/")
  dist_x     <- sqrt(rowSums(diff_x_std^2, na.rm = FALSE))

  # ---- Return ------------------------------------------------------
  structure(
    list(
      B_z               = B_z,
      B_x               = B_x,
      B_pred            = B_pred,
      dist_z            = dist_z,
      dist_x            = dist_x,
      Z_obs             = Z_obs,
      x_obs             = data_sub[, var_names, drop = FALSE],
      pred_obs          = pred_obs,
      class_obs         = class_obs,
      z_boundaries_list = z_boundaries_list,
      z_boundary_type   = z_boundary_type,
      nr_boundaries     = nr_boundaries,
      position_store    = position_store,
      bl_result         = bl_result
    ),
    class = "bl_boundary"
  )
}


# --------------------------------------------------------------------------
# S3 print method
# --------------------------------------------------------------------------

#' @export
print.bl_boundary <- function(x, ...) {
  cat("<bl_boundary>\n")
  cat(sprintf("  Observations   : %d\n", nrow(x$Z_obs)))
  cat(sprintf("  Boundaries used: %d\n", x$nr_boundaries))
  cat(sprintf("  Dist Z (mean)  : %.4f\n", mean(x$dist_z, na.rm = TRUE)))
  cat(sprintf("  Dist X (mean)  : %.4f\n", mean(x$dist_x, na.rm = TRUE)))
  cat(sprintf("  B_pred range   : [%.4f, %.4f]\n",
              min(x$B_pred, na.rm = TRUE), max(x$B_pred, na.rm = TRUE)))
  invisible(x)
}
