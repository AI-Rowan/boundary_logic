############################################################
# Phase 3: Local counterfactual analysis
# bl_select_target(), set_filters(), bl_find_boundary_local()
# Adapted from: scripts/1.3 Optimal Rotation.R
############################################################


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

#' Generate eigenvector index pairs in canonical order
#'
#' @param max_dim  Maximum dimension index (number of eigenvectors available).
#' @param max_pairs Maximum number of pairs to return.
#' @return Integer matrix with 2 columns; each row is one (i, j) pair.
#' @noRd
.generate_pairs <- function(max_dim, max_pairs) {
  pairs <- vector("list", max_pairs)
  count <- 0L
  for (j in 2L:max_dim) {
    for (i in 1L:(j - 1L)) {
      count <- count + 1L
      pairs[[count]] <- c(i, j)
      if (count >= max_pairs) {
        return(do.call(rbind, pairs[seq_len(count)]))
      }
    }
  }
  if (count == 0L) return(matrix(integer(0), nrow = 0L, ncol = 2L))
  do.call(rbind, pairs[seq_len(count)])
}


#' Rotate the biplot to align with the target observation
#'
#' Adapts \code{Biplot_rotation()} from the research script.
#'
#' @param x_target_st Named numeric vector (length p) of the standardised
#'   target point.
#' @param V p x p loading matrix.
#' @param proj_pair Integer vector of length 2 giving the eigenvector indices.
#' @return Named list with \code{Vr_rot} (p x 2) and \code{tVr_rot} (2 x p).
#' @noRd
.bl_rotate <- function(x_target_st, V, proj_pair) {
  p    <- length(x_target_st)
  Vr   <- V[, proj_pair, drop = FALSE]           # p x 2

  Y    <- rbind(-x_target_st, rep(0, p), x_target_st)   # 3 x p
  YV   <- Y %*% V                                        # 3 x p
  YVr  <- Y %*% Vr                                       # 3 x 2

  # Pad YVr with zeros so it is 3 x p
  YVr_padded <- cbind(YVr, matrix(0, nrow = 3L, ncol = p - 2L))

  svd_res <- svd(t(YV) %*% YVr_padded)           # p x p SVD
  A       <- svd_res$v %*% t(svd_res$u)           # p x p rotation

  Vrho  <- V %*% t(A)                             # rotated full loading matrix
  tVrho <- solve(Vrho)

  list(
    Vr_rot  = Vrho[, proj_pair, drop = FALSE],    # p x 2
    tVr_rot = tVrho[proj_pair, , drop = FALSE]    # 2 x p
  )
}


#' Apply actionability constraints to back-projected contour points
#'
#' @param Bx           n-row data frame of back-projected X-space contour
#'   points.
#' @param x_obs_named  Named numeric vector of observed feature values.
#' @param filters      A \code{"bl_filters"} object or NULL.
#' @return Logical vector of length \code{nrow(Bx)}.
#' @noRd
.apply_actionability <- function(Bx, x_obs_named, filters) {
  if (is.null(filters) || length(filters$constraints) == 0L) {
    return(rep(TRUE, nrow(Bx)))
  }
  keep <- rep(TRUE, nrow(Bx))
  for (nm in names(filters$constraints)) {
    v <- filters$constraints[[nm]]
    if (is.character(v)) {
      keep <- keep & switch(
        v,
        "decrease" = Bx[[nm]] <= x_obs_named[[nm]],
        "increase" = Bx[[nm]] >= x_obs_named[[nm]],
        "fixed"    = abs(Bx[[nm]] - x_obs_named[[nm]]) < 1e-6,
        stop(sprintf("Unknown constraint type '%s' for variable '%s'.", v, nm),
             call. = FALSE)
      )
    } else if (is.numeric(v) && length(v) == 2L) {
      keep <- keep & Bx[[nm]] >= v[1L] & Bx[[nm]] <= v[2L]
    } else {
      stop(sprintf(
        "Invalid constraint for '%s': must be 'decrease', 'increase', 'fixed', or c(min, max).",
        nm), call. = FALSE)
    }
  }
  keep
}


# ---------------------------------------------------------------------------
# Exported functions
# ---------------------------------------------------------------------------

#' Select a target observation for local counterfactual analysis
#'
#' Prepares a single observation as the target for
#' \code{\link{bl_find_boundary_local}}. The target can be specified either as
#' an integer row index into a data frame, or as a single-row data frame
#' containing all feature columns.
#'
#' @param bl_result A \code{"bl_result"} object produced by
#'   \code{\link{bl_assemble}}.
#' @param target Integer row index into \code{data}, or a single-row data frame
#'   with all feature columns.
#' @param data Data frame to index when \code{target} is an integer. Defaults
#'   to \code{bl_result$test_data}.
#'
#' @return A list of class \code{"bl_target"} with the following fields:
#'   \describe{
#'     \item{x_obs}{Single-row data frame of observed feature values.}
#'     \item{z_obs}{1 x 2 matrix of projected coordinates in Z-space.}
#'     \item{pred_prob}{Scalar predicted probability.}
#'     \item{pred_class}{Integer predicted class (0 or 1).}
#'     \item{row_id}{Integer row index, or \code{NA_integer_} when
#'       \code{target} is a data frame.}
#'   }
#'
#' @export
bl_select_target <- function(bl_result, target, data = NULL) {
  if (!inherits(bl_result, "bl_result"))
    stop("'bl_result' must be a 'bl_result' object.", call. = FALSE)

  var_names   <- bl_result$var_names
  X_center    <- bl_result$X_center
  X_sd        <- bl_result$X_sd
  standardise <- bl_result$standardise
  proj_dims   <- bl_result$proj_dims
  V           <- bl_result$V
  cutoff      <- bl_result$cutoff
  rounding    <- bl_result$rounding

  # Resolve x_obs and row_id
  if (is.data.frame(target)) {
    if (nrow(target) != 1L)
      stop("When 'target' is a data frame it must have exactly 1 row.", call. = FALSE)
    missing_vars <- setdiff(var_names, names(target))
    if (length(missing_vars) > 0L)
      stop(sprintf("'target' data frame is missing columns: %s",
                   paste(missing_vars, collapse = ", ")), call. = FALSE)
    x_obs  <- target[, var_names, drop = FALSE]
    row_id <- NA_integer_
  } else if (is.numeric(target) || is.integer(target)) {
    row_id <- as.integer(target[[1L]])
    if (is.null(data)) data <- bl_result$test_data
    if (is.null(data))
      stop("No 'data' supplied and 'bl_result$test_data' is NULL.", call. = FALSE)
    if (row_id < 1L || row_id > nrow(data))
      stop(sprintf("'target' row index %d is out of range [1, %d].",
                   row_id, nrow(data)), call. = FALSE)
    x_obs <- data[row_id, var_names, drop = FALSE]
  } else {
    stop("'target' must be an integer row index or a single-row data frame.",
         call. = FALSE)
  }

  # Project to Z-space
  sv    <- if (isTRUE(standardise)) X_sd else FALSE
  X_st  <- scale(as.matrix(x_obs), center = X_center, scale = sv)
  Z_obs <- X_st %*% V[, proj_dims, drop = FALSE]
  colnames(Z_obs) <- c("x", "y")

  # Score
  pred_prob  <- .pred_function(
    model_use  = bl_result$model,
    model_type = bl_result$model_type,
    rounding   = rounding,
    new_data   = x_obs
  )[[1L]]
  pred_class <- as.integer(pred_prob >= cutoff)

  structure(
    list(
      x_obs      = x_obs,
      z_obs      = Z_obs,
      pred_prob  = pred_prob,
      pred_class = pred_class,
      row_id     = row_id
    ),
    class = "bl_target"
  )
}


#' Set actionability constraints for local counterfactual search
#'
#' Defines per-variable constraints that restrict the counterfactual search in
#' \code{\link{bl_find_boundary_local}}. Constraints are applied in X-space
#' after back-projecting contour segments.
#'
#' @param bl_target A \code{"bl_target"} object from
#'   \code{\link{bl_select_target}}.
#' @param ... Named per-variable constraints. Each argument name must be a
#'   variable in the model. Valid constraint values are:
#'   \describe{
#'     \item{\code{"decrease"}}{The counterfactual value must be \eqn{\le} the
#'       observed value.}
#'     \item{\code{"increase"}}{The counterfactual value must be \eqn{\ge} the
#'       observed value.}
#'     \item{\code{"fixed"}}{The counterfactual value must equal the observed
#'       value (within \code{1e-6}).}
#'     \item{\code{c(min, max)}}{The counterfactual value must lie within the
#'       specified range.}
#'   }
#'
#' @return A list of class \code{"bl_filters"} with fields:
#'   \describe{
#'     \item{constraints}{Named list of per-variable constraints.}
#'     \item{bl_target}{The \code{"bl_target"} object passed in.}
#'   }
#'
#' @export
set_filters <- function(bl_target, ...) {
  if (!inherits(bl_target, "bl_target"))
    stop("'bl_target' must be a 'bl_target' object from bl_select_target().",
         call. = FALSE)

  constraints <- list(...)
  var_names   <- names(bl_target$x_obs)

  unknown <- setdiff(names(constraints), var_names)
  if (length(unknown) > 0L)
    stop(sprintf(
      "Unknown variable name(s) in constraints: %s. Available: %s",
      paste(unknown, collapse = ", "),
      paste(var_names, collapse = ", ")
    ), call. = FALSE)

  valid_chars <- c("decrease", "increase", "fixed")
  for (nm in names(constraints)) {
    v <- constraints[[nm]]
    if (is.character(v)) {
      match.arg(v, valid_chars)
    } else if (is.numeric(v) && length(v) == 2L) {
      # ok
    } else {
      stop(sprintf(
        "Constraint for '%s' must be one of %s or a numeric vector c(min, max).",
        nm, paste(sprintf("'%s'", valid_chars), collapse = ", ")
      ), call. = FALSE)
    }
  }

  structure(
    list(constraints = constraints, bl_target = bl_target),
    class = "bl_filters"
  )
}


#' Find the nearest decision boundary for a single target observation
#'
#' Searches across multiple eigenvector pairs for the closest decision boundary
#' to the target observation in the rotated biplot space. For each pair, the
#' loading matrix is rotated so that the target observation lies along the first
#' rotated axis, then a fine grid is scored and the boundary contour is
#' extracted.
#'
#' @section Rotation:
#' The rotation adapts \code{Biplot_rotation()} from the research script
#' (\file{scripts/1.3 Optimal Rotation.R}). An SVD-based orthogonal rotation
#' matrix \eqn{A} is computed such that the rotated loading matrix aligns the
#' target observation with the projection plane defined by \code{proj_pair}.
#' This gives each eigenvector pair the best possible chance of revealing a
#' nearby boundary.
#'
#' @section Constraints:
#' The Z-space training hull polygon is \strong{not} applied here because the
#' per-pair rotation invalidates the polygon from \code{bl_build_grid()}.
#' Two constraint types are applied instead:
#' \describe{
#'   \item{train_ranges}{Back-projected X-space values must lie within the
#'     per-variable ranges observed in the training data.}
#'   \item{set_filters}{Actionability constraints supplied via
#'     \code{\link{set_filters}}.}
#' }
#'
#' @param bl_result   A \code{"bl_result"} object.
#' @param bl_target   A \code{"bl_target"} object from
#'   \code{\link{bl_select_target}}.
#' @param set_filters A \code{"bl_filters"} object from
#'   \code{\link{set_filters}}, or \code{NULL} (no actionability constraints).
#' @param max_pairs   Integer; maximum number of eigenvector pairs to try.
#'   Default \code{10L}.
#' @param m           Integer; grid resolution along each axis. Default
#'   \code{200L}.
#' @param verbose     Logical; print progress messages. Default \code{TRUE}.
#'
#' @return A list of class \code{"bl_local_result"} with the following fields:
#'   \describe{
#'     \item{B_z}{1 x 2 matrix; best boundary point in rotated Z-space.}
#'     \item{B_x}{Data frame (1 row); back-projected boundary in X-space.}
#'     \item{B_pred}{Scalar predicted probability at the boundary.}
#'     \item{dist_z}{Euclidean distance from target to boundary in Z-space.}
#'     \item{Z_target}{1 x 2 matrix; target coordinates in rotated Z-space.}
#'     \item{best_pair}{Integer vector of length 2; winning eigenvector pair.}
#'     \item{Vr_rot}{p x 2 rotated loading matrix for the best pair.}
#'     \item{tVr_rot}{2 x p inverse loading matrix for the best pair.}
#'     \item{Z_train_rot}{n x 2 matrix; training data in rotated Z-space.}
#'     \item{Zgrid}{Grid coordinates in rotated Z-space.}
#'     \item{grid_prob}{Numeric vector of model scores at each grid point.}
#'     \item{col_value}{Colour vector for grid points.}
#'     \item{ct_local}{Raw contour list from \code{grDevices::contourLines}.}
#'     \item{xseq, yseq}{Grid axis sequences.}
#'     \item{min_val, max_val}{Grid axis limits.}
#'     \item{all_distances}{Named numeric vector; best Z-distance per pair.}
#'     \item{solution_found}{Logical.}
#'     \item{blocking_constraint}{Character message if no solution, else NULL.}
#'     \item{bl_target}{The \code{"bl_target"} object.}
#'     \item{bl_result}{The \code{"bl_result"} object.}
#'   }
#'
#' @importFrom grDevices contourLines colorRampPalette
#' @export
bl_find_boundary_local <- function(bl_result, bl_target,
                                   set_filters = NULL,
                                   max_pairs   = 10L,
                                   m           = 200L,
                                   verbose     = TRUE) {

  if (!inherits(bl_result, "bl_result"))
    stop("'bl_result' must be a 'bl_result' object.", call. = FALSE)
  if (!inherits(bl_target, "bl_target"))
    stop("'bl_target' must be a 'bl_target' object from bl_select_target().",
         call. = FALSE)
  if (!is.null(set_filters) && !inherits(set_filters, "bl_filters"))
    stop("'set_filters' must be a 'bl_filters' object from set_filters() or NULL.",
         call. = FALSE)

  V            <- bl_result$V
  X_center     <- bl_result$X_center
  X_sd         <- bl_result$X_sd
  standardise  <- bl_result$standardise
  cutoff       <- bl_result$cutoff
  rounding     <- bl_result$rounding
  var_names    <- bl_result$var_names
  p            <- bl_result$num_vars
  train_ranges <- bl_result$train_ranges
  b_margin     <- 1 / (10^rounding)

  x_obs <- as.numeric(bl_target$x_obs[, var_names])

  # Standardise target and training data
  sv           <- if (isTRUE(standardise)) X_sd else rep(1, p)
  x_target_st  <- (x_obs - X_center) / sv
  X_train_st   <- scale(
    as.matrix(bl_result$train_data[, var_names, drop = FALSE]),
    center = X_center,
    scale  = if (isTRUE(standardise)) X_sd else FALSE
  )

  pairs <- .generate_pairs(p, as.integer(max_pairs))
  if (nrow(pairs) == 0L)
    stop("Cannot generate any eigenvector pairs from the available dimensions.",
         call. = FALSE)

  pair_names    <- apply(pairs, 1L, function(r) paste0("(", r[1], ",", r[2], ")"))
  all_distances <- setNames(rep(NA_real_, nrow(pairs)), pair_names)
  best_dist     <- Inf
  best_result   <- NULL

  for (k in seq_len(nrow(pairs))) {
    pair <- pairs[k, ]
    if (isTRUE(verbose))
      message(sprintf("  Pair %d/%d: eigenvectors (%d, %d) ...",
                      k, nrow(pairs), pair[1L], pair[2L]))

    # Rotation
    rot     <- .bl_rotate(x_target_st, V, pair)
    Vr_rot  <- rot$Vr_rot    # p x 2
    tVr_rot <- rot$tVr_rot   # 2 x p

    # Project target and training data into rotated Z-space
    Z_target    <- matrix(x_target_st %*% Vr_rot, nrow = 1L)
    colnames(Z_target) <- c("x", "y")
    Z_train_rot <- X_train_st %*% Vr_rot

    # Grid bounds from training data range + 10% padding
    z_range <- range(Z_train_rot)
    z_pad   <- diff(z_range) * 0.10
    min_val <- z_range[1L] - z_pad
    max_val <- z_range[2L] + z_pad

    # Build m x m grid in rotated Z-space
    xseq  <- seq(min_val, max_val, length.out = m)
    yseq  <- xseq
    Zgrid <- as.matrix(expand.grid(xseq, yseq))
    colnames(Zgrid) <- c("x", "y")

    # Back-project grid to X-space
    Xgrid <- as.data.frame(Zgrid %*% tVr_rot)
    if (isTRUE(standardise)) Xgrid <- sweep(Xgrid, 2L, X_sd, "*")
    Xgrid <- sweep(Xgrid, 2L, X_center, "+")
    colnames(Xgrid) <- var_names

    # Score grid (chunked)
    chunk     <- 50000L
    grid_prob <- numeric(nrow(Xgrid))
    for (i0 in seq(1L, nrow(Xgrid), by = chunk)) {
      j0 <- min(i0 + chunk - 1L, nrow(Xgrid))
      grid_prob[i0:j0] <- .pred_function(
        model_use  = bl_result$model,
        model_type = bl_result$model_type,
        rounding   = rounding,
        new_data   = Xgrid[i0:j0, , drop = FALSE]
      )
    }

    # Extract contour lines
    grid_mat <- matrix(grid_prob, ncol = m, byrow = FALSE)
    ct_local <- grDevices::contourLines(
      xseq, yseq, grid_mat,
      levels = c(cutoff - b_margin, cutoff + b_margin)
    )

    if (length(ct_local) == 0L) {
      if (isTRUE(verbose)) message("    No boundary contours found.")
      next
    }

    # Filter and validate contour segments
    target_pred_class <- bl_target$pred_class
    x_obs_named       <- setNames(x_obs, var_names)
    valid_Mi          <- list()

    for (i in seq_along(ct_local)) {
      lev <- ct_local[[i]]$level
      # Keep only opposing-class contours
      seg_class <- as.integer(lev >= cutoff)
      if (seg_class == target_pred_class) next

      Mi  <- cbind(ct_local[[i]]$x, ct_local[[i]]$y)
      colnames(Mi) <- c("x", "y")

      # Back-project segment to X-space
      Bx <- as.data.frame(Mi %*% tVr_rot)
      if (isTRUE(standardise)) Bx <- sweep(Bx, 2L, X_sd, "*")
      Bx <- sweep(Bx, 2L, X_center, "+")
      colnames(Bx) <- var_names

      # train_ranges filter
      keep <- if (is.list(train_ranges)) {
        get_filter_logical_vector(Bx, train_ranges)
      } else {
        rep(TRUE, nrow(Bx))
      }

      # Actionability filter (Phase 3 only — no hull polygon)
      if (!is.null(set_filters)) {
        keep <- keep & .apply_actionability(Bx, x_obs_named, set_filters)
      }

      if (!any(keep)) next

      # Model consistency check (same as global bl_find_boundary)
      Mi_keep <- Mi[keep, , drop = FALSE]
      Bx_keep <- Bx[keep, , drop = FALSE]
      p_bnd   <- .pred_function(
        model_use  = bl_result$model,
        model_type = bl_result$model_type,
        rounding   = rounding,
        new_data   = Bx_keep
      )
      ok <- if (lev < cutoff) (p_bnd < cutoff) else (p_bnd >= cutoff)
      if (!any(ok)) next

      valid_Mi[[length(valid_Mi) + 1L]] <- Mi_keep[ok, , drop = FALSE]
    }

    if (length(valid_Mi) == 0L) {
      if (isTRUE(verbose)) message("    All segments eliminated by constraints.")
      next
    }

    all_Mi <- do.call(rbind, valid_Mi)
    idx    <- .nearest_idx_block(Z_target, all_Mi, block = 1L)
    if (is.na(idx[[1L]])) next

    B_z_local <- all_Mi[idx[[1L]], , drop = FALSE]
    dist_z    <- sqrt(sum((Z_target - B_z_local)^2))

    all_distances[pair_names[k]] <- dist_z
    if (isTRUE(verbose))
      message(sprintf("    Distance: %.4f", dist_z))

    if (dist_z < best_dist) {
      best_dist <- dist_z

      # Back-project counterfactual to X-space
      B_x_mat <- B_z_local %*% tVr_rot
      if (isTRUE(standardise)) B_x_mat <- sweep(B_x_mat, 2L, X_sd, "*")
      B_x_mat <- sweep(B_x_mat, 2L, X_center, "+")
      B_x     <- as.data.frame(B_x_mat)
      colnames(B_x) <- var_names

      B_pred <- .pred_function(
        model_use  = bl_result$model,
        model_type = bl_result$model_type,
        rounding   = rounding,
        new_data   = B_x
      )

      col_vec   <- grDevices::colorRampPalette(
        c("deepskyblue", "white", "lightsalmon")
      )(101L)
      col_value <- col_vec[floor(grid_prob * 100) + 1L]

      best_result <- list(
        B_z         = B_z_local,
        B_x         = B_x,
        B_pred      = B_pred[[1L]],
        dist_z      = dist_z,
        Z_target    = Z_target,
        best_pair   = pair,
        Vr_rot      = Vr_rot,
        tVr_rot     = tVr_rot,
        Z_train_rot = Z_train_rot,
        Zgrid       = Zgrid,
        grid_prob   = grid_prob,
        col_value   = col_value,
        ct_local    = ct_local,
        xseq        = xseq,
        yseq        = yseq,
        min_val     = min_val,
        max_val     = max_val
      )
    }
  }

  solution_found <- !is.null(best_result)

  if (!solution_found) {
    if (!is.null(set_filters) && length(set_filters$constraints) > 0L) {
      blocking <- paste(
        "No counterfactual satisfying all constraints was found across",
        nrow(pairs), "eigenvector pair(s).",
        "Try relaxing or removing actionability filters in set_filters()."
      )
    } else {
      blocking <- paste(
        "No counterfactual found across", nrow(pairs), "eigenvector pair(s).",
        "The decision boundary may not cross the projected feature space.",
        "Try increasing max_pairs."
      )
    }
    message("bl_find_boundary_local: ", blocking)
  } else {
    blocking <- NULL
    if (isTRUE(verbose))
      message(sprintf("Best pair: (%d, %d) | Distance: %.4f",
                      best_result$best_pair[1L],
                      best_result$best_pair[2L],
                      best_result$dist_z))
  }

  null_pair_fields <- list(
    B_z = NULL, B_x = NULL, B_pred = NA_real_, dist_z = NA_real_,
    Z_target = NULL, best_pair = NA_integer_, Vr_rot = NULL,
    tVr_rot = NULL, Z_train_rot = NULL, Zgrid = NULL,
    grid_prob = NULL, col_value = NULL, ct_local = NULL,
    xseq = NULL, yseq = NULL, min_val = NULL, max_val = NULL
  )

  structure(
    c(
      if (solution_found) best_result else null_pair_fields,
      list(
        all_distances       = all_distances,
        solution_found      = solution_found,
        blocking_constraint = blocking,
        bl_target           = bl_target,
        bl_result           = bl_result
      )
    ),
    class = "bl_local_result"
  )
}


# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Plot the local biplot for a bl_local_result
#'
#' Renders the rotated biplot for the best eigenvector pair found by
#' \code{\link{bl_find_boundary_local}}. Training points are coloured by
#' confusion type (TP = red, TN = blue, FP = purple, FN = orange). The target
#' observation is shown as a yellow filled circle and the nearest boundary
#' point as a cross with an arrow.
#'
#' @param x          A \code{"bl_local_result"} object.
#' @param show_grid  Logical; overlay the scored grid as coloured tiles.
#'   Default \code{TRUE}.
#' @param show_axes  Logical; draw loading arrows. Default \code{TRUE}.
#' @param cex_z      Point size for training observations. Default \code{0.7}.
#' @param contour_col Colour for contour lines. Default \code{"black"}.
#' @param ...        Additional arguments (currently ignored).
#'
#' @importFrom graphics plot points lines arrows text
#' @importFrom dplyr case_when
#' @export
plot.bl_local_result <- function(x, show_grid = TRUE, show_axes = TRUE,
                                 cex_z = 0.7, contour_col = "black", ...) {
  if (!x$solution_found) {
    message("No solution found — nothing to plot.")
    return(invisible(x))
  }

  bl_result   <- x$bl_result
  bl_target   <- x$bl_target
  Z_train_rot <- x$Z_train_rot
  Z_target    <- x$Z_target
  B_z         <- x$B_z
  min_val     <- x$min_val
  max_val     <- x$max_val
  Vr_rot      <- x$Vr_rot
  tVr_rot     <- x$tVr_rot
  best_pair   <- x$best_pair

  var_names  <- bl_result$var_names
  cutoff     <- bl_result$cutoff
  rounding   <- bl_result$rounding

  # Confusion colours for training points
  train_df <- bl_result$train_data
  train_pred <- .pred_function(
    model_use  = bl_result$model,
    model_type = bl_result$model_type,
    rounding   = rounding,
    new_data   = train_df[, var_names, drop = FALSE]
  )
  train_pred_class <- as.integer(train_pred >= cutoff)

  # Derive true labels: try to get from train_data response column
  resp_col <- setdiff(names(train_df), var_names)
  if (length(resp_col) >= 1L) {
    true_class <- as.integer(as.character(train_df[[resp_col[1L]]]))
  } else {
    true_class <- rep(NA_integer_, nrow(train_df))
  }

  train_col <- dplyr::case_when(
    !is.na(true_class) & true_class == 1L & train_pred_class == 1L ~ "red",
    !is.na(true_class) & true_class == 0L & train_pred_class == 0L ~ "blue",
    !is.na(true_class) & true_class == 0L & train_pred_class == 1L ~ "purple",
    !is.na(true_class) & true_class == 1L & train_pred_class == 0L ~ "orange",
    TRUE ~ "grey50"
  )

  # Labels for plot title
  row_id    <- bl_target$row_id
  row_label <- if (is.na(row_id)) "external" else as.character(row_id)
  pred_prob  <- bl_target$pred_prob
  pred_class <- bl_target$pred_class

  # Base plot
  graphics::plot(
    NA,
    xlim = c(min_val, max_val),
    ylim = c(min_val, max_val),
    xlab = paste0("Rotated axis (", best_pair[1L], ")"),
    ylab = paste0("Rotated axis (", best_pair[2L], ")"),
    main = sprintf(
      "Local biplot \u2013 target %s  [p=%.3f, class %d]",
      row_label, pred_prob, pred_class
    )
  )

  # Grid
  if (show_grid) {
    graphics::points(x$Zgrid[, 1L], x$Zgrid[, 2L],
                     col = x$col_value, pch = 15L, cex = 0.35)
  }

  # Training points
  graphics::points(Z_train_rot[, 1L], Z_train_rot[, 2L],
                   col = train_col, pch = 16L, cex = cex_z)

  # Contour lines
  for (cl in x$ct_local) {
    graphics::lines(cl$x, cl$y, col = contour_col, lwd = 1.2)
  }

  # Loading arrows
  if (show_axes) {
    p         <- nrow(Vr_rot)
    ax_dir    <- t(tVr_rot)          # p x 2
    half_range <- (max_val - min_val) / 2
    ax_norms  <- sqrt(rowSums(ax_dir^2))
    ax_unit   <- ax_dir / pmax(ax_norms, 1e-12)
    ax_scaled <- ax_unit * 0.72 * half_range

    for (j in seq_len(p)) {
      graphics::arrows(
        0, 0, ax_scaled[j, 1L], ax_scaled[j, 2L],
        length = 0.07, col = "grey60", lwd = 1.0
      )
      graphics::text(
        ax_scaled[j, 1L] * 1.1, ax_scaled[j, 2L] * 1.1,
        labels = var_names[j], cex = 0.65, col = "grey40"
      )
    }
  }

  # Target point
  graphics::points(Z_target[, 1L], Z_target[, 2L],
                   pch = 21L, bg = "yellow", col = "black", cex = 1.8)
  if (!is.na(row_id)) {
    graphics::text(Z_target[, 1L], Z_target[, 2L],
                   labels = row_label, pos = 3L, cex = 0.7)
  }

  # CF cross
  graphics::points(B_z[, 1L], B_z[, 2L],
                   pch = 4L, col = "grey30", cex = 1.2, lwd = 1.5)

  # CF arrow
  graphics::arrows(
    Z_target[1L, 1L], Z_target[1L, 2L],
    B_z[1L, 1L],      B_z[1L, 2L],
    length = 0.08, col = "grey30", lwd = 1.2
  )

  invisible(x)
}


#' Print method for bl_local_result
#'
#' @param x A \code{"bl_local_result"} object.
#' @param ... Currently ignored.
#' @export
print.bl_local_result <- function(x, ...) {
  cat("-- bl_local_result --\n")
  cat(sprintf("  Solution found : %s\n", x$solution_found))
  if (x$solution_found) {
    cat(sprintf("  Best pair      : (%d, %d)\n",
                x$best_pair[1L], x$best_pair[2L]))
    cat(sprintf("  Distance (Z)   : %.4f\n", x$dist_z))
    cat(sprintf("  Boundary pred  : %.4f\n", x$B_pred))
  } else {
    cat(sprintf("  Reason         : %s\n", x$blocking_constraint))
  }
  cat("\n  Target:\n")
  print(x$bl_target)
  cat("\n  All Z-distances by pair:\n")
  ad <- x$all_distances
  df <- data.frame(pair = names(ad), distance = ad, row.names = NULL,
                   stringsAsFactors = FALSE)
  print(df, row.names = FALSE)
  invisible(x)
}


#' Print method for bl_target
#'
#' @param x A \code{"bl_target"} object.
#' @param ... Currently ignored.
#' @export
print.bl_target <- function(x, ...) {
  cat("-- bl_target --\n")
  rid <- if (is.na(x$row_id)) "external" else as.character(x$row_id)
  cat(sprintf("  Row ID         : %s\n", rid))
  cat(sprintf("  Predicted prob : %.4f\n", x$pred_prob))
  cat(sprintf("  Predicted class: %d\n",  x$pred_class))
  cat("  Feature values :\n")
  vals <- as.numeric(x$x_obs)
  names(vals) <- names(x$x_obs)
  print(round(vals, 4L))
  invisible(x)
}


#' Print method for bl_filters
#'
#' @param x A \code{"bl_filters"} object.
#' @param ... Currently ignored.
#' @export
print.bl_filters <- function(x, ...) {
  cat("-- bl_filters --\n")
  if (length(x$constraints) == 0L) {
    cat("  (no constraints set)\n")
  } else {
    for (nm in names(x$constraints)) {
      v <- x$constraints[[nm]]
      if (is.character(v)) {
        cat(sprintf("  %-20s : %s\n", nm, v))
      } else {
        cat(sprintf("  %-20s : [%.4g, %.4g]\n", nm, v[1L], v[2L]))
      }
    }
  }
  invisible(x)
}
