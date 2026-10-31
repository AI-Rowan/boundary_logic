############################################################
# Phase 3: Local counterfactual analysis
# bl_select_target(), set_filters(), bl_find_local_cf()
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
    Vr_rot  = Vrho[, c(1L, 2L), drop = FALSE],    # p x 2 -- SVD always puts target info in cols 1-2
    tVr_rot = tVrho[c(1L, 2L), , drop = FALSE]    # 2 x p
  )
}


# Search tolerance for a "fixed" constraint, in model units. A vertex of the
# 2D boundary contour rarely holds a feature at exactly its observed value,
# so the search accepts a change of up to this much; .restore_fixed() then
# tries the observed value itself.
.fixed_tol <- 0.5


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
        "fixed"    = abs(Bx[[nm]] - x_obs_named[[nm]]) <= .fixed_tol,
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


#' Warn when a data-frame target lies outside the training ranges
#'
#' A target given in raw units to a model fitted on standardised features is
#' the usual cause: `bl_set_scaling()` is display-only. When a scaling is
#' recorded and converting the target would leave fewer features out of range,
#' the warning suggests `bl_apply_scaling()`.
#'
#' @param x_obs     One-row data frame of the target's feature values.
#' @param bl_result A `"bl_result"` object.
#' @return `NULL`, invisibly; warns as a side effect.
#' @noRd
.warn_target_range <- function(x_obs, bl_result) {
  rng <- bl_result$train_ranges
  if (!is.list(rng)) return(invisible(NULL))
  vars <- intersect(names(rng), names(x_obs))
  out_of_range <- function(x) {
    vars[vapply(vars, function(v) {
      val <- x[[v]][1L]
      !is.na(val) && (val < rng[[v]][1L] || val > rng[[v]][2L])
    }, logical(1L))]
  }

  out_now <- out_of_range(x_obs)
  if (length(out_now) == 0L) return(invisible(NULL))

  msg <- sprintf(paste(
    "The target lies outside the training range for: %s.",
    "The local search looks for counterfactuals within that range."),
    paste(out_now, collapse = ", "))
  if (!is.null(bl_result$scaling)) {
    out_conv <- tryCatch(out_of_range(bl_apply_scaling(x_obs, bl_result)),
                         error = function(e) out_now)
    if (length(out_conv) < length(out_now))
      msg <- paste(msg,
                   "It looks like raw units: bl_set_scaling() is display-only, so",
                   "convert the target first with bl_apply_scaling(target, bl_result).")
  }
  warning(msg, call. = FALSE)
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Exported functions
# ---------------------------------------------------------------------------

#' Select a target observation for local counterfactual analysis
#'
#' Prepares a single observation as the target for
#' \code{\link{bl_find_local_cf}}. The target can be specified either as
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
#' @details The target's true class is read from the same row that supplies
#'   its features: \code{data[target, ]} for an index, or \code{target} itself
#'   for a data frame. The class column is \code{"class"} when present,
#'   otherwise the first non-feature column. When there is no such column the
#'   true class is unknown (\code{NA}).
#'
#'   A data-frame target must be in the model's units. When
#'   \code{\link{bl_set_scaling}} recorded a scaling, the model works in the
#'   standardised units, so convert a raw-unit target first with
#'   \code{\link{bl_apply_scaling}}. A data-frame target with any feature
#'   outside \code{bl_result$train_ranges} gives a warning, which suggests
#'   \code{bl_apply_scaling()} when the target looks like raw units.
#'
#' @return A list of class \code{"bl_target"} with the following fields:
#'   \describe{
#'     \item{x_obs}{Single-row data frame of observed feature values.}
#'     \item{z_obs}{1 x 2 matrix of projected coordinates in Z-space.}
#'     \item{pred_prob}{Scalar predicted probability.}
#'     \item{pred_class}{Integer predicted class (0 or 1).}
#'     \item{row_id}{Integer row index, or \code{NA_integer_} when
#'       \code{target} is a data frame.}
#'     \item{class_obs}{Integer true class of the target (0 or 1), or
#'       \code{NA_integer_} when the source has no class column. Used by
#'       \code{\link{plot.bl_local_result}} to colour the target by its
#'       confusion group.}
#'     \item{scaling}{The display-only raw-unit scaling from
#'       \code{bl_result$scaling}, or \code{NULL}.}
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
    .warn_target_range(x_obs, bl_result)
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

  # True class, read from the same row that supplies x_obs. Prefer the
  # package's "class" column; otherwise the first non-feature column. A
  # non-numeric label is treated as unknown (NA).
  src  <- if (is.data.frame(target)) target else data[row_id, , drop = FALSE]
  resp <- setdiff(names(src), var_names)
  if ("class" %in% resp) resp <- "class"
  class_obs <- if (length(resp) >= 1L) {
    suppressWarnings(as.integer(as.character(src[[resp[1L]]])))
  } else {
    NA_integer_
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
    new_data   = x_obs
  )[[1L]]
  pred_class <- as.integer(pred_prob >= cutoff)

  structure(
    list(
      x_obs      = x_obs,
      z_obs      = Z_obs,
      pred_prob  = pred_prob,
      pred_class = pred_class,
      row_id     = row_id,
      class_obs  = class_obs,
      scaling    = bl_result$scaling   # display-only raw-unit scaling (or NULL)
    ),
    class = "bl_target"
  )
}


#' Set actionability constraints for local counterfactual search
#'
#' Defines per-variable constraints that restrict the counterfactual search in
#' \code{\link{bl_find_local_cf}}. Constraints are applied in X-space
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
#'     \item{\code{"fixed"}}{The counterfactual value must be within \eqn{\pm
#'       0.5} of the observed value, in model units, during the boundary search
#'       (0.5 training SD for a standardised feature). It is a tolerance, not an
#'       equality, because a point on the 2D decision boundary rarely holds a
#'       feature at exactly its observed value. After the search,
#'       \code{\link{bl_find_local_cf}} sets each fixed feature it moved back to
#'       its observed value and keeps that point if it still crosses the
#'       boundary (see \strong{Fixed constraints} there); \code{print()} of the
#'       result explains any move. When a sparse counterfactual is subsequently
#'       generated by \code{\link{bl_find_sparse_cf}}, any variable marked
#'       \code{"fixed"} always reverts to its observed value regardless of its
#'       Shapley classification.}
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
#' loading matrix is rotated so that the target observation is aligned with the
#' plane of that pair (see \strong{Rotation}), then a fine grid is scored and
#' the boundary contour is extracted.
#'
#' @section Rotation:
#' The rotation adapts \code{Biplot_rotation()} from the research script
#' (\file{scripts/1.3 Optimal Rotation.R} on the \code{original_PhD_code}
#' branch of the repository). An SVD-based orthogonal rotation
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
#' @section Fixed constraints:
#' During the search a \code{"fixed"} feature may move by up to 0.5 in model
#' units, because a vertex of the 2D boundary contour rarely holds a feature at
#' exactly its observed value. After the best pair is chosen, each fixed
#' feature the search moved is set back to its observed value and the point is
#' re-scored. If the restored point still crosses the decision boundary and
#' lies within \code{train_ranges}, it becomes the counterfactual
#' (\code{B_x}, \code{B_z}, \code{B_pred}); otherwise the search's point is
#' kept. Either way \code{fixed_check} records what happened, and
#' \code{print()} explains it. The restored point satisfies the fixed
#' constraint exactly but lies off the 2D boundary, so it need not be the
#' nearest valid point; the pair choice and the reported distances are those
#' of the search.
#'
#' @section Cross-pair selection criterion:
#' The best eigenvector pair is selected by minimum distance from the target
#' to the nearest valid boundary vertex. Two distance metrics are available:
#' \describe{
#'   \item{\code{"mahalanobis"} (default)}{Squared Mahalanobis distance in the
#'     original X-space using the within-class covariance matrix \eqn{W} stored
#'     on \code{bl_result$metric_inv}: \eqn{d_M^2 = v^T W^{-1} v} where
#'     \eqn{v = B_x - x_{obs}}. This is the PhD's cross-pair criterion.}
#'   \item{\code{"euclidean"}}{Plain Euclidean distance in the rotated 2D
#'     Z-space: \eqn{d_z = ||Z_{target} - B_z||_2}.}
#' }
#' Under CVA the two give the same pair: the CVA loadings satisfy
#' \eqn{V^T W V = I / (n - g)}, every rotated basis inherits this, and so
#' \eqn{d_M^2 = (n - g) d_z^2} for every pair. Under PCA they are different
#' metrics (standardised Euclidean versus the binary-class \eqn{W}) and can
#' pick different pairs; the thesis's rule for PCA is the Euclidean one.
#' Within-pair selection (the nearest valid contour vertex) always uses 2D
#' Euclidean distance in the rotated Z-space; only the cross-pair comparison
#' is affected by \code{distance}. See
#' \code{documentation/biplot-rotation/biplot_rotation_and_pair_selection.md}
#' section 3.4 in the package repository.
#'
#' @param bl_result   A \code{"bl_result"} object.
#' @param bl_target   A \code{"bl_target"} object from
#'   \code{\link{bl_select_target}}.
#' @param set_filters A \code{"bl_filters"} object from
#'   \code{\link{set_filters}}, or \code{NULL} (no actionability constraints).
#' @param max_pairs   Integer; maximum number of eigenvector pairs to try.
#'   Default \code{21L}: every pair for up to 7 features
#'   (\code{choose(7, 2) = 21}). With more features, the pairs with the
#'   highest components are left out.
#' @param m           Integer; grid resolution along each axis. Default
#'   \code{200L}.
#' @param distance    Character; cross-pair selection metric. One of
#'   \code{"mahalanobis"} (default) or \code{"euclidean"}. See \strong{Cross-pair
#'   selection criterion} section above.
#' @param verbose     Logical; print progress messages. Default \code{TRUE}.
#'
#' @return A list of class \code{"bl_local_result"} with the following fields:
#'   \describe{
#'     \item{B_z}{1 x 2 matrix; the counterfactual in rotated Z-space: the best
#'       boundary point, or the projection of the restored point when a
#'       \code{"fixed"} restore was applied (see \strong{Fixed constraints}).}
#'     \item{B_x}{Data frame (1 row); the counterfactual in X-space: the
#'       back-projected boundary point, or the restored point.}
#'     \item{B_pred}{Scalar predicted probability at the counterfactual.}
#'     \item{dist_z}{Euclidean distance from target to boundary in the rotated
#'       Z-space of \code{best_pair} only (per-pair values are in
#'       \code{all_distances}). It is the selection criterion only when
#'       \code{distance = "euclidean"}. \code{NA} when no solution is found.}
#'     \item{dist_mahalanobis}{Squared Mahalanobis distance from target to
#'       boundary in X-space, for \code{best_pair} only (per-pair values are in
#'       \code{all_distances_mahalanobis}). The selection criterion when
#'       \code{distance = "mahalanobis"}. \code{NA} when no solution is found.}
#'     \item{distance}{Character; the selection metric used (\code{"mahalanobis"}
#'       or \code{"euclidean"}).}
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
#'     \item{all_distances}{Named numeric vector; Z-distance per pair tried.}
#'     \item{all_distances_mahalanobis}{Named numeric vector; squared
#'       Mahalanobis distance per pair tried.}
#'     \item{solution_found}{Logical.}
#'     \item{blocking_constraint}{Character message if no solution, else NULL.}
#'     \item{bl_target}{The \code{"bl_target"} object (the observed point).}
#'     \item{bl_counterfactual}{A \code{"bl_counterfactual"} object (the boundary
#'       point), or \code{NULL} when no solution was found. Mirrors
#'       \code{bl_target}: fields \code{x_cf} (1-row data frame, model units),
#'       \code{z_cf}, \code{pred_prob} (boundary probability), \code{pred_class}
#'       (from \code{B_pred >= cutoff}), and \code{scaling}. Its print method
#'       shows the counterfactual in original (raw) units when a scaling was set
#'       via \code{\link{bl_set_scaling}}.}
#'     \item{bl_result}{The \code{"bl_result"} object.}
#'     \item{set_filters}{The \code{"bl_filters"} object passed in, or
#'       \code{NULL}. Stored so that \code{\link{bl_find_sparse_cf}} can identify
#'       \code{"fixed"} variables and revert them to their observed values.}
#'     \item{fixed_check}{\code{NULL} when no solution was found or the search
#'       moved no \code{"fixed"} feature. Otherwise a list: \code{vars} (the
#'       fixed features the search moved), \code{applied} (\code{TRUE} when the
#'       restored point is the counterfactual), \code{pred_restored} (model
#'       probability at the restored point), \code{in_range} (whether it lies
#'       within \code{train_ranges}), and \code{search_point} (the search's own
#'       \code{B_x}, \code{B_z}, \code{B_pred}).}
#'   }
#'
#' @importFrom grDevices contourLines colorRampPalette
#' @importFrom stats setNames
#' @export
bl_find_local_cf <- function(bl_result, bl_target,
                                   set_filters = NULL,
                                   max_pairs   = 21L,
                                   m           = 200L,
                                   distance    = c("mahalanobis", "euclidean"),
                                   verbose     = TRUE) {

  if (!inherits(bl_result, "bl_result"))
    stop("'bl_result' must be a 'bl_result' object.", call. = FALSE)
  if (!inherits(bl_target, "bl_target"))
    stop("'bl_target' must be a 'bl_target' object from bl_select_target().",
         call. = FALSE)
  if (!is.null(set_filters) && !inherits(set_filters, "bl_filters"))
    stop("'set_filters' must be a 'bl_filters' object from set_filters() or NULL.",
         call. = FALSE)
  distance <- match.arg(distance)

  V            <- bl_result$V
  X_center     <- bl_result$X_center
  X_sd         <- bl_result$X_sd
  standardise  <- bl_result$standardise
  cutoff       <- bl_result$cutoff
  b_margin     <- bl_result$b_margin
  var_names    <- bl_result$var_names
  p            <- bl_result$num_vars
  train_ranges <- bl_result$train_ranges
  metric_inv   <- bl_result$metric_inv

  # Fall back to Euclidean if metric_inv is missing (older bl_result objects)
  if (distance == "mahalanobis" && is.null(metric_inv)) {
    warning(
      "bl_result$metric_inv not found; falling back to distance = 'euclidean'. ",
      "Rebuild bl_result with the current package version to enable Mahalanobis.",
      call. = FALSE
    )
    distance <- "euclidean"
  }

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

  pair_names                <- apply(pairs, 1L, function(r) paste0("(", r[1], ",", r[2], ")"))
  all_distances             <- setNames(rep(NA_real_, nrow(pairs)), pair_names)
  all_distances_mahalanobis <- setNames(rep(NA_real_, nrow(pairs)), pair_names)
  best_selector             <- Inf
  best_result               <- NULL

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

      # Actionability filter first (Phase 3 only -- no hull polygon).
      # Typically more selective than train_ranges (especially with "fixed"
      # constraints), so running it first lets us short-circuit before the
      # train_ranges pass on the vertices it eliminates.
      keep <- if (!is.null(set_filters)) {
        .apply_actionability(Bx, x_obs_named, set_filters)
      } else {
        rep(TRUE, nrow(Bx))
      }
      if (!any(keep)) next

      # train_ranges feasibility filter -- only on survivors.
      if (is.list(train_ranges)) {
        in_range <- get_filter_logical_vector(Bx[keep, , drop = FALSE], train_ranges)
        keep[keep] <- in_range
      }

      if (!any(keep)) next

      # Model consistency check (same as global bl_find_boundary)
      Mi_keep <- Mi[keep, , drop = FALSE]
      Bx_keep <- Bx[keep, , drop = FALSE]
      p_bnd   <- .pred_function(
        model_use  = bl_result$model,
        model_type = bl_result$model_type,
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

    # Back-project counterfactual to X-space (needed for every pair when
    # distance = "mahalanobis", since the selector uses X-space coordinates).
    B_x_mat <- B_z_local %*% tVr_rot
    if (isTRUE(standardise)) B_x_mat <- sweep(B_x_mat, 2L, X_sd, "*")
    B_x_mat <- sweep(B_x_mat, 2L, X_center, "+")
    B_x     <- as.data.frame(B_x_mat)
    colnames(B_x) <- var_names

    # Squared Mahalanobis distance in X-space (using stored within-class
    # metric W^-1). NA when metric_inv is unavailable (fallback handled above).
    dist_mah <- if (!is.null(metric_inv)) {
      v <- as.numeric(B_x[1L, var_names]) - x_obs
      as.numeric(t(v) %*% metric_inv %*% v)
    } else {
      NA_real_
    }

    all_distances[pair_names[k]]             <- dist_z
    all_distances_mahalanobis[pair_names[k]] <- dist_mah

    if (isTRUE(verbose)) {
      if (distance == "mahalanobis") {
        message(sprintf("    Distance: %.4f Z  |  %.4f Mahalanobis (sqrt)",
                        dist_z, sqrt(dist_mah)))
      } else {
        message(sprintf("    Distance: %.4f Z", dist_z))
      }
    }

    selector <- if (distance == "mahalanobis") dist_mah else dist_z
    if (is.finite(selector) && selector < best_selector) {
      best_selector <- selector

      # Score the counterfactual through the model (only for the winning pair
      # at any point in the loop -- cheap, but no need to do it for losers).
      B_pred <- .pred_function(
        model_use  = bl_result$model,
        model_type = bl_result$model_type,
        new_data   = B_x
      )

      col_vec   <- grDevices::colorRampPalette(
        c("deepskyblue", "white", "lightsalmon")
      )(101L)
      col_value <- col_vec[floor(grid_prob * 100) + 1L]

      best_result <- list(
        B_z              = B_z_local,
        B_x              = B_x,
        B_pred           = B_pred[[1L]],
        dist_z           = dist_z,
        dist_mahalanobis = dist_mah,
        Z_target         = Z_target,
        best_pair        = pair,
        Vr_rot           = Vr_rot,
        tVr_rot          = tVr_rot,
        Z_train_rot      = Z_train_rot,
        Zgrid            = Zgrid,
        grid_prob        = grid_prob,
        col_value        = col_value,
        ct_local         = ct_local,
        xseq             = xseq,
        yseq             = yseq,
        min_val          = min_val,
        max_val          = max_val
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
    message("bl_find_local_cf: ", blocking)
  } else {
    blocking <- NULL
    if (isTRUE(verbose)) {
      if (distance == "mahalanobis") {
        message(sprintf(
          "Best pair: (%d, %d) | Selector: Mahalanobis (sqrt) = %.4f | Z = %.4f",
          best_result$best_pair[1L], best_result$best_pair[2L],
          sqrt(best_result$dist_mahalanobis), best_result$dist_z
        ))
      } else {
        message(sprintf("Best pair: (%d, %d) | Distance (Z): %.4f",
                        best_result$best_pair[1L],
                        best_result$best_pair[2L],
                        best_result$dist_z))
      }
    }
  }

  # Restore "fixed" features the search moved, when the restored point still
  # crosses the boundary (post-selection; the pair choice is unchanged).
  fixed_check <- NULL
  if (solution_found) {
    rf <- .restore_fixed(best_result, x_obs, var_names, set_filters,
                         bl_result, bl_target, sv)
    if (!is.null(rf)) {
      best_result <- rf$best_result
      fixed_check <- rf$fixed_check
      if (isTRUE(verbose)) {
        vars_txt <- paste(fixed_check$vars, collapse = ", ")
        if (fixed_check$applied) {
          message(sprintf(paste(
            "Fixed feature(s) %s restored to observed value(s); the restored",
            "point still crosses the boundary (pred %.3f) and is used as the",
            "counterfactual."), vars_txt, fixed_check$pred_restored))
        } else {
          message(sprintf(paste(
            "Fixed feature(s) %s: the restored point (pred %.3f) %s, so the",
            "search's point is kept."), vars_txt, fixed_check$pred_restored,
            if (fixed_check$in_range) "does not cross the boundary"
            else "lies outside the training ranges"))
        }
      }
    }
  }

  null_pair_fields <- list(
    B_z = NULL, B_x = NULL, B_pred = NA_real_,
    dist_z = NA_real_, dist_mahalanobis = NA_real_,
    Z_target = NULL, best_pair = NA_integer_, Vr_rot = NULL,
    tVr_rot = NULL, Z_train_rot = NULL, Zgrid = NULL,
    grid_prob = NULL, col_value = NULL, ct_local = NULL,
    xseq = NULL, yseq = NULL, min_val = NULL, max_val = NULL
  )

  structure(
    c(
      if (solution_found) best_result else null_pair_fields,
      list(
        all_distances             = all_distances,
        all_distances_mahalanobis = all_distances_mahalanobis,
        distance                  = distance,
        solution_found            = solution_found,
        blocking_constraint       = blocking,
        bl_target                 = bl_target,
        bl_counterfactual         = if (solution_found)
          .make_bl_counterfactual(best_result, bl_result) else NULL,
        bl_result                 = bl_result,
        set_filters               = set_filters,
        fixed_check               = fixed_check
      )
    ),
    class = "bl_local_result"
  )
}


# --------------------------------------------------------------------------
# Private: package the boundary counterfactual as a printable object
# --------------------------------------------------------------------------

#' Build a bl_counterfactual object from the winning local search result
#'
#' Mirrors the `bl_target` object: stores the counterfactual feature values in
#' model units plus the display-only `scaling`, so `print.bl_counterfactual()`
#' can show them in original (raw) units. No new computation -- just packages the
#' already-found boundary point (`B_x`/`B_z`/`B_pred`).
#'
#' @param best_result The winning per-pair result list (carries `B_x`, `B_z`,
#'   `B_pred`).
#' @param bl_result   The parent `bl_result` (for `cutoff` and `scaling`).
#' @return A `"bl_counterfactual"` object.
#' @noRd
.make_bl_counterfactual <- function(best_result, bl_result) {
  structure(
    list(
      x_cf       = best_result$B_x,
      z_cf       = best_result$B_z,
      pred_prob  = best_result$B_pred,
      pred_class = as.integer(best_result$B_pred >= bl_result$cutoff),
      scaling    = bl_result$scaling
    ),
    class = "bl_counterfactual"
  )
}


#' Restore "fixed" features after the local search
#'
#' The search lets a `"fixed"` feature move by up to `.fixed_tol` model units.
#' After the best pair is chosen, this sets each fixed feature the search moved
#' back to its observed value and re-scores the point. If the restored point
#' still crosses the decision boundary and lies within `train_ranges`, it
#' replaces the search's point as the counterfactual. Only the fixed features
#' change, so the other actionability constraints still hold.
#'
#' @param best_result The winning per-pair result list (`B_x`, `B_z`, `B_pred`,
#'   `Vr_rot`, ...).
#' @param x_obs       Numeric vector of observed values (model units), in
#'   `var_names` order.
#' @param var_names   Character vector of feature names.
#' @param set_filters A `"bl_filters"` object or `NULL`.
#' @param bl_result   The parent `bl_result` (model, cutoff, `X_center`,
#'   `train_ranges`).
#' @param bl_target   The `bl_target` (for its predicted class).
#' @param sv          Per-feature scale used by the search (`X_sd` when
#'   standardised, otherwise 1s).
#' @return `NULL` when the search moved no fixed feature; otherwise a list with
#'   `best_result` (updated when the restore is applied) and `fixed_check`.
#' @noRd
.restore_fixed <- function(best_result, x_obs, var_names, set_filters,
                           bl_result, bl_target, sv) {
  cons <- if (!is.null(set_filters)) set_filters$constraints else NULL
  if (length(cons) == 0L) return(NULL)
  is_fixed <- vapply(cons, function(v) is.character(v) && identical(v, "fixed"),
                     logical(1L))
  fixed <- names(cons)[is_fixed]
  if (length(fixed) == 0L) return(NULL)

  x_obs_named <- setNames(x_obs, var_names)
  B_x_search  <- best_result$B_x
  moved <- fixed[abs(as.numeric(unlist(B_x_search[1L, fixed])) -
                       x_obs_named[fixed]) > 1e-8]
  if (length(moved) == 0L) return(NULL)

  B_x_r <- B_x_search
  B_x_r[1L, moved] <- x_obs_named[moved]
  pred_r <- .pred_function(
    model_use  = bl_result$model,
    model_type = bl_result$model_type,
    new_data   = B_x_r
  )[[1L]]

  crosses  <- as.integer(pred_r >= bl_result$cutoff) != bl_target$pred_class
  in_range <- if (is.list(bl_result$train_ranges)) {
    get_filter_logical_vector(B_x_r, bl_result$train_ranges)[[1L]]
  } else {
    TRUE
  }
  applied <- isTRUE(crosses) && isTRUE(in_range)

  fixed_check <- list(
    vars          = moved,
    applied       = applied,
    pred_restored = pred_r,
    in_range      = in_range,
    search_point  = list(B_x    = B_x_search,
                         B_z    = best_result$B_z,
                         B_pred = best_result$B_pred)
  )

  if (applied) {
    # Same projection as bl_select_target(); a point on the plane re-projects
    # to itself, the restored point lands at its biplot projection.
    x_r_st <- (as.numeric(unlist(B_x_r[1L, var_names])) - bl_result$X_center) / sv
    B_z_r  <- matrix(x_r_st %*% best_result$Vr_rot, nrow = 1L)
    colnames(B_z_r) <- colnames(best_result$B_z)
    best_result$B_x    <- B_x_r
    best_result$B_z    <- B_z_r
    best_result$B_pred <- pred_r
  }

  list(best_result = best_result, fixed_check = fixed_check)
}


#' Constraint label per feature, for the comparison tables
#'
#' @param constraints Named list of constraints (`bl_filters$constraints`), or
#'   `NULL`.
#' @param var_names   Character vector of feature names.
#' @return Character vector, one label per feature: the constraint as given, a
#'   range as `"[lo, hi]"`, or `"none"`.
#' @noRd
.constraint_labels <- function(constraints, var_names) {
  out <- setNames(rep("none", length(var_names)), var_names)
  for (nm in intersect(names(constraints), var_names)) {
    v <- constraints[[nm]]
    out[[nm]] <- if (is.character(v)) v else sprintf("[%.4g, %.4g]", v[1L], v[2L])
  }
  unname(out)
}


#' Observed vs counterfactual comparison table
#'
#' Shared by `print.bl_local_result()` and `print.bl_sparse_result()` so the
#' two tables have the same layout. Levels convert as `value * scale + center`
#' and the change (a delta) as `value * scale`, at display time only.
#'
#' @param x_obs       Numeric vector of observed values (model units), in
#'   `var_names` order.
#' @param cols        Named list of numeric vectors (model units, `var_names`
#'   order), e.g. `list(counterfactual = B_x)`. `change` is measured from the
#'   last one.
#' @param var_names   Character vector of feature names (the row order).
#' @param scaling     `bl_result$scaling`, or `NULL`.
#' @param constraints Named list of constraints, or `NULL`.
#' @return A list with `table` (data frame, rows = features, unrounded),
#'   `units` (`"raw"` or `"model"`) and `digits` (decimals to print per
#'   feature row; see `.fmt_digits()`).
#' @noRd
.cf_compare_table <- function(x_obs, cols, var_names, scaling, constraints) {
  # The first conversion warns once if the scaling is incomplete; the other
  # columns then stay in model units without repeating the warning.
  observed <- .scale_to_raw(x_obs, scaling, var_names, "level")
  raw      <- !is.null(scaling) && !anyNA(scaling$center[var_names]) &&
    !anyNA(scaling$scale[var_names])
  sc       <- if (raw) scaling else NULL

  tab <- data.frame(observed = unname(observed), row.names = var_names)
  for (nm in names(cols))
    tab[[nm]] <- .scale_to_raw(unname(cols[[nm]]), sc, var_names, "level")
  last       <- unname(cols[[length(cols)]])
  tab$change <- .scale_to_raw(last - unname(x_obs), sc, var_names, "delta")
  tab$constraint <- .constraint_labels(constraints, var_names)

  # Decimals per feature row: 2, widened where a difference from observed
  # would otherwise print as 0.00 (e.g. a change of -0.0016).
  diffs <- as.matrix(tab[c(names(cols), "change")])
  diffs[, names(cols)] <- diffs[, names(cols)] - tab$observed
  d      <- .fmt_digits(diffs)
  dim(d) <- dim(diffs)

  list(table = tab, units = if (raw) "raw" else "model",
       digits = as.integer(apply(d, 1L, max)))
}


#' Decimals needed to print a value at 2 dp without hiding it
#'
#' @param x          Numeric vector (or matrix).
#' @param digits     Default number of decimals.
#' @param max_digits Most decimals used.
#' @return Integer vector, one per element of `x`: `digits`, or more when `x`
#'   is not zero but rounds to zero at `digits` decimals (2 significant
#'   digits, at most `max_digits`). A value that rounds to zero even at
#'   `max_digits` keeps `digits`.
#' @noRd
.fmt_digits <- function(x, digits = 2L, max_digits = 4L) {
  x      <- as.numeric(x)
  d      <- rep(as.integer(digits), length(x))
  hidden <- !is.na(x) & round(x, digits) == 0 & round(x, max_digits) != 0
  d[hidden] <- pmin(as.integer(max_digits),
                    as.integer(1 - floor(log10(abs(x[hidden])))))
  d
}


#' Format numbers in fixed notation
#'
#' @param x Numeric vector.
#' @param d Integer decimals, one per element (recycled).
#' @return Character vector: never scientific notation, never `"-0.00"`.
#' @noRd
.fmt_fixed <- function(x, d) {
  x <- as.numeric(x)
  d <- rep_len(as.integer(d), length(x))
  x <- round(x, d)
  x[!is.na(x) & x == 0] <- 0
  sprintf("%.*f", d, x)
}


#' Format a comparison table for printing
#'
#' @param cmp The list returned by `.cf_compare_table()`.
#' @return `cmp$table` with every numeric column formatted by `.fmt_fixed()`,
#'   using the row decimals `cmp$digits`.
#' @noRd
.format_cf_table <- function(cmp) {
  tab <- cmp$table
  for (nm in setdiff(names(tab), "constraint"))
    tab[[nm]] <- .fmt_fixed(tab[[nm]], cmp$digits)
  tab
}


#' Note explaining a "fixed" feature the local search moved
#'
#' @param fixed_check The `fixed_check` field of a `bl_local_result`, or `NULL`.
#' @param x_obs       Numeric vector of observed values (model units), in
#'   `var_names` order.
#' @param var_names   Character vector of feature names.
#' @param scaling     `bl_result$scaling`, or `NULL`.
#' @return Character vector of wrapped note lines; empty when there is nothing
#'   to explain.
#' @noRd
.fixed_notes <- function(fixed_check, x_obs, var_names, scaling) {
  if (is.null(fixed_check) || length(fixed_check$vars) == 0L) return(character(0))
  vars <- fixed_check$vars
  raw  <- !is.null(scaling) && !anyNA(scaling$center[vars]) &&
    !anyNA(scaling$scale[vars])
  sc   <- if (raw) scaling else NULL
  fmt  <- function(v, flag = "") {
    formatC(round(v, 4L), format = "f", digits = 4L, flag = flag,
            drop0trailing = TRUE)
  }

  x_obs_named <- setNames(x_obs, var_names)
  search_x    <- as.numeric(unlist(fixed_check$search_point$B_x[1L, vars]))
  search_raw  <- .scale_to_raw(search_x, sc, vars, "level")
  move_raw    <- .scale_to_raw(search_x - x_obs_named[vars], sc, vars, "delta")
  tol_raw     <- .scale_to_raw(rep(.fixed_tol, length(vars)), sc, vars, "delta")

  moves <- paste(sprintf('%s is "fixed", but the boundary search moved it to %s (%s).',
                         vars, fmt(search_raw), fmt(move_raw, "+")),
                 collapse = " ")
  tol_txt <- if (length(vars) == 1L) {
    sprintf("+/- %s here", fmt(tol_raw))
  } else {
    paste(sprintf("+/- %s for %s", fmt(tol_raw), vars), collapse = ", ")
  }
  why <- sprintf(paste(
    'The search lets a "fixed" feature move by up to %s model units (%s),',
    "because a point on the 2D decision boundary rarely holds a feature at",
    "exactly its observed value."), fmt(.fixed_tol), tol_txt)

  vars_txt <- paste(vars, collapse = ", ")
  value_s  <- if (length(vars) == 1L) "value" else "values"
  outcome <- if (isTRUE(fixed_check$applied)) {
    sprintf(paste(
      "Restoring %s to its observed %s still crosses the boundary (prediction",
      "%.4f), so the restored point is the selected counterfactual."),
      vars_txt, value_s, fixed_check$pred_restored)
  } else if (!isTRUE(fixed_check$in_range)) {
    sprintf(paste(
      "Restoring %s to its observed %s would leave the training feature",
      "ranges, so the search's point is kept."), vars_txt, value_s)
  } else {
    sprintf(paste(
      "Restoring %s to its observed %s gives prediction %.4f, which does not",
      "cross the boundary, so the search's point is kept."),
      vars_txt, value_s, fixed_check$pred_restored)
  }

  strwrap(paste(moves, why, outcome), width = 78L,
          initial = "  Note: ", prefix = "        ")
}


# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Plot the local biplot for a bl_local_result
#'
#' Renders the rotated biplot for the best eigenvector pair found by
#' \code{\link{bl_find_local_cf}} using the same biplotEZ pipeline as
#' \code{\link{plot_biplotEZ}}. The biplotEZ object is patched with the rotated
#' loading matrix (\code{Vr_rot}) and rotated training coordinates
#' (\code{Z_train_rot}) so that axes, grid, points, and contours are all drawn
#' in the locally-rotated Z-space. Training points are coloured by confusion
#' type (TP = red, TN = blue, FP = purple, FN = orange). The target observation
#' is shown as a filled circle coloured by the confusion scheme when its true
#' class is known (\code{bl_target$class_obs}, resolved by
#' \code{\link{bl_select_target}} from the target's own row), or red/blue by
#' predicted class when it is not. The nearest boundary point is shown as a cross
#' with an arrow.
#'
#' @param x               A \code{"bl_local_result"} object.
#' @param no_grid         Logical; if \code{TRUE} the prediction grid is hidden.
#'   Default \code{FALSE}.
#' @param plot_points     Logical; if \code{FALSE} (default) training points are
#'   hidden (show only the target and counterfactual). Set to \code{TRUE} to
#'   also show training data. To overlay test or holdout data instead, use
#'   \code{\link{plot_biplotEZ}} with its \code{points} argument.
#' @param no_contour      Logical; if \code{TRUE} boundary contour lines are
#'   hidden. Default \code{FALSE}.
#' @param cex_z           Numeric; size of training data points. Default
#'   \code{0.5}.
#' @param label_dir       Character; biplotEZ axis label direction.
#'   \code{"Paral"} (default) follows the plot border -- left/right borders are
#'   vertical, top/bottom are horizontal. \code{"Hor"} forces horizontal.
#'   \code{"Orthog"} draws labels perpendicular to the axis.
#' @param label_cex       Numeric; variable name (axis label) size. Default
#'   \code{1}.
#' @param tick_label_cex  Numeric; axis tick label size. Default \code{0.6}.
#' @param ticks_v         Integer; number of ticks per variable axis. Default
#'   \code{1L}.
#' @param ticks_var       Integer, integer vector, character, or character
#'   vector; variable(s) whose tick-mark count should override \code{ticks_v}.
#'   Supply variable names (e.g. \code{"Sepal.Length"}) or integer indices
#'   (e.g. \code{c(1L, 3L)}). \code{0} (default) means no override -- every
#'   axis uses \code{ticks_v}.
#' @param ticks_n         Integer or integer vector; tick-mark count(s) for
#'   the variable(s) named in \code{ticks_var}. A single value applies to all
#'   listed variables; a vector must match the length of \code{ticks_var} and
#'   pairs with it position by position. A name not in the model is dropped,
#'   with a warning, together with its value. Default \code{5L}.
#' @param which           Integer vector; indices of variables to draw axes for.
#'   Defaults to all variables.
#' @param X_names         Character vector; custom variable names for axis
#'   labels. Defaults to \code{bl_result$var_names}.
#' @param label_offset_var Integer, integer vector, character, or character
#'   vector; variable(s) whose axis labels should be shifted outward. Supply
#'   variable names or integer indices. \code{0} (default) means no offset.
#' @param label_offset_dist Numeric or numeric vector; outward offset distance(s)
#'   in margin lines. A single value applies to all listed variables; a vector
#'   must match the length of \code{label_offset_var} and pairs with it
#'   position by position. A name not in the model is dropped, with a warning,
#'   together with its value. Useful range 1--3. Default \code{1.5}.
#' @param show_arrows    Logical; if \code{TRUE} (default) an arrow is drawn
#'   from the target to its counterfactual.
#' @param arrow_col       Character; colour for the CF cross and arrow. Default
#'   \code{"grey30"}.
#' @param contour_col     Character; colour for contour lines. Default
#'   \code{"black"}.
#' @param contour_lwd     Numeric; line width for contour lines. Default
#'   \code{1.5}.
#' @param contour_lty     Integer; line type for contour lines. Default
#'   \code{1L}.
#' @param print_summary   Logical; print a local CF summary to the console.
#'   Default \code{TRUE}. Set to \code{FALSE} when called from
#'   \code{\link{plot.bl_sparse_result}}.
#' @param new_title       Character; overrides the auto-generated local biplot
#'   title (\code{"Local biplot -- target N [pair ...]"}). \code{NA} (default)
#'   keeps the auto-generated title.
#' @param ...             Additional arguments (currently ignored). There is no
#'   \code{rotate_deg}: the local biplot's orientation is fixed by the optimal
#'   rotation (\code{.bl_rotate()}), so passing \code{rotate_deg} is an error.
#'   \code{rotate_deg} is only for viewing the global biplot
#'   (\code{\link{plot_biplotEZ}}, \code{\link{plot.bl_surrogate}}).
#'
#' @importFrom graphics points lines arrows text
#' @importFrom grDevices dev.cur dev.new
#' @importFrom dplyr case_when
#' @export
plot.bl_local_result <- function(x,
                                 no_grid           = FALSE,
                                 plot_points       = FALSE,
                                 no_contour        = FALSE,
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
                                 show_arrows       = TRUE,
                                 arrow_col         = "grey30",
                                 contour_col       = "black",
                                 contour_lwd       = 1.5,
                                 contour_lty       = 1L,
                                 print_summary     = TRUE,
                                 new_title         = NA,
                                 ...) {
  # rotate_deg is for viewing the global biplot only. Without this check it
  # would fall into ... and be silently ignored.
  if ("rotate_deg" %in% ...names())
    stop("'rotate_deg' is not available for the local biplot: its orientation ",
         "is fixed by the optimal rotation (.bl_rotate()). Use rotate_deg only ",
         "to view the global biplot (plot_biplotEZ(), plot(<bl_surrogate>)).",
         call. = FALSE)

  if (!x$solution_found) {
    message("No solution found -- nothing to plot.")
    return(invisible(x))
  }

  label_dir <- match.arg(label_dir, c("Paral", "Hor", "Orthog"))

  if (grDevices::dev.cur() == 1L) grDevices::dev.new()

  bl_result   <- x$bl_result
  bl_target   <- x$bl_target
  Vr_rot      <- x$Vr_rot
  tVr_rot     <- x$tVr_rot
  Z_train_rot <- x$Z_train_rot
  Z_target    <- x$Z_target
  B_z         <- x$B_z
  best_pair   <- x$best_pair

  var_names  <- bl_result$var_names
  num_vars   <- bl_result$num_vars
  proj_dims  <- bl_result$proj_dims
  cutoff     <- bl_result$cutoff

  # Defaults for axis arguments (mirrors plot_biplotEZ)
  if (is.null(which))   which   <- seq_len(num_vars)
  if (is.null(X_names)) X_names <- var_names

  label_line_vec <- .make_label_line_vec(label_offset_var, label_offset_dist,
                                         num_vars, var_names)
  ticks_vec <- .make_ticks_vec(ticks_v, ticks_var, ticks_n, num_vars, var_names)

  # Confusion colours for training points
  train_df         <- bl_result$train_data
  train_pred       <- .pred_function(
    model_use  = bl_result$model,
    model_type = bl_result$model_type,
    new_data   = train_df[, var_names, drop = FALSE]
  )
  train_pred_class <- as.integer(train_pred >= cutoff)

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

  row_id    <- bl_target$row_id
  row_label <- if (is.na(row_id)) "external" else as.character(row_id)
  pred_prob  <- bl_target$pred_prob
  pred_class <- bl_target$pred_class

  # ---- Target point colour ---------------------------------------------
  # Known true class (bl_target$class_obs, resolved by bl_select_target()):
  # confusion colour. Unknown true class (no class column): red/blue by
  # prediction.
  # A bl_target built before class_obs existed has no such field.
  target_true_class <- bl_target$class_obs
  if (is.null(target_true_class)) target_true_class <- NA_integer_
  target_bg <- if (is.na(target_true_class)) {
    if (pred_class == 1L) "red" else "blue"
  } else {
    dplyr::case_when(
      target_true_class == 1L & pred_class == 1L ~ "red",
      target_true_class == 0L & pred_class == 0L ~ "blue",
      target_true_class == 0L & pred_class == 1L ~ "purple",
      target_true_class == 1L & pred_class == 0L ~ "orange",
      TRUE ~ "grey50"
    )
  }

  # ---- Patch biplotEZ object for rotated local space -------------------
  # Replace the proj_dims columns of Z and Lmat with the rotated local
  # counterparts, then recompute ax.one.unit. Same approach as the
  # rotate_deg path in plot_biplotEZ().
  biplot_plot <- bl_result$biplot_obj
  biplot_plot$Z[, proj_dims]    <- Z_train_rot
  biplot_plot$Lmat[, proj_dims] <- Vr_rot
  biplot_plot$ax.one.unit       <-
    (1 / diag(t(tVr_rot) %*% tVr_rot)) * t(tVr_rot)
  biplot_plot$Title <- sprintf(
    "Local biplot \u2013 target %s  [pair (%d,%d) | p=%.3f, class %d]",
    row_label, best_pair[1L], best_pair[2L], pred_prob, pred_class
  )
  if (!is.na(new_title)) biplot_plot$Title <- new_title

  x_grid   <- x$Zgrid
  ct_local <- x$ct_local

  # ---- Optional raw-unit axis relabel (display only) -------------------
  # means/sd are unaffected by the local rotation above, so the affine
  # relabel is applied to the rotated object directly. No-op when scaling is
  # NULL.
  biplot_plot <- .bl_rescale_biplot_axes(biplot_plot, bl_result$scaling, var_names)

  # ---- Step 1: biplotEZ base plot (axes canvas) ------------------------
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

  # ---- Step 2: prediction grid -----------------------------------------
  if (!isTRUE(no_grid)) {
    graphics::points(x_grid[, 1L], x_grid[, 2L],
                     col = x$col_value, pch = 15L, cex = 0.5)
  }

  # ---- Step 3: training points (confusion colours) ---------------------
  if (isTRUE(plot_points)) {
    graphics::points(Z_train_rot[, 1L], Z_train_rot[, 2L],
                     col = train_col, pch = 16L, cex = cex_z)
  }

  # ---- Step 4: axes redrawn on top -------------------------------------
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
  if (!isTRUE(no_contour)) {
    for (cl in ct_local) {
      graphics::lines(cl$x, cl$y,
                      col = contour_col, lwd = contour_lwd, lty = contour_lty)
    }
  }

  # ---- Step 6: target observation (coloured circle) --------------------
  graphics::points(Z_target[, 1L], Z_target[, 2L],
                   pch = 21L, bg = target_bg, col = "black", cex = 1.8)
  if (!is.na(row_id)) {
    graphics::text(Z_target[, 1L], Z_target[, 2L],
                   labels = row_label, pos = 3L, cex = 0.7)
  }

  # ---- Step 7: CF cross + arrow ----------------------------------------
  graphics::points(B_z[, 1L], B_z[, 2L],
                   pch = 4L, col = arrow_col, cex = 1.2, lwd = 1.5)
  if (isTRUE(show_arrows)) {
    graphics::arrows(
      Z_target[1L, 1L], Z_target[1L, 2L],
      B_z[1L, 1L],      B_z[1L, 2L],
      length = 0.08, col = arrow_col, lwd = 1.2
    )
  }

  # ---- Step 8: local CF summary ----------------------------------------
  if (isTRUE(print_summary)) {
    cat("\n--- Local CF summary ---\n")
    cat(sprintf("  Target         : %s\n",         row_label))
    cat(sprintf("  Pred prob      : %.4f  (class %d)\n", pred_prob, pred_class))
    cat(sprintf("  Best pair      : (%d, %d)\n",   best_pair[1L], best_pair[2L]))
    cat(sprintf("  Selector       : %s\n",
                if (!is.null(x$distance)) x$distance else "euclidean"))
    cat(sprintf("  Distance (Z)   : %.4f\n",       x$dist_z))
    if (!is.na(x$dist_mahalanobis))
      cat(sprintf("  Mahalanobis    : %.4f (sqrt)\n", sqrt(x$dist_mahalanobis)))
    cat(sprintf("  Boundary pred  : %.4f\n",       x$B_pred))
    cat("------------------------\n")
  }

  invisible(x)
}


#' Print method for bl_local_result
#'
#' Prints, in order: the search summary (best pair, selector, distances,
#' prediction at the counterfactual); the Target and Counterfactual blocks; the
#' per-pair distances; and, when a solution was found, a table comparing the
#' observed and counterfactual value of each feature with its change and its
#' \code{\link{set_filters}} constraint. Values are in original (raw) units
#' when a scaling was set via \code{\link{bl_set_scaling}}, otherwise in model
#' units. When the search moved a \code{"fixed"} feature, a note under the
#' table explains why and whether the restored point was used (see
#' \strong{Fixed constraints} in \code{\link{bl_find_local_cf}}).
#'
#' Numbers print in fixed notation with 2 decimals; a feature whose change
#' would print as 0.00 gets more decimals (up to 4) in its row, and the Target
#' and Counterfactual blocks use the same decimals as the table.
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
    cat(sprintf("  Selector       : %s\n",
                if (!is.null(x$distance)) x$distance else "euclidean"))
    cat(sprintf("  Distance (Z)   : %.4f\n", x$dist_z))
    if (!is.na(x$dist_mahalanobis))
      cat(sprintf("  Mahalanobis    : %.4f (sqrt of squared d_M^2 = %.4f)\n",
                  sqrt(x$dist_mahalanobis), x$dist_mahalanobis))
    cat(sprintf("  Boundary pred  : %.4f\n", x$B_pred))
    if (isTRUE(x$fixed_check$applied))
      cat(sprintf("  Fixed restored : %s\n",
                  paste(x$fixed_check$vars, collapse = ", ")))
  } else {
    cat(sprintf("  Reason         : %s\n", x$blocking_constraint))
  }
  # The comparison table sets the decimals that the Target and
  # Counterfactual blocks then share.
  cmp    <- NULL
  digits <- NULL
  if (isTRUE(x$solution_found) && !is.null(x$B_x)) {
    var_names <- x$bl_result$var_names
    scaling   <- x$bl_result$scaling
    x_obs     <- as.numeric(unlist(x$bl_target$x_obs[1L, var_names]))
    cf        <- as.numeric(unlist(x$B_x[1L, var_names]))
    cmp <- .cf_compare_table(x_obs, list(counterfactual = cf), var_names,
                             scaling, x$set_filters$constraints)
    digits <- setNames(cmp$digits, var_names)
  }

  cat("\n  Target:\n")
  print(x$bl_target, digits = digits)
  if (!is.null(x$bl_counterfactual)) {
    cat("\n  Counterfactual:\n")
    print(x$bl_counterfactual, digits = digits)
  }
  cat("\n  Per-pair distances:\n")
  df <- data.frame(
    pair             = names(x$all_distances),
    dist_z           = x$all_distances,
    dist_mahalanobis = x$all_distances_mahalanobis,
    row.names        = NULL,
    stringsAsFactors = FALSE
  )
  print(df, row.names = FALSE)

  if (!is.null(cmp)) {
    cat(sprintf("\n  Observed vs counterfactual (%s units):\n", cmp$units))
    print(.format_cf_table(cmp))
    notes <- .fixed_notes(x$fixed_check, x_obs, var_names, scaling)
    if (length(notes) > 0L) cat("\n", paste(notes, collapse = "\n"), "\n", sep = "")
  }
  invisible(x)
}


#' Print method for bl_target
#'
#' @param x A \code{"bl_target"} object.
#' @param digits Optional named integer vector of decimals per feature; by
#'   default 2, or more where a value would print as 0.00.
#' @param ... Currently ignored.
#' @export
print.bl_target <- function(x, digits = NULL, ...) {
  cat("-- bl_target --\n")
  rid <- if (is.na(x$row_id)) "external" else as.character(x$row_id)
  cat(sprintf("  Row ID         : %s\n", rid))
  cat(sprintf("  Predicted prob : %.4f\n", x$pred_prob))
  cat(sprintf("  Predicted class: %d\n",  x$pred_class))
  unit_note <- if (!is.null(x$scaling)) " (raw units)" else ""
  cat(sprintf("  Feature values%s :\n", unit_note))
  vals  <- as.numeric(x$x_obs)
  nms   <- names(x$x_obs)
  vals  <- .scale_to_raw(vals, x$scaling, nms, "level")
  .print_values(vals, nms, digits)
  invisible(x)
}


#' Print feature values in fixed notation
#'
#' @param vals   Numeric vector of feature values.
#' @param nms    Feature names.
#' @param digits Optional named integer vector of decimals per feature; a
#'   feature it does not name gets `.fmt_digits()`.
#' @return `NULL`, invisibly; prints a named character vector.
#' @noRd
.print_values <- function(vals, nms, digits = NULL) {
  d <- .fmt_digits(vals)
  if (!is.null(digits)) {
    given <- unname(digits[nms])
    d[!is.na(given)] <- as.integer(given[!is.na(given)])
  }
  out <- .fmt_fixed(vals, d)
  names(out) <- nms
  print(out, quote = FALSE, right = TRUE)
  invisible(NULL)
}


#' Print method for bl_counterfactual
#'
#' Companion to [print.bl_target()]: shows the boundary counterfactual found by
#' [bl_find_local_cf()] in original (raw) feature units when a scaling was
#' recorded via [bl_set_scaling()], otherwise in model units.
#'
#' @param x A \code{"bl_counterfactual"} object (the \code{bl_counterfactual}
#'   field of a \code{"bl_local_result"}).
#' @param digits Optional named integer vector of decimals per feature; by
#'   default 2, or more where a value would print as 0.00.
#' @param ... Currently ignored.
#' @export
print.bl_counterfactual <- function(x, digits = NULL, ...) {
  cat("-- bl_counterfactual --\n")
  cat(sprintf("  Boundary prob  : %.4f\n", x$pred_prob))
  cat(sprintf("  Predicted class: %d\n",   x$pred_class))
  unit_note <- if (!is.null(x$scaling)) " (raw units)" else ""
  cat(sprintf("  Feature values%s :\n", unit_note))
  vals <- as.numeric(x$x_cf)
  nms  <- names(x$x_cf)
  vals <- .scale_to_raw(vals, x$scaling, nms, "level")
  .print_values(vals, nms, digits)
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
