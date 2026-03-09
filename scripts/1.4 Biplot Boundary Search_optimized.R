############################################################
# Script: Biplot Boundary Search (Optimized & Documented)
# Author: Adriaan Rowan
# Date: 2025-09-13
# Description:
#   Memory-efficient and faster boundary search for biplot-based models.
#   
#   Inputs are expected to match those produced by `biplot_plane()` and
#   the rest of your pipeline in “3 Biplot Search – STB”.
#   Uses blockwise nearest-neighbour search to avoid n×m allocations.
############################################################

# --------------------------- helpers --------------------------------------

# Fast nearest-neighbour indices without building a full n×m distance matrix.
# Z: n×2 matrix of points; B: m×2 boundary polyline.
# Returns integer vector of length n giving the row of B that is nearest.
.nearest_idx_block <- function(Z, B, block = 5000L) {
  Z <- as.matrix(Z); B <- as.matrix(B)
  if (ncol(Z) != 2L || ncol(B) != 2L) stop(".nearest_idx_block: Z and B must have 2 columns (x,y).")
  n <- nrow(Z); m <- nrow(B)
  if (m == 0L || n == 0L) return(rep(NA_integer_, n))
  best_idx <- rep.int(NA_integer_, n)
  best_d2  <- rep.int(Inf, n)
  
  from <- 1L
  while (from <= n) {
    to  <- min(from + block - 1L, n)
    Zb  <- Z[from:to, , drop = FALSE]
    
    # Iterate over B in moderate chunks to limit memory.
    b_from <- 1L
    stepB  <- max(1000L, as.integer(1e7 / max(1L, (to - from + 1L))))  # heuristic
    while (b_from <= m) {
      b_to <- min(b_from + stepB - 1L, m)
      Bb   <- B[b_from:b_to, , drop = FALSE]
      if (nrow(Bb) == 0L) { b_from <- b_to + 1L; next }
      
      # Use outer() to build nb×kb difference matrices; base R has no broadcasting.
      dx <- outer(Zb[, 1], Bb[, 1], "-")
      dy <- outer(Zb[, 2], Bb[, 2], "-")
      d2 <- dx * dx + dy * dy
      
      # For each row in the block, get the argmin in this Bb slice.
      min_loc <- max.col(-d2)                       # index of min per row
      d2_min  <- d2[cbind(seq_len(nrow(Zb)), min_loc)]
      
      upd <- d2_min < best_d2[from:to]
      if (any(upd)) {
        best_d2[from:to][upd] <- d2_min[upd]
        best_idx[from:to][upd] <- (b_from - 1L) + min_loc[upd]
      }
      b_from <- b_to + 1L
    }
    from <- to + 1L
  }
  best_idx
}

# Filter polyline points to those inside polygon (matrix in, matrix out)
.poly_clip <- function(Mxy, poly) {
  if (is.null(Mxy) || length(Mxy) == 0L) return(matrix(numeric(0), ncol = 2))
  if (!inherits(poly, "SpatialPolygons")) return(as.matrix(Mxy))
  inside <- sp::point.in.polygon(Mxy[,1], Mxy[,2],
                                 poly@polygons[[1]]@Polygons[[1]]@coords[,1],
                                 poly@polygons[[1]]@Polygons[[1]]@coords[,2])
  Mxy[inside > 0L, , drop = FALSE]
}

# Detect closed polylines (first == last)
.is_closed <- function(Mxy) {
  Mxy <- as.matrix(Mxy)
  if (nrow(Mxy) < 2L) return(FALSE)
  (Mxy[1,1] == Mxy[nrow(Mxy),1]) && (Mxy[1,2] == Mxy[nrow(Mxy),2])
}

# ---------------------------- main ----------------------------------------

Biplot_boundary_search <- function(
    X = biplot_grid,
    tdp = NA,
    model_select = model_select,
    standardise = FALSE,             # default FALSE (CVA-compatible)
    Xcenter = X_center, Xsd = X_sd,
    CVA_method = FALSE,
    proj = c(1, 2),
    V = V,
    tV = tV,
    set_filters = NA,
    train_ranges = NA
) {
  # Unpack grid inputs
  if (is.na(tdp[1])) tdp <- seq_len(nrow(X$Z.st))
  ct            <- X$ct
  data_polygon  <- X$polygon
  Z.st          <- as.matrix(X$Z.st[tdp, , drop = FALSE])
  biplot_pred_use <- X$pred_use[tdp]
  
  nr_boundaries <- length(ct)
  if (nr_boundaries == 0L) {
    return(list(
      B_boundary_z_value = matrix(0, nrow = nrow(Z.st), ncol = 2),
      B_boundary_x_value = as.data.frame(matrix(NA_real_, nrow = nrow(Z.st),
                                                ncol = length(var_names),
                                                dimnames = list(NULL, var_names))),
      B_boundary_x_pred = rep(NA_real_, nrow(Z.st)),
      z_boundaries_list = list(),
      z_boundary_type = numeric(0),
      nr_boundaries = 0L,
      bb_position = list(),
      position_store = matrix(0, nrow = nrow(Z.st), ncol = 2)
    ))
  }
  
  # Build boundary lists & types
  z_boundaries_list <- vector("list", nr_boundaries)
  z_boundary_type   <- numeric(nr_boundaries)
  for (i in seq_len(nr_boundaries)) {
    Mi <- cbind(ct[[i]]$x, ct[[i]]$y)
    colnames(Mi) <- c("x", "y")
    if (inherits(data_polygon, "SpatialPolygons")) {
      Mi <- .poly_clip(Mi, data_polygon)
    }
    z_boundaries_list[[i]] <- Mi
    z_boundary_type[i]     <- ct[[i]]$level
  }
  
  # Retain relevant contours only
  boundary_used <- integer(0)
  for (i in seq_len(nr_boundaries)) {
    Mi <- z_boundaries_list[[i]]
    if (is.null(Mi) || nrow(Mi) == 0L) next
    is_circle <- .is_closed(Mi)
    keep <- TRUE
    if (is_circle && inherits(data_polygon, "SpatialPolygons")) {
      inside_ct <- sp::point.in.polygon(Z.st[,1], Z.st[,2], Mi[,1], Mi[,2])
      keep <- any(inside_ct > 0L)
    }
    if (keep) boundary_used <- c(boundary_used, i)
  }
  if (length(boundary_used) == 0L) boundary_used <- seq_len(nr_boundaries)
  
  # Rebuild lists with only the used boundaries
  z_boundary_type <- z_boundary_type[boundary_used]
  z_boundaries_list <- lapply(boundary_used, function(j){
    Mi <- cbind(ct[[j]]$x, ct[[j]]$y)
    colnames(Mi) <- c("x", "y")
    if (inherits(data_polygon, "SpatialPolygons")) Mi <- .poly_clip(Mi, data_polygon)
    Mi
  })
  nr_boundaries <- length(z_boundaries_list)
  
  # Prune by model decision consistency (keep points on the correct side)
  if (isTRUE(CVA_method)) standardise <- FALSE
  for (i in seq_len(nr_boundaries)) {
    Mi <- z_boundaries_list[[i]]
    if (is.null(Mi) || nrow(Mi) == 0L) next
    
    Bx <- Mi %*% tV[proj, , drop = FALSE]
    if (isTRUE(standardise)) Bx <- sweep(Bx, 2, Xsd, "*")
    Bx <- sweep(Bx, 2, Xcenter, "+")
    Bx <- as.data.frame(Bx)
    colnames(Bx) <- var_names
    
    rows1 <- if (is.list(set_filters)) get_filter_logical_vector(Bx, set_filters) else rep(TRUE, nrow(Bx))
    rows2 <- if (is.list(train_ranges)) get_filter_logical_vector(Bx, train_ranges) else rep(TRUE, nrow(Bx))
    keep_rows <- (rows1 & rows2)
    
    if (!any(keep_rows)) { z_boundaries_list[[i]] <- Mi[0, , drop = FALSE]; next }
    
    Bx_keep <- Bx[keep_rows, , drop = FALSE]
    Mi_keep <- Mi[keep_rows, , drop = FALSE]
    
    p <- pred_function(model_use = model_use,
                       model_select = model_select,
                       rounding = rounding,
                       new_data = Bx_keep)
    
    ok <- if (z_boundary_type[i] < cut_off) (p < cut_off) else (p >= cut_off)
    if (!any(ok)) {
      z_boundaries_list[[i]] <- Mi[0, , drop = FALSE]
    } else {
      z_boundaries_list[[i]] <- Mi_keep[ok, , drop = FALSE]
    }
  }
  
  # Drop emptied contours
  non_empty <- vapply(z_boundaries_list, function(M) !is.null(M) && nrow(M) > 0L, logical(1))
  if (!all(non_empty)) {
    z_boundaries_list <- z_boundaries_list[non_empty]
    z_boundary_type   <- z_boundary_type[non_empty]
    nr_boundaries     <- length(z_boundaries_list)
  }
  if (nr_boundaries == 0L) {
    return(list(
      B_boundary_z_value = matrix(0, nrow = nrow(Z.st), ncol = 2),
      B_boundary_x_value = as.data.frame(matrix(NA_real_, nrow = nrow(Z.st),
                                                ncol = length(var_names),
                                                dimnames = list(NULL, var_names))),
      B_boundary_x_pred = rep(NA_real_, nrow(Z.st)),
      z_boundaries_list = list(),
      z_boundary_type = numeric(0),
      nr_boundaries = 0L,
      bb_position = list(),
      position_store = matrix(0, nrow = nrow(Z.st), ncol = 2)
    ))
  }
  
  # Nearest points per boundary (streamed, with outer())
  B_boundary_z <- vector("list", nr_boundaries)
  bb_position  <- vector("list", nr_boundaries)
  for (i in seq_len(nr_boundaries)) {
    Mi <- z_boundaries_list[[i]]
    idx <- .nearest_idx_block(Z.st, Mi)       # <- uses outer(), blockwise
    bb_position[[i]] <- idx
    B_boundary_z[[i]] <- Mi[idx, , drop = FALSE]
  }
  
  # Choose best opposing-side contour per point
  boundary_class_0 <- z_boundary_type <  cut_off
  boundary_class_1 <- z_boundary_type >= cut_off
  
  z_boundary_row <- matrix(0, ncol = 2, nrow = nr_boundaries)
  B0 <- matrix(0, nrow = nrow(Z.st), ncol = 2)
  B1 <- matrix(0, nrow = nrow(Z.st), ncol = 2)
  position_store <- matrix(0L, nrow = nrow(Z.st), ncol = 2)
  
  for (j in seq_len(nrow(Z.st))) {
    for (i in seq_len(nr_boundaries)) {
      Bij <- B_boundary_z[[i]][j, ]
      if (length(Bij) == 0L || any(!is.finite(Bij))) {
        z_boundary_row[i, ] <- c(1e9, 1e9)
      } else {
        z_boundary_row[i, ] <- Bij
      }
    }
    d2_all <- rowSums((z_boundary_row - matrix(Z.st[j, ], nrow = nr_boundaries, ncol = 2, byrow = TRUE))^2)
    
    d2_1 <- d2_all[boundary_class_1]
    d2_0 <- d2_all[boundary_class_0]
    
    pos1 <- if (length(d2_1)) which.min(d2_1) else NA_integer_
    pos0 <- if (length(d2_0)) which.min(d2_0) else NA_integer_
    
    if (is.na(pos1)) {
      position_store[j, 1] <- 0L; B0[j, ] <- c(NA_real_, NA_real_)
    } else {
      k1 <- which(boundary_class_1)[pos1]
      position_store[j, 1] <- k1
      B0[j, ] <- B_boundary_z[[k1]][j, ]
    }
    
    if (is.na(pos0)) {
      position_store[j, 2] <- 0L; B1[j, ] <- c(NA_real_, NA_real_)
    } else {
      k0 <- which(boundary_class_0)[pos0]
      position_store[j, 2] <- k0
      B1[j, ] <- B_boundary_z[[k0]][j, ]
    }
  }
  
  group1 <- (biplot_pred_use == 1)
  group0 <- !group1
  B0[group1, ] <- 0
  B1[group0, ] <- 0
  B_boundary_z_value <- B0 + B1
  
  # Back-project to X space & score
  B_boundary_x_value <- B_boundary_z_value %*% tV[proj, , drop = FALSE]
  if (isTRUE(standardise)) B_boundary_x_value <- sweep(B_boundary_x_value, 2, Xsd, "*")
  B_boundary_x_value <- sweep(B_boundary_x_value, 2, Xcenter, "+")
  B_boundary_x_value <- as.data.frame(B_boundary_x_value)
  colnames(B_boundary_x_value) <- var_names
  
  B_boundary_x_pred <- pred_function(model_use = model_use,
                                     model_select = model_select,
                                     rounding = rounding,
                                     new_data = B_boundary_x_value)
  
  return(list(
    B_boundary_z_value = B_boundary_z_value,
    B_boundary_x_value = B_boundary_x_value,
    B_boundary_x_pred = B_boundary_x_pred,
    z_boundaries_list = z_boundaries_list,
    z_boundary_type = z_boundary_type,
    nr_boundaries = nr_boundaries,
    bb_position = bb_position,
    position_store = position_store
  ))
}