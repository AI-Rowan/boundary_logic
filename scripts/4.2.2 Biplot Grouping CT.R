############################################################
# Script: 4.2.2 Biplot grouping with contour
# Description:
#   Assign each projected point (rows of Z.st) to the nearest contour
#   (ct level) inside the training-data polygon. This is a grouping
#   helper used by the surrogate workflow
#
# Inputs expected in the environment:
#   - Z.st          : n x 2 matrix of biplot coordinates (from biplot_grouping$Z.st)
#   - ct            : list of contour sets as returned by biplot_plane(... )$ct
#   - data_polygon  : SpatialPolygons object (from biplot_grouping$polygon)
#
# Outputs placed in the environment:
#   - boundary_set_nr : integer vector (length nrow(Z.st)) indicating the
#                       closest contour set index for each point
#   - boundary_used   : integer vector of indices of contour sets retained
#
# Notes vs. boundary-search code:
#   * No removal of "incorrect" predictions (no CF enforcing).
#   * No reordering to force counterfactual pairing.
############################################################

# ----------------------------- 1) Build boundaries ------------------------
nr_boundaries    <- length(ct)
z_boundaries_list <- vector("list", nr_boundaries)
z_boundary_type   <- numeric(nr_boundaries)   # Stores ct level (e.g., 0.495 or 0.505)

for (i in seq_len(nr_boundaries)) {
  boundary_vals <- cbind(ct[[i]]$x, ct[[i]]$y)      # two-column (x,y)
  colnames(boundary_vals) <- c("x", "y")
  z_boundaries_list[[i]] <- boundary_vals
  z_boundary_type[i]     <- ct[[i]]$level
}

# ---------------- 1.2) Retain boundaries that cover/contain data ---------
boundary_used <- integer(0)
counter <- 0L

for (i in seq_len(nr_boundaries)) {
  # Keep only contour points inside polygon
  data_points <- as.data.frame(z_boundaries_list[[i]])
  coordinates(data_points) <- ~ x + y
  data_points$inside_polygon <- !is.na(over(data_points, data_polygon))
  z_boundaries_list[[i]] <- as.matrix(z_boundaries_list[[i]][data_points$inside_polygon, , drop = FALSE])

  # Detect closed contour
  circle_contour <- 0L
  if (length(z_boundaries_list[[i]]) > 0L) {
    circle_contour <- ifelse(
      z_boundaries_list[[i]][1, 1] == z_boundaries_list[[i]][nrow(z_boundaries_list[[i]]), 1] &
      z_boundaries_list[[i]][1, 2] == z_boundaries_list[[i]][nrow(z_boundaries_list[[i]]), 2],
      1L, 0L
    )
  }

  # If closed, check whether it surrounds (any) data
  interior_data <- 0L
  if (length(z_boundaries_list[[i]]) > 0L && circle_contour == 1L) {
    interior_data <- sum(point.in.polygon(Z.st[, 1], Z.st[, 2],
                                          z_boundaries_list[[i]][, 1], z_boundaries_list[[i]][, 2]))
  }

  # Keep if not a circle, or a circle that encloses data
  if ((circle_contour == 0L || interior_data > 0L) && length(z_boundaries_list[[i]]) > 0L) {
    counter <- counter + 1L
    boundary_used[counter] <- i
  }
}

# ---------------- 1.3) Rebuild with retained boundaries only -------------
if (!is.null(boundary_used)) {
  nr_boundaries <- length(boundary_used)
} else {
  nr_boundaries <- length(ct)
  boundary_used <- seq_len(nr_boundaries)
}

z_boundaries_list <- vector("list", nr_boundaries)
z_boundary_type   <- numeric(nr_boundaries)

for (i in seq_len(nr_boundaries)) {
  j <- boundary_used[i]
  boundary_vals <- cbind(ct[[j]]$x, ct[[j]]$y)
  colnames(boundary_vals) <- c("x", "y")

  data_points <- as.data.frame(boundary_vals)
  coordinates(data_points) <- ~ x + y
  data_points$inside_polygon <- !is.na(over(data_points, data_polygon))
  z_boundaries_list[[i]] <- as.matrix(boundary_vals[data_points$inside_polygon, , drop = FALSE])

  z_boundary_type[i] <- ct[[j]]$level
}

# ---------------- 2) Distances & nearest-per-boundary --------------------
# For each boundary set i, compute distances from each Z.st row to all
# contour points in that set; then select the nearest point index.
point.dist.matrix <- vector("list", nr_boundaries)
for (i in seq_len(nr_boundaries)) {
  point.dist.matrix[[i]] <- as.matrix(pdist(Z.st, z_boundaries_list[[i]]))
}

nearest_idx <- vector("list", nr_boundaries)
for (i in seq_len(nr_boundaries)) {
  if (nrow(point.dist.matrix[[i]]) > 0L) {
    nearest_idx[[i]] <- apply(point.dist.matrix[[i]], 1, which.min)
  } else {
    nearest_idx[[i]] <- integer(0)
  }
}

# ---------------- 3) Select the closest boundary set per point -----------
z_boundary_row   <- matrix(0, ncol = 2, nrow = nr_boundaries)
boundary_set_nr  <- integer(nrow(Z.st))

for (j in seq_len(nrow(Z.st))) {
  for (i in seq_len(nr_boundaries)) {
    if (length(nearest_idx[[i]]) == 0L) {
      z_boundary_row[i, ] <- c(1e8, 1e8)  # ignore empty sets
    } else {
      z_boundary_row[i, ] <- z_boundaries_list[[i]][nearest_idx[[i]][j], , drop = FALSE]
    }
  }
  # Choose boundary set whose nearest point is closest to Z.st[j, ]
  boundary_set_nr[j] <- apply(as.matrix(pdist(Z.st[j, , drop = FALSE], z_boundary_row)), 1, which.min)
}

# (Outputs in environment): boundary_set_nr, boundary_used
