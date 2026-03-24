############################################################
# Contour inspection — extract, back-project, and score
# decision boundary contour points
#
# Purpose : Investigates the contour lines used as potential
#           counterfactuals by:
#             1. Extracting Z-space contour coordinates
#             2. Back-projecting to original X-space
#             3. Scoring through the fitted model
#             4. Flagging points inside / outside the
#                training data convex hull polygon
#             5. Overlaying back onto the biplot
#             6. Matching selected counterfactuals from
#                bl_find_boundary() back to contour points
#
# Requires: a bl_result object named bl_results and a
#           bl_boundary object named bl_bnd already built
#           in the session (e.g. from
#           01_iris_SVM_PCA_biplot.R)
#
# Run interactively: Ctrl+Enter line by line
############################################################


# ---- 1. Inspect the raw contour structure --------------------------------
ct <- bl_results$biplot_grid$ct

cat("Contour segments :", length(ct), "\n")
cat("Contour levels   :", unique(sapply(ct, `[[`, "level")), "\n")
# levels will be e.g. 0.49 / 0.51 for rounding = 2L

# Points per segment
sapply(ct, function(cl) length(cl$x))


# ---- 2. Combine all contour points into one Z-space matrix ---------------
ct_z <- do.call(rbind, lapply(ct, function(cl) {
  cbind(x     = cl$x,
        y     = cl$y,
        level = cl$level)
}))
ct_z <- as.data.frame(ct_z)


# ---- 3. Back-project from Z-space to X-space -----------------------------
# Z -> X_standardised via the inverse projection matrix (tV)
tVr  <- bl_results$tV[bl_results$proj_dims, ]          # 2 x p
X_st <- as.matrix(ct_z[, c("x", "y")]) %*% tVr        # n x p

# Unscale back to original feature space
if (bl_results$standardise) {
  X_orig <- sweep(X_st, 2, bl_results$X_sd,     "*")
}
X_orig <- sweep(X_orig, 2, bl_results$X_center, "+")
X_orig <- as.data.frame(X_orig)
names(X_orig) <- bl_results$var_names


# ---- 4. Score each contour point through the model -----------------------
ct_preds <- bl_predict(bl_results, data = X_orig)

# Attach Z-space coordinates and contour level
ct_preds$z1            <- ct_z$x
ct_preds$z2            <- ct_z$y
ct_preds$contour_level <- ct_z$level


# ---- 5. Polygon membership -----------------------------------------------
if (!is.null(bl_results$polygon)) {
  Z_df <- data.frame(x = ct_z$x, y = ct_z$y)
  sp::coordinates(Z_df) <- ~x + y
  ct_preds$inside_polygon <- !is.na(sp::over(Z_df, bl_results$polygon))
} else {
  ct_preds$inside_polygon <- TRUE   # no polygon — treat all as inside
}


# ---- 6. Inspect ----------------------------------------------------------
ct_preds$dist_to_cutoff <- abs(ct_preds$pred_prob - bl_results$cutoff)

head(ct_preds)
summary(ct_preds$pred_prob)

# Points inside the polygon only
ct_preds[ct_preds$inside_polygon, ]


# ---- 7. Match bl_find_boundary() counterfactuals to contour points ------
# B_z values are interpolated along contour segments so they are not exact
# coordinate matches. For each valid B_z row, find the nearest contour point
# in ct_preds by Z-space Euclidean distance, then tag it.

B_z <- bl_bnd$B_z   # n x 2 matrix; NA where no boundary was found

# Initialise columns
ct_preds$is_counterfactual <- FALSE
ct_preds$cf_count          <- 0L    # how many observations map to this point
ct_preds$cf_obs_rows       <- NA_character_  # which observation row indices

ct_mat <- as.matrix(ct_preds[, c("z1", "z2")])

valid_obs <- which(is.finite(B_z[, 1L]) & is.finite(B_z[, 2L]))

for (j in valid_obs) {
  bz <- B_z[j, , drop = FALSE]

  # Euclidean distance from this B_z to every contour point
  dists <- sqrt(rowSums(sweep(ct_mat, 2, bz, "-")^2))
  nearest <- which.min(dists)

  ct_preds$is_counterfactual[nearest] <- TRUE
  ct_preds$cf_count[nearest]          <- ct_preds$cf_count[nearest] + 1L

  # Accumulate observation row indices as a comma-separated string
  existing <- ct_preds$cf_obs_rows[nearest]
  ct_preds$cf_obs_rows[nearest] <- if (is.na(existing)) {
    as.character(j)
  } else {
    paste(existing, j, sep = ",")
  }
}

# Summary
cat("Contour points selected as counterfactuals :",
    sum(ct_preds$is_counterfactual), "\n")
cat("Observations with no boundary found        :",
    sum(!is.finite(B_z[, 1L])), "\n")

# Inspect the matched counterfactual points
ct_preds[ct_preds$is_counterfactual, ]

# Points used by more than one observation (shared boundary region)
ct_preds[ct_preds$cf_count > 1L, ]


# ---- 8. Overlay back onto the biplot -------------------------------------
plot_biplotEZ(bl_results)

# All contour points: green = inside polygon, grey = outside
points(ct_preds$z1, ct_preds$z2,
       pch = 16, cex = 0.4,
       col = ifelse(ct_preds$inside_polygon, "green3", "grey70"))

# Highlight the points actually selected as counterfactuals
cf_pts <- ct_preds[ct_preds$is_counterfactual, ]
points(cf_pts$z1, cf_pts$z2,
       pch = 4, cex = 0.8, lwd = 1.5, col = "red")

cf_pts <- ct_preds[ct_preds$inside_polygon, ]
points(cf_pts$z1, cf_pts$z2,
       pch = 4, cex = 0.8, lwd = 1.5, col = "blue")


cf_pts <- ct_preds[ct_preds$row==380, ]
points(cf_pts$z1, cf_pts$z2,
       pch = 4, cex = 0.8, lwd = 1.5, col = "blue")

# ===========================================================
# BIPLOT POINT PICKER
# Click on the active biplot to find the nearest data point.
# Press Escape to stop clicking.
# ===========================================================

# ---- Usage ---------------------------------------------------------------
# First render the biplot, then call the picker on the same plot device.

plot_biplotEZ(bl_results)
picked <- bl_pick_point(bl_results)

# Inspect all picked points
picked

# Pick from test data instead
# plot_biplotEZ(bl_results, points = bl_project_points(bl_results$test_data, bl_results))
# picked <- bl_pick_point(bl_results, data = bl_results$test_data)


cf_pts
ct_preds[360:590,]

getOption("max.print" = 100000)




bp_Z <- bl_results$biplot_obj$Z[, bl_results$proj_dims]

sv <- if (bl_results$standardise) bl_results$X_sd else FALSE
X_st <- scale(as.matrix(bl_results$train_data[, bl_results$var_names]),
              center = bl_results$X_center, scale = sv)
manual_Z <- X_st %*% bl_results$V[, bl_results$proj_dims]

all.equal(bp_Z, manual_Z, check.attributes = FALSE)
