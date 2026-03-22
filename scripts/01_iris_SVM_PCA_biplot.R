############################################################
# Iris — SVM, PCA biplot with decision boundary and counterfactuals
#
# Dataset : datasets::iris (binary: versicolor vs. rest)
# Model   : SVM (RBF kernel, via kernlab)
# Biplot  : PCA
# Grid    : 50 x 50
#
# Purpose : Demonstrates Phase 1 and Phase 2 (Step 7) of the
#           boundary logic workflow on the iris dataset:
#
#           Phase 1 (Steps 1-6) : data prep, SVM fit, PCA biplot
#           Phase 2 (Step 7)    : nearest boundary point per
#                                 observation with arrows overlaid
#
# Run interactively: Ctrl+Enter line by line
############################################################

# ---- Load the package (development workflow) ---------------------------
rm(list = ls())
devtools::load_all()


# ===========================================================
# PHASE 1 — Build the pipeline
# ===========================================================

# ---- Step 1: Prepare data ---------------------------------------------
# target_class = "versicolor" creates a binary outcome:
#   class 1 = versicolor, class 0 = all other species
bl_dat <- bl_prepare_data(
  data           = datasets::iris[,c(1:3,5)],
  class_col      = "Species",
  target_class   = "versicolor",
  train_fraction = 0.8,
  seed           = 42L
)

print(bl_dat)


# ---- Step 2: Filter outliers ------------------------------------------
bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 1)

plot_biplotEZ(bl_dat)


# ---- Step 3: Fit SVM --------------------------------------------------
# SVM (RBF kernel) fitted via parsnip + kernlab.
# probability = TRUE is handled internally so that .pred_function()
# can extract class-1 probabilities.

bl_mod <- bl_fit_model(
  train_data = bl_filt$train_data,
  var_names  = bl_filt$var_names,
  model_type = "GLM"
)

print(bl_mod)

# ---- Steps 4-6: Build projection, grid, and assemble ------------------
# method = "PCA" — variance-based projection.
# m = 50L        — 50 x 50 prediction grid.
# rounding = 2L  — contour band = cutoff ± 0.01

bl_results <- bl_build_result(
  bl_data  = bl_filt,
  bl_model = bl_mod,
  method   = "PCA",
  m        = 300,
  calc_hull = T,
  #title    = "Iris — SVM, PCA biplot",
  rounding = 2L
)

print(bl_results)

# Reference biplot — training data, confusion colours
plot_biplotEZ(bl_results, rotate_deg = 180)

# Predicted-class colours only (red = class 1, blue = class 0)
plot_biplotEZ(bl_results, confusion_cols = F, cex = 2, rotate_deg = 180)


# ===========================================================
# PHASE 2 — Step 7: Counterfactuals with boundary arrows
# ===========================================================

# ---- Step 7: Find nearest boundary point for every observation --------
# bl_find_boundary() searches the decision boundary contour for each row
# and back-projects the nearest boundary point to X-space.
#
# Key outputs:
#   B_z    — Z-space boundary coordinate per observation
#   B_x    — X-space counterfactual per observation
#   dist_z — Euclidean distance to boundary in Z-space
#   dist_x — Standardised distance to boundary in X-space
bl_find_boundary
bl_bnd <- bl_find_boundary(bl_results, data = bl_results$train_data)
print(bl_bnd)

# Biplot with per-observation boundary arrows overlaid.
# Each arrow runs from the observation to its nearest boundary point.
# Confusion colours (TP/TN/FP/FN):
plot_biplotEZ(bl_results, boundary = bl_bnd)

# Predicted-class colours only (red / blue):
plot_biplotEZ(bl_results, boundary = bl_bnd, confusion_cols = FALSE, 
              rotate_deg = 180, cex = 2)



