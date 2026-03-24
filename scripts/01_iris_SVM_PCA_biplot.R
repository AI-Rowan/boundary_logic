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

# irisd <- iris
# irisd[,c(1:4)] <- scale(irisd[,c(1:4)],T,T)

# ---- Step 1: Prepare data ---------------------------------------------
# target_class = "versicolor" creates a binary outcome:
#   class 1 = versicolor, class 0 = all other species
bl_dat <- bl_prepare_data(
  data           = datasets::iris[,c(1:4,5)], #irisd[,c(1:3,5)],#
  class_col      = "Species",
  target_class   = "versicolor",
  train_fraction = 0.8,
  seed           = 121L
)




print(bl_dat)


# ---- Step 2: Filter outliers ------------------------------------------
bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 1)



# ---- Step 3: Fit SVM --------------------------------------------------
# SVM (RBF kernel) fitted via parsnip + kernlab.
# probability = TRUE is handled internally so that .pred_function()
# can extract class-1 probabilities.

bl_mod <- bl_fit_model(
  train_data = bl_filt$train_data,
  var_names  = bl_filt$var_names,
  model_type = "SVM"
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
  standardise = T,
  m        = 200,
  calc_hull = T,
  outlie = 1,
  #title    = "Iris — SVM, PCA biplot",
  rounding = 2L
)


# Predicted-class colours only (red = class 1, blue = class 0)
plot_biplotEZ(bl_results, cex = 2, rotate_deg = 0, confusion_cols = T,
              no_grid = F,
              label_cex = 1.5, tick_label_cex = 1, ticks_v = 2)


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

bl_bnd <- bl_find_boundary(bl_results, data = bl_results$train_data)

# Biplot with per-observation boundary arrows overlaid.
# Each arrow runs from the observation to its nearest boundary point.
# Confusion colours (TP/TN/FP/FN):
plot_biplotEZ(bl_results, boundary = bl_bnd)

# Predicted-class colours only (red / blue):
plot_biplotEZ(bl_results, boundary = bl_bnd, confusion_cols = FALSE, cex = 2,
              rotate_deg = 0, label_cex = 1.5, tick_label_cex = 1, ticks_v = 2)

plot(bl_bnd, cex = 2)
plot(bl_bnd, type = "boxplot")


bl_surr <- bl_surrogate(bl_results)
print(bl_surr, )
plot(bl_surr, label_cex = 1.5, cex = 2)





pred_summary <- bl_predict(bl_results)
print(pred_summary)

tdp <- 15
tgt <- bl_select_target(bl_results, target = tdp)
print(tgt)

# Highlight the target on the main biplot
test_unlabeled <- bl_results$test_data[ tdp, !(names(bl_results$test_data) %in% "class")]
holdout_pts <- bl_project_points(test_unlabeled, bl_results)

bl_bnd <- bl_find_boundary(bl_results, data = bl_results$test_data[15,])

plot_biplotEZ(bl_results, points = holdout_pts)


plot_biplotEZ(

  bl_results, 
  points = holdout_pts,
  boundary = bl_bnd,
  
  cex = 1, rotate_deg = 0, confusion_cols = T,
  no_grid = F,
  label_cex = 1.5, tick_label_cex = 1, ticks_v = 2,


  target_point = unlist(tgt$x_obs),
  target_label = tdp
)

bl_local <- bl_find_local_cf(
  bl_result   = bl_results,
  bl_target   = tgt
)

bl_shapley_values <- bl_shapley(bl_local)
print(bl_shapley_values)

plot(bl_shapley_values)


plot(bl_local)
