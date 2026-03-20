############################################################
# Pima diabetes — Phases 1, 2 & 3: Global and local interpretation
#
# Dataset : scripts/pima diabetes.csv
#           768 observations, 8 numeric features, binary outcome
# Model   : GAM with spline terms (mgcv::gam)
# Biplot  : CVA
#
# Purpose : Demonstrates the full boundary logic workflow:
#
#           Phase 1 (Steps 1-6): data prep, model fitting, biplot
#           Phase 2 (Steps 7-9): global interpretations
#             7. Find nearest boundary point for each observation
#             8. Distance-to-boundary plot (jitter + boxplot)
#                  — also prints robustness score and per-variable totals
#             9. Surrogate model
#           Phase 3 (Steps 10-18): local interpretation of one target
#            10. Inspect predictions and select target
#            11. (removed — merged into Step 10)
#            12. Set actionability constraints
#            13. Find local counterfactual via biplot rotation
#            14. Local biplot plot
#            15. Shapley contribution plot
#            16. Sparse counterfactual
#            17. (Optional) Unconstrained local search
#            18. External data point
#
# Run interactively: Ctrl+Enter line by line
############################################################

# ---- Load the package (development workflow) ---------------------------
rm(list = ls())
devtools::load_all()
library(dplyr)


# ===========================================================
# PHASE 1 — Build the pipeline
# ===========================================================

# ---- Step 1: Prepare data ---------------------------------------------
pima_raw <- read.csv(
  file.path(dirname(rstudioapi::getActiveDocumentContext()$path),
            "pima diabetes.csv")
)

pima_raw <- pima_raw %>% filter(Insulin > 0, SkinThickness > 0)

bl_dat <- bl_prepare_data(
  data           = pima_raw,
  class_col      = "Outcome",
  target_class   = NULL,       # already 0/1
  train_fraction = 0.8,
  seed           = 121L
)

print(bl_dat)


# ---- Step 2: Filter outliers ------------------------------------------
bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 0.9)

#print(bl_filt)


# ---- Step 3: Fit GAM and wrap -----------------------------------------
gam_formula <- class ~ s(Glucose) + s(BloodPressure) +
  s(Pregnancies) + s(SkinThickness) + Insulin +
  s(BMI) + DiabetesPedigreeFunction + Age

custom_gam <- mgcv::gam(
  formula = gam_formula,
  data    = bl_filt$train_data,
  family  = binomial(link = "logit")
)

summary(custom_gam)

bl_mod <- bl_wrap_model(
  model      = custom_gam,
  model_type = "custom",
  var_names  = bl_filt$var_names,
  predict_fn = function(m, new_data) {
    as.numeric(mgcv::predict.gam(m, newdata = new_data, type = "response"))
  },
  train_data = bl_filt$train_data
)

print(bl_mod)


# ---- Steps 4-6: Build projection, grid, and assemble ------------------
# bl_build_result() wraps bl_build_projection() + bl_build_grid() +
# bl_assemble() into a single call.
#   - Pass bl_filt when outlier filtering was applied; bl_dat otherwise.
#   - Set method = "PCA" for a variance-based projection.
#   - Omit bl_mod to get an exploratory biplot without a prediction surface:
#       bl_proj <- bl_build_result(bl_filt, method = "CVA"); plot(bl_proj)
#
# rounding controls only the decision boundary contour band width:
#   rounding = 2  ->  band = cutoff ± 0.01  (contours at 0.49 / 0.51)
#   rounding = 3  ->  band = cutoff ± 0.001 (contours at 0.499 / 0.501)
# All model predictions are always rounded to 3 decimal places regardless.
#
# Note: train_fraction = 1 means all data is used for training (no holdout
# test set). bl_results$test_data is automatically set to train_data so
# that bl_predict(), bl_project_points(), and bl_find_boundary() all work.

bl_results <- bl_build_result(
  bl_data  = bl_filt,
  bl_model = bl_mod,
  method   = "CVA",
  title    = "Pima diabetes — GAM, CVA biplot",
  rounding = 2L   # contour band = cutoff ± 0.001
)

print(bl_results)

# Reference biplot — training data (confusion colours)
plot_biplotEZ(bl_results)

# Project all observations and overlay on the biplot
test_pts <- bl_project_points(bl_results$test_data, bl_results)
plot_biplotEZ(bl_results, points = test_pts)


# ===========================================================
# PHASE 2 — Global interpretations (population level)
# ===========================================================

# ---- Step 7: Find nearest boundary point for every observation --------
# bl_find_boundary() searches the decision boundary contour for each row
# of data and back-projects the nearest boundary point to X-space.
#
# Default data = bl_results$test_data (= train_data when train_fraction = 1).
# To run on a custom subset use the tdp argument:
#   bl_find_boundary(bl_results, tdp = 1:50)
#
# Key outputs:
#   B_z    — Z-space boundary coordinate per observation
#   B_x    — X-space counterfactual per observation
#   dist_z — Euclidean distance to boundary in Z-space
#   dist_x — Standardised distance to boundary in X-space

bl_bnd <- bl_find_boundary(bl_results)
print(bl_bnd)

# Biplot with per-observation boundary arrows overlaid.
# Each arrow runs from the observation to its nearest boundary point.
plot_biplotEZ(bl_results, points = test_pts, boundary = bl_bnd)


# ---- Step 8: Distance-to-boundary plot --------------------------------
# Decomposes each observation's distance to the boundary into a signed,
# standardised contribution per variable. Variables are sorted from least
# to most important (bottom to top).
#
# Y-axis: "variable : total absolute standardised distance" (importance proxy)
# X-axis: signed standardised distance — positive = observation is above
#         the boundary in that variable's direction.
#
# Colour scheme (when true labels are present):
#   TP = red, TN = blue, FP = purple, FN = orange
#
# Two chart types:
#   plot(bl_bnd)                   — jitter (default; one point per obs)
#   plot(bl_bnd, type = "boxplot") — box-and-whisker by confusion group

plot(bl_bnd)
plot(bl_bnd, type = "boxplot")
# Both calls also print the overall robustness scalar and per-variable totals
# to the console. Use bl_robustness(bl_bnd) to access these values
# programmatically (e.g. to compare robustness across models).


# ---- Step 9: Surrogate model -----------------------------------------
# Assigns each observation a class based solely on its position in Z-space
# relative to the hull-clipped surrogate contour lines. This provides a
# purely spatial approximation of the classifier confined to the training hull.
#
# Default data = bl_results$train_data.
# Observations that project outside the training hull receive surrogate_pred = NA.
#
# Accuracy reported against:
#   accuracy_vs_model  — agreement with the fitted model's predictions
#   accuracy_vs_labels — accuracy vs true class labels (if available)
#
# plot() shows the prediction grid clipped to the hull, with observations
# coloured by their surrogate-assigned class (red = 1, blue = 0).

bl_surr <- bl_surrogate(bl_results)
print(bl_surr)
plot(bl_surr)


# ===========================================================
# PHASE 3 — Local interpretation of a single target point
# ===========================================================

# ---- Step 10: Inspect predictions and select target -------------------
# bl_predict() scores every row and returns a tidy data frame with
# row, pred_prob, pred_class, true_class, confusion, and all features.
# Default data = bl_results$test_data.
# Review pred_summary to choose the row you want to explain locally,
# then assign it to tdp. bl_select_target() wraps that observation into
# a bl_target object, projecting it into Z-space and scoring it through
# the model.

pred_summary <- bl_predict(bl_results)
print(pred_summary)

tdp <- 1
tgt <- bl_select_target(bl_results, target = tdp)
print(tgt)

# Highlight the target on the main biplot
plot_biplotEZ(
  bl_results,
  points       = test_pts,
  target_point = unlist(tgt$x_obs),
  target_label = tdp
)


# ---- Step 12: Set actionability constraints ---------------------------
# Specify which variable directions are clinically plausible for this
# patient. Omit set_filters() entirely if no constraints are needed.
#
# Valid constraint types:
#   "decrease"  — counterfactual value must be <= observed
#   "increase"  — counterfactual value must be >= observed
#   "fixed"     — counterfactual search allows +/- 0.5 of the observed
#                 value. When a sparse CF is generated, fixed variables
#                 always revert to the original observed value.
#   c(min, max) — counterfactual must lie within this absolute range

flt <- set_filters(
  tgt,
  Glucose     = "decrease",  # reducing Glucose is the clinical lever
  Age         = "fixed",     # Age is not actionable
  Pregnancies = "fixed"      # Pregnancies is not actionable
)

print(flt)


# ---- Step 13: Find local counterfactual via biplot rotation ----------
# For each eigenvector pair (1,2), (1,3), (2,3), (1,4), ...:
#   - Rotates the biplot so the target lies on the projection plane
#   - Builds a local prediction grid in the rotated Z-space
#   - Applies train_ranges (feasibility) + set_filters (actionability)
#   - NOTE: the Z-space convex hull polygon is NOT applied here —
#     rotation invalidates the original polygon coordinates
#   - Finds the nearest opposing-class boundary point
# The pair with the overall minimum Z-space distance is selected.

bl_local <- bl_find_local_cf(
  bl_result   = bl_results,
  bl_target   = tgt,
  set_filters = flt
)

print(bl_local)


# ---- Step 14: Local biplot plot ---------------------------------------
# Rendered using the biplotEZ pipeline (same style as plot_biplotEZ).
# The biplotEZ object is patched with the rotated loading matrix (Vr_rot)
# and rotated training coordinates (Z_train_rot) for the best pair.
#
# Default (no_points = TRUE): shows only the target + counterfactual.
#   - Filled circle = target observation (confusion colour or red/blue)
#   - Grey cross + arrow = nearest boundary (counterfactual)
#
# To also show training data:
#   plot(bl_local, no_points = FALSE)
#   Training colours: TP=red, TN=blue, FP=purple, FN=orange
#
# Supported parameters (same as plot_biplotEZ):
#   no_grid, no_points, no_contour, cex_z, label_dir, tick_label_cex,
#   ticks_v, which, X_names, label_offset_var, label_offset_dist,
#   show_arrows, arrow_col, contour_col, contour_lwd, contour_lty

plot(bl_local)


# ---- Step 15: Shapley contribution plot --------------------------------
# Decomposes the path from target to counterfactual into per-variable
# contributions using the exact Shapley formula (p <= 14 variables).
#
# Each variable is labelled as:
#   Supports    (darkblue) — moving in this direction helps cross the boundary
#   Contradicts (grey)     — moving in this direction works against crossing
#
# Y-axis labels: "variable: observed_value -> required_change"
# X-axis: Shapley value = contribution to the probability shift

bl_shapley_values <- bl_shapley(bl_local)
print(bl_shapley_values)

plot(bl_shapley_values)


# ---- Step 16: Sparse counterfactual -----------------------------------
# Automatically zeroes out Contradicts variables — they revert to their
# observed values. Only Supports variables take their counterfactual value.
# Variables constrained as "fixed" in set_filters() always revert to their
# observed value in the sparse CF, regardless of Shapley class.
#
# round_to = NULL  keeps exact CF values
# round_to = 0.5   rounds to the nearest 0.5 (recommended default)
# round_to = 1     rounds to nearest whole number
#
# solution_valid = TRUE  means the sparse CF still crosses the boundary.
# solution_valid = FALSE means the zeroed/rounded CF is not enough;
#   try round_to = 1, remove round_to, or relax actionability filters.

bl_sparse <- bl_find_sparse_cf(bl_shapley_values, round_to = NULL)

# print() shows three prediction values to verify classification changes:
#   Observed  — model score for the original data point
#   Full CF   — model score at the boundary counterfactual
#   Sparse CF — model score for the sparse (zeroed) counterfactual
print(bl_sparse)

# Local biplot with both CFs overlaid (training points hidden by default).
#   Grey cross   = full counterfactual (all variables changed)
#   Green cross  = sparse CF (valid — crosses boundary)
#   Yellow cross = sparse CF (invalid — does not cross boundary)
# To also show training data: plot(bl_sparse, no_points = FALSE)
# All plot(bl_local) parameters are supported via ...
plot(bl_sparse)


# ===========================================================
# ALTERNATIVE: run without actionability constraints
# Useful as a baseline before adding constraints, or when
# bl_local$solution_found == FALSE to confirm a boundary exists.
# ===========================================================

# ---- Step 17 (optional): Unconstrained local search ------------------
bl_local_free <- bl_find_local_cf(bl_results, tgt)
print(bl_local_free)

plot(bl_local_free)

bl_shap_free   <- bl_shapley(bl_local_free)
print(bl_shap_free)
bl_sparse_free <- bl_find_sparse_cf(bl_shap_free, round_to = NULL)
print(bl_sparse_free)

plot(bl_shap_free)
plot(bl_sparse_free)


# ===========================================================
# EXTERNAL DATA POINT
# Analyse a new patient who does not appear in the train/test data.
# Pass a single-row data frame directly as the target.
# No class column is required.
# ===========================================================

# ---- Step 18: External patient ----------------------------------------
new_patient <- data.frame(
  Pregnancies              = 2,
  Glucose                  = 148,
  BloodPressure            = 72,
  SkinThickness            = 35,
  Insulin                  = 100,
  BMI                      = 33.6,
  DiabetesPedigreeFunction = 0.627,
  Age                      = 50
)

tgt_ext <- bl_select_target(bl_results, target = new_patient)
print(tgt_ext)

# Highlight on main biplot
plot_biplotEZ(
  bl_results,
  points       = test_pts,
  target_point = unlist(tgt_ext$x_obs),
  target_label = "new"
)

# Local search — no actionability constraints for this new patient
bl_local_ext <- bl_find_local_cf(bl_results, tgt_ext)
print(bl_local_ext)

plot(bl_local_ext)

bl_shap_ext   <- bl_shapley(bl_local_ext)
bl_sparse_ext <- bl_find_sparse_cf(bl_shap_ext)

plot(bl_shap_ext)
plot(bl_sparse_ext)
