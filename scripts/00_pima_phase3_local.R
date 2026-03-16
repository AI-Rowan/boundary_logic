############################################################
# Pima diabetes — Phase 3: Local counterfactual interpretation
#
# Dataset : scripts/pima diabetes.csv
#           768 observations, 8 numeric features, binary outcome
# Model   : GAM with spline terms (mgcv::gam)
# Biplot  : CVA
#
# Purpose : Demonstrates the Phase 3 local interpretation workflow
#           for a single target observation:
#             1. Select a target data point
#             2. Optionally constrain actionable variable directions
#             3. Find the nearest decision boundary via biplot rotation
#                (up to 10 eigenvector pairs tested)
#             4. Shapley plot: per-variable contribution to the move
#             5. Sparse counterfactual: only Supports variables changed
#             6. Visualise the final result on the rotated local biplot
#
#           Also shows how to pass an external (new) data point that is
#           not in the training or test set.
#
# Run interactively: Ctrl+Enter line by line
############################################################

# ---- Load the package (development workflow) ---------------------------
rm(list = ls())
devtools::load_all()
library(dplyr)


# ===========================================================
# PHASE 1 — Build the pipeline (same setup as 00_pima_GAM_CVA.R)
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
bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 0.8)

print(bl_filt)


# ---- Step 3: Fit GAM and wrap -----------------------------------------
gam_formula <- class ~ s(Glucose) + s(BloodPressure) +
  s(Pregnancies) + SkinThickness + Insulin +
  s(BMI) + DiabetesPedigreeFunction + s(Age)

custom_gam <- mgcv::gam(
  formula = gam_formula,
  data    = bl_filt$train_data,
  family  = binomial(link = "logit")
)

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


# ---- Step 4: CVA projection -------------------------------------------
bl_proj <- bl_build_projection(
  train_data = bl_filt$train_data,
  var_names  = bl_filt$var_names,
  method     = "CVA",
  proj_dims  = c(1L, 2L),
  bl_model   = bl_mod,
  title      = "Pima diabetes — GAM, CVA biplot (Phase 3 local)"
)

print(bl_proj)


# ---- Step 5: Build prediction grid ------------------------------------
bl_grid <- bl_build_grid(
  train_data    = bl_filt$train_data,
  bl_projection = bl_proj,
  bl_model      = bl_mod,
  m             = 200L
)

print(bl_grid)


# ---- Step 6: Assemble result object -----------------------------------
result <- bl_assemble(
  bl_data          = bl_dat,
  bl_filter_result = bl_filt,
  bl_model         = bl_mod,
  bl_projection    = bl_proj,
  bl_grid          = bl_grid
)

print(result)

# Reference biplot — training data
plot_biplotEZ(result)

# Reference biplot — test data (four-colour confusion scheme)
test_pts <- bl_project_points(result$test_data, result)
plot_biplotEZ(result, points = test_pts)


# ===========================================================
# PHASE 3 — Local interpretation of a single target point
# ===========================================================

# ---- Step 7: Inspect test predictions to choose a target --------------
# Print a summary to identify a suitable class-1 observation (diabetic)
# whose counterfactual we want to explain.

pred_summary <- data.frame(
  row       = seq_len(nrow(result$test_data)),
  Glucose   = result$test_data$Glucose,
  BMI       = result$test_data$BMI,
  Age       = result$test_data$Age,
  pred_prob = round(
    sapply(seq_len(nrow(result$test_data)), function(i) {
      bl_select_target(result, target = i)$pred_prob
    }), 3
  ),
  true_class = result$test_data$class
)

print(pred_summary)

# Choose the target index — a class-1 observation with a probability
# well above 0.5 gives a clear and illustrative counterfactual.
tdp <- which(result$test_data$class == 1 &
               pred_summary$pred_prob >= 0.70)[1]

cat(sprintf("\nSelected target: row %d  |  pred = %.3f  |  true class = %d\n",
            tdp, pred_summary$pred_prob[tdp],
            result$test_data$class[tdp]))


# ---- Step 8: Select the target ----------------------------------------
# bl_select_target() wraps the observation into a bl_target object,
# projecting it into Z-space and scoring it through the model.
# Default data = result$test_data.

tgt <- bl_select_target(result, target = tdp)
print(tgt)

# Highlight the target on the main biplot
plot_biplotEZ(
  result,
  points       = test_pts,
  target_point = unlist(tgt$x_obs),
  target_label = tdp
)


# ---- Step 9: Set actionability constraints ----------------------------
# Specify which variable directions are clinically plausible for this
# patient. Omit set_filters() entirely if no constraints are needed.
#
# Valid constraint types:
#   "decrease"  — counterfactual value must be <= observed
#   "increase"  — counterfactual value must be >= observed
#   "fixed"     — variable cannot change
#   c(min, max) — counterfactual must lie within this absolute range

flt <- set_filters(
  tgt,
  Glucose = "decrease",   # reducing Glucose is the clinical lever
  Age     = "fixed"       # Age is not actionable
)

print(flt)


# ---- Step 10: Find local counterfactual via biplot rotation -----------
# For each eigenvector pair (1,2), (1,3), (2,3), (1,4), ...:
#   - Rotates the biplot so the target lies on the projection plane
#   - Builds a local prediction grid in the rotated Z-space
#   - Applies train_ranges (feasibility) + set_filters (actionability)
#   - NOTE: the Z-space convex hull polygon is NOT applied here —
#     rotation invalidates the original polygon coordinates
#   - Finds the nearest opposing-class boundary point
# The pair with the overall minimum Z-space distance is selected.

bl_local <- bl_find_boundary_local(
  bl_result   = result,
  bl_target   = tgt,
  set_filters = flt
)

print(bl_local)


# ---- Step 11: Local biplot plot ---------------------------------------
# Shows the rotated biplot for the best eigenvector pair only.
# Training points: TP=red, TN=blue, FP=purple, FN=orange.
# Yellow filled circle = target observation.
# Grey cross           = nearest boundary (counterfactual).
# Arrow                = direction of minimum change.

plot(bl_local)


# ---- Step 12: Shapley contribution plot --------------------------------
# Decomposes the path from target to counterfactual into per-variable
# contributions using the exact Shapley formula (p <= 14 variables).
#
# Each variable is labelled as:
#   Supports    (darkblue) — moving in this direction helps cross the boundary
#   Contradicts (grey)     — moving in this direction works against crossing
#
# Y-axis labels: "variable: observed_value -> required_change"
# X-axis: Shapley value = contribution to the probability shift

bl_shap <- bl_shapley(bl_local)
print(bl_shap)

plot(bl_shap)


# ---- Step 13: Sparse counterfactual -----------------------------------
# Automatically zeroes out Contradicts variables — they revert to their
# observed values. Only Supports variables take their counterfactual
# value. round_to = 0 rounds remaining changes to whole numbers.
#
# solution_valid = TRUE  means the sparse CF still crosses the boundary.
# solution_valid = FALSE means the rounded/zeroed CF is no longer enough;
#   try round_to = 1, remove round_to, or relax actionability filters.

bl_sparse <- bl_sparse_cf(bl_shap, round_to = 0)
print(bl_sparse)

# Local biplot with both the full CF and the sparse CF overlaid.
#   Grey cross  = full counterfactual (all variables changed)
#   Blue cross  = sparse CF (valid — crosses boundary)
#   Red cross   = sparse CF (invalid — does not cross boundary)
plot(bl_sparse)


# ===========================================================
# ALTERNATIVE: run without actionability constraints
# Useful as a baseline before adding constraints, or when
# bl_local$solution_found == FALSE to confirm a boundary exists.
# ===========================================================

# ---- Step 14 (optional): Unconstrained local search -------------------
bl_local_free <- bl_find_boundary_local(result, tgt)
print(bl_local_free)

plot(bl_local_free)

bl_shap_free   <- bl_shapley(bl_local_free)
bl_sparse_free <- bl_sparse_cf(bl_shap_free, round_to = 0)

plot(bl_shap_free)
plot(bl_sparse_free)


# ===========================================================
# EXTERNAL DATA POINT
# Analyse a new patient who does not appear in the train/test data.
# Pass a single-row data frame directly as the target.
# No class column is required.
# ===========================================================

# ---- Step 15: External patient ----------------------------------------
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

tgt_ext <- bl_select_target(result, target = new_patient)
print(tgt_ext)

# Highlight on main biplot
plot_biplotEZ(
  result,
  points       = test_pts,
  target_point = unlist(tgt_ext$x_obs),
  target_label = "new"
)

# Local search — no actionability constraints for this new patient
bl_local_ext <- bl_find_boundary_local(result, tgt_ext)
print(bl_local_ext)

plot(bl_local_ext)

bl_shap_ext   <- bl_shapley(bl_local_ext)
bl_sparse_ext <- bl_sparse_cf(bl_shap_ext)

plot(bl_shap_ext)
plot(bl_sparse_ext)
