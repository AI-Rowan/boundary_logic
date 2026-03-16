############################################################
# Pima diabetes — custom GAM + CVA biplot
#
# Dataset : scripts/pima diabetes.csv
#           768 observations, 8 numeric features, binary outcome
# Model   : GAM fitted with mgcv::gam()
#             s(Glucose)       — cubic regression spline
#             s(BloodPressure) — cubic regression spline
#             all other variables — linear terms
# Biplot  : CVA (uses TP/TN/FP/FN confusion labels via bl_model)
#
# Run interactively: Ctrl+Enter line by line
############################################################

# ---- Load the package (development workflow) ---------------------------
rm()
devtools::load_all()


# ---- Step 1: Prepare data ---------------------------------------------
# Outcome is already numeric 0/1, so target_class = NULL
pima_raw <- read.csv(
  file.path(dirname(rstudioapi::getActiveDocumentContext()$path),
            "pima diabetes.csv")
)

library(dplyr)
#str(pima_raw)
#pima_raw <- pima_raw[pima_raw$Insulin >0 & pima_raw$SkinThickness>0,]
pima_raw <- pima_raw %>% filter(Insulin >0, SkinThickness>0)

bl_dat <- bl_prepare_data(
  data           = pima_raw,
  class_col      = "Outcome",
  target_class   = NULL,          # already 0/1
  train_fraction = 0.8,
  seed           = 121L
)

print(bl_dat)


# ---- Step 2 (optional): Filter outliers --------------------------------
bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 0.8)

print(bl_filt)


# ---- Step 3: Fit custom GAM -------------------------------------------
# Spline terms for Glucose and BloodPressure; linear for all others.
# mgcv::gam() is used directly for full formula control.
# The model is wrapped via bl_wrap_model() so it integrates with the
# rest of the boundary logic pipeline.

gam_formula <- class ~ s(Glucose) + s(BloodPressure) +
  s(Pregnancies) + SkinThickness + Insulin +
  s(BMI) + DiabetesPedigreeFunction + s(Age)

custom_gam <- mgcv::gam(
  formula = gam_formula,
  data    = bl_filt$train_data,
  family  = binomial(link = "logit")
)

summary(custom_gam)

# Wrap the fitted model so boundary logic functions can use it
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
# bl_model supplies confusion labels (TP/TN/FP/FN) automatically,
# satisfying CVA's requirement for >= 3 class levels.
bl_proj <- bl_build_projection(
  train_data = bl_filt$train_data,
  var_names  = bl_filt$var_names,
  method     = "CVA",
  proj_dims  = c(1L, 2L),
  bl_model   = bl_mod,
  title      = "Pima diabetes — custom GAM, CVA biplot"
)

print(bl_proj)


# ---- Step 5: Build prediction grid ------------------------------------
# calc_ct_in_hull = TRUE is required for CVA biplots
bl_grid <- bl_build_grid(
  train_data      = bl_filt$train_data,
  bl_projection   = bl_proj,
  bl_model        = bl_mod,
  m               = 200L
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
summary(result)


# ---- Step 7a: Plot CVA biplot — training data -------------------------
# Default: bl_project_points(result$train_data, result) used automatically
plot_biplotEZ(result)


# ---- Step 7b: Plot CVA biplot — test data -----------------------------
# Project test data into the biplot space and pass as the points argument.
# Confusion colours (TP/TN/FP/FN) are assigned because test_data has a
# "class" column.
test_pts <- bl_project_points(result$test_data, result)
#print(test_pts)

plot_biplotEZ(result, points = test_pts)


# ---- Step 8: Find counterfactuals for test data -----------------------
# bl_find_boundary() defaults to result$test_data.
# train_ranges are always applied to ensure CFs stay within training space.
bl_bnd <- bl_find_boundary(result)

print(bl_bnd)

# Quick look at the counterfactual X-space values vs the originals
head(bl_bnd$B_x)
head(bl_bnd$x_obs)

plot_biplotEZ(result, points = test_pts, boundary = bl_bnd) #, show_arrows = FALSE)

# ---- Step 9: Distance-to-boundary plots ------------------------------
# Jitter plot (default): one point per observation per variable
plot(bl_bnd)

# Box-whisker plot: distribution of distances by class per variable
plot(bl_bnd, type = "boxplot")

# Standalone robustness score (no plot)
bl_robustness(bl_bnd)


# ---- Step 10: Surrogate model -----------------------------------------
# bl_surrogate() uses ct_surrogate (hull-clipped contours pre-computed in
# bl_build_grid()) to assign each training observation to the nearest
# contour region. Accuracy vs the fitted model and vs true labels is
# reported.

bl_surr <- bl_surrogate(result, data = result$test_data)
print(bl_surr)

# Biplot with surrogate region colouring and accuracy output
plot(bl_surr)


############################################################
# Phase 3: Local interpretation of a single target point
############################################################

# ---- Step 11: Select a target observation -----------------------------
# Choose a test data point predicted as diabetic (class 1) to explain
# what changes would move it across the decision boundary.
# Inspect predictions to find a suitable candidate:
head(data.frame(
  row      = seq_len(nrow(result$test_data)),
  Glucose  = result$test_data$Glucose,
  pred     = round(bl_bnd$pred_obs, 3),
  class    = result$test_data$class
))

# Select the target by row index (default data = result$test_data).
# Change the index to any class-1 observation of interest.
tdp <- which(result$test_data$class == 1)[1]

tgt <- bl_select_target(result, target = tdp)
print(tgt)

# Highlight the target on the main biplot before proceeding
plot_biplotEZ(
  result,
  points       = test_pts,
  target_point = unlist(tgt$x_obs),
  target_label = tdp
)


# ---- Step 12: Set actionability constraints (optional) ----------------
# Constrain which variable directions are clinically plausible.
# Omit set_filters() entirely to search with no constraints beyond
# the training data variable ranges.
flt <- set_filters(
  tgt,
  Glucose  = "decrease",   # Glucose can only be reduced
  Age      = "fixed"       # Age cannot change
)

print(flt)


# ---- Step 13: Find local counterfactual with biplot rotation ----------
# Tests up to 10 eigenvector pair combinations (1,2),(1,3),(2,3),...
# For each pair: rotates the biplot to align the target, builds a local
# prediction grid, applies train_ranges + actionability constraints,
# finds the nearest opposing-class boundary point.
# Hull polygon is NOT applied — rotation invalidates the Z-space hull.

bl_local <- bl_find_boundary_local(
  bl_result   = result,
  bl_target   = tgt,
  set_filters = flt
)

print(bl_local)


# ---- Step 14: Local biplot plot ---------------------------------------
# Shows the rotated biplot for the best eigenvector pair.
# Training points: TP=red, TN=blue, FP=purple, FN=orange.
# Yellow circle = target; grey cross = counterfactual; arrow shows direction.
plot(bl_local)


# ---- Step 15: Shapley contribution plot --------------------------------
# Decomposes the move from target to counterfactual into per-variable
# contributions. Variables are classified as:
#   Supports   (darkblue): moving in this direction helps reach the boundary
#   Contradicts (grey)   : moving in this direction works against it
# Y-axis labels show observed value -> required change.

bl_shap <- bl_shapley(bl_local)
print(bl_shap)

plot(bl_shap)


# ---- Step 16: Sparse counterfactual -----------------------------------
# Automatically zeroes out Contradicts variables (no change from observed).
# round_to = 0 rounds remaining changes to whole numbers for a sparser,
# more interpretable solution.
# solution_valid = TRUE means the sparse CF still crosses the boundary.

bl_sparse <- bl_sparse_cf(bl_shap, round_to = 0)
print(bl_sparse)

# Biplot showing both full CF (grey cross) and sparse CF:
#   blue cross = valid (crosses boundary)
#   red cross  = invalid (does not cross boundary — relax round_to or filters)
plot(bl_sparse)


# ---- Step 17: Re-run without actionability constraints ----------------
# If no solution was found (bl_local$solution_found == FALSE), remove
# the filters to confirm a boundary exists before tightening constraints.

# bl_local_unconstrained <- bl_find_boundary_local(result, tgt)
# plot(bl_local_unconstrained)


# ---- Step 18: External data point (not in train/test) -----------------
# Pass a single-row data frame directly as the target.
# All feature columns must be present; no class column required.

# new_patient <- data.frame(
#   Pregnancies            = 2,
#   Glucose                = 148,
#   BloodPressure          = 72,
#   SkinThickness          = 35,
#   Insulin                = 100,
#   BMI                    = 33.6,
#   DiabetesPedigreeFunction = 0.627,
#   Age                    = 50
# )
# tgt_ext   <- bl_select_target(result, target = new_patient)
# local_ext <- bl_find_boundary_local(result, tgt_ext)
# plot(local_ext)
