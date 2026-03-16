############################################################
# Pima diabetes — Phase 2 global interpretations
#            with an UNLABELED holdout sample
#
# Dataset : scripts/pima diabetes.csv
#           768 observations, 8 numeric features, binary outcome
# Model   : GAM with spline terms (mgcv::gam)
# Biplot  : CVA
#
# Purpose : Demonstrates the full Phase 2 workflow on a test set
#           where the true outcome labels have been withheld.
#           When labels are absent:
#             - biplot colours = predicted class (red/blue, 2 colours)
#             - distance plot  = predicted class (red/blue, 2 colours)
#             - surrogate      = accuracy vs model only (no vs-labels score)
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

# Exclude zero-imputed Insulin and SkinThickness (physiologically implausible)
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


# ---- Step 4: CVA projection -------------------------------------------
bl_proj <- bl_build_projection(
  train_data = bl_filt$train_data,
  var_names  = bl_filt$var_names,
  method     = "CVA",
  proj_dims  = c(1L, 2L),
  bl_model   = bl_mod,
  title      = "Pima diabetes — GAM, CVA biplot (unlabeled test)"
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
summary(result)


# ===========================================================
# CREATE UNLABELED HOLDOUT
# Remove the "class" column from the test set to simulate a
# real deployment scenario where ground truth is unknown.
# ===========================================================

test_unlabeled <- result$test_data[ , !(names(result$test_data) %in% "class")]

cat("\nUnlabeled holdout: ", nrow(test_unlabeled), "observations,",
    ncol(test_unlabeled), "features — no 'class' column\n")


# ===========================================================
# PHASE 2 — Global interpretations on the UNLABELED holdout
# ===========================================================

# ---- Step 7a: Training biplot (reference) -----------------------------
# Shown first so the analyst can see the full decision surface.
# Four-colour TP/TN/FP/FN confusion scheme (training labels are known).
plot_biplotEZ(result)


# ---- Step 7b: Unlabeled holdout projected onto the biplot -------------
# Because test_unlabeled has no "class" column, bl_project_points()
# falls back to predicted-class colouring: red = class 1, blue = class 0.
holdout_pts <- bl_project_points(test_unlabeled, result)

plot_biplotEZ(result, points = holdout_pts)


# ---- Step 8: Counterfactuals for the unlabeled holdout ----------------
# bl_find_boundary() accepts any data frame of raw feature values via the
# data argument. The hull polygon and train_ranges constraints always apply
# to ensure counterfactuals remain within the training space.
bl_bnd <- bl_find_boundary(result, data = test_unlabeled)

print(bl_bnd)

# X-space view: original values vs counterfactual values
head(bl_bnd$x_obs)
head(bl_bnd$B_x)

# Biplot with CF arrows overlaid
# Points are predicted-class coloured (red/blue); crosses mark counterfactuals.
plot_biplotEZ(result, points = holdout_pts, boundary = bl_bnd)

# Suppress arrows if only counterfactual positions are needed
# plot_biplotEZ(result, points = holdout_pts, boundary = bl_bnd,
#               show_arrows = FALSE)


# ---- Step 9: Distance-to-boundary plots -------------------------------
# Because labels are absent, both plot types use the two-colour scheme:
# red = predicted class 1, blue = predicted class 0.

# Jitter plot: one point per observation per variable
plot(bl_bnd)

# Box-whisker: distribution of distances by predicted class
plot(bl_bnd, type = "boxplot")

# Standalone robustness score (no plot)
bl_robustness(bl_bnd)


# ---- Step 10: Surrogate model on the unlabeled holdout ----------------
# Without true labels, only accuracy_vs_model is reported.
# Observations projected outside the training-data polygon are excluded
# (surrogate_pred = NA); accuracy is computed on in-hull observations only.
bl_surr <- bl_surrogate(result, data = test_unlabeled)

print(bl_surr)

# Surrogate biplot: hull-clipped background, ct_surrogate contours,
# points coloured by surrogate assignment (red = class 1, blue = class 0).
plot(bl_surr)


# ===========================================================
# COMPARISON: re-run Step 7b–10 WITH labels to see the
# difference in colour scheme
# ===========================================================

# ---- Step 7b (labeled): Confusion colours on the same holdout ---------
# result$test_data retains the true "class" column.
test_pts_labeled <- bl_project_points(result$test_data, result)
plot_biplotEZ(result, points = test_pts_labeled)


# ---- Steps 8-9 (labeled): Four-colour TP/TN/FP/FN scheme -------------
bl_bnd_labeled <- bl_find_boundary(result)   # defaults to result$test_data

plot_biplotEZ(result, points = test_pts_labeled, boundary = bl_bnd_labeled)
plot(bl_bnd_labeled)
plot(bl_bnd_labeled, type = "boxplot")


# ---- Step 10 (labeled): Surrogate with accuracy vs labels -------------
bl_surr_labeled <- bl_surrogate(result, data = result$test_data)
print(bl_surr_labeled)
plot(bl_surr_labeled)
