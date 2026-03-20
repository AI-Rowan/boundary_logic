############################################################
# Pima diabetes — SHAP analysis (global + local)
#
# Dataset : scripts/pima diabetes.csv  (same as 00_pima_phase3_local.R)
# Model   : GAM with spline terms (mgcv::gam)
# Packages: fastshap, shapviz
#
# Purpose : Model-agnostic SHAP interpretation of the same GAM fitted
#           in the boundary logic workflow. Provides a reference for
#           comparing SHAP-based feature attribution against the
#           boundary logic Shapley decomposition (Phase 3).
#
#           Global (Steps 1-5):
#             1. Data prep + model fit (mirrors Phase 1 exactly)
#             2. Compute SHAP values for all training observations
#             3. Importance bar chart — mean |SHAP| per variable
#             4. Beeswarm summary — SHAP distribution vs feature value
#             5. Dependence plots for top variables
#
#           Local (Steps 6-8):
#             6. Select a target observation (same criteria as Phase 3)
#             7. Compute local SHAP values for that observation
#             8. Waterfall and force plots
#
# Run interactively: Ctrl+Enter line by line
############################################################

rm(list = ls())
devtools::load_all()   # boundary logic package — for data prep only
library(dplyr)
library(mgcv)
library(fastshap)
library(shapviz)


# ===========================================================
# STEP 1: Data preparation and model fit
# Mirrors 00_pima_phase3_local.R Steps 1-3 exactly.
# ===========================================================

pima_raw <- read.csv(
  file.path(dirname(rstudioapi::getActiveDocumentContext()$path),
            "pima diabetes.csv")
)

pima_raw <- pima_raw %>%
  filter(Insulin > 0, SkinThickness > 0) #%>%  mutate(Age = 90 - Age)

# Use bl_prepare_data() and bl_filter_outliers() so the training split
# and hull filter are identical to the boundary logic workflow.
bl_dat <- bl_prepare_data(
  data           = pima_raw,
  class_col      = "Outcome",
  target_class   = NULL,
  train_fraction = 0.8,
  seed           = 121L
)

bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 0.9)

var_names  <- bl_filt$var_names
train_data <- bl_filt$train_data
test_data  <- bl_filt$test_data

X_train <- train_data[, var_names]
X_test  <- test_data[,  var_names]

# Fit the same GAM as in 00_pima_phase3_local.R
gam_formula <- class ~ s(Glucose) + s(BloodPressure) +
  s(Pregnancies) + s(SkinThickness) + Insulin +
  s(BMI) + DiabetesPedigreeFunction + Age

fit_gam <- mgcv::gam(
  formula = gam_formula,
  data    = train_data,
  family  = binomial(link = "logit")
)

summary(fit_gam)

# Prediction wrapper — fastshap requires f(object, newdata) -> numeric vector
pred_fn <- function(object, newdata) {
  as.numeric(mgcv::predict.gam(object, newdata = newdata, type = "response"))
}


# ===========================================================
# GLOBAL SHAP (Steps 2-5)
# Explains the full training set. Uses Monte Carlo permutation
# sampling (nsim = 100) to approximate marginal Shapley values.
# Increase nsim for more stable estimates at the cost of speed.
# ===========================================================

# ---- Step 2: Compute SHAP values for all training observations --------
# X provides the background distribution for marginal expectations.
# fastshap approximates E[f(X) | X_S = x_S] by randomly imputing
# the complement features from X for each permutation.

set.seed(42L)
shap_global <- fastshap::explain(
  object       = fit_gam,
  X            = X_train,
  pred_wrapper = pred_fn,
  nsim         = 100L
)

# Build a shapviz object — attaches feature values for dependence plots
sv_global <- shapviz::shapviz(shap_global, X = X_train)

print(sv_global)


# ---- Step 3: Importance bar chart — mean |SHAP| per variable ----------
# Variables are ranked by their mean absolute SHAP value across all
# training observations. This is the standard SHAP feature importance.

shapviz::sv_importance(sv_global, kind = "bar") +
  ggplot2::labs(title = "Global SHAP importance — Pima GAM (training data)")


# ---- Step 4: Beeswarm summary plot ------------------------------------
# Each point is one observation. The x-axis shows the SHAP value
# (contribution to predicted log-odds) for that variable. Points are
# coloured by the feature value (blue = low, red = high), showing
# whether high or low values of each variable push the prediction up
# or down.

shapviz::sv_importance(sv_global, kind = "beeswarm") +
  ggplot2::labs(title = "Global SHAP beeswarm — Pima GAM (training data)")


# ---- Step 5: Dependence plots for top variables -----------------------
# Shows how the SHAP value for a variable changes with its feature
# value. The colour is the SHAP value of the variable with the
# strongest interaction (chosen automatically).
# Replace the variable names below with those that rank highest in Step 3.

shapviz::sv_dependence(sv_global, v = "Glucose") +
  ggplot2::labs(title = "SHAP dependence — Glucose")

shapviz::sv_dependence(sv_global, v = "BMI") +
  ggplot2::labs(title = "SHAP dependence — BMI")

shapviz::sv_dependence(sv_global, v = "Age") +
  ggplot2::labs(title = "SHAP dependence — Age")


# ===========================================================
# LOCAL SHAP (Steps 6-8)
# Explains a single observation — same target selection as Phase 3
# (class-1, predicted probability >= 0.70 on the test set).
# ===========================================================

# ---- Step 6: Select target observation --------------------------------
pred_prob_test <- pred_fn(fit_gam, X_test)

tdp <- which(test_data$class == 1 & pred_prob_test >= 0.70)[1]

cat(sprintf("\nSelected target: row %d  |  pred = %.3f  |  true class = %d\n",
            tdp, pred_prob_test[tdp], test_data$class[tdp]))

cat("\nFeature values:\n")
print(X_test[tdp, ])


# ---- Step 7: Compute local SHAP values --------------------------------
# newdata = the single observation to explain.
# X       = training data as the background distribution (same as global).
# Use more permutations (nsim = 500) for a stable single-observation
# estimate.

set.seed(42L)
shap_local <- fastshap::explain(
  object       = fit_gam,
  X            = X_train,
  pred_wrapper = pred_fn,
  nsim         = 500L,
  newdata      = X_test[tdp, , drop = FALSE]
)

sv_local <- shapviz::shapviz(shap_local, X = X_test[tdp, , drop = FALSE])


# ---- Step 8: Local SHAP plots -----------------------------------------
# Waterfall: shows each variable's contribution as a step from the
# baseline (mean training prediction) to the final predicted value.
# Bars to the right push the prediction toward class 1 (diabetic);
# bars to the left push toward class 0.

shapviz::sv_waterfall(sv_local) +
  ggplot2::labs(title = sprintf(
    "Local SHAP waterfall — observation %d  (pred = %.3f, true class = %d)",
    tdp, pred_prob_test[tdp], test_data$class[tdp]
  ))

# Force plot: compact alternative to the waterfall showing the same
# push/pull contributions as overlapping bars on a single axis.

shapviz::sv_force(sv_local) +
  ggplot2::labs(title = sprintf(
    "Local SHAP force plot — observation %d", tdp
  ))
