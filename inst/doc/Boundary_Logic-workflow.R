## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 6,
  fig.height = 5
)

## ----step1--------------------------------------------------------------------
library(boundarylogic)

bl_dat <- bl_prepare_data(
  data           = iris,
  class_col      = "Species",
  target_class   = "versicolor",
  train_fraction = 0.8,
  seed           = 121L
)

print(bl_dat)

## ----step1-check--------------------------------------------------------------
# Class balance in the training set
table(bl_dat$train_data$class)

# First few rows
head(bl_dat$train_data)

## ----step2--------------------------------------------------------------------
bl_filt <- bl_filter_outliers(
  bl_data       = bl_dat,
  hull_fraction = 0.91,
  verbose       = TRUE
)

print(bl_filt)

## ----step3--------------------------------------------------------------------
bl_mod <- bl_fit_model(
  train_data = bl_filt$train_data,
  var_names  = bl_filt$var_names,
  model_type = "GLM",   # logistic regression
  cutoff     = 0.5
)

print(bl_mod)

## ----steps4-6-----------------------------------------------------------------
bl_results <- bl_build_result(
  bl_data  = bl_filt,     # pass bl_dat here if Step 2 was skipped
  bl_model = bl_mod,
  method   = "PCA",       # or "CVA" — see note below
  m        = 100L,        # grid resolution; use 200L for final analysis
  rounding = 3L           # contour band width: cutoff ± 0.001
)

print(bl_results)

## ----steps4-6-summary---------------------------------------------------------
summary(bl_results)

## ----step7--------------------------------------------------------------------
plot_biplotEZ(bl_results)

## ----step7-test---------------------------------------------------------------
test_pts <- bl_project_points(bl_results$test_data, bl_results)
plot_biplotEZ(bl_results, points = test_pts)

## ----phase2-boundary----------------------------------------------------------
# Find boundary for all test observations
bl_bnd <- bl_find_boundary(bl_results, data = bl_results$test_data)
print(bl_bnd)

# Overlay counterfactuals on the biplot
plot_biplotEZ(bl_results, points = test_pts, boundary = bl_bnd)

## ----phase2-boundary-plot-----------------------------------------------------
plot(bl_bnd)               # jitter plot (default)
plot(bl_bnd, type = "box") # boxplot

## ----phase2-robustness--------------------------------------------------------
bl_rob <- bl_robustness(bl_bnd)
print(bl_rob)

## ----phase2-surrogate---------------------------------------------------------
bl_sur <- bl_surrogate(bl_results)
print(bl_sur)

# Visualise the surrogate regions
plot(bl_sur)

## ----phase3-select------------------------------------------------------------
pred_summary <- bl_predict(bl_results)
print(pred_summary)

# Inspect pred_summary to choose a row. For a clear counterfactual,
# pick a class-1 observation with a predicted probability well above 0.5.
# Here we select the first predicted class-1 row automatically.
tdp <- which(pred_summary$pred_class == 1)[1]
tgt <- bl_select_target(bl_results, target = tdp)
print(tgt)

# Highlight the target on the main biplot
plot_biplotEZ(
  bl_results,
  points       = test_pts,
  target_point = unlist(tgt$x_obs),
  target_label = tdp
)

## ----phase3-external----------------------------------------------------------
new_obs <- data.frame(
  Sepal.Length = 5.8, Sepal.Width = 2.7,
  Petal.Length = 4.1, Petal.Width = 1.0
)
tgt_ext <- bl_select_target(bl_results, target = new_obs)
print(tgt_ext)

## ----phase3-filters-----------------------------------------------------------
# Filters are illustrated here but not applied in the local search below,
# because the constraints that are clinically meaningful depend on the
# specific dataset and target. In interactive use, pass flt to
# bl_find_local_cf() via the set_filters argument.
flt <- set_filters(
  tgt,
  Petal.Length = "decrease",   # only consider reductions
  Sepal.Length = "fixed"       # treat as non-actionable
)

print(flt)

## ----phase3-local-cf----------------------------------------------------------
# Unconstrained search — finds the nearest boundary point across all
# eigenvector pairs. Pass set_filters = flt to add actionability
# constraints in interactive use.
bl_local <- bl_find_local_cf(
  bl_result = bl_results,
  bl_target = tgt
)

print(bl_local)

# Local biplot: target (filled circle) + counterfactual (grey cross)
# Default: training points hidden. Pass no_points = FALSE to show them.
plot(bl_local)
plot(bl_local, no_points = FALSE)   # also show training data

## ----phase3-shapley-----------------------------------------------------------
bl_shapley_values <- bl_shapley(bl_local)
print(bl_shapley_values)

plot(bl_shapley_values)

## ----phase3-sparse------------------------------------------------------------
bl_sparse <- bl_find_sparse_cf(bl_shapley_values, round_to = NULL)

# print() shows three prediction values to verify the class flip:
#   Observed  — model score for the original data point
#   Full CF   — model score at the full counterfactual
#   Sparse CF — model score for the sparse (zeroed) counterfactual
print(bl_sparse)

# Local biplot with both CFs overlaid:
#   Grey cross   = full counterfactual (all variables changed)
#   Green cross  = sparse CF (valid — crosses boundary)
#   Yellow cross = sparse CF (invalid — does not cross boundary)
# To also show training data: plot(bl_sparse, no_points = FALSE)
plot(bl_sparse)

