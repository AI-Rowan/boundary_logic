## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----load-pkg-----------------------------------------------------------------
library(boundarylogic)
library(dplyr)

## ----step1--------------------------------------------------------------------
pima_csv <- system.file("extdata", "pima_diabetes.csv",
                        package = "boundarylogic")

pima_raw <- read.csv(pima_csv)

pima_raw <- pima_raw %>%
  filter(Insulin > 0, SkinThickness > 0)

bl_dat <- bl_prepare_data(
  data           = pima_raw,
  class_col      = "Outcome",
  target_class   = NULL,       # already 0/1
  train_fraction = 0.8,
  seed           = 121L
)

print(bl_dat)

## ----step1-check--------------------------------------------------------------
table(bl_dat$train_data$class)

## ----step2--------------------------------------------------------------------
bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 0.9)
print(bl_filt)

## ----step3--------------------------------------------------------------------
gam_formula <- class ~ s(Glucose) + s(BloodPressure) +
  s(Pregnancies) + s(SkinThickness) + Insulin +
  s(BMI) + DiabetesPedigreeFunction + Age

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

## ----steps4-6-----------------------------------------------------------------
bl_results <- bl_build_result(
  bl_data  = bl_filt,
  bl_model = bl_mod,
  method   = "CVA",
  title    = "Pima diabetes — GAM, CVA biplot",
  rounding = 2L
)

print(bl_results)

## ----steps4-6-summary---------------------------------------------------------
summary(bl_results)

## ----biplot-train-------------------------------------------------------------
plot_biplotEZ(bl_results)

## ----biplot-test--------------------------------------------------------------
test_pts <- bl_project_points(bl_results$test_data, bl_results)
plot_biplotEZ(bl_results, points = test_pts)

## ----step7-boundary-----------------------------------------------------------
bl_bnd <- bl_find_boundary(bl_results)
print(bl_bnd)

## ----step7-biplot-------------------------------------------------------------
plot_biplotEZ(bl_results, points = test_pts, boundary = bl_bnd)

## ----step8-jitter-------------------------------------------------------------
plot(bl_bnd)

## ----step8-boxplot------------------------------------------------------------
plot(bl_bnd, type = "boxplot")

## ----step9-surrogate----------------------------------------------------------
bl_surr <- bl_surrogate(bl_results)
print(bl_surr)
plot(bl_surr)

## ----step10-select------------------------------------------------------------
pred_summary <- bl_predict(bl_results)
print(pred_summary)

# Select the first predicted diabetic patient (class 1) for local explanation.
# With seed = 121L this is row 1 of test_data; the dynamic selection below
# ensures the vignette knits correctly regardless of the data split.
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

## ----step12-filters-----------------------------------------------------------
flt <- set_filters(
  tgt,
  Glucose     = "decrease",   # reducing Glucose is the clinical lever
  Age         = "fixed",      # Age is not actionable
  Pregnancies = "fixed"       # Pregnancies is not actionable
)

print(flt)

## ----step13-local-cf----------------------------------------------------------
bl_local <- bl_find_local_cf(
  bl_result   = bl_results,
  bl_target   = tgt,
  set_filters = flt
)

print(bl_local)

## ----step14-local-biplot------------------------------------------------------
plot(bl_local)

## ----step14-with-training-----------------------------------------------------
plot(bl_local, no_points = FALSE)

## ----step15-shapley-----------------------------------------------------------
bl_shapley_values <- bl_shapley(bl_local)
print(bl_shapley_values)

plot(bl_shapley_values)

## ----step16-sparse------------------------------------------------------------
bl_sparse <- bl_find_sparse_cf(bl_shapley_values, round_to = NULL)

# print() shows three predictions to verify the class flip:
#   Observed  — model score for the original patient
#   Full CF   — model score at the full counterfactual
#   Sparse CF — model score for the sparse (zeroed) counterfactual
print(bl_sparse)

# Grey cross   = full counterfactual (all variables changed)
# Green cross  = sparse CF (valid — crosses boundary)
# Yellow cross = sparse CF (invalid — does not cross boundary)
plot(bl_sparse)

## ----step17-unconstrained-----------------------------------------------------
bl_local_free   <- bl_find_local_cf(bl_results, tgt)
print(bl_local_free)

plot(bl_local_free)

bl_shap_free   <- bl_shapley(bl_local_free)
bl_sparse_free <- bl_find_sparse_cf(bl_shap_free, round_to = NULL)

plot(bl_shap_free)
plot(bl_sparse_free)

## ----step18-external----------------------------------------------------------
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

# Highlight on the main biplot
plot_biplotEZ(
  bl_results,
  points       = test_pts,
  target_point = unlist(tgt_ext$x_obs),
  target_label = "new"
)

## ----step18-local-------------------------------------------------------------
bl_local_ext <- bl_find_local_cf(bl_results, tgt_ext)
print(bl_local_ext)

plot(bl_local_ext)

## ----step18-shapley-----------------------------------------------------------
bl_shap_ext   <- bl_shapley(bl_local_ext)
bl_sparse_ext <- bl_find_sparse_cf(bl_shap_ext)

plot(bl_shap_ext)
plot(bl_sparse_ext)

