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
devtools::load_all()


# ---- Step 1: Prepare data ---------------------------------------------
# Outcome is already numeric 0/1, so target_class = NULL
pima_raw <- read.csv(
  file.path(dirname(rstudioapi::getActiveDocumentContext()$path),
            "pima diabetes.csv")
)

bl_dat <- bl_prepare_data(
  data           = pima_raw,
  class_col      = "Outcome",
  target_class   = NULL,          # already 0/1
  train_fraction = 0.8,
  seed           = 121L
)

print(bl_dat)


# ---- Step 2 (optional): Filter outliers --------------------------------
bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 0.9)

print(bl_filt)


# ---- Step 3: Fit custom GAM -------------------------------------------
# Spline terms for Glucose and BloodPressure; linear for all others.
# mgcv::gam() is used directly for full formula control.
# The model is wrapped via bl_wrap_model() so it integrates with the
# rest of the boundary logic pipeline.

gam_formula <- class ~ s(Glucose) + s(BloodPressure) +
  Pregnancies + SkinThickness + Insulin +
  BMI + DiabetesPedigreeFunction + Age

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
  #calc_ct_in_hull = TRUE
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
print(test_pts)

plot_biplotEZ(result, points = test_pts)
