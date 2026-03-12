############################################################
# Pima diabetes — external test data workflow
#
# Scenario : The test data was NOT part of the data passed
#            to bl_prepare_data(). True labels are unknown.
#            Points are coloured by predicted class only:
#              Red  — predicted class 1 (diabetic)
#              Blue — predicted class 0 (non-diabetic)
#
# The ONLY workflow difference from the standard script is:
#   1. Split BEFORE bl_prepare_data(); hold out the test rows.
#   2. Call bl_prepare_data() with train_fraction = 1 so the
#      entire labelled set becomes training data.
#   3. Drop the Outcome column from the held-out test set so
#      bl_project_points() sees no "class" column and falls
#      back to prediction-only colouring automatically.
#
# No function changes are needed.
############################################################

# ---- Load the package (development workflow) ---------------------------
rm()
devtools::load_all()



# ---- Step 1: Load and split data BEFORE the pipeline ------------------
# We hold out 20 % as "external" test data whose labels we will not use.

pima_raw <- read.csv(
  file.path(dirname(rstudioapi::getActiveDocumentContext()$path),
            "pima diabetes.csv")
)

set.seed(42L)
n          <- nrow(pima_raw)
test_idx   <- sample(n, size = floor(0.2 * n), replace = FALSE)

train_raw  <- pima_raw[-test_idx, ]   # 80 %  — labelled, fed to pipeline
test_raw   <- pima_raw[ test_idx, ]   # 20 %  — external, labels withheld


# ---- Step 2: Prepare data — training rows only -------------------------
# train_fraction = 1 means bl_prepare_data() allocates all rows to
# training; result$test_data will be empty (that is intentional here).

bl_dat <- bl_prepare_data(
  data           = train_raw,
  class_col      = "Outcome",
  target_class   = NULL,      # already 0/1
  train_fraction = 0.8,         # <-- all labelled data becomes training
  seed           = 121L
)

print(bl_dat)


# ---- Step 3 (optional): Filter outliers --------------------------------
bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 0.8)

print(bl_filt)


# ---- Step 4: Fit custom GAM -------------------------------------------
gam_formula <- class ~ s(Glucose) + s(BloodPressure) +
  Pregnancies + SkinThickness + Insulin +
  BMI + DiabetesPedigreeFunction + Age

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


# ---- Step 5: CVA projection -------------------------------------------
bl_proj <- bl_build_projection(
  train_data = bl_filt$train_data,
  var_names  = bl_filt$var_names,
  method     = "CVA",
  proj_dims  = c(1L, 2L),
  bl_model   = bl_mod,
  title      = "Pima diabetes — external test, prediction colours"
)

print(bl_proj)


# ---- Step 6: Build prediction grid ------------------------------------
bl_grid <- bl_build_grid(
  train_data    = bl_filt$train_data,
  bl_projection = bl_proj,
  bl_model      = bl_mod,
  m             = 200L
)

print(bl_grid)


# ---- Step 7: Assemble result object -----------------------------------
result <- bl_assemble(
  bl_data          = bl_dat,
  bl_filter_result = bl_filt,
  bl_model         = bl_mod,
  bl_projection    = bl_proj,
  bl_grid          = bl_grid
)

print(result)
summary(result)


# ---- Step 8a: Reference biplot — training data (confusion colours) ----
# TP = red, TN = blue, FP = purple, FN = orange
plot_biplotEZ(result)


# ---- Step 8b: External test data — prediction colours only ------------
# Key step: drop the Outcome column so bl_project_points() has no "class"
# column and automatically colours by predicted class (red / blue).

external_test <- test_raw[, bl_filt$var_names]   # features only, no Outcome

ext_pts <- bl_project_points(external_test, result)
print(ext_pts)

# Red  = predicted diabetic (class 1)
# Blue = predicted non-diabetic (class 0)
plot_biplotEZ(result, points = ext_pts)


# ---- (Optional) Step 8c: Verify — what if we DID know the labels? -----
# Provide the true labels for a sanity check on prediction quality.
# This step is purely diagnostic; in the real scenario you would skip it.

external_test_labelled           <- test_raw[, c(bl_filt$var_names, "Outcome")]
names(external_test_labelled)[ncol(external_test_labelled)] <- "class"

labelled_pts <- bl_project_points(external_test_labelled, result)
print(labelled_pts)

# Now coloured as TP/TN/FP/FN so you can see where the model errs
plot_biplotEZ(result, points = labelled_pts)
