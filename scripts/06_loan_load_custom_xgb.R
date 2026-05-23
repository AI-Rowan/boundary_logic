############################################################
# Loan default — Load a pre-saved custom XGBoost model
#
# Dataset : inst/extdata/loan_data.csv
# Model   : XGBoost booster loaded from disk via xgboost::xgb.load()
#           or readRDS(), then registered with bl_wrap_model()
#
# Use this script when the XGBoost model has already been fitted and
# saved externally (e.g. from script 05 or a separate training pipeline).
# The data preparation steps are identical to scripts 03 and 05 so that
# bl_filt matches the data the model was originally trained on.
#
# Run interactively: place cursor inside a {} block and press Ctrl+Enter
############################################################

# ---- Load the package (development workflow) ---------------------------
{
  rm(list = ls())
  devtools::load_all()
  library(xgboost)
}


# ===========================================================
# PHASE 1a — Reproduce the data preparation pipeline
# (must match the preparation used when the model was trained)
# ===========================================================

# ---- Step 1: Load and encode loan data ---------------------------------
{
  loan_raw <- read.csv(
    system.file("extdata", "loan_data.csv", package = "boundarylogic")
  )

  loan_encoded <- loan_raw

  gender_map   <- c(female = 0, male = 1)
  defaults_map <- c("No" = 0, "Yes" = 1)

  loan_encoded$person_gender                  <- gender_map[loan_encoded$person_gender]
  loan_encoded$previous_loan_defaults_on_file <- defaults_map[loan_encoded$previous_loan_defaults_on_file]

  loan_encoded[] <- lapply(loan_encoded, as.numeric)

  vars_to_remove <- c("loan_intent", "person_education", "person_home_ownership")
  loan_encoded   <- loan_encoded[, setdiff(names(loan_encoded), vars_to_remove),
                                  drop = FALSE]
}


# ---- Step 2: Domain filter — keep applicants without prior defaults ----
{
  loan_filtered <- loan_encoded[loan_encoded$previous_loan_defaults_on_file == 0, ]

  feature_cols <- setdiff(
    names(loan_filtered),
    c("loan_status", "previous_loan_defaults_on_file")
  )
}


# ---- Step 3: Prepare data and filter outliers --------------------------
# Use the same seed and hull_fraction as the original training run so
# bl_filt$train_data matches what the model was fitted on.
{
  bl_dat <- bl_prepare_data(
    data           = loan_filtered,
    class_col      = "loan_status",
    feature_cols   = feature_cols,
    train_fraction = 0.8,
    seed           = 121L
  )

  bl_filt   <- bl_filter_outliers(bl_dat, hull_fraction = 0.9)
  var_names <- bl_filt$var_names

  cat("Training rows:", nrow(bl_filt$train_data), "\n")
  cat("Features     :", paste(var_names, collapse = ", "), "\n")
}


# ===========================================================
# LOAD THE SAVED MODEL AND REGISTER WITH BOUNDARYLOGIC
# ===========================================================

# ---- Step 4: Load XGBoost model from disk ------------------------------
{
  # ---- Option A: XGBoost binary format ----------------------------------
  # Saved previously with: xgboost::xgb.save(xgb_fit, "path/to/model.xgb")
  # xgb_fit <- xgboost::xgb.load("path/to/model.xgb")

  # ---- Option B: R object format ----------------------------------------
  # Saved previously with: saveRDS(xgb_fit, "path/to/model.rds")
  # xgb_fit <- readRDS("path/to/model.rds")

  # For demonstration, fit a quick model here so the script is self-contained.
  # Replace these three lines with one of the xgb.load / readRDS calls above.
  dtrain  <- xgboost::xgb.DMatrix(
    data  = as.matrix(bl_filt$train_data[, var_names]),
    label = bl_filt$train_data[["class"]]
  )
  xgb_fit <- xgboost::xgb.train(
    params  = list(objective = "binary:logistic", eval_metric = "auc",
                   eta = 0.05, max_depth = 6, nthread = 1L),
    data    = dtrain,
    nrounds = 100L,
    verbose = 0L
  )
}


# ---- Step 5: Wrap the loaded model for boundarylogic -------------------
# model must be list(model = <xgb.Booster>, features = <character vector>).
# var_names must match the column order the booster was trained on.
{
  bl_mod <- bl_wrap_model(
    model      = list(model = xgb_fit, features = var_names),
    model_type = "XGB",
    var_names  = var_names,
    train_data = bl_filt$train_data   # supplies accuracy and Gini; omit if unavailable
  )
  print(bl_mod)
}


# ---- Step 6: Build result object and biplot ----------------------------
{
  bl_results <- bl_build_result(
    bl_data  = bl_filt,
    bl_model = bl_mod,
    method   = "CVA",
    title    = "Loan default -- loaded XGB model, CVA biplot",
    b_margin = 0.001
  )
  print(bl_results)

  plot_biplotEZ(bl_results, label_dir = "Hor")

  test_pts <- bl_project_points(bl_results$test_data, bl_results,
                                 filter_to_polygon = TRUE)
  plot_biplotEZ(bl_results, points = test_pts)
}


# ===========================================================
# PHASE 2 — Global interpretations
# ===========================================================

# ---- Step 7: Nearest boundary point ------------------------------------
{
  bl_bnd <- bl_find_boundary(bl_results)
  print(bl_bnd)
  plot_biplotEZ(bl_results, points = test_pts)
  # To inspect individual counterfactuals: bl_pick_point(bl_results, bl_boundary = bl_bnd)
}


# ---- Step 8: Distance-to-boundary plot ---------------------------------
{
  plot(bl_bnd)
  plot(bl_bnd, type = "boxplot")
}


# ---- Step 9: Surrogate model -------------------------------------------
{
  bl_surr <- bl_surrogate(bl_results)
  print(bl_surr)
  plot(bl_surr)
}


# ===========================================================
# PHASE 3 — Local interpretation
# ===========================================================

# ---- Step 10: Select a target observation ------------------------------
{
  pred_summary <- bl_predict(bl_results)
  print(pred_summary)

  tdp <- 1
  tgt <- bl_select_target(bl_results, target = tdp)
  print(tgt)

  plot_biplotEZ(
    bl_results,
    points       = test_pts,
    target_point = unlist(tgt$x_obs),
    target_label = tdp
  )
}


# ---- Step 11: Actionability constraints --------------------------------
{
  flt <- set_filters(
    tgt,
    loan_amnt     = "decrease",
    loan_int_rate = "fixed"
  )
  print(flt)
}


# ---- Step 12: Local counterfactual via biplot rotation -----------------
{
  bl_local <- bl_find_local_cf(
    bl_result   = bl_results,
    bl_target   = tgt,
    set_filters = flt
  )
  print(bl_local)
  plot(bl_local)
}


# ---- Step 13: Shapley contributions ------------------------------------
{
  bl_shapley_values <- bl_shapley(bl_local)
  print(bl_shapley_values)
  plot(bl_shapley_values)
}


# ---- Step 14: Sparse counterfactual ------------------------------------
{
  bl_sparse <- bl_find_sparse_cf(bl_shapley_values, round_to = NULL)
  print(bl_sparse)
  plot(bl_sparse)
}
