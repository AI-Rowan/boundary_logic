############################################################
# Loan default — Iterative Phases 1, 2 & 3
# Custom XGBoost variant
#
# Dataset : inst/extdata/loan_data.csv
#           45,000 observations, binary outcome
# Model   : xgboost::xgb.train() with custom hyperparameters,
#           wrapped via bl_wrap_model()
# Biplot  : CVA
#
# Key difference from script 03:
#   Steps 4-6 and 4b-6b fit XGBoost directly with xgboost::xgb.train()
#   using a train/validation watchlist and early stopping, then register
#   the fitted booster via bl_wrap_model() instead of bl_fit_model().
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
# PHASE 1a — Load and explore
# ===========================================================

# ---- Step 1: Load and integer-encode categorical variables -------------
{
  loan_raw <- read.csv(
    system.file("extdata", "loan_data.csv", package = "boundarylogic")
  )

  loan_encoded <- loan_raw

  gender_map     <- c(female = 0, male = 1)
  education_map  <- c("High School" = 0, "Associate" = 1,
                      "Bachelor" = 2, "Master" = 3, "Doctorate" = 4)
  ownership_map  <- c("MORTGAGE" = 0, "OTHER" = 1, "OWN" = 2, "RENT" = 3)
  intent_map     <- c("DEBTCONSOLIDATION" = 0, "EDUCATION" = 1,
                      "HOMEIMPROVEMENT" = 2, "MEDICAL" = 3,
                      "PERSONAL" = 4, "VENTURE" = 5)
  defaults_map   <- c("No" = 0, "Yes" = 1)

  loan_encoded$person_gender                  <- gender_map[loan_encoded$person_gender]
  loan_encoded$person_education               <- education_map[loan_encoded$person_education]
  loan_encoded$person_home_ownership          <- ownership_map[loan_encoded$person_home_ownership]
  loan_encoded$loan_intent                    <- intent_map[loan_encoded$loan_intent]
  loan_encoded$previous_loan_defaults_on_file <- defaults_map[loan_encoded$previous_loan_defaults_on_file]

  loan_encoded[] <- lapply(loan_encoded, as.numeric)

  vars_to_remove <- c("loan_intent", "person_education", "person_home_ownership")
  loan_encoded <- loan_encoded[, setdiff(names(loan_encoded), vars_to_remove), drop = FALSE]

  str(loan_encoded)
}


# ---- Step 2: Exploratory biplot — full data, no model -----------------
{
  bl_dat_exp  <- bl_prepare_data(
    data           = loan_encoded,
    class_col      = "loan_status",
    train_fraction = 1,
    seed           = 121L
  )

  bl_filt_exp <- bl_filter_outliers(bl_dat_exp, hull_fraction = 1)

  bl_proj <- bl_build_result(
    bl_data = bl_filt_exp,
    method  = "PCA",
    title   = "Loan default — exploratory PCA biplot (all data)"
  )

  plot_biplotEZ(bl_proj, label_dir = "Hor")
}


# ---- Step 3: Domain filter — keep applicants without prior defaults ----
{
  loan_filtered <- loan_encoded[loan_encoded$previous_loan_defaults_on_file == 0, ]

  feature_cols <- setdiff(
    names(loan_filtered),
    c("loan_status", "previous_loan_defaults_on_file")
  )

  cat("Retained rows:", nrow(loan_filtered), "\n")
  cat("Features:", paste(feature_cols, collapse = ", "), "\n")
}


# ---- Steps 4-6: Prepare data, fit custom XGBoost, build biplot --------
{
  bl_dat  <- bl_prepare_data(
    data           = loan_filtered,
    class_col      = "loan_status",
    feature_cols   = feature_cols,
    train_fraction = 0.8,
    seed           = 121L
  )
  print(bl_dat)

  bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 0.9)

  # ---- Fit custom XGBoost -----------------------------------------------
  # Use 90 % of the filtered training data to train and 10 % as a validation
  # watchlist for early stopping.  Set seed before the split so results
  # are reproducible.
  train_df  <- bl_filt$train_data
  var_names <- bl_filt$var_names

  set.seed(42L)
  val_idx   <- sample(nrow(train_df), size = floor(0.1 * nrow(train_df)))
  xgb_train <- train_df[-val_idx, ]
  xgb_val   <- train_df[ val_idx, ]

  dtrain <- xgboost::xgb.DMatrix(
    data  = as.matrix(xgb_train[, var_names]),
    label = xgb_train[["class"]]
  )
  dval <- xgboost::xgb.DMatrix(
    data  = as.matrix(xgb_val[, var_names]),
    label = xgb_val[["class"]]
  )

  xgb_params <- list(
    objective         = "binary:logistic",
    eval_metric       = "auc",
    eta               = 0.05,        # learning rate
    max_depth         = 6,           # tree depth
    subsample         = 0.8,         # row subsampling per tree
    colsample_bytree  = 0.8,         # column subsampling per tree
    min_child_weight  = 5,           # minimum leaf weight (regularises small splits)
    gamma             = 0.1,         # minimum loss reduction to split
    nthread           = 1L           # single thread for reproducibility
  )

  xgb_fit <- xgboost::xgb.train(
    params                = xgb_params,
    data                  = dtrain,
    nrounds               = 2000L,
    watchlist             = list(train = dtrain, val = dval),
    early_stopping_rounds = 50L,
    verbose               = 1L,
    print_every_n         = 100L
  )

  cat("\nBest iteration:", xgb_fit$best_iteration, "\n")
  cat("Best val AUC:  ", xgb_fit$best_score, "\n")

  # ---- Wrap for boundarylogic -------------------------------------------
  # bl_wrap_model() with model_type = "XGB" expects:
  #   model = list(model = <xgb.Booster>, features = <character vector>)
  bl_mod <- bl_wrap_model(
    model      = list(model = xgb_fit, features = var_names),
    model_type = "XGB",
    var_names  = var_names,
    train_data = train_df
  )
  print(bl_mod)

  bl_results <- bl_build_result(
    bl_data  = bl_filt,
    bl_model = bl_mod,
    method   = "CVA",
    title    = "Loan default — custom XGB, CVA biplot",
    rounding = 3L
  )
  print(bl_results)

  plot_biplotEZ(bl_results, label_dir = "Hor")

  test_pts <- bl_project_points(bl_results$test_data, bl_results,
                                 filter_to_polygon = TRUE)
  plot_biplotEZ(bl_results, points = test_pts)
}


# ===========================================================
# PHASE 2 (first pass) — Global interpretations
# ===========================================================

# ---- Step 7: Find nearest boundary point for each observation ----------
{
  bl_bnd <- bl_find_boundary(bl_results)
  print(bl_bnd)
  hist(bl_bnd$B_pred, main = "Predicted probability at counterfactual", xlab = "p")
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
# STEP 10 — Prune least-important variables
# ===========================================================
{
  rob     <- bl_robustness(bl_bnd)
  var_imp <- sort(rob$sum_of_distance)
  print(round(var_imp, 2))

  # Update this list after reviewing var_imp output above
  vars_to_drop    <- c("person_gender", "person_emp_exp",
                       "cb_person_cred_hist_length", "person_income")
  feature_cols_v2 <- setdiff(bl_filt$var_names, vars_to_drop)
  cat("Retained features:", paste(feature_cols_v2, collapse = ", "), "\n")
}


# ===========================================================
# PHASE 1 (refit) — Custom XGBoost on reduced feature set
# ===========================================================

# ---- Steps 4b-6b: Rebuild pipeline with pruned features ----------------
{
  bl_dat_v2  <- bl_prepare_data(
    data           = loan_filtered,
    class_col      = "loan_status",
    feature_cols   = feature_cols_v2,
    train_fraction = 0.8,
    seed           = 121L
  )

  bl_filt_v2 <- bl_filter_outliers(bl_dat_v2, hull_fraction = 0.9)

  # ---- Fit custom XGBoost on reduced features ---------------------------
  train_df_v2  <- bl_filt_v2$train_data
  var_names_v2 <- bl_filt_v2$var_names

  set.seed(42L)
  val_idx_v2   <- sample(nrow(train_df_v2), size = floor(0.1 * nrow(train_df_v2)))
  xgb_train_v2 <- train_df_v2[-val_idx_v2, ]
  xgb_val_v2   <- train_df_v2[ val_idx_v2, ]

  dtrain_v2 <- xgboost::xgb.DMatrix(
    data  = as.matrix(xgb_train_v2[, var_names_v2]),
    label = xgb_train_v2[["class"]]
  )
  dval_v2 <- xgboost::xgb.DMatrix(
    data  = as.matrix(xgb_val_v2[, var_names_v2]),
    label = xgb_val_v2[["class"]]
  )

  xgb_fit_v2 <- xgboost::xgb.train(
    params                = xgb_params,   # same hyperparameters as Phase 1
    data                  = dtrain_v2,
    nrounds               = 2000L,
    watchlist             = list(train = dtrain_v2, val = dval_v2),
    early_stopping_rounds = 50L,
    verbose               = 1L,
    print_every_n         = 100L
  )

  cat("\nBest iteration:", xgb_fit_v2$best_iteration, "\n")
  cat("Best val AUC:  ", xgb_fit_v2$best_score, "\n")

  bl_mod_v2 <- bl_wrap_model(
    model      = list(model = xgb_fit_v2, features = var_names_v2),
    model_type = "XGB",
    var_names  = var_names_v2,
    train_data = train_df_v2
  )
  print(bl_mod_v2)

  bl_results_v2 <- bl_build_result(
    bl_data  = bl_filt_v2,
    bl_model = bl_mod_v2,
    method   = "CVA",
    title    = "Loan default — custom XGB, reduced features, CVA biplot",
    rounding = 2L
  )
  print(bl_results_v2)

  plot_biplotEZ(bl_results_v2)
  test_pts_v2 <- bl_project_points(bl_results_v2$test_data, bl_results_v2)
  plot_biplotEZ(bl_results_v2, points = test_pts_v2)
}


# ===========================================================
# PHASE 2 (second pass) — Global interpretations, reduced model
# ===========================================================

# ---- Step 7b: Boundary on reduced model --------------------------------
{
  bl_bnd_v2 <- bl_find_boundary(bl_results_v2)
  print(bl_bnd_v2)
  plot_biplotEZ(bl_results_v2, points = test_pts_v2)
  # To inspect individual counterfactuals: bl_pick_point(bl_results_v2, bl_boundary = bl_bnd_v2)
}


# ---- Step 8b: Distance-to-boundary plot --------------------------------
{
  plot(bl_bnd_v2)
  plot(bl_bnd_v2, type = "boxplot")
}


# ---- Step 9b: Surrogate model ------------------------------------------
{
  bl_surr_v2 <- bl_surrogate(bl_results_v2)
  print(bl_surr_v2)
  plot(bl_surr_v2)
}


# ===========================================================
# PHASE 3 — Local interpretation (reduced model)
# ===========================================================

# ---- Step 11: Inspect predictions and select target -------------------
{
  pred_summary <- bl_predict(bl_results_v2)
  print(pred_summary)

  tdp <- 1
  tgt <- bl_select_target(bl_results_v2, target = tdp)
  print(tgt)

  plot_biplotEZ(
    bl_results_v2,
    points       = test_pts_v2,
    target_point = unlist(tgt$x_obs),
    target_label = tdp
  )
}


# ---- Step 12: Set actionability constraints ----------------------------
{
  flt <- set_filters(
    tgt,
    loan_amnt     = "decrease",
    loan_int_rate = "fixed"
  )
  print(flt)
}


# ---- Step 13: Find local counterfactual via biplot rotation -----------
{
  bl_local <- bl_find_local_cf(
    bl_result   = bl_results_v2,
    bl_target   = tgt,
    set_filters = flt
  )
  print(bl_local)
}


# ---- Step 14: Local biplot plot ----------------------------------------
{
  plot(bl_local)
}


# ---- Step 15: Shapley contribution plot --------------------------------
{
  bl_shapley_values <- bl_shapley(bl_local)
  print(bl_shapley_values)
  plot(bl_shapley_values)
}


# ---- Step 16: Sparse counterfactual ------------------------------------
{
  bl_sparse <- bl_find_sparse_cf(bl_shapley_values, round_to = NULL)
  print(bl_sparse)
  plot(bl_sparse)
}


# ===========================================================
# ALTERNATIVE: unconstrained local search
# ===========================================================

# ---- Step 17 (optional): Unconstrained local search ------------------
{
  bl_local_free  <- bl_find_local_cf(bl_results_v2, tgt)
  print(bl_local_free)
  plot(bl_local_free)

  bl_shap_free   <- bl_shapley(bl_local_free)
  bl_sparse_free <- bl_find_sparse_cf(bl_shap_free, round_to = NULL)

  plot(bl_shap_free)
  plot(bl_sparse_free)
}


# ===========================================================
# EXTERNAL APPLICANT
# ===========================================================

# ---- Step 18: External applicant --------------------------------------
{
  new_applicant <- data.frame(
    person_age                 = 32,
    person_gender              = 1,      # male
    person_income              = 248000,
    person_emp_exp             = 5,
    person_home_ownership      = 3,      # RENT
    loan_amnt                  = 12000,
    loan_int_rate              = 13.5,
    loan_percent_income        = 12000 / 248000,
    cb_person_cred_hist_length = 4,
    credit_score               = 420
  )
  new_applicant <- new_applicant[, intersect(names(new_applicant), feature_cols_v2),
                                  drop = FALSE]

  tgt_ext <- bl_select_target(bl_results_v2, target = new_applicant)
  print(tgt_ext)

  plot_biplotEZ(
    bl_results_v2,
    points       = test_pts_v2,
    target_point = unlist(tgt_ext$x_obs),
    target_label = "new"
  )

  bl_local_ext  <- bl_find_local_cf(bl_results_v2, tgt_ext)
  print(bl_local_ext)
  plot(bl_local_ext)

  bl_shap_ext   <- bl_shapley(bl_local_ext)
  bl_sparse_ext <- bl_find_sparse_cf(bl_shap_ext)

  plot(bl_shap_ext)
  plot(bl_sparse_ext)
}