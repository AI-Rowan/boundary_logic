############################################################
# Loan default — Iterative Phases 1, 2 & 3
#
# Dataset : inst/extdata/loan_data.csv
#           45,000 observations, 14 columns (5 categorical), binary outcome
#           Source: https://github.com/TSMathi/loan_approval_analysis/tree/main
# Model   : XGBoost (xgboost)
# Biplot  : PCA
#
# Iterative workflow:
#
#   Step 1   : Load and integer-encode categorical variables
#   Step 2   : Exploratory biplot (no model) — inspect cluster structure
#   Step 3   : Domain filter: keep only prior-defaulter sub-population
#   Steps 4-6: Fit XGBoost model, build first PCA biplot
#
#   Phase 2 (Steps 7-9): Global interpretations
#     7. Find nearest boundary point for each observation
#     8. Distance-to-boundary plot (jitter + boxplot)
#     9. Surrogate model
#
#   Step 10  : Extract per-variable importance, prune weak variables
#   Steps 4b-6b: Refit XGBoost on reduced feature set, rebuild biplot
#   Phase 2 second pass (Steps 7b-9b): re-inspect with reduced model
#
#   Phase 3 (Steps 11-18): Local interpretation of one target
#    11. Inspect predictions and select target
#    12. Set actionability constraints
#    13. Find local counterfactual via biplot rotation
#    14. Local biplot plot
#    15. Shapley contribution plot
#    16. Sparse counterfactual
#    17. (Optional) Unconstrained local search
#    18. External applicant
#
# Run interactively: place cursor inside a {} block and press Ctrl+Enter
############################################################

# ---- Load the package (development workflow) ---------------------------
{
  rm(list = ls())
  devtools::load_all()
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

  # Coerce all columns to double so the training data and prediction grid
  # share the same type — read.csv() reads integer-valued columns as integer,
  # but bl_build_grid() generates sequences as double.
  loan_encoded[] <- lapply(loan_encoded, as.numeric)

  str(loan_encoded)

  vars_to_remove <- c("loan_intent", "person_education", "person_home_ownership")
  loan_encoded <- loan_encoded[, setdiff(names(loan_encoded), vars_to_remove), drop = FALSE]
}


# ---- Step 2: Exploratory biplot — full data, no model -----------------
# Build a PCA biplot coloured by loan_status before fitting any model.
# Purpose: inspect class separation and variable loading directions
# to guide modelling and feature decisions.
#
# bl_build_result() with no bl_model argument returns a bl_projection
# object (not a full bl_result). plot_biplotEZ() renders it in the
# same style as a model-backed biplot.
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

  plot_biplotEZ(
    bl_proj,
    label_dir         = "Hor",  # "Hor" = horizontal, "Orthog" = orthogonal to axis
    label_offset_var  = 0L,     # variable index/indices to shift, e.g. c(1L, 3L)
    label_offset_dist = 0.5     # outward distance per shifted label
  )
}


# ---- Step 3: Domain filter — keep applicants without prior defaults ----
# Retain only applicants with no prior default on file
# (previous_loan_defaults_on_file == 0). This sub-population is the
# focus of the credit-risk model.
#
# Because previous_loan_defaults_on_file is now constant (all 0) within
# the filtered set, it is excluded from feature_cols to avoid zero-variance
# issues in the model and projection.
{
  loan_filtered <- loan_encoded[loan_encoded$previous_loan_defaults_on_file == 0, ]

  feature_cols <- setdiff(
    names(loan_filtered),
    c("loan_status", "previous_loan_defaults_on_file")
  )

  cat("Retained rows:", nrow(loan_filtered), "\n")
  cat("Features:", paste(feature_cols, collapse = ", "), "\n")
}


# ---- Steps 4-6: Prepare data, fit XGBoost, build biplot ---------------
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

  # NOTE: GAM via bl_fit_model() is currently broken — use XGB or GLM.
  # To fit a custom GAM, use bl_wrap_model() with mgcv::gam() directly
  # (see scripts/00_pima_Boundary_Logic.R Step 3 for the pattern).
  bl_mod <- bl_fit_model(
    train_data = bl_filt$train_data,
    var_names  = bl_filt$var_names,
    model_type = "XGB"
  )
  print(bl_mod)

  # CVA is the default biplot method for this workflow — it maximises
  # class separation in the projection plane.
  bl_results <- bl_build_result(
    bl_data  = bl_filt,
    bl_model = bl_mod,
    method   = "CVA",
    title    = "Loan default (prior defaulters) — XGB, CVA biplot",
    rounding = 3L
  )
  print(bl_results)
  bl_results$test_data
  # Reference biplot — training data coloured by confusion category
  plot_biplotEZ(
    bl_results,
    label_dir         = "Hor",  # "Hor" = horizontal, "Orthog" = orthogonal to axis
    label_offset_var  = 0L,     # variable index/indices to shift, e.g. c(1L, 3L)
    label_offset_dist = 0.5     # outward distance per shifted label
  )

  # Project all observations and overlay on the biplot
  test_pts <- bl_project_points(bl_results$test_data, bl_results)
  sum(test_pts$inside_polygon==F)
  plot_biplotEZ(bl_results, points = test_pts)
}


# ===========================================================
# PHASE 2 (first pass) — Global interpretations
# ===========================================================


# ---- Step 7: Find nearest boundary point for each observation ----------
{
  bl_bnd <- bl_find_boundary(bl_results)
  print(bl_bnd)
  hist(bl_bnd$B_pred)
  bl_bnd$B_pred
  bl_bnd$x_obs
  bl_bnd$pred_obs
  plot_biplotEZ(bl_results, points = test_pts, boundary = bl_bnd)
}


# ---- Step 8: Distance-to-boundary plot ---------------------------------
# Y-axis label shows "variable : total absolute standardised distance" —
# a variable-level importance proxy. Sorted ascending (least important first).
# Colour scheme: TP = red, TN = blue, FP = purple, FN = orange.
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
# bl_robustness() prints sum_of_distance: the total absolute standardised
# distance contribution per variable (ascending = weakest first).
# Review var_imp, then explicitly list the variables to drop in vars_to_drop.
{
  rob     <- bl_robustness(bl_bnd)
  var_imp <- sort(rob$sum_of_distance)     # ascending: weakest variable first
  print(round(var_imp, 2))

  # List variables to remove based on the var_imp output above.
  vars_to_drop <- c("person_gender", "person_emp_exp", "cb_person_cred_hist_length", "person_income")
  
  feature_cols_v2 <- setdiff(bl_filt$var_names, vars_to_drop)
  cat("Retained features:", paste(feature_cols_v2, collapse = ", "), "\n")
}


# ===========================================================
# PHASE 1 (refit) — XGBoost on reduced feature set
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

  bl_mod_v2 <- bl_fit_model(
    train_data = bl_filt_v2$train_data,
    var_names  = bl_filt_v2$var_names,
    model_type = "XGB"
  )
  print(bl_mod_v2)

  bl_results_v2 <- bl_build_result(
    bl_data  = bl_filt_v2,
    bl_model = bl_mod_v2,
    method   = "CVA",
    title    = "Loan default — XGB, reduced features, CVA biplot",
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
  plot_biplotEZ(bl_results_v2, points = test_pts_v2, boundary = bl_bnd_v2)
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
# PHASE 3 — Local interpretation (using reduced model bl_results_v2)
# ===========================================================

# ---- Step 11: Inspect predictions and select target -------------------
# Review pred_summary to choose a target. False Negatives (predicted 0,
# true 1) are most actionable: the applicant will default but the model
# clears them — what would need to change to flag them correctly?
{
  pred_summary <- bl_predict(bl_results_v2)
  print(pred_summary)

  tdp <- 2
  tgt <- bl_select_target(bl_results_v2, target = tdp)
  print(tgt)

  # Highlight the target on the main biplot
  plot_biplotEZ(
    bl_results_v2,
    points       = test_pts_v2,
    target_point = unlist(tgt$x_obs),
    target_label = tdp
  )
}


# ---- Step 12: Set actionability constraints ----------------------------
# Omit any pruned variables from set_filters() — they are no longer in
# the model and passing them will raise an error.
#
# Valid constraint types:
#   "decrease"  — counterfactual value must be <= observed
#   "increase"  — counterfactual value must be >= observed
#   "fixed"     — constrained to ± 0.5 of observed; always reverts in sparse CF
#   c(min, max) — counterfactual must lie within this absolute range
{
  flt <- set_filters(
    tgt,
  #  person_age     = "fixed",       # not actionable
   # person_gender  = "fixed",       # not actionable
  #  person_emp_exp = "increase",    # can only grow over time
    loan_amnt      = "decrease",    # borrow less to reduce repayment risk
    loan_int_rate  = "fixed"        # set by the lender, not the applicant
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
# Analyse a new loan application not in the training data.
# Only include variables in feature_cols_v2 (remove any pruned variables).
# No class column is required.
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
  # Drop any columns that were pruned in Step 10
  new_applicant <- new_applicant[, intersect(names(new_applicant), feature_cols_v2),
                                  drop = FALSE]

  tgt_ext <- bl_select_target(bl_results_v2, target = new_applicant)
  print(tgt_ext)

  # Highlight on main biplot
  plot_biplotEZ(
    bl_results_v2,
    points       = test_pts_v2,
    target_point = unlist(tgt_ext$x_obs),
    target_label = "new"
  )

  # Local search — no actionability constraints for this new applicant
  bl_local_ext  <- bl_find_local_cf(bl_results_v2, tgt_ext)
  print(bl_local_ext)
  plot(bl_local_ext)

  bl_shap_ext   <- bl_shapley(bl_local_ext)
  bl_sparse_ext <- bl_find_sparse_cf(bl_shap_ext)

  plot(bl_shap_ext)
  plot(bl_sparse_ext)
}

