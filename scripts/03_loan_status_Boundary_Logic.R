############################################################
# Loan default — Iterative Phases 1, 2 & 3
#
# Dataset : inst/extdata/loan_data.csv
#           45,000 observations, 14 columns (5 categorical), binary outcome
#           Source: https://github.com/TSMathi/loan_approval_analysis/tree/main
# Model   : SVM baseline (bl_fit_model) + XGBoost full analysis (bl_wrap_model)
# Biplot  : PCA
#
# Iterative workflow:
#
#   Step 1   : Load and integer-encode categorical variables
#   Step 2   : Exploratory biplot (no model) — inspect cluster structure
#   Step 3   : Domain filter: keep only prior-defaulter sub-population
#   Steps 4-6: Fit XGBoost model, build first PCA biplot
#
#   Phase 2 (Steps 7-8): Global interpretations
#     7. Find nearest boundary point + distance-to-boundary plots
#     8. Surrogate model
#
#   Prune least-important variables
#   Steps 4b-6b: Refit XGBoost on reduced feature set, rebuild biplot
#   Phase 2 second pass (Steps 7b-8b): re-inspect with reduced model
#
#   Phase 3 (Steps 9-16): Local interpretation of one target
#     9. Inspect predictions and select target
#    10. Set actionability constraints
#    11. Find local counterfactual via biplot rotation
#    12. Local biplot plot
#    13. Shapley contribution plot
#    14. Sparse counterfactual
#    15. (Optional) Unconstrained local search
#    16. External applicant
#
#   Comparison with SHAP: validate global/local attributions against
#                          fastshap + shapviz (optional Suggests packages)
#
# Run interactively: place cursor inside a {} block and press Ctrl+Enter
############################################################

# ---- Load the package (development workflow) ---------------------------
{
  rm(list = ls())
  devtools::load_all()
  # library(boundary_logic) # alternative if downloaded the package from github
}
#devtools::check()

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
  bl_dat_exp <- bl_prepare_data(
    data           = loan_encoded,
    class_col      = "loan_status",
    train_fraction = 1,
    seed           = 121L,
    hull_fraction  = 1
  )

  bl_proj <- bl_build_result(
    bl_data = bl_dat_exp,
    method  = "PCA",
    title   = "Loan default — exploratory PCA biplot (all data)"
  )

  plot(
    bl_proj,
    label_dir         = "Paral",  # "Hor" = horizontal, "Orthog" = orthogonal to axis
    label_offset_var  = c("person_age",
                          "person_gender",
                          "person_income",                 
                          "person_emp_exp",
                          "loan_amnt",
                          "loan_int_rate",                 
                          "loan_percent_income",
                          "cb_person_cred_hist_length",
                          "credit_score",                  
                          "previous_loan_defaults_on_file"),     # variable index/indices to shift, e.g. c(1L, 3L)
    label_offset_dist = c(0,0,0.5,1,0.5,1,0,1.5,0,0 )    # outward distance per shifted label
  )
}

names(loan_encoded)
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
# Standard path: bl_prepare_data() handles feature selection, class encoding,
# seeded train/test split, and outlier filtering in one call.
#
# Alternative path using bl_wrap_data() — use this when your data is already
# split (e.g. from a cross-validation framework or external pipeline):
#
#   data_clean <- loan_filtered[, c(feature_cols, "loan_status")]
#   names(data_clean)[names(data_clean) == "loan_status"] <- "class"
#   set.seed(121L)
#   n          <- nrow(data_clean)
#   train_idx  <- sample(n, size = floor(0.8 * n), replace = FALSE)
#   train_data <- data_clean[ train_idx, , drop = FALSE]
#   test_data  <- data_clean[-train_idx, , drop = FALSE]
#   rownames(train_data) <- NULL; rownames(test_data) <- NULL
#   

#   bl_dat <- bl_wrap_data(
#     train_data   = train_data,
#     test_data    = test_data,
#     var_names    = feature_cols,
#     target_class = NULL       # loan_status already encoded as 0/1
#   )
#   # bl_wrap_data() does not filter outliers; call bl_filter_outliers() if needed:
#   # bl_dat <- bl_filter_outliers(bl_dat, hull_fraction = 0.9)

{
  bl_dat <- bl_prepare_data(
    data           = loan_filtered,
    class_col      = "loan_status",
    feature_cols   = feature_cols,
    train_fraction = 0.8,
    seed           = 121L,
    hull_fraction  = 0.9
  )
  print(bl_dat)

  # ---- Step 5a: SVM via bl_fit_model() (quick baseline) ------------------
  # bl_mod_svm <- bl_fit_model(
  #   train_data = bl_dat$train_data,
  #   var_names  = bl_dat$var_names,
  #   model_type = "SVM"
  # )
  # print(bl_mod_svm)

  # ---- Step 5b: XGB via bl_wrap_model() with explicit predict_fn ---------
  
  # ---- Option A: XGBoost binary format ----------------------------------
  # Saved previously with: xgboost::xgb.save(xgb_fit, "path/to/model.xgb")
  # xgb_fit <- xgboost::xgb.load("path/to/model.xgb")
  
  # ---- Option B: R object format ----------------------------------------
  # Saved previously with: saveRDS(xgb_fit, "path/to/model.rds")
  # xgb_fit <- readRDS("path/to/model.rds")

  # For demonstration, fit a quick model here so the script is self-contained.
  # Replace these lines with one of the xgb.load / readRDS calls above.
  
  # Using model_type = "custom" makes the prediction contract explicit.

  

  
  xgb_data <- xgboost::xgb.DMatrix(
    data  = as.matrix(bl_dat$train_data[, bl_dat$var_names]),
    label = bl_dat$train_data$class
  )
  xgb_fit <- xgboost::xgb.train(
    params  = list(objective   = "binary:logistic",
                   eval_metric = "logloss",
                   max_depth   = 3,
                   learning_rate = 0.1),
    data    = xgb_data,
    nrounds = 200,
    verbose = 0
  )
  
  # wrap the custom model into a "bl_model" object (bl_mod), as required by bl_build_result()'s bl_model argument
  bl_mod <- bl_wrap_model(
    model      = xgb_fit,
    model_type = "custom",
    var_names  = bl_dat$var_names,
    predict_fn = function(m, new_data) {
      mat <- xgboost::xgb.DMatrix(as.matrix(new_data))
      as.numeric(predict(m, newdata = mat))
    },
    train_data = bl_dat$train_data
  )
  print(bl_mod)

  # CVA is the default biplot method for this workflow — it maximises
  # class separation in the projection plane.
  bl_results <- bl_build_result(
    bl_data  = bl_dat,
    bl_model = bl_mod,
    method   = "CVA",
    title    = "Loan default (prior defaulters) -- XGB, CVA biplot",
    b_margin = 0.001
  )
  print(bl_results)
  bl_results$test_data
  # Reference biplot — training data coloured by confusion category
  plot(
    bl_results,
    label_dir         = "Hor",  # "Hor" = horizontal, "Orthog" = orthogonal to axis
    label_offset_var  = 0L,     # variable index/indices to shift, e.g. c(1L, 3L)
    label_offset_dist = 0.5     # outward distance per shifted label
  )

  # Project all observations and overlay on the biplot
  # Create test data points to plot instead of the default training data
  # bl_results$test_data can be any data.frame with the necessary column names.
  # Convert to the format required by plot() using the bl_project_points() function. 
  # filter_to_polygon = TRUE removes new data points that are outside the training data biplot polygon
  
  test_pts <- bl_project_points(bl_results$test_data, bl_results, filter_to_polygon = TRUE )   # removes out-of-polygon points before plotting)
  plot(bl_results, points = test_pts)
}


# ===========================================================
# PHASE 2 (first pass) — Global interpretations
# ===========================================================


# ---- Step 7: Find boundary + distance-to-boundary plots ---------------
# Y-axis label shows "variable : total absolute standardised distance" —
# a variable-level importance proxy. Sorted ascending (least important first).
# Colour scheme: TP = red, TN = blue, FP = purple, FN = orange.
{
  bl_bnd <- bl_find_boundary(bl_results)
  hist(bl_bnd$B_pred)

  plot(bl_results, points = test_pts)
  # To inspect individual counterfactuals: bl_pick_point(bl_results, bl_boundary = bl_bnd)

  plot(bl_bnd)
  plot(bl_bnd, type = "boxplot")
}


# ---- Step 8: Surrogate model -------------------------------------------
{
  bl_surr <- bl_surrogate(bl_results)
  print(bl_surr)
  plot(bl_surr)
}


# ===========================================================
# Prune least-important variables
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
  
  feature_cols_v2 <- setdiff(bl_dat$var_names, vars_to_drop)
  cat("Retained features:", paste(feature_cols_v2, collapse = ", "), "\n")
}


# ===========================================================
# PHASE 1 (refit) — XGBoost on reduced feature set
# ===========================================================

# ---- Steps 4b-6b: Rebuild pipeline with pruned features ----------------
{
  bl_dat_v2 <- bl_prepare_data(
    data           = loan_filtered,
    class_col      = "loan_status",
    feature_cols   = feature_cols_v2,
    train_fraction = 0.8,
    seed           = 121L,
    hull_fraction  = 0.9
  )

  # Direct XGB path — no predict_fn required; see Step 5b for the
  # custom/predict_fn alternative demonstrated on the full feature set.
  xgb_data_v2 <- xgboost::xgb.DMatrix(
    data  = as.matrix(bl_dat_v2$train_data[, bl_dat_v2$var_names]),
    label = bl_dat_v2$train_data$class
  )
  xgb_fit_v2 <- xgboost::xgb.train(
    params  = list(objective   = "binary:logistic",
                   eval_metric = "logloss",
                   max_depth   = 3,
                   learning_rate = 0.1),
    data    = xgb_data_v2,
    nrounds = 200,
    verbose = 0
  )
  bl_mod_v2 <- bl_wrap_model(
    model      = list(model = xgb_fit_v2, features = bl_dat_v2$var_names),
    model_type = "XGB",
    var_names  = bl_dat_v2$var_names,
    train_data = bl_dat_v2$train_data
  )
  print(bl_mod_v2)

  bl_results_v2 <- bl_build_result(
    bl_data  = bl_dat_v2,
    bl_model = bl_mod_v2,
    method   = "CVA",
    title    = "Loan default -- XGB, reduced features, CVA biplot",
    b_margin = 0.01
  )
  print(bl_results_v2)

  plot(bl_results_v2)
  test_pts_v2 <- bl_project_points(bl_results_v2$test_data, bl_results_v2)
  plot(bl_results_v2, points = test_pts_v2)
}


# ===========================================================
# PHASE 2 (second pass) — Global interpretations, reduced model
# ===========================================================

# ---- Step 7b: Boundary + distance-to-boundary plots (reduced model) ---
{
  bl_bnd_v2 <- bl_find_boundary(bl_results_v2)
  print(bl_bnd_v2)
  plot(bl_results_v2, points = test_pts_v2)
  # To inspect individual counterfactuals: bl_pick_point(bl_results_v2, bl_boundary = bl_bnd_v2)

  plot(bl_bnd_v2)
  plot(bl_bnd_v2, type = "boxplot")
}


# ---- Step 8b: Surrogate model ------------------------------------------
{
  bl_surr_v2 <- bl_surrogate(bl_results_v2)
  print(bl_surr_v2)
  plot(bl_surr_v2,
       label_dir         = "Paral",   # default; "Hor" and "Orthog" also accepted
       label_offset_var  = c("person_age",
                             "loan_amnt",
                             "loan_int_rate",
                             "loan_percent_income",
                             "credit_score"),         # or a character/integer vector of variable names/indices
       label_offset_dist = c(0,0,0.5,0,0),
       ticks_v = 5)
}


# ===========================================================
# PHASE 3 — Local interpretation (using reduced model bl_results_v2)
# ===========================================================

# ---- Step 9: Inspect predictions and select target --------------------
# Review pred_summary to choose a target. False Negatives (predicted 0,
# true 1) are most actionable: the applicant will default but the model
# clears them — what would need to change to flag them correctly?
{
  pred_summary <- bl_predict(bl_results_v2)

  tdp <- 1
  
  # Generate a subset of data to plot on the biplot. Use bl_project_points to convert a data.frame bl_results_v2$test_data[tdp,] into the format required for plotting on biplot
  # Can be empty, or only the target data point to avoid unnecessary information on the biplot
  test_points <- bl_project_points(bl_results_v2$test_data[tdp:10,], bl_results_v2)

  target_value <- bl_results_v2$test_data[tdp, ]
  # Highlight the target on the main biplot

  plot(
    bl_results_v2,
    points       = test_points, # specify the points to plot. If emply, will plot the training data
    plot_points  = FALSE, # removes the points from the plot
    target_point = target_value,
    target_label = tdp,
    label_dir         = "Paral",   # default; "Hor" and "Orthog" also accepted
    label_offset_var  = c("person_age",
                          "loan_amnt",
                          "loan_int_rate",
                          "loan_percent_income",
                          "credit_score"),         # or a character/integer vector of variable names/indices
    label_offset_dist = c(0,0,0.5,0,0),
    ticks_v = 3
  )
}





# ---- Step 10: Set actionability constraints ---------------------------
# Omit any pruned variables from set_filters() — they are no longer in
# the model and passing them will raise an error.
#
# Valid constraint types:
#   "decrease"  — counterfactual value must be <= observed
#   "increase"  — counterfactual value must be >= observed
#   "fixed"     — constrained to ± 0.5 of observed; always reverts in sparse CF
#   c(min, max) — counterfactual must lie within this absolute range
{
  tgt <- bl_select_target(bl_results_v2, target = tdp)

  flt <- set_filters(
    tgt,
    person_age     = "fixed",       # not actionable
   # loan_amnt      = "fixed",    # borrow less to reduce repayment risk
    loan_int_rate  = "increase",        # set by the lender, not the applicant
    credit_score = "increase"
  )
  print(flt)
}


# ---- Step 11: Find local counterfactual via biplot rotation ----------
{
  bl_local <- bl_find_local_cf(
    bl_result   = bl_results_v2,
    set_filters = flt,
    bl_target   = tgt
    )
  print(bl_local)
}


# ---- Step 12: Local biplot plot ----------------------------------------
{
  plot(bl_local,
       label_dir         = "Paral",   # default; "Hor" and "Orthog" also accepted
       label_offset_var  = c("person_age",
                             "loan_amnt",
                             "loan_int_rate",
                             "loan_percent_income",
                             "credit_score"),         # or a character/integer vector of variable names/indices
       label_offset_dist = c(0,0.5,0,0,0),
       ticks_var = c("person_age",
                      "loan_amnt",
                      "loan_int_rate",
                      "loan_percent_income",
                      "credit_score"),
       ticks_n = c(2,4,200,200,20)
  )
}


# ---- Step 13: Shapley contribution plot --------------------------------
{
  bl_shapley_values <- bl_shapley(bl_local)
  print(bl_shapley_values)
  plot(bl_shapley_values)
}


# ---- Step 14: Sparse counterfactual ------------------------------------
{
  bl_sparse <- bl_find_sparse_cf(bl_shapley_values, round_to = NULL)
  print(bl_sparse)
  plot(bl_sparse,
       label_dir         = "Paral",   # default; "Hor" and "Orthog" also accepted
       label_offset_var  = c("person_age",
                             "loan_amnt",
                             "loan_int_rate",
                             "loan_percent_income",
                             "credit_score"),         # or a character/integer vector of variable names/indices
       label_offset_dist = c(0,0.5,0,0,0),
       ticks_var = c("person_age",
                     "loan_amnt",
                     "loan_int_rate",
                     "loan_percent_income",
                     "credit_score"),
       ticks_n = c(2,4,200,200,20))
}


# ===========================================================
# ALTERNATIVE: unconstrained local search
# ===========================================================

# ---- Step 15 (optional): Unconstrained local search ------------------
{
  bl_local_free  <- bl_find_local_cf(bl_results_v2, tgt)
  print(bl_local_free)
  plot(bl_local_free,
       label_dir         = "Paral",   # default; "Hor" and "Orthog" also accepted
       label_offset_var  = c("person_age",
                             "loan_amnt",
                             "loan_int_rate",
                             "loan_percent_income",
                             "credit_score"),         # or a character/integer vector of variable names/indices
       label_offset_dist = c(0,0,0,0,1),
       ticks_var = c("person_age",
                     "loan_amnt",
                     "loan_int_rate",
                     "loan_percent_income",
                     "credit_score"),
       ticks_n = c(2,4,200,200,20)
  )

  bl_shap_free   <- bl_shapley(bl_local_free)
  bl_sparse_free <- bl_find_sparse_cf(bl_shap_free, round_to = NULL)

  plot(bl_shap_free)
  plot(bl_sparse_free,
       label_dir         = "Paral",   # default; "Hor" and "Orthog" also accepted
       label_offset_var  = c("person_age",
                             "loan_amnt",
                             "loan_int_rate",
                             "loan_percent_income",
                             "credit_score"),         # or a character/integer vector of variable names/indices
       label_offset_dist = c(0,0,0,0,1),
       ticks_var = c("person_age",
                     "loan_amnt",
                     "loan_int_rate",
                     "loan_percent_income",
                     "credit_score"),
       ticks_n = c(2,4,200,200,20)
  )
}


# ===========================================================
# EXTERNAL APPLICANT
# Analyse a new loan application not in the training data.
# Only include variables in feature_cols_v2 (remove any pruned variables).
# No class column is required.
# ===========================================================

# ---- Step 16: External applicant --------------------------------------
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
  # Drop any columns that were pruned in Step 9
  new_applicant <- new_applicant[, intersect(names(new_applicant), feature_cols_v2),
                                  drop = FALSE]


  # Local search — no actionability constraints for this new applicant
  tgt_ext       <- bl_select_target(bl_results_v2, target = new_applicant)
  bl_local_ext  <- bl_find_local_cf(bl_results_v2, tgt_ext)
  print(bl_local_ext)
  plot(bl_local_ext)

  bl_shap_ext   <- bl_shapley(bl_local_ext)
  bl_sparse_ext <- bl_find_sparse_cf(bl_shap_ext)
  print(bl_sparse_ext)
  plot(bl_shap_ext)
  plot(bl_sparse_ext)
}


# ===========================================================
# COMPARISON WITH SHAP
# Validate Boundary Logic's global (Phase 2) and local (Step 13)
# attributions against real SHAP values for the same XGBoost models,
# using fastshap + shapviz (optional Suggests packages).
# ===========================================================

# ---- SHAP setup ----------------------------------------------------------
{
  has_shap <- requireNamespace("fastshap", quietly = TRUE) &&
              requireNamespace("shapviz",  quietly = TRUE)

  if (!has_shap) {
    cat("SHAP packages not available. Install fastshap and shapviz to run",
        "this section:\n  install.packages(c(\"fastshap\", \"shapviz\"))\n")
  }

  if (has_shap) {
    # fastshap needs a prediction wrapper f(object, newdata) -> numeric vector.
    # For XGBoost the new data must be wrapped in an xgb.DMatrix.
    pred_fn <- function(object, newdata) {
      as.numeric(predict(object, xgboost::xgb.DMatrix(as.matrix(newdata))))
    }

    # Background / explanation data for each model
    X_train_full <- bl_dat$train_data[,    bl_dat$var_names,    drop = FALSE]
    X_train_v2   <- bl_dat_v2$train_data[, bl_dat_v2$var_names, drop = FALSE]
    X_test_v2    <- bl_dat_v2$test_data[,  bl_dat_v2$var_names, drop = FALSE]

    # SHAP baseline = mean predicted probability over the (reduced-model)
    # training background. fastshap does not attach one, so shapviz needs it.
    shap_baseline <- mean(pred_fn(xgb_fit_v2, X_train_v2))
  }
}


# ---- Global: SHAP vs Boundary Logic variable selection --------------------
# Phase 2's distance-to-boundary importance (rob$sum_of_distance) picked the
# five variables in feature_cols_v2. Ask SHAP the same question on the full
# nine-feature model: rank by mean absolute SHAP value and compare its top
# five against the distance-to-boundary top five.
{
  if (has_shap) {
    set.seed(42L)
    shap_full <- fastshap::explain(
      object       = xgb_fit,            # full-feature model
      X            = X_train_full,
      pred_wrapper = pred_fn,
      nsim         = 50L
    )

    shap_imp  <- names(sort(colMeans(abs(as.matrix(shap_full))), decreasing = TRUE))
    shap_top5 <- shap_imp[1:5]

    # Distance-to-boundary ranking (strongest first) from Phase 2
    dist_imp <- names(sort(rob$sum_of_distance, decreasing = TRUE))
    dist_top5 <- names(sort(rob$sum_of_distance, decreasing = TRUE))[1:5]

    print(data.frame(
     # rank          = 1:5,
      SHAP          = shap_imp,
      dist_to_bound = dist_imp
    ))
  }
}


# ---- Global SHAP of the retained five variables ----------------------------
# Run SHAP on the reduced model -- the five variables the workflow proceeded
# with -- and view the standard global summaries. Compare against
# plot(bl_bnd_v2) and bl_surrogate(bl_results_v2) above.
{
  if (has_shap) {
    set.seed(42L)
    shap_v2 <- fastshap::explain(
      object       = xgb_fit_v2,         # reduced model
      X            = X_train_v2,
      pred_wrapper = pred_fn,
      nsim         = 50L
    )
    sv_v2 <- shapviz::shapviz(as.matrix(shap_v2), X = X_train_v2,
                              baseline = shap_baseline)

    print(shapviz::sv_importance(sv_v2, kind = "bar") +
            ggplot2::labs(title = "Global SHAP importance - reduced XGB model"))
  }
}

plot(bl_bnd_v2)

{
  if (has_shap) {
    print(shapviz::sv_importance(sv_v2, kind = "beeswarm") +
            ggplot2::labs(title = "Global SHAP beeswarm - reduced XGB model"))
  }
}


# ---- Local: SHAP vs bl_shapley attribution ----------------------------------
# For the same target applicant explained in Step 13, compute local SHAP and
# view the waterfall and force plots. Compare against plot(bl_shapley_values).
{
  if (has_shap) {
    set.seed(42L)
    shap_local <- fastshap::explain(
      object       = xgb_fit_v2,
      X            = X_train_v2,
      pred_wrapper = pred_fn,
      nsim         = 200L,
      newdata      = X_test_v2[tdp, , drop = FALSE],
      adjust       = TRUE          # force additivity: sum(shap) = f(x) - baseline
    )
    sv_local <- shapviz::shapviz(as.matrix(shap_local),
                                 X = X_test_v2[tdp, , drop = FALSE],
                                 baseline = shap_baseline)

    print(shapviz::sv_waterfall(sv_local) +
            ggplot2::labs(title = sprintf("Local SHAP waterfall - applicant %d", tdp)))
  }
}

{
  if (has_shap) {
    print(shapviz::sv_force(sv_local) +
            ggplot2::labs(title = sprintf("Local SHAP force plot - applicant %d", tdp)))
  }
}


