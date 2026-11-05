############################################################
# Loan default — Iterative Phases 1, 2 & 3
#
# Dataset : inst/extdata/loan_data.csv
#           45,000 observations, 14 columns (5 categorical), binary outcome
#           Source: https://github.com/TSMathi/loan_approval_analysis/tree/main
# Model   : XGBoost (bl_wrap_model); an SVM baseline via bl_fit_model() is
#           left commented out in Step 5a
# Biplot  : PCA (exploratory, Step 2) and CVA (model biplots, Steps 4-16)
#
# Iterative workflow:
#
#   Step 1   : Load and integer-encode categorical variables
#   Step 2   : Exploratory biplot (no model) — inspect cluster structure
#   Step 3   : Domain filter: keep only applicants with no prior default on file
#   Steps 4-6: Fit XGBoost model, build first CVA biplot
#
#   Phase 2 (Step 7): Global interpretations
#     7. Find nearest boundary point + distance-to-boundary plots
#
#   Prune least-important variables
#   Steps 4b-6b: Refit XGBoost on reduced feature set, rebuild biplot
#   Phase 2 second pass (Step 7b): re-inspect with reduced model
#   Steps 8 and 8b: Surrogate models, full and reduced (not used directly
#                   in this workflow)
#
#   Phase 3 (Steps 9-16): Local interpretation of one target
#     9. Inspect predictions, choose the target row, and view its global
#        counterfactual (interactive bl_pick_point())
#    10. Select the target and set actionability constraints
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
# bl_build_result() with no bl_model argument returns a model-free
# bl_result (no prediction grid or contours). plot_biplotEZ() renders it
# coloured by true class.
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
# bl_prepare_data() (outlier hull) and bl_build_result() (hull + grid bounds)
# compute by drawing on a temporary off-screen device. They close only that
# device and re-activate yours, so plots you have open stay the active ones.
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
    title    = "Loan default (no prior default) -- XGB, CVA biplot",
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
  # To inspect individual counterfactuals (bl_bnd rows are test_data rows;
  # Step 9 runs this for the target):
  # bl_pick_point(bl_results, data = bl_results$test_data, bl_boundary = bl_bnd)

  plot(bl_bnd)
  plot(bl_bnd, type = "boxplot")
}


# ===========================================================
# Prune least-important variables
# ===========================================================
# bl_robustness() returns sum_of_distance: the total absolute standardised
# distance contribution per variable (sorted ascending = weakest first).
# Review var_imp, then explicitly list the variables to drop in vars_to_drop.
{
  rob     <- bl_robustness(bl_bnd)
  var_imp <- sort(rob$sum_of_distance)     # ascending: weakest variable first
  print(round(var_imp, 2))

  # The overall score is the sum of the per-variable totals. It is printed
  # for reference only; the pruning uses the per-variable totals above.
  cat(sprintf("Overall robustness (not used here): %.2f\n", rob$robustness))

  # List variables to remove based on the var_imp output above.
  # person_income goes although it outranks credit_score: loan_percent_income
  # is loan_amnt / person_income, so income is redundant once both loan
  # variables are kept. credit_score stays because it is actionable (Step 10).
  vars_to_drop <- c("person_gender", "person_emp_exp", "cb_person_cred_hist_length", "person_income", "person_age")
  
  feature_cols_v2 <- setdiff(bl_dat$var_names, vars_to_drop)
  cat("Retained features:", paste(feature_cols_v2, collapse = ", "), "\n")
}


# ===========================================================
# PHASE 1 (refit) — XGBoost on reduced feature set
# ===========================================================

# ---- Steps 4b-6b: Rebuild pipeline with pruned features ----------------
{
  # ---- Standardise the retained features (raw-unit biplot axes demo) ----
  # Scale-sensitive models (SVM, NNET, logistic regression) often fit better on
  # standardised inputs, so a common workflow standardises features before
  # modelling. XGBoost is scale-invariant per feature -- its predictions here are
  # unchanged -- but the PCA/CVA projection geometry does depend on scale, and
  # without help the biplot axes would then read in standardised units.
  # bl_set_scaling() records the transform standardised = (raw - center) / scale
  # so every biplot axis can be relabelled back into the original loan units.
  # The scaling is display-only: the model and grid still see standardised data.
  #
  # The standardisation is fit on the TRAINING rows only (no test leakage) and
  # the same training centre/scale are then applied to the test split. So we
  # split first on raw data -- the hull filter standardises internally, so
  # filtering on raw vs standardised is equivalent -- then scale.
  bl_dat_v2 <- bl_prepare_data(
    data           = loan_filtered,
    class_col      = "loan_status",
    feature_cols   = feature_cols_v2,
    train_fraction = 0.8,
    seed           = 121L,
    hull_fraction  = 0.9
  )

  # Scaling derived from the (filtered) training features only.
  train_scl  <- scale(bl_dat_v2$train_data[, feature_cols_v2])
  scl_center <- attr(train_scl, "scaled:center")
  scl_scale  <- attr(train_scl, "scaled:scale")

  # Apply the training transform to both splits (test uses the training stats).
  bl_dat_v2$train_data[, feature_cols_v2] <- train_scl
  bl_dat_v2$test_data[,  feature_cols_v2] <-
    scale(bl_dat_v2$test_data[, feature_cols_v2],
          center = scl_center, scale = scl_scale)

  # Record the pre-standardisation transform for raw-unit biplot axes.
  bl_dat_v2 <- bl_set_scaling(bl_dat_v2,
                              center = scl_center,
                              scale  = scl_scale,
                              method = "z-score")

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
  print(bl_results_v2)   # prints "Axis units : raw (z-score scaling stored)"

  # Axes now read in original loan units (years, currency, rate, score) even
  # though the model and projection operate on standardised features.
  plot(bl_results_v2)
  test_pts_v2 <- bl_project_points(bl_results_v2$test_data, bl_results_v2)
  plot(bl_results_v2, points = test_pts_v2)
}


# ===========================================================
# PHASE 2 (second pass) — Global interpretations, reduced model
# ===========================================================
# person_age was pruned in "Prune least-important variables" (vars_to_drop),
# so the label and tick arguments of the later plots (Steps 8b, 9, 12, 14 and
# 15) name only the four retained features. Each label_offset_dist / ticks_n
# vector pairs position by position with its label_offset_var / ticks_var
# vector.

# ---- Step 7b: Boundary + distance-to-boundary plots (reduced model) ---
{
  bl_bnd_v2 <- bl_find_boundary(bl_results_v2)
  print(bl_bnd_v2)
  plot(bl_results_v2, points = test_pts_v2)
  # To inspect individual counterfactuals (bl_bnd_v2 rows are test_data rows;
  # Step 9 runs this for the target):
  # bl_pick_point(bl_results_v2, data = bl_results_v2$test_data, bl_boundary = bl_bnd_v2)

  plot(bl_bnd_v2)
  plot(bl_bnd_v2, type = "boxplot")
}


# ---- Surrogate models (Steps 8 and 8b) ---------------------------------
# The surrogate models are not used directly in this workflow: no later
# step reads bl_surr or bl_surr_v2. bl_surrogate() classifies each training
# row by its nearest boundary contour in the biplot, not by the model, and
# reports how often the two agree -- a check of how faithfully the 2D biplot
# represents the model. Step 8 uses the full model (bl_results), Step 8b
# the reduced model (bl_results_v2).

# ---- Step 8: Surrogate model -------------------------------------------
{
  bl_surr <- bl_surrogate(bl_results)
  print(bl_surr)
  plot(bl_surr)
}


# ---- Step 8b: Surrogate model ------------------------------------------
{
  bl_surr_v2 <- bl_surrogate(bl_results_v2)
  print(bl_surr_v2)
  plot(bl_surr_v2,
       label_dir         = "Paral",   # default; "Hor" and "Orthog" also accepted
       label_offset_var  = c("loan_amnt",
                             "loan_int_rate",
                             "loan_percent_income",
                             "credit_score"),         # or a character/integer vector of variable names/indices
       label_offset_dist = c(0,0.5,0,0),
       ticks_v = 5)
}


# ===========================================================
# PHASE 3 — Local interpretation (using reduced model bl_results_v2)
# ===========================================================

# ---- Step 9: Inspect predictions and select target --------------------
# Review pred_summary to choose a target. False Negatives (predicted 0,
# true 1) are often the most interesting: the applicant will default but the
# model clears them -- what would need to change to flag them correctly?
# Here tdp = 13 is a true negative (predicted 0.111, did not default): what
# would make the model predict default?
{
  pred_summary <- bl_predict(bl_results_v2)

  tdp <- 13
  
  # A subset of points for the biplot: bl_project_points() converts rows of a
  # data frame (here test_data rows 13 down to 10) into the format plot()
  # takes. With plot_points = FALSE below they are not drawn; their prediction
  # summary is still printed.
  test_points <- bl_project_points(bl_results_v2$test_data[tdp:10,], bl_results_v2)

  target_value <- bl_results_v2$test_data[tdp, ]
  # Highlight the target on the main biplot

  plot(
    bl_results_v2,
    points       = test_points, # the points to plot; if empty, the training data
    plot_points  = FALSE, # removes the points from the plot
    target_point = target_value,
    target_label = tdp,
    label_dir         = "Paral",   # default; "Hor" and "Orthog" also accepted
    label_offset_var  = c("loan_amnt",
                          "loan_int_rate",
                          "loan_percent_income",
                          "credit_score"),         # or a character/integer vector of variable names/indices
    label_offset_dist = c(0,0.5,0,0),
    ticks_v = 3
  )

  # Interactive: click the target (test_data row 13) to draw the arrow from it
  # to its global counterfactual, the nearest boundary point of this biplot
  # (bl_bnd_v2, Step 7b), and print the row. Press Esc to finish. bl_bnd_v2
  # holds one row per test_data row, so pick from test_data too. Skipped when
  # the script is source()d non-interactively.
  if (interactive())
    bl_pick_point(bl_results_v2, data = bl_results_v2$test_data,
                  bl_boundary = bl_bnd_v2)
}





# ---- Step 10: Set actionability constraints ---------------------------
# Omit any pruned variables from set_filters() — they are no longer in
# the model and passing them will raise an error.
#
# Valid constraint types:
#   "decrease"  — counterfactual value must be <= observed
#   "increase"  — counterfactual value must be >= observed
#   "fixed"     — the search allows +/- 0.5 of observed, then restores the
#                 observed value; if the restored point no longer crosses the
#                 boundary or leaves the training ranges, no counterfactual is
#                 reported
#   c(min, max) — counterfactual must lie within this range, given in model
#                 units (z-scores for bl_results_v2)
#   "integer"   — the value must be a whole number in raw units (a count or a
#                 score); combine it with a direction, c("increase", "integer"),
#                 or a range, list(c(min, max), "integer"). Without it a
#                 feature is continuous.
{
  # bl_select_target() reads the target's features AND its true class
  # (tgt$class_obs) from the same row -- here test_data[tdp, ]. The local
  # biplot (Step 12) colours the target by its confusion group, class_obs
  # against the model's predicted class: TP red, TN blue, FP purple, FN orange.
  tgt <- bl_select_target(bl_results_v2, target = tdp)

  # To explain a training row instead, index train_data; class_obs then comes
  # from that same training row:
  # tgt <- bl_select_target(bl_results_v2, target = 5L,
  #                         data = bl_results_v2$train_data)

  flt <- set_filters(
    tgt,
    #loan_percent_income     = "fixed",       #
    loan_amnt      = "fixed",    # make not actionable
    loan_int_rate  = "increase",        # set by the lender, not the applicant
    credit_score = c("increase","integer")   # a whole-number score
  )
  print(flt)
}


# ---- Step 11: Find local counterfactual via biplot rotation ----------
{
  # hull = TRUE (the default) keeps only boundary points inside the convex
  # hull of the training data and the target on each rotated plane, the
  # thesis's plausibility rule; hull = FALSE drops it.
  bl_local <- bl_find_local_cf(
    bl_result   = bl_results_v2,
    set_filters = flt,
    bl_target   = tgt,
    hull        = TRUE
    )
  # print(bl_local) ends with the observed vs counterfactual table (one row
  # per feature, with its Step 10 constraint, in raw loan units) and notes
  # explaining a "fixed" feature the search moved and restored, and an
  # "integer" feature rounded to a whole value (credit_score here). When the
  # restore or the rounding fails, no counterfactual is reported and the
  # Reason line says why.
  print(bl_local)
}


# ---- Step 12: Local biplot plot ----------------------------------------
{
  plot(bl_local,
       label_dir         = "Paral",   # default; "Hor" and "Orthog" also accepted
       label_offset_var  = c("loan_amnt",
                             "loan_int_rate",
                             "loan_percent_income",
                             "credit_score"),         # or a character/integer vector of variable names/indices
       label_offset_dist = c(0,0,1,0),
       ticks_var = c("loan_amnt",
                     "loan_int_rate",
                     "loan_percent_income",
                     "credit_score"),
       ticks_n = c(4,5,5,5)
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
  bl_sparse <- bl_find_sparse_cf(bl_shapley_values)
  print(bl_sparse)
  plot(bl_sparse,
       new_title = "Sparse counterfactual, test row 13",
       label_offset_var  = c("loan_amnt",
                             "loan_int_rate",
                             "loan_percent_income",
                             "credit_score"),         # or a character/integer vector of variable names/indices
       label_offset_dist = c(0,0,1,0),
       ticks_var = c("loan_amnt",
                     "loan_int_rate",
                     "loan_percent_income",
                     "credit_score"),
       ticks_n = c(4,5,5,5)
  )
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
       label_offset_var  = c("loan_amnt",
                             "loan_int_rate",
                             "loan_percent_income",
                             "credit_score"),         # or a character/integer vector of variable names/indices
       label_offset_dist = c(0,0,0,1),
       ticks_var = c("loan_amnt",
                     "loan_int_rate",
                     "loan_percent_income",
                     "credit_score"),
       ticks_n = c(4,10,10,10)
  )

  bl_shap_free   <- bl_shapley(bl_local_free)
  bl_sparse_free <- bl_find_sparse_cf(bl_shap_free)
  print(bl_sparse_free)
  plot(bl_shap_free)
  plot(bl_sparse_free,
       label_dir         = "Paral",   # default; "Hor" and "Orthog" also accepted
       label_offset_var  = c("loan_amnt",
                             "loan_int_rate",
                             "loan_percent_income",
                             "credit_score"),         # or a character/integer vector of variable names/indices
       label_offset_dist = c(0,0,0,1),
       ticks_var = c("loan_amnt",
                     "loan_int_rate",
                     "loan_percent_income",
                     "credit_score"),
       ticks_n = c(4,10,10,10)
  )
}


# ===========================================================
# EXTERNAL APPLICANT
# Analyse a new loan application not in the training data.
# Only include variables in feature_cols_v2 (remove any pruned variables).
# Give the values in raw loan units; bl_apply_scaling() converts them.
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
  # Drop any columns pruned in "Prune least-important variables"
  new_applicant <- new_applicant[, intersect(names(new_applicant), feature_cols_v2),
                                  drop = FALSE]
  # bl_set_scaling() is display-only: bl_results_v2's model works in z-scores
  # (Steps 4b-6b), so convert the raw application to the model's units first.
  new_applicant <- bl_apply_scaling(new_applicant, bl_results_v2)


  # Local search — no actionability constraints for this new applicant.
  # new_applicant has no class column, so tgt_ext$class_obs is NA and the
  # local biplot colours the target by its predicted class only (red = 1,
  # blue = 0). A one-row frame that does carry a "class" column gets its
  # confusion colour (TP/TN/FP/FN) instead.
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
# four variables in feature_cols_v2. Ask SHAP the same question on the full
# nine-feature model: rank by mean absolute SHAP value and print the two full
# rankings side by side (SHAP against distance-to-boundary).
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


# ---- Global SHAP of the retained four variables ----------------------------
# Run SHAP on the reduced model -- the four variables the workflow proceeded
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

plot(bl_bnd_v2, type = "boxplot")

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


