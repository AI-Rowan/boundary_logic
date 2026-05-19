# =============================================================================
# 04_loan_wrap_data_demo.R
# Phase 1 demonstration — using bl_wrap_data() instead of bl_prepare_data()
#
# PURPOSE
# -------
# bl_prepare_data() handles three tasks in one call:
#   1. Feature column selection
#   2. Binary class encoding
#   3. Seeded train/test split
#
# bl_wrap_data() is the alternative entry point for users who have already
# split their data (e.g. from a cross-validation framework, an external
# pipeline, or a pre-defined holdout set). It skips the split and simply
# validates and packages whatever train/test frames you supply.
#
# This script performs Steps 1–7 of the loan-default Phase 1 workflow using
# bl_wrap_data(). The seed, fraction, and feature set match script 03 exactly
# so the biplot outputs are directly comparable.
#
# Run: devtools::load_all() then source this file, or open it in RStudio
#      and run sections interactively.
# =============================================================================

# ---- Load the package (development workflow) ---------------------------
{
  rm(list = ls())
  devtools::load_all()
}

# ---- Step 1: Load and integer-encode categorical variables -----------------
# Identical to script 03. Raw CSV contains character columns; these must be
# integer-encoded before any boundarylogic function is called.
{
  loan_raw <- read.csv(
    system.file("extdata", "loan_data.csv", package = "boundarylogic")
  )

  loan_encoded <- loan_raw

  gender_map    <- c(female = 0, male = 1)
  education_map <- c("High School" = 0, "Associate" = 1,
                     "Bachelor" = 2, "Master" = 3, "Doctorate" = 4)
  ownership_map <- c("MORTGAGE" = 0, "OTHER" = 1, "OWN" = 2, "RENT" = 3)
  intent_map    <- c("DEBTCONSOLIDATION" = 0, "EDUCATION" = 1,
                     "HOMEIMPROVEMENT" = 2, "MEDICAL" = 3,
                     "PERSONAL" = 4, "VENTURE" = 5)
  defaults_map  <- c("No" = 0, "Yes" = 1)

  loan_encoded$person_gender                  <- gender_map[loan_encoded$person_gender]
  loan_encoded$person_education               <- education_map[loan_encoded$person_education]
  loan_encoded$person_home_ownership          <- ownership_map[loan_encoded$person_home_ownership]
  loan_encoded$loan_intent                    <- intent_map[loan_encoded$loan_intent]
  loan_encoded$previous_loan_defaults_on_file <- defaults_map[loan_encoded$previous_loan_defaults_on_file]

  # Coerce all columns to double so the training data and prediction grid
  # share the same type — read.csv() reads integer-valued columns as integer,
  # but bl_build_grid() generates sequences as double.
  loan_encoded[] <- lapply(loan_encoded, as.numeric)

  vars_to_remove <- c("loan_intent", "person_education", "person_home_ownership")
  loan_encoded   <- loan_encoded[, setdiff(names(loan_encoded), vars_to_remove), drop = FALSE]
}


# ---- Step 2: Domain filter — keep applicants without prior defaults --------
# Identical to script 03. Retain only applicants with no prior default on
# file. previous_loan_defaults_on_file is then constant (all 0) within the
# filtered set, so it is excluded from feature_cols.
{
  loan_filtered <- loan_encoded[loan_encoded$previous_loan_defaults_on_file == 0, ]

  feature_cols <- setdiff(
    names(loan_filtered),
    c("loan_status", "previous_loan_defaults_on_file")
  )

  cat("Retained rows:", nrow(loan_filtered), "\n")
  cat("Features:", paste(feature_cols, collapse = ", "), "\n")
}


# ---- Step 3: Manual train/test split (replaces bl_prepare_data) -----------
# bl_prepare_data() does three things internally:
#   1. Selects feature columns
#   2. Renames the target column to "class" (required by bl_wrap_data)
#   3. Performs a seeded random split at train_fraction
#
# Here we do those steps explicitly so bl_wrap_data() can accept the result.
# Use the same seed (121L) and fraction (0.8) as script 03 for comparability.
{
  # Select features + target; rename target to "class"
  data_clean <- loan_filtered[, c(feature_cols, "loan_status")]
  names(data_clean)[names(data_clean) == "loan_status"] <- "class"
  # "class" must be numeric 0/1 — loan_status is already encoded that way

  # Seeded 80/20 split
  set.seed(121L)
  n         <- nrow(data_clean)
  train_idx <- sample(n, size = floor(0.8 * n), replace = FALSE)
  train_data <- data_clean[ train_idx, , drop = FALSE]
  test_data  <- data_clean[-train_idx, , drop = FALSE]
  rownames(train_data) <- NULL
  rownames(test_data)  <- NULL

  cat("Train rows:", nrow(train_data), "  Test rows:", nrow(test_data), "\n")
}


# ---- Step 4: Wrap pre-split data with bl_wrap_data() ----------------------
# bl_wrap_data() validates that:
#   - train_data (and test_data if supplied) are data frames
#   - both contain a column named "class" with only numeric 0/1 values
#   - var_names (if supplied) exist in both frames
# It returns a bl_data object — structurally identical to bl_prepare_data()
# output and compatible with all downstream Phase 1 functions.
#
# target_class = NULL because loan_status was already encoded as 0/1 before
# the split; no further class conversion is needed.
{
  bl_dat <- bl_wrap_data(
    train_data   = train_data,
    test_data    = test_data,
    var_names    = feature_cols,
    target_class = NULL
  )
  print(bl_dat)
}


# ---- Step 5: Outlier filter ------------------------------------------------
# bl_wrap_data() only validates and packages the pre-split data — it does no
# cleaning. bl_filter_outliers() is still needed because:
#
#   - A manual split does not remove outliers; it only divides rows.
#   - Outliers in the training set distort the biplot projection and shift the
#     decision boundary away from the bulk of the data.
#   - The convex hull polygon produced here flows into bl_build_grid(), where
#     it clips the prediction contours to the observed-data region.
#   - hull_fraction = 0.9 trims the outermost ~10 % of training points;
#     hull_fraction = 1 retains all points but still builds the polygon.
#
# bl_filter_outliers() accepts either a bl_data object (from bl_prepare_data
# or bl_wrap_data) or a bl_filter_result — both have identical structure.
# The call below is therefore identical to script 03.
{
  bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 0.9)
}


# ---- Step 6: Fit XGBoost model --------------------------------------------
{
  bl_mod <- bl_fit_model(
    train_data = bl_filt$train_data,
    var_names  = bl_filt$var_names,
    model_type = "XGB"
  )
  print(bl_mod)
}


# ---- Step 7: Build biplot result and plot ----------------------------------
{
  bl_results <- bl_build_result(
    bl_data  = bl_filt,
    bl_model = bl_mod,
    method   = "CVA",
    title    = "Loan default — bl_wrap_data() demo, XGB + CVA biplot",
    rounding = 3L
  )
  print(bl_results)

  # Reference biplot — training data coloured by confusion category
  plot_biplotEZ(
    bl_results,
    label_dir         = "Hor",
    label_offset_var  = 0L,
    label_offset_dist = 0.5
  )

  # Project test data and overlay on biplot
  test_pts <- bl_project_points(bl_results$test_data, bl_results)
  cat("Test points outside polygon:", sum(!test_pts$inside_polygon), "\n")
  plot_biplotEZ(bl_results, points = test_pts)
}
