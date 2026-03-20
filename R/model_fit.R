############################################################
# bl_fit_model(): fit a classification model and report accuracy/Gini
# Refactored from: scripts/2.2 Model_use fitting v2.R
############################################################

#' Fit a classification model for boundary logic analysis
#'
#' Fits one of the supported ML classifiers on binary 0/1 training data
#' and returns the fitted model together with training accuracy and Gini
#' coefficient.
#'
#' @param train_data   Data frame with feature columns named by `var_names`
#'   plus a column named `"class"` (numeric 0/1).
#' @param var_names    Character vector of feature column names.
#' @param model_type   Character scalar specifying the model family. One of:
#'   `"GLM"` (logistic regression), `"GAM"` (generalised additive model),
#'   `"GBM"` (gradient boosting), `"LDA"` (linear discriminant analysis),
#'   `"SVM"` (support vector machine), `"NNET"` (neural network),
#'   `"RForrest"` (rpart decision tree), `"XGB"` (XGBoost).
#'   Default `"GLM"`.
#' @param cutoff       Numeric decision threshold for computing accuracy.
#'   Default `0.5`.
#' @param rounding     Integer; floor-based rounding precision for predicted
#'   probabilities. Default `2`.
#' @param model_params Named list of hyperparameter overrides. Unrecognised
#'   keys are silently ignored. Parameter names follow **parsnip** conventions
#'   for all types except `"GBM"` (which uses the gbm package directly).
#'   Common overrides:
#'   - GBM: `list(n.trees = 1000, shrinkage = 0.05)`
#'   - NNET: `list(hidden_units = 10, penalty = 0.01, epochs = 500)`
#'   - XGB: `list(trees = 300, tree_depth = 5, learn_rate = 0.05)`
#'   - RForrest: `list(min_n = 10)`
#'
#' @return A list of class `"bl_model"` with components:
#' \describe{
#'   \item{`model`}{The fitted model object.}
#'   \item{`model_type`}{Character; the model family used.}
#'   \item{`var_names`}{Character vector of feature names.}
#'   \item{`cutoff`}{The decision threshold.}
#'   \item{`rounding`}{The rounding precision.}
#'   \item{`accuracy`}{Numeric; classification accuracy on `train_data`.}
#'   \item{`gini`}{Numeric; Gini coefficient on `train_data`.}
#' }
#'
#' @examples
#' bl_dat <- bl_prepare_data(datasets::iris,
#'                            class_col    = "Species",
#'                            target_class = "versicolor")
#' bl_mod <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
#'                        model_type = "GLM")
#' print(bl_mod)
#'
#' @export
bl_fit_model <- function(train_data,
                         var_names,
                         model_type    = "GLM",
                         cutoff        = 0.5,
                         model_params  = list()) {

  # ---- Validation -------------------------------------------------------
  stop_if_not_data_frame(train_data, "train_data")
  stop_if_not_character(var_names, "var_names")
  if (!"class" %in% names(train_data))
    stop("'train_data' must contain a column named 'class'.", call. = FALSE)
  stop_if_not_scalar_numeric(cutoff, "cutoff")

  valid_types <- c("GLM", "GAM", "GBM", "LDA", "SVM", "NNET", "RForrest", "XGB")
  if (!model_type %in% valid_types)
    stop(sprintf("model_type '%s' is not supported. Choose from: %s.",
                 model_type, paste(valid_types, collapse = ", ")),
         call. = FALSE)

  # ---- Fit the model ----------------------------------------------------
  fit_result <- .fit_model(train_data, var_names, model_type, model_params)

  # ---- Compute training accuracy and Gini --------------------------------
  pred_prob <- .pred_function(
    model_use  = fit_result$model,
    model_type = model_type,
    new_data   = train_data[, var_names, drop = FALSE]
  )
  pred_class <- as.numeric(pred_prob >= cutoff)
  actual     <- train_data[["class"]]

  accuracy <- mean(pred_class == actual, na.rm = TRUE)
  gini     <- calc_gini(actual, pred_prob)

  # ---- Return -----------------------------------------------------------
  structure(
    list(
      model      = fit_result$model,
      model_type = model_type,
      var_names  = var_names,
      cutoff     = cutoff,
      accuracy   = accuracy,
      gini       = gini
    ),
    class = "bl_model"
  )
}


# --------------------------------------------------------------------------
# S3 print method
# --------------------------------------------------------------------------

#' @export
print.bl_model <- function(x, ...) {
  cat("<bl_model>\n")
  cat(sprintf("  Model type    : %s\n", x$model_type))
  cat(sprintf("  Cutoff        : %g\n", x$cutoff))
  if (!is.na(x$accuracy)) {
    cat(sprintf("  Train Accuracy: %.4f  |  Train Gini: %.4f\n",
                x$accuracy, x$gini))
  } else {
    cat("  Train Accuracy: (not computed)\n")
  }
  invisible(x)
}


# --------------------------------------------------------------------------
# bl_wrap_model(): register an externally fitted model as a bl_model
# --------------------------------------------------------------------------

#' Wrap an externally fitted model for use in boundary logic analysis
#'
#' Use this function when you have fitted your own classification model
#' outside of `bl_fit_model()` — for example, using custom hyperparameters,
#' a caret workflow, or a model type not natively supported. The wrapped
#' object is a `"bl_model"` compatible with `bl_build_grid()` and
#' `bl_assemble()`.
#'
#' @section Supported model types:
#' For standard types (`"GLM"`, `"GAM"`, `"GBM"`, `"LDA"`, `"SVM"`,
#' `"NNET"`, `"RForrest"`, `"XGB"`), the existing prediction dispatch is
#' used. Your model object must be compatible with the prediction call for
#' that type (e.g., a `glm` object for `"GLM"`).
#'
#' For `"custom"` model types, supply a `predict_fn` — a function with
#' signature `function(model, new_data)` that returns a numeric vector of
#' probabilities (one per row of `new_data`).
#'
#' @note For `"XGB"` models, `model` must be a list of the form
#'   `list(model = <xgb.Booster>, features = <character vector>)`, matching
#'   the format produced internally by `bl_fit_model()`.
#'
#' @param model        The fitted model object (or for `"XGB"`, the wrapper
#'   list described above).
#' @param model_type   Character scalar identifying the model family. One of
#'   the standard types or `"custom"`.
#' @param var_names    Character vector of feature column names used during
#'   training.
#' @param predict_fn   Function with signature `function(model, new_data)`
#'   returning a numeric probability vector. Required when
#'   `model_type = "custom"`; ignored otherwise.
#' @param train_data   Optional data frame (feature columns + `"class"`) used
#'   to compute training accuracy and Gini. If `NULL`, these are returned as
#'   `NA`.
#' @param cutoff       Numeric decision threshold. Default `0.5`.
#'
#'   **Note**: only `cutoff = 0.5` gives valid results in the current
#'   methodology. Other values are accepted for exploration but should be
#'   treated with caution.
#' @param rounding     Integer; floor-based rounding precision. Default `2`.
#'
#' @return A list of class `"bl_model"` compatible with `bl_build_grid()`
#'   and `bl_assemble()`.
#'
#' @examples
#' # Fit a GLM with custom settings, then wrap it
#' bl_dat <- bl_prepare_data(datasets::iris, class_col = "Species",
#'                            target_class = "versicolor")
#' my_glm <- glm(class ~ ., data = bl_dat$train_data,
#'               family = binomial(link = "logit"))
#' bl_mod <- bl_wrap_model(my_glm, model_type = "GLM",
#'                          var_names  = bl_dat$var_names,
#'                          train_data = bl_dat$train_data)
#'
#' # Custom model with predict_fn
#' my_fn  <- function(m, new_data) predict(m, new_data, type = "response")
#' bl_mod <- bl_wrap_model(my_glm, model_type = "custom",
#'                          var_names  = bl_dat$var_names,
#'                          predict_fn = my_fn,
#'                          train_data = bl_dat$train_data)
#'
#' @export
bl_wrap_model <- function(model,
                           model_type,
                           var_names,
                           predict_fn  = NULL,
                           train_data  = NULL,
                           cutoff      = 0.5) {

  stop_if_not_character(var_names, "var_names")
  stop_if_not_scalar_numeric(cutoff, "cutoff")

  valid_types <- c("GLM", "GAM", "GBM", "LDA", "SVM", "NNET",
                   "RForrest", "XGB", "custom")
  if (!model_type %in% valid_types)
    stop(sprintf(
      "model_type '%s' is not supported. Choose from: %s.",
      model_type, paste(valid_types, collapse = ", ")
    ), call. = FALSE)

  if (model_type == "custom" && is.null(predict_fn))
    stop("'predict_fn' is required when model_type = 'custom'.", call. = FALSE)

  if (model_type == "custom" && !is.function(predict_fn))
    stop("'predict_fn' must be a function.", call. = FALSE)

  # For custom type, bundle model and predict_fn together so that
  # .pred_function() can locate both.
  model_use <- if (model_type == "custom") {
    list(model = model, predict_fn = predict_fn)
  } else {
    model
  }

  # ---- Optionally compute accuracy and Gini from training data ----------
  if (!is.null(train_data)) {
    stop_if_not_data_frame(train_data, "train_data")
    if (!"class" %in% names(train_data))
      stop("'train_data' must contain a column named 'class'.", call. = FALSE)

    pred_prob  <- .pred_function(
      model_use  = model_use,
      model_type = model_type,
      new_data   = train_data[, var_names, drop = FALSE]
    )
    pred_class <- as.numeric(pred_prob >= cutoff)
    actual     <- train_data[["class"]]
    accuracy   <- mean(pred_class == actual, na.rm = TRUE)
    gini       <- calc_gini(actual, pred_prob)
  } else {
    accuracy <- NA_real_
    gini     <- NA_real_
  }

  structure(
    list(
      model      = model_use,
      model_type = model_type,
      var_names  = var_names,
      cutoff     = cutoff,
      accuracy   = accuracy,
      gini       = gini
    ),
    class = "bl_model"
  )
}
