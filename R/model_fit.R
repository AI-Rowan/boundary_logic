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
#'   keys are silently ignored. Common overrides:
#'   - GBM: `list(n.trees = 1000, shrinkage = 0.05)`
#'   - NNET: `list(size = 10, decay = 0.01)`
#'   - XGB: `list(max_depth = 5, nrounds = 300)`
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
#' bl_dat <- bl_prepare_data(iris,
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
                         rounding      = 2L,
                         model_params  = list()) {

  # ---- Validation -------------------------------------------------------
  stop_if_not_data_frame(train_data, "train_data")
  stop_if_not_character(var_names, "var_names")
  if (!"class" %in% names(train_data))
    stop("'train_data' must contain a column named 'class'.", call. = FALSE)
  stop_if_not_scalar_numeric(cutoff, "cutoff")
  stop_if_not_positive_integer(rounding, "rounding")

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
    rounding   = rounding,
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
      rounding   = rounding,
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
  cat(sprintf("  Train Accuracy: %.4f  |  Train Gini: %.4f\n",
              x$accuracy, x$gini))
  invisible(x)
}
