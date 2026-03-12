############################################################
# Internal prediction utilities — tidymodels edition
# Refactored from: scripts/2.2 Model_use fitting v2.R
#
# Primary path: models fitted by bl_fit_model() (except "GBM") are
# stored as workflows::workflow objects; prediction uses the unified
# predict(<workflow>, new_data, type = "prob") interface.
#
# Legacy / fallback path: raw model objects (GBM, or models passed
# via bl_wrap_model()) are handled by a model_type switch.
############################################################


# --------------------------------------------------------------------------
# SVM probability extraction — legacy fallback only
# Used when a raw e1071 SVM object is passed via bl_wrap_model().
# --------------------------------------------------------------------------

#' Extract class-1 probability from a raw e1071 SVM model
#'
#' @param model_use An e1071 SVM model trained with `probability = TRUE`.
#' @param new_data  Data frame for prediction.
#'
#' @return Numeric vector of class-1 probabilities.
#' @keywords internal
.svm_pred <- function(model_use, new_data) {
  preds_attr <- attr(
    predict(model_use, newdata = new_data, probability = TRUE),
    "probabilities"
  )
  prob_heading <- colnames(preds_attr)
  prob_nr      <- which(prob_heading == "1")[1L]
  if (is.na(prob_nr)) {
    stop("SVM model does not have a class labelled '1'. ",
         "Ensure the model was trained with a binary 0/1 response.",
         call. = FALSE)
  }
  preds_attr[, prob_nr]
}


# --------------------------------------------------------------------------
# Generic prediction wrapper
# --------------------------------------------------------------------------

#' Standardise prediction calls across all supported model families
#'
#' Dispatches to one of three paths:
#'
#' 1. **tidymodels workflow** (`inherits(model_use, "workflow")`): all model
#'    types fitted by `bl_fit_model()` except `"GBM"`. Uses
#'    `predict(<workflow>, new_data, type = "prob")` and extracts `.pred_1`.
#'
#' 2. **Legacy / raw model**: `"GBM"` (raw `gbm` object) and any model passed
#'    via `bl_wrap_model()` using a raw fitted object. Model-type-specific
#'    predict calls are used.
#'
#' 3. **Custom**: `model_type = "custom"` from `bl_wrap_model()`. Calls the
#'    user-supplied `predict_fn`.
#'
#' @param model_use    Fitted model. A `workflows::workflow` for tidymodels
#'   types; a raw `gbm` object for `"GBM"`; a list with `model` and
#'   `predict_fn` for `"custom"`.
#' @param model_type   Character scalar identifying the model family.
#' @param rounding     Integer; floor-based rounding decimal places. Default `2`.
#' @param new_data     Data frame of observations to score.
#'
#' @importFrom xgboost xgb.DMatrix
#' @return Numeric vector of predicted probabilities, length `nrow(new_data)`,
#'   rounded via `floor(x * 10^rounding) / 10^rounding`.
#' @keywords internal
.pred_function <- function(model_use, model_type, rounding = 2L, new_data) {

  floor_round <- function(x) floor(x * (10^rounding)) / (10^rounding)

  pred_value <- if (model_type == "custom") {
    # User-supplied predict_fn from bl_wrap_model()
    floor_round(model_use$predict_fn(model_use$model, new_data))

  } else if (inherits(model_use, "workflow")) {
    # tidymodels workflow: unified predict interface
    # predict() returns a tibble with .pred_0 and .pred_1
    preds <- predict(model_use, new_data = new_data, type = "prob")
    floor_round(preds$.pred_1)

  } else {
    # Legacy: raw model object (GBM fitted by bl_fit_model(), or a raw
    # model object passed via bl_wrap_model())
    switch(
      model_type,

      "GBM" = floor_round(
        predict(model_use, newdata = new_data,
                n.trees = model_use$n.trees, type = "response")
      ),

      "GLM" = ,
      "GAM" = floor_round(
        predict(model_use, newdata = new_data, type = "response")
      ),

      "LDA" = floor_round(
        predict(model_use, newdata = new_data)$posterior[, 2]
      ),

      "SVM" = floor_round(.svm_pred(model_use, new_data)),

      "NNET" = floor_round(
        predict(model_use, newdata = new_data)
      ),

      "RForrest" = floor_round(
        predict(model_use, newdata = new_data)[, 2]
      ),

      "XGB" = {
        new_mat <- as.matrix(new_data[, model_use$features, drop = FALSE])
        floor_round(predict(model_use$model,
                            newdata = xgboost::xgb.DMatrix(new_mat)))
      },

      stop(sprintf(
        "Unknown model_type '%s'. Must be one of: GLM, GAM, GBM, LDA, SVM, NNET, RForrest, XGB, custom.",
        model_type
      ), call. = FALSE)
    )
  }

  as.numeric(pred_value)
}
