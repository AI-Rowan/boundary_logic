############################################################
# Internal prediction utilities
# Refactored from: scripts/2.2 Model_use fitting v2.R
# All global references replaced with explicit arguments.
############################################################

# --------------------------------------------------------------------------
# SVM probability extraction
# --------------------------------------------------------------------------

#' Extract class-1 probability from an e1071 SVM model
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
#' Applies floor-based rounding and dispatches to the correct prediction
#' method for each model type.
#'
#' @param model_use    Fitted model object. For XGB this must be a list with
#'   fields `model` (an `xgb.Booster`) and `features` (character vector of
#'   column names used during training).
#' @param model_type   Character scalar identifying the model family. One of
#'   `"GLM"`, `"GAM"`, `"GBM"`, `"LDA"`, `"SVM"`, `"NNET"`, `"RForrest"`,
#'   `"XGB"`.
#' @param rounding     Integer number of decimal places for floor-based
#'   rounding. Default `2`.
#' @param new_data     Data frame of observations to score.
#'
#' @return Numeric vector of predicted probabilities / scores, length
#'   `nrow(new_data)`, rounded via `floor(x * 10^rounding) / 10^rounding`.
#' @keywords internal
.pred_function <- function(model_use, model_type, rounding = 2L, new_data) {

  floor_round <- function(x) floor(x * (10^rounding)) / (10^rounding)

  pred_value <- switch(
    model_type,

    # GLM / GAM: type = "response" gives probabilities directly
    "GLM" = ,
    "GAM" = {
      floor_round(predict(model_use, newdata = new_data, type = "response"))
    },

    # GBM: requires n.trees; use the number the model was trained with
    "GBM" = {
      floor_round(
        predict(model_use, newdata = new_data,
                n.trees = model_use$n.trees, type = "response")
      )
    },

    # LDA: extract class-2 posterior
    "LDA" = {
      floor_round(predict(model_use, newdata = new_data)$posterior[, 2])
    },

    # SVM: use the internal probability helper
    "SVM" = {
      floor_round(.svm_pred(model_use, new_data))
    },

    # NNET: predict() returns probabilities directly (single-output net)
    "NNET" = {
      floor_round(predict(model_use, newdata = new_data))
    },

    # rpart (labelled RForrest in original code)
    "RForrest" = {
      floor_round(predict(model_use, newdata = new_data)[, 2])
    },

    # XGBoost: model_use is list(model = <booster>, features = <char vec>)
    "XGB" = {
      new_mat <- as.matrix(new_data[, model_use$features, drop = FALSE])
      preds   <- predict(model_use$model,
                         newdata = xgboost::xgb.DMatrix(new_mat))
      floor_round(preds)
    },

    # Unknown model type
    stop(sprintf(
      "Unknown model_type '%s'. Must be one of: GLM, GAM, GBM, LDA, SVM, NNET, RForrest, XGB.",
      model_type
    ), call. = FALSE)
  )

  as.numeric(pred_value)
}
