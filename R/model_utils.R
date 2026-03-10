############################################################
# Internal model fitting dispatcher
# Refactored from: scripts/2.2 Model_use fitting v2.R
# Replaces the global boolean-flag pattern with a single model_type argument.
############################################################

# Default hyperparameters — match the original research scripts exactly
.default_model_params <- list(
  GBM  = list(shrinkage = 0.15, n.minobsinnode = 10L,
              n.trees = 500L, interaction.depth = 3L,
              distribution = "bernoulli"),
  NNET = list(size = 20L, decay = 0.001, maxit = 1000L,
              linout = FALSE, trace = FALSE),
  XGB  = list(max_depth = 3L, eta = 0.1, subsample = 0.9,
              colsample_bytree = 0.9, nrounds = 200L,
              objective = "binary:logistic", eval_metric = "logloss"),
  SVM  = list(probability = TRUE, type = "C-classification"),
  LDA  = list(),
  GLM  = list(),
  GAM  = list(),
  RForrest = list(method = "class", minsplit = 5L)
)


#' Fit a classification model (internal dispatcher)
#'
#' @param train_data  Data frame with feature columns + a column named
#'   `"class"` (numeric 0/1).
#' @param var_names   Character vector of feature column names.
#' @param model_type  One of `"GLM"`, `"GAM"`, `"GBM"`, `"LDA"`, `"SVM"`,
#'   `"NNET"`, `"RForrest"`, `"XGB"`.
#' @param model_params Named list of hyperparameter overrides. Unknown keys
#'   are silently ignored per model.
#'
#' @return A named list: `list(model = <fitted object>, model_type = <char>)`.
#' @keywords internal
.fit_model <- function(train_data, var_names, model_type, model_params = list()) {

  # Merge user params over defaults (user values win)
  defaults <- .default_model_params[[model_type]]
  if (is.null(defaults)) defaults <- list()
  params   <- utils::modifyList(defaults, model_params)

  # Build class ~ x1 + x2 + ... formula
  formula_str  <- paste("class ~", paste(var_names, collapse = " + "))
  formula_data <- stats::as.formula(formula_str)

  model <- switch(
    model_type,

    "GLM" = {
      stats::glm(formula_data, data = train_data,
                 family = stats::binomial(link = "logit"))
    },

    "GAM" = {
      gam_terms <- paste(sprintf("s(%s)", var_names), collapse = " + ")
      gam_form  <- stats::as.formula(paste("class ~", gam_terms))
      mgcv::gam(gam_form, data = train_data,
                family = stats::binomial(link = "logit"))
    },

    "GBM" = {
      gbm::gbm(
        formula_data,
        data              = train_data,
        distribution      = params$distribution,
        shrinkage         = params$shrinkage,
        n.minobsinnode    = params$n.minobsinnode,
        n.trees           = params$n.trees,
        interaction.depth = params$interaction.depth,
        verbose           = FALSE
      )
    },

    "LDA" = {
      MASS::lda(formula_data, data = train_data)
    },

    "SVM" = {
      e1071::svm(formula_data, data = train_data,
                 probability = params$probability,
                 type        = params$type)
    },

    "NNET" = {
      nnet::nnet(
        formula_data,
        data   = train_data,
        size   = params$size,
        decay  = params$decay,
        maxit  = params$maxit,
        linout = params$linout,
        trace  = params$trace
      )
    },

    "RForrest" = {
      train_rpart      <- train_data
      class_col_idx    <- which(names(train_rpart) == "class")
      train_rpart[, class_col_idx] <- as.factor(train_rpart[, class_col_idx])
      rpart::rpart(formula_data, data = train_data,
                   method   = params$method,
                   minsplit = params$minsplit,
                   model    = TRUE)
    },

    "XGB" = {
      X_train   <- as.matrix(train_data[, var_names, drop = FALSE])
      y_train   <- as.numeric(train_data[["class"]])
      dtrain    <- xgboost::xgb.DMatrix(data = X_train, label = y_train)
      xgb_param <- list(
        objective        = params$objective,
        eval_metric      = params$eval_metric,
        max_depth        = params$max_depth,
        eta              = params$eta,
        subsample        = params$subsample,
        colsample_bytree = params$colsample_bytree
      )
      booster <- xgboost::xgb.train(
        params  = xgb_param,
        data    = dtrain,
        nrounds = params$nrounds,
        verbose = 0L
      )
      # Return a wrapper so pred_function can align features at prediction time
      list(model = booster, features = colnames(X_train))
    },

    stop(sprintf(
      "Unknown model_type '%s'. Must be one of: GLM, GAM, GBM, LDA, SVM, NNET, RForrest, XGB.",
      model_type
    ), call. = FALSE)
  )

  list(model = model, model_type = model_type)
}
