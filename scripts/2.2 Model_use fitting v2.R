
#############################
# 1.1) SVM Probability Wrapper
#############################
#' Predict class-1 probabilities using an e1071 SVM model
#'
#' @param model_use Model object (e1071::svm) trained with probability=TRUE
#' @param data Newdata for prediction (data.frame)
#'
#' @return Numeric vector of class-1 probabilities
svm_pred <- function(model_use = model_use, data){
  prob_heading <- colnames(attr(predict(model_use, newdata = data, probability = TRUE), "probabilities"))
  prob_nr <- which(prob_heading == "1")[1]
  preds <- attr(predict(model_use, newdata = data, probability = TRUE), "probabilities")[, prob_nr]
  return(preds)
}




########################################
# 1.2) Generic Prediction Function Wrapper: Added XGBoost
########################################

suppressPackageStartupMessages({
  library(xgboost)   # NEW
  # ... your other libraries ...
})

#' Standardize prediction calls across model families
#'
#' Applies rounding and supports different model families.
#'
#' @param model_use Model object
#' @param model_select Character flag identifying model family
#'   (e.g., "NNET", "RForrest", "LDA", "SVM", "XGB").
#' @param rounding Integer number of decimal places (floor-based rounding)
#' @param new_data Data frame for prediction
#'
#' @return Numeric vector of predicted probabilities/scores
pred_function <-  function(model_use,
                           model_select,
                           rounding = 2,
                           new_data) {
  
  # Include XGB in the options so default GLM-style branch is skipped
  m_options <- c("NNET", "RForrest", "LDA", "SVM", "XGB")
  
  # Default (e.g., GLM, GAM): type="response" is probability
  if(!(model_select %in% m_options))
    pred_value <- floor(predict(model_use, newdata = new_data, type = "response") * (10^rounding)) / (10^rounding)
  
  if(model_select == "NNET")
    pred_value <- floor(predict(model_use, newdata = new_data) * (10^rounding)) / (10^rounding)
  
  if(model_select == "RForrest")
    pred_value <- floor(predict(model_use, newdata = new_data)[, 2] * (10^rounding)) / (10^rounding)
  
  if(model_select == "LDA")
    pred_value <- floor(predict(model_use, newdata = new_data)$posterior[, 2] * (10^rounding)) / (10^rounding)
  
  if(model_select == "SVM")
    pred_value <- floor(svm_pred(model_use, new_data) * (10^rounding)) / (10^rounding)
  
  # --- NEW: XGBoost (binary:logistic) returns probabilities directly
  if(model_select == "XGB") {
    # model_use is a small wrapper: list(model=<xgb.Booster>, features=<char>)
    new_mat <- as.matrix(new_data[, model_use$features, drop = FALSE])
    preds <- predict(model_use$model, newdata = xgb.DMatrix(new_mat))
    pred_value <- floor(preds * (10^rounding)) / (10^rounding)
  }
  
  return(pred_value)
}


#############################
# 2.1) Modeling Hints (Optional)
#############################
# You can specify the model directly and skip the generic parts by setting:
#   model_select <- "<MODEL>"  # e.g., "GLM", "LDA", "SVM", etc.
#   model_use    <- model_<model>



#################################
# 2.2) Model Fitting Wrapper
#################################
#' Fit a classification model according to selected flags
#'
#' Builds a formula class ~ x1 + x2 + ... and fits the chosen model(s).
#' The last enabled model flag determines the returned model.
#'
#' @param train_data Data frame containing features in columns var_names and
#'   a numeric binary 'class' column.
#'
#' @return list(model_use = <model object>, model_select = <character label>)
model_fitting <- function(train_data = train_data) {
  
  # Build formula: class ~ x1 + x2 + ...
  formula_data <- paste("class ~ ", var_names[1])
  for (i in 2:num_vars) {
    formula_data <- paste(formula_data, " + ", var_names[i])
  }
  formula_data <- as.formula(formula_data)
  
  # Base model: GLM (logistic)
  if (glm_use == TRUE) {  
    model_select <- "GLM"
    model_glm <- glm(formula_data, data = train_data, family = binomial(link = "logit"))
    model_use <- model_glm
  }
  
  # LDA (if enabled)
  if (lda == TRUE) {
    model_lda <- lda(formula_data, data = train_data)
    model_select <- "LDA"
    # model_qda <- qda(formula_data, data = train_data)
  }
  
  # Other models when LDA flag is FALSE
  if (lda == FALSE) {
    
    # SVM (if enabled)
    if (svm_u == TRUE)  {
      model_svm <- svm(formula_data, data = train_data, probability = TRUE, type = "C-classification")
      model_select <- "SVM"
    }
    
    # GBM (if enabled)
    if (gbm_used == TRUE) {
      model_gbm <- gbm(
        formula_data,
        data = train_data,
        distribution = "bernoulli",
        shrinkage = 0.15,
        n.minobsinnode = 10,
        n.trees = 500,
        interaction.depth = 3
      )
      model_select <- "GBM"
    }
    
    # rpart (if enabled)
    if (rpart_used == TRUE) {
      train_data_rpart <- train_data
      train_data_rpart[, num_vars + 1] <- as.factor(train_data_rpart[, num_vars + 1])
      model_rpart <- rpart(formula_data, data = train_data, method = "class", model = TRUE, minsplit = 5)
      model_select <- "RForrest"
    }
    
    # GAM (if enabled)
    if (gam_used == TRUE) {
      formula_data_gam <- paste("class ~ s(", var_names[1], ")")
      for (i in 2:num_vars) {
        formula_data_gam <- paste(formula_data_gam, " + s(", var_names[i], ")")
      }
      formula_data_gam <- as.formula(formula_data_gam)
      model_gam <- gam(formula_data_gam, data = train_data, family = binomial(link = "logit"))
      model_select <- "GAM"
    }
    
    # NNET (if enabled)
    if (nnet_used == TRUE) {
      library(nnet)
      model_nn <- nnet(
        formula_data,
        data = train_data,
        size = c(20),
        decay = 0.001,
        maxit = 1000,
        linout = FALSE,
        trace = FALSE
      )
      model_select <- "NNET"
    }
    
    # -----------------------
    # NEW: XGBoost (if enabled)
    # -----------------------
    if (xgb_used == TRUE) {
      # Build numeric design matrix from selected features
      X_train <- as.matrix(train_data[, var_names, drop = FALSE])
      y_train <- as.numeric(train_data[["class"]])  # expects 0/1
      
      dtrain <- xgb.DMatrix(data = X_train, label = y_train)
      
      xgb_params <- list(
        objective = "binary:logistic",
        eval_metric = "logloss",
        max_depth = 3,
        eta = 0.1,
        subsample = 0.9,
        colsample_bytree = 0.9
      )
      
      model_xgb_booster <- xgb.train(
        params  = xgb_params,
        data    = dtrain,
        nrounds = 200,
        verbose = 0
      )
      
      # Small wrapper to carry feature names for prediction-time alignment
      model_xgb <- list(
        model    = model_xgb_booster,
        features = colnames(X_train)
      )
      
      model_select <- "XGB"
    }
  }
  
  # Select final model object according to flags (last one wins)
  if (lda == TRUE)        model_use <- model_lda
  if (gam_used == TRUE)   model_use <- model_gam
  if (gbm_used == TRUE)   model_use <- model_gbm
  if (rpart_used == TRUE) model_use <- model_rpart
  if (svm_u == TRUE)      model_use <- model_svm
  if (nnet_used == TRUE)  model_use <- model_nn
  if (exists("model_xgb", inherits = FALSE) && xgb_used == TRUE) {
    model_use <- model_xgb
    model_select <- "XGB"
  }
  
  return(list(
    model_use    = model_use,
    model_select = model_select
  ))
}
