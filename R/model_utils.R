############################################################
# Internal model fitting dispatcher — tidymodels edition
# Refactored from: scripts/2.2 Model_use fitting v2.R
#
# All model types except "GBM" are fitted as parsnip + workflows
# objects. "GBM" uses the gbm package directly because gbm is not
# a registered engine in parsnip.
############################################################

# Default hyperparameters.
# NNET, XGB, and RForrest use parsnip parameter names.
.default_model_params <- list(
  GBM      = list(shrinkage = 0.15, n.minobsinnode = 10L,
                  n.trees = 500L, interaction.depth = 3L,
                  distribution = "bernoulli"),
  NNET     = list(hidden_units = 20L, penalty = 0.001, epochs = 1000L),
  XGB      = list(trees = 200L, tree_depth = 3L, learn_rate = 0.1,
                  sample_size = 0.9, colsample_bytree = 0.9,
                  objective = "binary:logistic", eval_metric = "logloss"),
  SVM      = list(),
  LDA      = list(),
  GLM      = list(),
  GAM      = list(),
  RForrest = list(min_n = 5L)
)


#' Fit a classification model (internal dispatcher)
#'
#' @param train_data  Data frame with feature columns + a column named
#'   `"class"` (numeric 0/1).
#' @param var_names   Character vector of feature column names.
#' @param model_type  One of `"GLM"`, `"GAM"`, `"GBM"`, `"LDA"`, `"SVM"`,
#'   `"NNET"`, `"RForrest"`, `"XGB"`.
#' @param model_params Named list of hyperparameter overrides.
#'
#' @importFrom parsnip logistic_reg gen_additive_mod mlp svm_rbf decision_tree boost_tree set_engine set_mode fit
#' @importFrom workflows workflow add_formula add_model
#' @importFrom discrim discrim_linear
#' @importFrom gbm gbm
#' @importFrom stats as.formula
#' @importFrom utils modifyList
#' @return A named list: `list(model = <fitted object>, model_type = <char>)`.
#'   For all types except `"GBM"`, `model` is a fitted `workflows::workflow`.
#'   For `"GBM"`, `model` is a raw `gbm` object.
#' @keywords internal
.fit_model <- function(train_data, var_names, model_type, model_params = list()) {

  # Merge user params over defaults (user values win)
  defaults <- .default_model_params[[model_type]]
  if (is.null(defaults)) defaults <- list()
  params <- utils::modifyList(defaults, model_params)

  # ---- GBM: not parsnip-supported; use direct gbm call ------------------
  if (model_type == "GBM") {
    formula_data <- stats::as.formula(
      paste("class ~", paste(var_names, collapse = " + "))
    )
    model <- gbm::gbm(
      formula_data,
      data              = train_data,
      distribution      = params$distribution,
      shrinkage         = params$shrinkage,
      n.minobsinnode    = params$n.minobsinnode,
      n.trees           = params$n.trees,
      interaction.depth = params$interaction.depth,
      verbose           = FALSE
    )
    return(list(model = model, model_type = model_type))
  }

  # ---- tidymodels: factor outcome required for classification ------------
  train_tm       <- train_data
  train_tm$class <- factor(train_tm$class, levels = c("0", "1"))

  # GAM uses smooth terms in the formula; all others use linear terms
  if (model_type == "GAM") {
    gam_terms   <- paste(sprintf("s(%s)", var_names), collapse = " + ")
    formula_str <- paste("class ~", gam_terms)
  } else {
    formula_str <- paste("class ~", paste(var_names, collapse = " + "))
  }
  form <- stats::as.formula(formula_str)

  # ---- Build parsnip model specification --------------------------------
  model_spec <- switch(
    model_type,

    "GLM" = parsnip::logistic_reg() |>
      parsnip::set_engine("glm") |>
      parsnip::set_mode("classification"),

    "GAM" = parsnip::gen_additive_mod() |>
      parsnip::set_engine("mgcv") |>
      parsnip::set_mode("classification"),

    "LDA" = discrim::discrim_linear() |>
      parsnip::set_engine("MASS") |>
      parsnip::set_mode("classification"),

    "SVM" = parsnip::svm_rbf() |>
      parsnip::set_engine("kernlab") |>
      parsnip::set_mode("classification"),

    "NNET" = parsnip::mlp(
      hidden_units = params$hidden_units,
      penalty      = params$penalty,
      epochs       = params$epochs
    ) |>
      parsnip::set_engine("nnet", trace = FALSE) |>
      parsnip::set_mode("classification"),

    "RForrest" = parsnip::decision_tree(
      min_n = params$min_n
    ) |>
      parsnip::set_engine("rpart") |>
      parsnip::set_mode("classification"),

    "XGB" = parsnip::boost_tree(
      trees       = params$trees,
      tree_depth  = params$tree_depth,
      learn_rate  = params$learn_rate,
      sample_size = params$sample_size
    ) |>
      parsnip::set_engine(
        "xgboost",
        colsample_bytree = params$colsample_bytree,
        objective        = params$objective,
        eval_metric      = params$eval_metric,
        verbose          = 0L
      ) |>
      parsnip::set_mode("classification"),

    stop(sprintf(
      "Unknown model_type '%s'. Must be one of: GLM, GAM, GBM, LDA, SVM, NNET, RForrest, XGB.",
      model_type
    ), call. = FALSE)
  )

  # ---- Build workflow and fit --------------------------------------------
  wf <- workflows::workflow() |>
    workflows::add_formula(form) |>
    workflows::add_model(model_spec)

  fitted_wf <- parsnip::fit(wf, data = train_tm)

  list(model = fitted_wf, model_type = model_type)
}
