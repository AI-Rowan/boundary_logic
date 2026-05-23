############################################################
# Internal model fitting dispatcher — tidymodels edition
# Refactored from: scripts/2.2 Model_use fitting v2.R
#
# Supported types: GLM, SVM, NNET, RForrest (all parsnip/workflows).
# Other types (GBM, GAM, LDA, XGB) must be fitted externally and
# registered via bl_wrap_model().
############################################################

# Default hyperparameters (parsnip parameter names).
.default_model_params <- list(
  NNET     = list(hidden_units = 20L, penalty = 0.001, epochs = 1000L),
  SVM      = list(),
  GLM      = list(),
  RForrest = list(min_n = 5L)
)


#' Fit a classification model (internal dispatcher)
#'
#' @param train_data  Data frame with feature columns + a column named
#'   `"class"` (numeric 0/1).
#' @param var_names   Character vector of feature column names.
#' @param model_type  One of `"GLM"`, `"SVM"`, `"NNET"`, `"RForrest"`.
#' @param model_params Named list of hyperparameter overrides.
#'
#' @importFrom parsnip logistic_reg mlp svm_rbf decision_tree set_engine set_mode fit
#' @importFrom workflows workflow add_formula add_model
#' @importFrom stats as.formula
#' @importFrom utils modifyList
#' @return A named list: `list(model = <fitted object>, model_type = <char>)`.
#'   `model` is a fitted `workflows::workflow` for all supported types.
#' @keywords internal
.fit_model <- function(train_data, var_names, model_type, model_params = list()) {

  # Merge user params over defaults (user values win)
  defaults <- .default_model_params[[model_type]]
  if (is.null(defaults)) defaults <- list()
  params <- utils::modifyList(defaults, model_params)

  # ---- tidymodels: factor outcome required for classification ------------
  train_tm       <- train_data
  train_tm$class <- factor(train_tm$class, levels = c("0", "1"))

  formula_str <- paste("class ~", paste(var_names, collapse = " + "))
  form        <- stats::as.formula(formula_str)

  # ---- Build parsnip model specification --------------------------------
  model_spec <- switch(
    model_type,

    "GLM" = parsnip::logistic_reg() |>
      parsnip::set_engine("glm") |>
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

    stop(sprintf(
      "Unknown model_type '%s'. Must be one of: GLM, SVM, NNET, RForrest.",
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
