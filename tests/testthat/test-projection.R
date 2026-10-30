bl_dat <- bl_prepare_data(iris, class_col = "Species",
                          target_class = "versicolor")

test_that("bl_build_projection returns bl_projection class for PCA", {
  result <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                                method = "PCA")
  expect_s3_class(result, "bl_projection")
})

test_that("V is a square matrix with nrow = ncol = num_vars", {
  result <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                                method = "PCA")
  p <- bl_dat$num_vars
  expect_equal(dim(result$V), c(p, p))
})

test_that("tV equals solve(V) to numerical tolerance", {
  result <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                                method = "PCA")
  expect_equal(result$tV, solve(result$V), tolerance = 1e-10)
})

test_that("X_sd matches apply(X, 2, sd) for PCA with standardise=TRUE", {
  X      <- as.matrix(bl_dat$train_data[, bl_dat$var_names])
  result <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                                method = "PCA", standardise = TRUE)
  expect_equal(result$X_sd, apply(X, 2, sd), tolerance = 1e-10)
})

test_that("CVA forces standardise = FALSE", {
  pred_col <- as.factor(ifelse(bl_dat$train_data$class == 1, "red", "blue"))
  result <- bl_build_projection(
    bl_dat$train_data, bl_dat$var_names,
    method = "CVA", standardise = TRUE,  # user passes TRUE
    cva_classes = pred_col
  )
  expect_false(result$standardise)
})

test_that("CVA with cva_classes=NULL throws informative error", {
  expect_error(
    bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                        method = "CVA", cva_classes = NULL),
    "method = 'CVA' requires"
  )
})

test_that("PCA round-trip projection stays in correct space", {
  result <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                                method = "PCA", standardise = FALSE)
  X      <- as.matrix(bl_dat$train_data[1L, bl_dat$var_names, drop = FALSE])
  X_c    <- sweep(X, 2L, result$X_center, "-")
  Z      <- X_c %*% result$V[, result$proj_dims, drop = FALSE]
  X_back <- sweep(Z %*% result$tV[result$proj_dims, , drop = FALSE],
                  2L, result$X_center, "+")
  # With only 2 of p=4 dims, this is an approximation — just check shape
  expect_equal(length(X_back), bl_dat$num_vars)
})

test_that("an empty confusion group is named, dropped, and keeps confusion colours", {
  # Predicts class 1 only for Petal.Length in [3.5, 4.4], where every iris is
  # versicolor: no false positives on any training split.
  fp_free <- bl_wrap_model(
    model = NULL, model_type = "custom", var_names = bl_dat$var_names,
    predict_fn = function(m, nd)
      as.numeric(nd$Petal.Length >= 3.5 & nd$Petal.Length <= 4.4))
  expect_message(
    res <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                               method = "CVA", bl_model = fp_free),
    "FP empty", fixed = TRUE)
  expect_identical(levels(res$cva_classes), c("FN", "TN", "TP"))
  expect_setequal(unique(res$point_col), c("green3", "steelblue", "red"))
})

test_that("an unused level in user-supplied cva_classes is dropped before W", {
  pl  <- bl_dat$train_data$Petal.Length
  cls <- factor(as.character(cut(pl, c(-Inf, 2.5, 4.8, Inf),
                                 labels = c("low", "mid", "high"))),
                levels = c("low", "mid", "high", "none"))
  res <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                             method = "CVA", cva_classes = cls)
  expect_identical(levels(res$cva_classes), c("low", "mid", "high"))
  expect_identical(res$metric_type, "W_cva")

  X <- as.matrix(bl_dat$train_data[, bl_dat$var_names])
  W <- Reduce(`+`, lapply(split(as.data.frame(X), droplevels(cls)),
                          function(d) crossprod(scale(as.matrix(d),
                                                      scale = FALSE))))
  expect_equal(res$metric, W / (nrow(X) - 3L), ignore_attr = TRUE)
})

test_that("confusion levels stay alphabetical, as as.factor() gives", {
  # The level order feeds biplotEZ::CVA(); changing it perturbs V.
  glm_mod <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
                          model_type = "GLM")
  res <- suppressMessages(
    bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                        method = "CVA", bl_model = glm_mod))
  expect_identical(levels(res$cva_classes),
                   sort(levels(res$cva_classes)))
})

test_that("four non-empty confusion groups give no message", {
  glm_mod <- bl_fit_model(bl_dat$train_data, bl_dat$var_names,
                          model_type = "GLM")
  pred   <- as.numeric(boundarylogic:::.pred_function(
    glm_mod$model, "GLM", bl_dat$train_data[, bl_dat$var_names]) >= 0.5)
  actual <- bl_dat$train_data$class
  groups <- c(sum(actual == 1 & pred == 1), sum(actual == 0 & pred == 0),
              sum(actual == 0 & pred == 1), sum(actual == 1 & pred == 0))
  skip_if(any(groups == 0L), "the GLM leaves a confusion group empty")
  expect_no_message(
    bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                        method = "CVA", bl_model = glm_mod))
})

test_that("print.bl_projection runs without error", {
  result <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                                method = "PCA")
  expect_output(print(result), "<bl_projection>")
})
