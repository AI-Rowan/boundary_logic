# Tests for bl_set_scaling(), the raw-unit biplot axis relabel
# (.bl_rescale_biplot_axes()), and raw-unit feature-value display
# (.scale_to_raw() in the Shapley / sparse / target paths).

make_std_data <- function() {
  df <- iris
  df$class <- as.numeric(df$Species == "versicolor")
  df$Species <- NULL
  feats <- setdiff(names(df), "class")
  z <- scale(df[, feats])
  std <- df
  std[, feats] <- z
  list(raw = df, std = std, feats = feats,
       center = attr(z, "scaled:center"),
       scale  = attr(z, "scaled:scale"))
}

test_that("bl_set_scaling stores center/scale/method ordered to var_names", {
  d   <- make_std_data()
  bl  <- bl_wrap_data(d$std, d$std, var_names = d$feats)
  bl  <- bl_set_scaling(bl, center = d$center, scale = d$scale)
  expect_type(bl$scaling, "list")
  expect_named(bl$scaling, c("center", "scale", "method"))
  expect_equal(names(bl$scaling$center), d$feats)
  expect_equal(names(bl$scaling$scale),  d$feats)
  expect_equal(bl$scaling$method, "z-score")
})

test_that("bl_set_scaling reorders a named vector to var_names order", {
  d  <- make_std_data()
  bl <- bl_wrap_data(d$std, d$std, var_names = d$feats)
  shuffled <- d$center[rev(d$feats)]
  bl <- bl_set_scaling(bl, center = shuffled, scale = d$scale)
  expect_equal(bl$scaling$center, d$center[d$feats])
})

test_that("bl_set_scaling rejects bad inputs", {
  d  <- make_std_data()
  bl <- bl_wrap_data(d$std, d$std, var_names = d$feats)
  expect_error(bl_set_scaling(bl, center = d$center, scale = rep(0, 4)),
               "must not contain zero")
  expect_error(bl_set_scaling(bl, center = unname(d$center)[1:2], scale = d$scale),
               "length")
  bad <- d$center; names(bad)[1] <- "wrong"
  expect_error(bl_set_scaling(bl, center = bad, scale = d$scale),
               "missing entries")
  expect_error(bl_set_scaling(list(), center = d$center, scale = d$scale),
               "bl_data")
})

test_that("scaling propagates into the model-free bl_result", {
  d  <- make_std_data()
  bl <- bl_wrap_data(d$std, d$std, var_names = d$feats)
  bl <- bl_set_scaling(bl, center = d$center, scale = d$scale)
  proj <- bl_build_result(bl, method = "PCA")   # no model -> model-free bl_result
  expect_false(is.null(proj$scaling))
  expect_equal(proj$scaling$center, d$center[d$feats])
})

test_that("scaling works on the bl_prepare_data (filter) path", {
  d   <- make_std_data()
  df  <- cbind(d$std, Species = iris$Species)  # multiclass col for target_class
  bl  <- bl_prepare_data(df, class_col = "Species", target_class = "versicolor",
                         feature_cols = d$feats, hull_fraction = 0.9,
                         verbose = FALSE)
  expect_s3_class(bl, "bl_filter_result")
  bl  <- bl_set_scaling(bl, center = d$center, scale = d$scale)
  expect_equal(bl$scaling$scale, d$scale[d$feats])
  proj <- bl_build_result(bl, method = "PCA")
  expect_equal(proj$scaling$center, d$center[d$feats])
})

test_that("bl_filter_outliers preserves a scaling set beforehand", {
  d  <- make_std_data()
  bl <- bl_wrap_data(d$std, d$std, var_names = d$feats)
  bl <- bl_set_scaling(bl, center = d$center, scale = d$scale)
  filt <- bl_filter_outliers(bl, hull_fraction = 0.9, verbose = FALSE)
  expect_false(is.null(filt$scaling))
  expect_equal(filt$scaling$center, d$center[d$feats])
})

test_that(".bl_rescale_biplot_axes is a no-op when scaling is NULL", {
  d  <- make_std_data()
  bl <- bl_wrap_data(d$std, d$std, var_names = d$feats)
  proj <- bl_build_result(bl, method = "PCA")
  bp   <- proj$biplot_obj
  out  <- boundarylogic:::.bl_rescale_biplot_axes(bp, NULL, d$feats)
  expect_identical(out, bp)
})

test_that("raw-unit relabel keeps geometry and shifts labels to raw scale", {
  d  <- make_std_data()
  bl <- bl_wrap_data(d$std, d$std, var_names = d$feats)
  bl <- bl_set_scaling(bl, center = d$center, scale = d$scale)

  for (method in c("PCA", "CVA")) {
    # CVA on 2-class data emits known biplotEZ notices; not under test here.
    proj <- suppressWarnings(
      bl_build_result(bl, method = method, standardise = (method == "PCA")))
    bp     <- proj$biplot_obj
    bp_raw <- boundarylogic:::.bl_rescale_biplot_axes(bp, proj$scaling, d$feats)

    # Geometry untouched
    expect_equal(bp_raw$Z, bp$Z)
    expect_equal(bp_raw$Lmat, bp$Lmat)
    expect_equal(bp_raw$ax.one.unit, bp$ax.one.unit)

    # means/sd transformed by the affine rule
    expect_equal(unname(bp_raw$means),
                 unname(bp$means * d$scale[d$feats] + d$center[d$feats]))
    expect_equal(unname(bp_raw$sd),
                 unname(bp$sd * d$scale[d$feats]))

    # Axis labels now sit on the raw scale: midpoint near the raw column mean,
    # far from the standardised mean (~0).
    ac <- biplotEZ::axes_coordinates(bp_raw)
    for (j in seq_along(d$feats)) {
      lab_mid <- mean(range(ac[[j]][, 3]))
      expect_lt(abs(lab_mid - d$center[d$feats][j]),
                abs(lab_mid - 0) + 1e-8)
    }
  }
})


# ---- Raw-unit feature-value display (.scale_to_raw + Shapley/sparse/target) ----

test_that(".scale_to_raw converts levels and deltas correctly", {
  d  <- make_std_data()
  bl <- bl_wrap_data(d$std, d$std, var_names = d$feats)
  bl <- bl_set_scaling(bl, center = d$center, scale = d$scale)
  sc <- bl$scaling

  v <- c(1, 0, -1, 2)
  # level: v * scale + center
  expect_equal(
    boundarylogic:::.scale_to_raw(v, sc, d$feats, "level"),
    v * d$scale[d$feats] + d$center[d$feats], ignore_attr = TRUE)
  # delta: v * scale (no centre)
  expect_equal(
    boundarylogic:::.scale_to_raw(v, sc, d$feats, "delta"),
    v * d$scale[d$feats], ignore_attr = TRUE)
  # NULL scaling -> unchanged
  expect_identical(boundarylogic:::.scale_to_raw(v, NULL, d$feats, "level"), v)
  # incomplete scaling -> warns and returns unchanged
  sc_bad <- sc; sc_bad$center <- sc_bad$center[1:2]; sc_bad$scale <- sc_bad$scale[1:2]
  expect_warning(
    out <- boundarylogic:::.scale_to_raw(v, sc_bad, d$feats, "level"),
    "does not cover")
  expect_identical(out, v)
})

# Build a small local-CF -> Shapley -> sparse pipeline on standardised iris.
build_shap_pipeline <- function(with_scaling = TRUE) {
  d   <- make_std_data()
  bl  <- bl_wrap_data(d$std, d$std, var_names = d$feats)
  if (with_scaling)
    bl <- bl_set_scaling(bl, center = d$center, scale = d$scale)
  mod <- bl_fit_model(bl$train_data, bl$var_names, model_type = "GLM")
  res <- suppressWarnings(
    bl_build_result(bl, mod, method = "PCA", m = 60L, b_margin = 0.01))
  tgt <- bl_select_target(res, target = 1L)
  loc <- bl_find_local_cf(res, tgt)
  list(d = d, res = res, tgt = tgt, loc = loc)
}

test_that("print.bl_target reports raw-unit feature values when scaling is set", {
  p   <- build_shap_pipeline(with_scaling = TRUE)
  txt <- paste(capture.output(print(p$tgt)), collapse = "\n")
  expect_match(txt, "raw units")
  # First feature's raw observed value should be within the raw data range.
  raw_obs1 <- as.numeric(p$tgt$x_obs[[1L]]) * p$d$scale[1L] + p$d$center[1L]
  expect_gte(raw_obs1, min(p$d$raw[[p$d$feats[1L]]]) - 1e-6)
  expect_lte(raw_obs1, max(p$d$raw[[p$d$feats[1L]]]) + 1e-6)
})

test_that("Shapley plot label and prints convert to raw units", {
  skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
  p <- build_shap_pipeline(with_scaling = TRUE)
  skip_if_not(p$loc$solution_found, "no local CF found for this seed")

  shp    <- bl_shapley(p$loc, seed = 1L)
  sparse <- bl_find_sparse_cf(shp)

  # Plot: y-axis title flips to the raw observed->counterfactual wording, and
  # the label factor levels carry raw-scale numbers (not ~[-2,2] standardised).
  g       <- plot(shp)
  expect_match(g$labels$y, "raw units")
  lvls    <- levels(g$data$varnames_p)
  raw_obs <- as.numeric(shp$shapley_df$pred_data) *
    p$d$scale[as.character(shp$shapley_df$varnames)] +
    p$d$center[as.character(shp$shapley_df$varnames)]
  expect_true(any(grepl(as.character(round(raw_obs[1L], 3)), lvls, fixed = TRUE)))

  # Prints carry the raw-units note.
  expect_match(paste(capture.output(print(shp)),    collapse = "\n"), "raw units")
  expect_match(paste(capture.output(print(sparse)), collapse = "\n"), "raw units")
})

test_that("no scaling -> Shapley displays unchanged (model units)", {
  p <- build_shap_pipeline(with_scaling = FALSE)
  skip_if_not(p$loc$solution_found, "no local CF found for this seed")
  shp <- bl_shapley(p$loc, seed = 1L)
  g   <- plot(shp)
  expect_match(g$labels$y, "change required")          # original wording
  expect_no_match(paste(capture.output(print(shp)), collapse = "\n"), "raw units")
})


# ---- bl_local$bl_counterfactual (mirror of bl_target, raw-unit print) ----

test_that("bl_local carries a bl_counterfactual mirroring B_x", {
  p <- build_shap_pipeline(with_scaling = TRUE)
  skip_if_not(p$loc$solution_found, "no local CF found for this seed")
  cf <- p$loc$bl_counterfactual
  expect_s3_class(cf, "bl_counterfactual")
  # x_cf is the model-unit B_x (storage stays in model units, like bl_target).
  expect_equal(as.numeric(cf$x_cf), as.numeric(p$loc$B_x))
  expect_equal(cf$pred_prob, p$loc$B_pred)
})

test_that("print.bl_counterfactual shows raw units matching B_x*scale+center", {
  p <- build_shap_pipeline(with_scaling = TRUE)
  skip_if_not(p$loc$solution_found, "no local CF found for this seed")
  cf  <- p$loc$bl_counterfactual
  txt <- paste(capture.output(print(cf)), collapse = "\n")
  expect_match(txt, "raw units")
  vn      <- p$d$feats
  raw_cf  <- as.numeric(p$loc$B_x[, vn]) * p$d$scale[vn] + p$d$center[vn]
  # The first feature's raw value should appear (at 2 dp) in the output.
  expect_match(txt, formatC(round(raw_cf[1L], 2L), format = "f", digits = 2L),
               fixed = TRUE)
})

test_that("print(bl_local) shows both Target and Counterfactual blocks", {
  p   <- build_shap_pipeline(with_scaling = TRUE)
  skip_if_not(p$loc$solution_found, "no local CF found for this seed")
  txt <- paste(capture.output(print(p$loc)), collapse = "\n")
  expect_match(txt, "Target:")
  expect_match(txt, "Counterfactual:")
})

test_that("bl_counterfactual is NULL when no scaling -> model-unit print", {
  p <- build_shap_pipeline(with_scaling = FALSE)
  skip_if_not(p$loc$solution_found, "no local CF found for this seed")
  cf  <- p$loc$bl_counterfactual
  expect_s3_class(cf, "bl_counterfactual")
  txt <- paste(capture.output(print(cf)), collapse = "\n")
  expect_no_match(txt, "raw units")
})


# ---- bl_apply_scaling(): new raw-unit data into the model's units ----------

build_scaled_result <- function(with_scaling = TRUE) {
  d  <- make_std_data()
  bl <- bl_wrap_data(d$std, d$std, var_names = d$feats)
  if (with_scaling)
    bl <- bl_set_scaling(bl, center = d$center, scale = d$scale)
  mod <- bl_fit_model(bl$train_data, bl$var_names, model_type = "GLM")
  res <- suppressWarnings(suppressMessages(
    bl_build_result(bl, mod, method = "PCA", m = 40L, b_margin = 0.01)))
  list(d = d, bl = bl, res = res)
}

test_that("bl_apply_scaling converts raw features and leaves other columns alone", {
  s   <- build_scaled_result()
  out <- bl_apply_scaling(s$d$raw, s$bl)
  expect_equal(as.matrix(out[s$d$feats]), as.matrix(s$d$std[s$d$feats]),
               ignore_attr = TRUE)
  expect_identical(out$class, s$d$raw$class)

  # Column order does not matter, and a bl_result works like the bl_data.
  rev_raw <- s$d$raw[, rev(names(s$d$raw))]
  expect_equal(bl_apply_scaling(rev_raw, s$res)[s$d$feats], out[s$d$feats])
})

test_that("bl_apply_scaling rejects bad inputs", {
  s <- build_scaled_result()
  expect_error(bl_apply_scaling(as.matrix(s$d$raw), s$bl), "data.frame")
  expect_error(bl_apply_scaling(s$d$raw, list()), "must be a 'bl_result'")
  expect_error(bl_apply_scaling(s$d$raw[, -1L], s$bl), "missing feature columns")
  chr <- s$d$raw
  chr$Sepal.Length <- as.character(chr$Sepal.Length)
  expect_error(bl_apply_scaling(chr, s$bl), "must be numeric")
  unscaled <- bl_wrap_data(s$d$std, s$d$std, var_names = s$d$feats)
  expect_error(bl_apply_scaling(s$d$raw, unscaled), "No scaling")
})

test_that("a raw-unit data-frame target warns and suggests bl_apply_scaling", {
  s       <- build_scaled_result()
  raw_row <- s$d$raw[1L, s$d$feats]
  expect_warning(bl_select_target(s$res, target = raw_row), "bl_apply_scaling")
  expect_no_warning(bl_select_target(s$res,
                                     target = bl_apply_scaling(raw_row, s$res)))
  expect_no_warning(bl_select_target(s$res, target = 1L))
})

test_that("an out-of-range target without a scaling warns without the hint", {
  s   <- build_scaled_result(with_scaling = FALSE)
  row <- s$d$std[1L, s$d$feats]
  row$Sepal.Length <- 10
  expect_warning(bl_select_target(s$res, target = row), "Sepal.Length")
  expect_no_match(
    tryCatch(bl_select_target(s$res, target = row),
             warning = function(w) conditionMessage(w)),
    "bl_apply_scaling")
})
