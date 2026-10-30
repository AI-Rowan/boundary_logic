# Phase 3: bl_select_target() resolves the target's true class (class_obs)
# from the same row that supplies its features, and plot.bl_local_result()
# colours the target from it.

lcf_bl_dat <- bl_prepare_data(iris, class_col = "Species",
                              target_class = "versicolor", verbose = FALSE)
lcf_bl_mod <- bl_fit_model(lcf_bl_dat$train_data, lcf_bl_dat$var_names,
                           model_type = "GLM")
lcf_bl_results <- suppressWarnings(
  bl_build_result(lcf_bl_dat, lcf_bl_mod, method = "PCA",
                  m = 60L, b_margin = 0.01))
lcf_train <- lcf_bl_results$train_data
lcf_test  <- lcf_bl_results$test_data

test_that("fixture: train_data[1] and test_data[1] differ in class", {
  expect_false(identical(lcf_train$class[1L], lcf_test$class[1L]))
})

test_that("an index into train_data gives train_data[1]'s class", {
  tgt <- bl_select_target(lcf_bl_results, target = 1L, data = lcf_train)
  expect_identical(tgt$class_obs, as.integer(lcf_train$class[1L]))
})

test_that("a one-row data-frame target gives the same class_obs", {
  tgt_idx <- bl_select_target(lcf_bl_results, target = 1L, data = lcf_train)
  tgt_df  <- bl_select_target(lcf_bl_results, target = lcf_train[1L, ])
  expect_identical(tgt_df$class_obs, tgt_idx$class_obs)
  expect_true(is.na(tgt_df$row_id))
})

test_that("the default (test_data) is unchanged", {
  tgt <- bl_select_target(lcf_bl_results, target = 1L)
  expect_identical(tgt$class_obs, as.integer(lcf_test$class[1L]))
  expect_identical(tgt$row_id, 1L)
  expect_equal(tgt$x_obs, lcf_test[1L, lcf_bl_results$var_names, drop = FALSE])
})

test_that("a data-frame target without a class column has class_obs NA", {
  x_only <- lcf_train[1L, lcf_bl_results$var_names, drop = FALSE]
  tgt    <- bl_select_target(lcf_bl_results, target = x_only)
  expect_identical(tgt$class_obs, NA_integer_)
})

test_that("plot.bl_local_result draws a train_data target and an old bl_target", {
  on.exit(grDevices::graphics.off(), add = TRUE)
  grDevices::pdf(NULL)
  tgt <- bl_select_target(lcf_bl_results, target = 1L, data = lcf_train)
  loc <- bl_find_local_cf(lcf_bl_results, tgt, verbose = FALSE)
  skip_if_not(isTRUE(loc$solution_found))
  expect_no_error(plot(loc, print_summary = FALSE))

  # A bl_target built before class_obs existed plots by prediction colour.
  loc_old <- loc
  loc_old$bl_target$class_obs <- NULL
  expect_no_error(plot(loc_old, print_summary = FALSE))
})


# ---- "fixed" restore after the search (.restore_fixed) ----------------------

# A minimal bl_result-like setup with a custom model whose probability depends
# on one feature only. The search point B_x has moved the fixed feature `b`
# from -0.3 to 0.1 (inside the 0.5 tolerance) and crosses the boundary.
rf_fixture <- function(depends_on = "a", b_range = c(-3, 3)) {
  predict_fn <- function(m, nd) stats::plogis(4 * nd[[depends_on]])
  list(
    bl_result = list(
      model        = list(model = NULL, predict_fn = predict_fn),
      model_type   = "custom",
      cutoff       = 0.5,
      X_center     = c(a = 0, b = 0),
      train_ranges = list(a = c(-3, 3), b = b_range)
    ),
    bl_target   = list(pred_class = 0L),
    x_obs       = c(-0.2, -0.3),
    var_names   = c("a", "b"),
    best_result = list(
      B_x    = data.frame(a = 0.1, b = 0.1),
      B_z    = matrix(c(0.1, 0.1), nrow = 1L, dimnames = list(NULL, c("x", "y"))),
      B_pred = 0.599,
      Vr_rot = diag(2)
    ),
    filters = structure(list(constraints = list(b = "fixed")), class = "bl_filters")
  )
}

run_restore <- function(f) {
  boundarylogic:::.restore_fixed(f$best_result, f$x_obs, f$var_names, f$filters,
                                 f$bl_result, f$bl_target, sv = c(1, 1))
}

squish <- function(x) gsub("\\s+", " ", paste(x, collapse = " "))

test_that("a restored point that still crosses becomes the counterfactual", {
  res <- run_restore(rf_fixture("a"))
  fc  <- res$fixed_check
  expect_true(fc$applied)
  expect_identical(fc$vars, "b")
  expect_equal(res$best_result$B_x$b, -0.3)
  expect_equal(res$best_result$B_x$a, 0.1)
  expect_equal(res$best_result$B_pred, floor(stats::plogis(0.4) * 1000) / 1000)
  expect_equal(as.numeric(res$best_result$B_z), c(0.1, -0.3))
  expect_equal(fc$search_point$B_x$b, 0.1)       # the search's point is kept
})

test_that("a restored point that does not cross is not used", {
  f   <- rf_fixture("b")
  res <- run_restore(f)
  expect_false(res$fixed_check$applied)
  expect_true(res$fixed_check$in_range)
  expect_identical(res$best_result, f$best_result)
})

test_that("a restored point outside train_ranges is not used", {
  res <- run_restore(rf_fixture("a", b_range = c(0, 3)))
  expect_false(res$fixed_check$applied)
  expect_false(res$fixed_check$in_range)
})

test_that("no moved fixed feature, or no fixed constraint, gives NULL", {
  f <- rf_fixture("a")
  f$best_result$B_x$b <- f$x_obs[2L]
  expect_null(run_restore(f))

  f <- rf_fixture("a")
  f$filters$constraints <- list(b = "increase")
  expect_null(run_restore(f))
})


# ---- Comparison tables and the fixed-move note --------------------------------

test_that(".cf_compare_table converts levels and the change, in var_names order", {
  sc  <- list(center = c(b = 10, a = 100), scale = c(b = 2, a = 5))
  cmp <- boundarylogic:::.cf_compare_table(
    x_obs = c(1, 2), cols = list(counterfactual = c(1.5, 2)),
    var_names = c("a", "b"), scaling = sc,
    constraints = list(b = "fixed", a = c(-1, 1)))
  expect_identical(cmp$units, "raw")
  expect_identical(rownames(cmp$table), c("a", "b"))
  expect_equal(cmp$table$observed,       c(1 * 5 + 100, 2 * 2 + 10))
  expect_equal(cmp$table$counterfactual, c(1.5 * 5 + 100, 2 * 2 + 10))
  expect_equal(cmp$table$change,         c(0.5 * 5, 0))
  expect_identical(cmp$table$constraint, c("[-1, 1]", "fixed"))
})

test_that(".cf_compare_table stays in model units and measures change from the last column", {
  cmp <- boundarylogic:::.cf_compare_table(
    x_obs = c(1, 2), cols = list(counterfactual = c(1.5, 2.5), sparse = c(1, 2.5)),
    var_names = c("a", "b"), scaling = NULL, constraints = NULL)
  expect_identical(cmp$units, "model")
  expect_named(cmp$table, c("observed", "counterfactual", "sparse", "change", "constraint"))
  expect_equal(cmp$table$change, c(0, 0.5))
  expect_identical(cmp$table$constraint, c("none", "none"))
})

test_that(".fixed_notes explains each outcome", {
  sc <- list(center = c(a = 0, b = 50), scale = c(a = 1, b = 10))
  f  <- rf_fixture("a")

  applied <- squish(boundarylogic:::.fixed_notes(
    run_restore(f)$fixed_check, f$x_obs, f$var_names, sc))
  expect_match(applied, 'b is "fixed", but the boundary search moved it to 51 (+4).',
               fixed = TRUE)
  expect_match(applied, "+/- 5 here", fixed = TRUE)
  expect_match(applied, "restored point is the selected counterfactual", fixed = TRUE)

  kept <- squish(boundarylogic:::.fixed_notes(
    run_restore(rf_fixture("b"))$fixed_check, f$x_obs, f$var_names, NULL))
  expect_match(kept, "+/- 0.5 here", fixed = TRUE)
  expect_match(kept, "does not cross the boundary", fixed = TRUE)

  out_of_range <- squish(boundarylogic:::.fixed_notes(
    run_restore(rf_fixture("a", b_range = c(0, 3)))$fixed_check,
    f$x_obs, f$var_names, NULL))
  expect_match(out_of_range, "training feature ranges", fixed = TRUE)

  expect_length(boundarylogic:::.fixed_notes(NULL, f$x_obs, f$var_names, NULL), 0L)
})


# ---- print(bl_local) and print(bl_sparse) on the iris pipeline --------------------

test_that("print(bl_local) ends with the observed vs counterfactual table", {
  tgt <- bl_select_target(lcf_bl_results, target = 1L)
  loc <- bl_find_local_cf(lcf_bl_results, tgt, verbose = FALSE)
  skip_if_not(isTRUE(loc$solution_found))
  expect_null(loc$fixed_check)

  txt       <- paste(capture.output(print(loc)), collapse = "\n")
  pos_pairs <- regexpr("Per-pair distances", txt, fixed = TRUE)
  pos_table <- regexpr("Observed vs counterfactual (model units)", txt, fixed = TRUE)
  expect_gt(pos_pairs, 0)
  expect_gt(pos_table, pos_pairs)
})

test_that("print(bl_local) shows no table when no solution is found", {
  tgt <- bl_select_target(lcf_bl_results, target = 1L)
  flt <- set_filters(tgt, Sepal.Length = c(1000, 1001))
  loc <- suppressMessages(
    bl_find_local_cf(lcf_bl_results, tgt, set_filters = flt, verbose = FALSE))
  expect_false(loc$solution_found)
  expect_null(loc$fixed_check)
  txt <- paste(capture.output(print(loc)), collapse = "\n")
  expect_match(txt, "Reason", fixed = TRUE)
  expect_no_match(txt, "Observed vs counterfactual", fixed = TRUE)
})

test_that("a fixed feature is restored or explained, and never 'used' in the sparse CF", {
  tgt       <- bl_select_target(lcf_bl_results, target = 1L)
  fixed_var <- lcf_bl_results$var_names[1L]
  flt <- do.call(set_filters, c(list(tgt), setNames(list("fixed"), fixed_var)))
  loc <- suppressMessages(
    bl_find_local_cf(lcf_bl_results, tgt, set_filters = flt, verbose = FALSE))
  skip_if_not(isTRUE(loc$solution_found))

  if (isTRUE(loc$fixed_check$applied))
    expect_equal(loc$B_x[[fixed_var]], tgt$x_obs[[fixed_var]])
  if (!is.null(loc$fixed_check))
    expect_match(paste(capture.output(print(loc)), collapse = "\n"), "Note:", fixed = TRUE)

  sp   <- bl_find_sparse_cf(bl_shapley(loc))
  used <- sp$shapley_df$used_in_sparse[as.character(sp$shapley_df$varnames) == fixed_var]
  expect_false(used)

  txt <- paste(capture.output(print(sp)), collapse = "\n")
  expect_match(txt, "Observed vs sparse counterfactual (model units)", fixed = TRUE)
  for (h in c("observed", "counterfactual", "sparse", "change", "constraint", "used"))
    expect_match(txt, h, fixed = TRUE)
})


# ---- Shapley "No effect" class and the greedy sparse CF -------------------------

# A minimal bl_local_result for features a and c, scored by the custom model
# f(a, c), with observed point x_obs and counterfactual B_x.
sp_local <- function(f, x_obs, B_x) {
  fl <- function(v) floor(v * 1000) / 1000
  p0 <- fl(f(x_obs[1L], x_obs[2L]))
  structure(list(
    solution_found = TRUE,
    bl_result   = list(var_names  = c("a", "c"),
                       model      = list(model = NULL,
                                         predict_fn = function(m, nd) f(nd$a, nd$c)),
                       model_type = "custom", cutoff = 0.5, scaling = NULL),
    bl_target   = list(x_obs = data.frame(a = x_obs[1L], c = x_obs[2L]),
                       pred_prob = p0, pred_class = as.integer(p0 >= 0.5),
                       row_id = 1L),
    B_x         = data.frame(a = B_x[1L], c = B_x[2L]),
    B_pred      = fl(f(B_x[1L], B_x[2L])),
    set_filters = NULL
  ), class = "bl_local_result")
}

test_that("a zero Shapley value is 'No effect' and stays at observed", {
  loc <- sp_local(function(a, c) stats::plogis(4 * a), c(-0.2, -0.3), c(0.1, 0.5))
  shp <- bl_shapley(loc)
  cls <- setNames(shp$shapley_df$Contribute, as.character(shp$shapley_df$varnames))
  expect_identical(cls[["c"]], "No effect")
  expect_identical(cls[["a"]], "Supports")

  sp <- bl_find_sparse_cf(shp)
  expect_equal(sp$x_sparse$c, -0.3)
  expect_equal(sp$x_sparse$a, 0.1)
  expect_false(sp$shapley_df$used_in_sparse[as.character(sp$shapley_df$varnames) == "c"])

  skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
  expect_no_error(ggplot2::ggplot_build(plot(shp)))
})

test_that("the greedy sparse CF stops once the largest change crosses", {
  loc <- sp_local(function(a, c) stats::plogis(4 * a + c), c(-0.2, -0.3), c(0.3, 0.2))
  sp  <- bl_find_sparse_cf(bl_shapley(loc))
  expect_identical(sp$sparse_path$variable, "a")
  expect_true(sp$sparse_path$crosses)
  expect_equal(sp$pred_sparse, 0.710)
  expect_true(sp$solution_valid)
  expect_equal(sp$x_sparse$c, -0.3)
  used <- setNames(sp$shapley_df$used_in_sparse, as.character(sp$shapley_df$varnames))
  expect_true(used[["a"]])
  expect_false(used[["c"]])
})

test_that("the greedy sparse CF adds features until the boundary is crossed", {
  loc <- sp_local(function(a, c) stats::plogis(a + c), c(-0.4, -0.4), c(0.35, 0.25))
  sp  <- bl_find_sparse_cf(bl_shapley(loc))
  expect_identical(sp$sparse_path$variable, c("a", "c"))
  expect_identical(sp$sparse_path$crosses, c(FALSE, TRUE))
  expect_equal(sp$sparse_path$pred, c(0.487, 0.645))
  expect_true(sp$solution_valid)

  txt <- paste(capture.output(print(sp)), collapse = "\n")
  expect_match(txt, "Path (Supports features, largest Shapley value first)", fixed = TRUE)
})

test_that("print(bl_shapley) uses the comparison-table layout", {
  loc <- sp_local(function(a, c) stats::plogis(a + c), c(-0.4, -0.4), c(0.35, 0.25))
  shp <- bl_shapley(loc)
  txt <- paste(capture.output(print(shp)), collapse = "\n")
  expect_match(txt, "Observed vs counterfactual, with Shapley values (model units)",
               fixed = TRUE)
  for (h in c("observed", "counterfactual", "change", "constraint", "shapley", "contribute"))
    expect_match(txt, h, fixed = TRUE)
  expect_match(txt, sprintf("Shapley sum    : %.4f", sum(shp$shapley_df$shapley_cause)),
               fixed = TRUE)
})


# ---- rotate_deg is for the global biplot only ------------------------------------

test_that("the local plots refuse rotate_deg before any device opens", {
  tgt <- bl_select_target(lcf_bl_results, target = 1L)
  loc <- bl_find_local_cf(lcf_bl_results, tgt, verbose = FALSE)
  skip_if_not(isTRUE(loc$solution_found))
  sp  <- bl_find_sparse_cf(bl_shapley(loc))

  dev_before <- grDevices::dev.cur()
  expect_error(plot(loc, rotate_deg = 30), "rotate_deg", fixed = TRUE)
  expect_error(plot(sp, rotate_deg = 30), "rotate_deg", fixed = TRUE)
  expect_identical(grDevices::dev.cur(), dev_before)
})

test_that("the local plots still draw, and the global biplot keeps rotate_deg", {
  on.exit(grDevices::graphics.off(), add = TRUE)
  grDevices::pdf(NULL)
  tgt <- bl_select_target(lcf_bl_results, target = 1L)
  loc <- bl_find_local_cf(lcf_bl_results, tgt, verbose = FALSE)
  skip_if_not(isTRUE(loc$solution_found))
  sp  <- bl_find_sparse_cf(bl_shapley(loc))

  expect_no_error(plot(loc, print_summary = FALSE))
  expect_no_error(capture.output(plot(sp)))
  expect_no_error(plot_biplotEZ(lcf_bl_results, rotate_deg = 30))
})
