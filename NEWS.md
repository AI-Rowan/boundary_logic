# boundarylogic (development version)

## Breaking changes

* `rounding` parameter replaced by `b_margin` across `bl_build_grid()`,
  `bl_build_result()`, and `bl_find_local_cf()`. `b_margin` sets the decision
  boundary contour band half-width *directly* (default `0.001`, equivalent to
  the old `rounding = 3L`), instead of indirectly via `b_margin = 1/10^rounding`.
  Model predictions are unaffected — they are always floor-rounded to 3 d.p.
* `no_points` renamed to `plot_points` in `plot_biplotEZ()`,
  `plot.bl_local_result()`, and `plot.bl_surrogate()`, with positive semantics
  (`TRUE` = show points). Defaults preserve prior visual behaviour for each
  function.
* `RForrest` model type renamed to `RForest` in `bl_fit_model()`,
  `bl_wrap_model()`, `.fit_model()`, and `.pred_function()`. The old spelling
  is still accepted at all four entry points but issues a warning and is
  redirected to `RForest`.
* `bl_fit_model()` now supports only four model types — GLM, SVM, NNET,
  RForest — all fitted via the parsnip/tidymodels backend. GAM, GBM, LDA, and
  XGB were removed from `bl_fit_model()`; use `bl_wrap_model()` for these (see
  its `predict_fn` contract for custom/GAM models, and the mandated
  `list(model = <xgb.Booster>, features = <character vector>)` form for XGB).
* `bl_prepare_data()` now performs outlier filtering internally via a new
  `hull_fraction` parameter (default `0.9`) and always returns a
  `"bl_filter_result"` object — the previous two-step
  `bl_prepare_data() |> bl_filter_outliers()` pattern is collapsed into one
  call. `bl_filter_outliers()` remains exported for advanced/iterative use
  (e.g. after `bl_wrap_data()`), but is no longer part of the standard
  workflow.

## New features

* Added Mahalanobis-aware distance measures, replacing per-feature Euclidean
  (divided by total SD) with squared Mahalanobis distance using the
  within-class covariance matrix `W`:
  - `bl_find_local_cf()` gains a `distance = c("mahalanobis", "euclidean")`
    parameter governing the *cross-pair* eigenvector-pair selector. New result
    fields `dist_mahalanobis`, `all_distances_mahalanobis`, `distance`.
  - `bl_robustness()` and `plot.bl_boundary()` gain the same `distance`
    parameter governing the per-feature standardisation denominator
    (`sqrt(diag(W))` instead of total SD `X_sd`), with a printed
    cross-correlation diagnostic and warning when the diagonal approximation
    captures less than ~75% of the full Mahalanobis distance.
  - `bl_build_projection()` now computes and stores the metric matrix and its
    Cholesky-based inverse (`metric`, `metric_inv`, `metric_type`), propagated
    through `bl_assemble()` into `bl_result`. See
    `documentation/mahalanobis_technical_note.md` for the full derivation.
  - Both new measures default to `"mahalanobis"`; pass `distance = "euclidean"`
    to reproduce the legacy behaviour exactly. Older `bl_result` objects
    without `metric_inv` fall back to Euclidean with a warning.
* Added a unified `plot()` interface for every biplot-producing object —
  `bl_result`, `bl_local_result`, `bl_sparse_result`, `bl_surrogate` — via a
  new `plot.bl_result()` S3 method and shared private helpers
  `.make_label_line_vec()` and `.apply_biplot_rotation()`.
* Added `ticks_var`/`ticks_n` parameters to `plot_biplotEZ()`,
  `plot.bl_local_result()`, and `plot.bl_surrogate()` for per-variable
  axis tick-mark count overrides (mirrors the existing
  `label_offset_var`/`label_offset_dist` pattern).
* `plot_biplotEZ()` gains `rotate_deg` (visual rotation), `label_cex`, and a
  `label_dir = "Paral"` default (border-adaptive label direction, replacing
  `"Hor"`); `label_offset_var` now accepts character variable names as well as
  integer indices.
* `bl_find_boundary()`'s out-of-range message reworded to clarify that
  identified counterfactuals outside the training feature ranges are
  *retained*, not discarded.

## Bug fixes and other improvements

* Fixed `.bl_rotate()`: when the selected eigenvector pair is not `(1, 2)`,
  the target point no longer collapses to the biplot origin. `Vrho`/`tVrho`
  are now sliced to columns/rows `(1, 2)` unconditionally, matching the SVD
  construction guarantee.
* Removed the unsuppressable zero-length-arrow warning from
  `plot_biplotEZ()`'s boundary overlay by moving counterfactual-arrow
  rendering exclusively into `bl_pick_point(bl_result, bl_boundary = ...)`,
  where it fires at most once per interactive click rather than once per
  observation at plot-flush time. The `boundary`, `show_arrows`, and
  `arrow_col` parameters were removed from `plot_biplotEZ()` accordingly.
* Swapped the filter order inside `bl_find_local_cf()`'s per-pair search:
  `set_filters` (actionability) now runs before `train_ranges` (feasibility),
  enabling an early per-segment exit on the more selective constraint.
* Added a condition-number check (`kappa(V, exact = FALSE)`) before
  `solve(V)` in `bl_build_projection()`, with a new `condition_number` field
  on the returned object and a warning above `1e10`.
* Clamped grid prediction probabilities to `[0, 1]` in `bl_build_grid()`,
  with a warning reporting how many cells were clamped.
* `bl_surrogate()` now reports `hull_coverage`, `n_in_hull`, and `n_total`,
  surfaced by `print.bl_surrogate()` alongside accuracy.
* `bl_shapley()` / `.shapley_perm_one()` default changed from `seed = 1L` to
  `seed = NULL` (only seeds the RNG when explicitly requested), so repeated
  calls are no longer silently deterministic.
* Guarded `bl_assemble()` against a `bl_model = NULL` crash; replaced two
  silent `invisible(NULL)` error paths in `bl_build_result()` with explicit
  `stop()` calls; guarded `print.bl_points()` against missing model
  predictions with a `"(no model)"` fallback; `.fit_model()` now warns on
  unknown hyperparameters passed via `model_params`.
* Fixed the `xgboost::xgboost()` → `xgboost::xgb.train()` API migration
  (`data` → `x`, `eta` → `learning_rate`, separate required `y`) in example
  scripts after an upstream xgboost release changed the interface.
* Numerous roxygen2 documentation corrections (e.g. `bl_predict()` rounding
  documented as 3 d.p., not 4; `bl_surrogate()`/`bl_build_projection()`
  `@return` fields synchronised with the objects they describe).

# boundarylogic 0.1.0 (2026-03-25)

* Initial implementation of the full three-phase Boundary Logic workflow:
  - **Phase 1** — data preparation (`bl_prepare_data()`, `bl_wrap_data()`,
    `bl_filter_outliers()`), model fitting/wrapping (`bl_fit_model()`,
    `bl_wrap_model()`), PCA/CVA biplot projection (`bl_build_projection()`),
    prediction-grid and decision-boundary contour construction
    (`bl_build_grid()`), and assembly into the central `bl_result` anchor
    object (`bl_assemble()`, `bl_build_result()`), with `plot_biplotEZ()` for
    visualisation.
  - **Phase 2** (global) — nearest-boundary counterfactual search
    (`bl_find_boundary()`), distance-to-boundary plots and per-variable
    robustness summaries (`plot.bl_boundary()`, `bl_robustness()`), and a
    spatial surrogate model (`bl_surrogate()`).
  - **Phase 3** (local) — single-observation target selection
    (`bl_select_target()`), actionability constraints (`set_filters()`),
    local counterfactual search via biplot rotation (`bl_find_local_cf()`),
    Shapley attribution of the counterfactual path (`bl_shapley()`), and
    sparse counterfactual generation (`bl_find_sparse_cf()`).
* Interactive biplot point picker `bl_pick_point()`.
* `bl_predict()` and `bl_project_points()` for scoring and projecting
  arbitrary data frames into the biplot's Z-space.
* Bundled the Pima Indians Diabetes and loan-default datasets, with
  accompanying vignettes covering the complete workflow end to end.
