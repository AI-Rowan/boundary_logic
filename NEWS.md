# boundarylogic 0.2.0 (2026-11-07)

Changes since 0.1.0. The measured effect of each change to the method, on the
loan harness and on the Pima example of Paper 1, is in
`documentation/methodology-changes/methodology-changes.md`.

## Breaking changes

### Data and models

* `bl_prepare_data()` removes outliers itself, through the new argument
  `hull_fraction` (default `1`: no rows removed, since the thesis makes
  outlier removal optional), and returns a `"bl_filter_result"`. The two-step
  `bl_prepare_data() |> bl_filter_outliers()` is one call. Pass
  `hull_fraction = 0.9` to remove the outermost 10% of the training rows, as
  the thesis's loan example, the loan harness and the vignette do.
  `bl_filter_outliers()` keeps its default of `0.9` (calling it means
  filtering is wanted) and stays exported for use after `bl_wrap_data()`.
* `bl_fit_model()` fits four model types, all through parsnip: GLM, SVM, NNET
  and RForest (a decision tree, `rpart`). It no longer fits GAM, GBM, LDA or
  XGB: fit these yourself and pass them to `bl_wrap_model()` (see its
  `predict_fn` contract for custom and GAM models; an XGB model is passed as
  `list(model = <xgboost model>, features = <character vector>)`).
* The model type `RForrest` is renamed `RForest` in `bl_fit_model()`,
  `bl_wrap_model()`, `.fit_model()` and `.pred_function()`. The old spelling
  still works at all four, with a warning.
* `bl_build_grid()` and `bl_build_result()` replace `rounding` with
  `b_margin`, the half-width of the decision-boundary contour band, set
  directly (default `0.001`, which `rounding = 3L` gave).
  `bl_find_local_cf()` reads it from `bl_result$b_margin`. Model predictions
  are unaffected: they are always floor-rounded to 3 decimals.

### Plots

* `plot_biplotEZ()` loses `boundary`, `show_arrows` and `arrow_col`.
  Counterfactual arrows are drawn by
  `bl_pick_point(bl_result, bl_boundary = ...)`, one per click: drawing every
  arrow at plot time raised zero-length-arrow warnings that could not be
  suppressed.
* `no_points` is renamed `plot_points` in `plot_biplotEZ()` and
  `plot.bl_local_result()`, with the positive meaning (`TRUE` shows the
  points). The defaults keep each plot's look. (`plot.bl_surrogate()` has
  neither argument; it always draws its points.)
* `plot.bl_local_result()` and `plot.bl_sparse_result()` stop with an error
  when given `rotate_deg` (0.1.0 ignored it). The local biplot's orientation
  comes from the optimal rotation (`.bl_rotate()`), which is part of the
  method, and is never turned further. `rotate_deg` is for viewing the global
  biplot: `plot_biplotEZ()` and `plot.bl_surrogate()`.
* In `plot_biplotEZ()`, `plot.bl_local_result()` and `plot.bl_surrogate()`
  the axis labels default to `label_dir = "Paral"` (the direction adapts to
  the border; it was `"Hor"`), and `label_offset_dist` defaults to `1.5` (it
  was `0.5`).

### Search results

* `bl_find_local_cf()` applies the thesis's plausibility rule (section
  4.3.1): for each component pair it builds the convex hull of the training
  data on the rotated plane, and keeps only the boundary points inside it,
  before the actionability constraints, the training ranges and the model
  re-score. The global biplot's hull is not used, since the rotation
  invalidates its coordinates. The hull is built from the training data
  alone (Appendix A of the thesis also includes the target), so a target
  outside the training data gets a counterfactual inside it. The result
  gains `hull` and `hull_polygon` (the winning plane's hull, stored also with
  the hull switched off, for the plot), and the no-solution message names the
  hull. On the loan harness the new applicant's counterfactual has a credit
  score of 534.92; without the hull it would be 489.44. See
  `documentation/local-convex-hull.md`.
* `set_filters()` holds all three of the thesis's constraints on a
  counterfactual (section 4.3.1, p.116): actionability (the per-variable
  constraints), plausibility (`hull = TRUE`, above) and feasibility
  (`train_ranges = TRUE`, every feature within its training range). Both
  switches are on by default, and `set_filters = NULL` in
  `bl_find_local_cf()` keeps them on. `set_filters(tgt, hull = FALSE,
  train_ranges = FALSE)` with no per-variable constraint gives the
  unconstrained counterfactual (thesis p.115). With `train_ranges = FALSE`
  the restored `"fixed"` and rounded `"integer"` points are not checked
  against the ranges either. The result gains the logical field
  `train_ranges`; a search with a rule switched off prints `Not applied`, and
  `print(<bl_filters>)` shows both switches. A feature named `hull` or
  `train_ranges` cannot be constrained through `set_filters()`.
* The empty-region rule (thesis section 3.10, step 4: a boundary area that
  contains no data is discarded) is judged on the training data, in both
  searches, through the private helper `.encloses_data()`. In
  `bl_find_boundary()` a contour that is still closed after the hull clip is
  kept only if it encloses training rows; in 0.1.0 the rows being searched
  decided, so a row's counterfactual could depend on the other rows passed in
  `data` or `tdp`. `bl_find_local_cf()` applies the same rule on each rotated
  plane when the hull is on. On the loan harness the global search keeps 18
  contours (with the searched rows deciding it would keep 12, and 346 of the
  4,429 test counterfactuals would differ), and without the rule the
  constrained counterfactual of test row 13 would lie on a closed contour
  with no training data inside.
* The distance follows the thesis rule by projection (thesis section 3.9.1):
  Euclidean for a PCA biplot, Mahalanobis for a CVA biplot. The new argument
  `distance = c("auto", "mahalanobis", "euclidean")`, whose default `"auto"`
  the private `.resolve_distance()` turns into that rule, sets it in
  `bl_find_local_cf()` (the comparison of component pairs and the
  `"integer"` rounding) and in `plot.bl_boundary()` / `bl_robustness()` (the
  per-variable standardisation). In 0.1.0 the local search compared the
  pairs by the distance on the rotated plane, and the importance divided each
  variable by its total standard deviation (`X_sd`).
  - Under CVA the pair choice is unchanged: the CVA loadings satisfy
    `V'WV = I/(n-g)`, so the squared Mahalanobis distance is `n - g` times
    the squared plane distance for every pair. The importance divides by the
    within-class standard deviation, `sqrt(diag(W))`, which raises features
    that separate the classes well; pass `distance = "euclidean"` for the
    0.1.0 importance.
  - Under PCA the pairs are compared by the Euclidean distance in X-space,
    each feature divided by its total standard deviation. For a standardised
    PCA this equals the plane distance at the search's point (the rotated
    loadings are orthonormal), and it also measures the `"integer"`
    candidates, which lie off the plane. The importance is as in 0.1.0.
* `bl_find_local_cf()`'s `max_pairs` default is `21L` (was `10L`), so the
  search tries every component pair for up to 7 features, as the thesis does.
  Results change for models with 6 or more features: more pairs are tried,
  the selected pair can change, and run time grows with the number of pairs.
  Pass `max_pairs = 10L` for the 0.1.0 search.
* A `"fixed"` constraint holds exactly. The search still accepts a fixed
  feature within 0.5 model units, because a point on the 2D boundary rarely
  holds a feature at exactly its observed value (0.1.0 reported that point).
  `bl_find_local_cf()` then sets each fixed feature the search moved back to
  its observed value. When the restored point still crosses the boundary and
  lies within `train_ranges`, it is the counterfactual. Otherwise there is
  none: `solution_found` is `FALSE`, and `blocking_constraint` (also given as
  a message) says where the search moved the fixed feature, why the restored
  point fails, and that changing the constraints on the actionable features
  may make a counterfactual possible. The new field `fixed_check` records the
  outcome and the search's own point. See
  `documentation/fixed-constraint-restore.md`.
* `bl_find_sparse_cf()` selects features greedily: it changes the
  `"Supports"` features one at a time, largest Shapley value first, re-scores
  the model after each, and stops once the boundary is crossed. Features not
  needed stay at their observed values. The new field `sparse_path` records
  each step, and `print()` shows it. In 0.1.0 every `"Supports"` feature was
  changed. See `documentation/sparse-counterfactual.md`.
* `bl_find_sparse_cf()`'s `round_to` is deprecated: a non-NULL value stops
  with a message. It rounded in model units (z-scores for a standardised
  model), one value for every feature, without re-checking the constraints,
  and could reverse a move and invalidate the sparse counterfactual. On the
  loan harness (target 13, `credit_score = "increase"`), `round_to = 1` moved
  `loan_percent_income` up from 0.05 to 0.0666 instead of down to 0.0484, and
  the sparse path never crossed the boundary (0.483, against 0.547 without
  rounding). Mark a whole-number feature `"integer"` in `set_filters()`
  instead (see New features). `round_to = NULL` still runs.
* `bl_shapley()` adds a `"No effect"` class for a Shapley value of 0 (up to
  `1e-8`): the feature's change never moves the prediction. In 0.1.0 such a
  feature was `"Supports"` for a class-0 target, and so was changed in the
  sparse counterfactual, and `"Contradicts"` for a class-1 target.
* `bl_shapley()`'s `seed` default is `NULL` (was `1L`): it seeds the random
  number generator only when asked, so repeated calls on the permutation path
  (more than `exact_max_vars` features) are no longer silently identical. The
  exact path (the default for up to 14 features) is deterministic and ignores
  `seed`.

## New features

* Raw units. `bl_set_scaling()` records the transform
  `standardised = (raw - center) / scale` applied to the features before the
  model was fitted, on a `bl_data` or `bl_filter_result`; it flows into
  `bl_result$scaling`. It is display-only: the projection, the grid and what
  the model sees do not change. With it, the biplot plot methods
  (`plot_biplotEZ()`, `plot.bl_projection()`, `plot.bl_local_result()`,
  `plot.bl_sparse_result()`, `plot.bl_surrogate()`) label their axes in raw
  units, `print.bl_result()` says so, and the prints and plots that show
  feature values (targets, counterfactuals, Shapley values, the sparse
  counterfactual) show them in raw units. `bl_apply_scaling(data, bl_result)`
  converts new raw-unit data, such as an external applicant, to the model's
  units, which every function that takes new data expects
  (`bl_select_target()`, `bl_project_points()`, `plot(..., target_point =)`).
  `bl_select_target()` warns when a data-frame target lies outside the
  training ranges, and suggests `bl_apply_scaling()` when converting would
  bring it inside.
* `bl_find_local_cf()` returns the counterfactual as a `"bl_counterfactual"`
  object (field `bl_counterfactual`), parallel to `bl_target`: the boundary point's
  values in model units, its prediction and class, and the scaling.
  `print.bl_counterfactual()` shows the values in raw units when a scaling is
  set. `print.bl_target()` and `print.bl_counterfactual()` take a `digits`
  argument, which `print.bl_local_result()` uses to keep its blocks
  consistent.
* Every biplot plot method takes `new_title` for a custom title:
  `plot.bl_local_result()` and `plot.bl_surrogate()` gain it, and
  `plot.bl_result()` and `plot.bl_sparse_result()` pass it on through `...`
  (`plot_biplotEZ()` had it). `plot.bl_projection()` does not take it.
* Mahalanobis distance. `bl_build_projection()` computes the within-class
  metric `W`, its inverse (by Cholesky, with a ridge fallback) and its type
  (`metric`, `metric_inv`, `metric_type`), which flow through
  `bl_assemble()` into `bl_result`. Under CVA it is the `W` of the CVA
  eigenproblem, so distances measured with it agree with the biplot axes.
  - `bl_find_local_cf()` returns `dist_mahalanobis`,
    `all_distances_mahalanobis`, `dist_x`, `all_distances_x` and `distance`,
    and its print shows them.
  - `bl_find_boundary()` returns `dist_mahalanobis`, the squared Mahalanobis
    distance of each observation's move to its boundary point, `v' W^-1 v`
    with `v = B_x - x_obs`: the same quantity as the local search's, so global
    and local distances can be compared. `print()` shows the mean of its
    square root.
  - `plot.bl_boundary()` prints a Mahalanobis check, which compares the
    squared length of the plotted per-feature distances,
    `sum(vec_to_boundary_sd^2)`, with the full `v' W^-1 v`, notes when the
    difference is at least 25% of the full distance, and returns the numbers
    as `mahalanobis_check`.
  - A `bl_result` without `metric_inv` (built by 0.1.0) falls back to
    Euclidean with a warning (`dist_mahalanobis` is then `NA`). See
    `documentation/boundary-distance-and-robustness.md`.
* `set_filters()` accepts an `"integer"` constraint for a whole-number
  feature (a count or a score): `"integer"` alone, with a direction,
  `c("increase", "integer")`, or with a range, `list(c(min, max),
  "integer")`. Without it a feature is continuous. `bl_find_local_cf()`
  searches as before, then rounds the counterfactual: for each integer
  feature it tries the whole values just below and above, in raw units when
  a scaling is set, keeps the combinations that still cross the boundary and
  meet every constraint and `train_ranges`, and takes the nearest by the
  search's selector. When none is valid it reports no counterfactual, with
  the reason. The Shapley values and the sparse counterfactual use the
  rounded point. The new field `integer_check` records the outcome, and
  `print()` explains it.
* `plot.bl_local_result()` (and `plot.bl_sparse_result()` through it)
  colours the prediction grid only inside the plane's convex hull of the
  training data, as the global biplot colours only inside its hull: the
  coloured area is where the search could choose boundary points. The new
  argument `trim_to_hull` (`NULL`, the default, follows the search; `TRUE`
  also shows the hull of a search with the hull switched off; `FALSE`
  colours the whole grid). Display only. The plot also widens its window
  when the target or the counterfactual lies outside the training points, so
  a new applicant beyond the training data is drawn.
* `plot()` works on a `bl_result` (new `plot.bl_result()`, which calls
  `plot_biplotEZ()`), as on a `bl_local_result`, `bl_sparse_result` and
  `bl_surrogate`.
* `plot_biplotEZ()`, `plot.bl_local_result()` and `plot.bl_surrogate()` gain
  `ticks_var`/`ticks_n` for per-variable tick counts (the pattern of
  `label_offset_var`/`label_offset_dist`), and `label_offset_var` accepts
  variable names as well as indices. `plot.bl_local_result()` gains
  `label_cex`, and `plot.bl_surrogate()` gains `rotate_deg`.
* `print()` of a `bl_local_result` ends with an observed / counterfactual /
  change / constraint table, plus a note on any `"fixed"` feature the search
  moved and any `"integer"` feature it rounded. `print()` of a `bl_shapley`
  uses the same layout with `shapley` and `contribute` columns and the
  Shapley sum, and `print()` of a `bl_sparse_result` with `sparse` and `used`
  columns and the sparse path. The tables, the Target and Counterfactual
  blocks and the sparse path print in fixed notation with 2 decimals, in raw
  units when a scaling is set; a feature whose change would print as 0.00
  gets more decimals (up to 4) in its row. The display rules are written up
  in `documentation/scaling_axis_relabel_note.md`.
* `plot.bl_shapley()` labels every variable `"variable: observed ->
  counterfactual"`, with the decimals of the tables, and shows the Shapley
  value on each bar, unclipped. The x-axis reads "Impact on prediction
  (Shapley value)".
* The `"XGB"` model type accepts models from `xgboost::xgboost()` (class
  `"xgboost"`), which are scored from the data frame because their predict
  method refuses an `xgb.DMatrix`. `bl_wrap_model()` stops unless such a model
  is a `binary:logistic` classifier on levels `"0"`/`"1"`: fitted with a
  numeric response, `xgboost()` gives a regression. `xgb.train()` boosters
  are scored as in 0.1.0.
* `bl_pick_point()` gains `bl_boundary`: given a `bl_find_boundary()`
  result, it draws the arrow from each clicked point to its boundary point.
* `bl_surrogate()` reports `hull_coverage`, `n_in_hull` and `n_total`, which
  `print.bl_surrogate()` shows beside the accuracy.
* `bl_build_projection()` checks the condition number of `V`
  (`kappa(V, exact = FALSE)`) before inverting it, returns it as
  `condition_number`, and warns above `1e10`.
* `bl_find_boundary()` says how many of its counterfactuals lie outside the
  training feature ranges; they are retained.
* `bl_find_local_cf()` warns when the target lies at the centre of the
  biplot: the Procrustes cross-product is then zero, the rotation is
  undefined, and every component pair gives the plane of components (1,2).

## Bug fixes and other improvements

* `.bl_rotate()`: when the selected component pair is not `(1, 2)`, the
  target no longer collapses to the biplot origin. `Vrho`/`tVrho` are sliced
  to columns/rows `(1, 2)`, which the SVD construction guarantees.
* `bl_find_boundary()`, `bl_surrogate()` and `bl_find_local_cf()` are
  deterministic. Their shared helper `.nearest_idx_block()` picked the
  nearest contour vertex with `max.col(-d2)`, whose default
  `ties.method = "random"` treats every vertex within 1e-5 of the *farthest*
  vertex's squared distance as tied and draws one from the session's
  random-number stream, so results could change between repeated calls. It
  now takes the exact nearest vertex (`ties.method = "first"`). The
  behaviour was inherited from the PhD code. In the loan harness, 1-3% of
  global counterfactuals could change under the old rule; taking the exact
  vertex moved the v1 robustness total from 9203.63 to 9203.29 and the
  surrogate accuracies from 0.8334/0.8843 to 0.8337/0.8845. The old rule also
  caused the repeat-call differences in the rotation scripts (`17_1`, `17_3`)
  that had been put down to the plane's under-determination. See
  `documentation/nearest-vertex-tie-break.md`.
* `bl_build_grid()` takes its `cutoff` from the model
  (`cutoff = bl_model$cutoff`) and stores it as `bl_grid$cutoff`;
  `bl_assemble()` stops if the grid's and the model's cutoffs differ. In
  0.1.0 the grid always drew its contours at `0.5`, while the searches used
  `bl_result$cutoff`, so a model with another cutoff gave a global boundary at
  one threshold and a local one at another. Results at the default `0.5` are
  unchanged.
* `bl_build_grid()` clamps the grid's predictions to `[0, 1]`, with a warning
  that reports how many cells were clamped, and `bl_find_local_cf()` clamps
  its grid scores the same way, with one warning per call. A custom
  `predict_fn` returning values outside `[0, 1]` gave missing colours,
  silently dropped cells, or stopped the search with an error.
* `plot_biplotEZ(grid_col =)` splits the two-colour grid at the model's
  cutoff, not at a fixed 0.5. `bl_build_grid()` no longer computes the
  unused `inside_data`.
* Every biplot uses the thesis colours for the four prediction outcomes: TP
  red, TN blue, FP purple, FN orange. `plot.bl_projection()` (through
  `bl_build_projection()`'s `point_col`) used green, steel blue, orange and
  red, so red meant FN there and TP everywhere else. The colours are held in
  one private constant, `.confusion_cols` (`R/utils.R`).
* `bl_build_projection()` with CVA: if a TP/TN/FP/FN group is empty, a
  message names it and the empty level is dropped, as is any empty level of a
  user-supplied `cva_classes`, so CVA and the metric `W` count only groups
  with rows. `plot.bl_projection()` uses the confusion colours for any subset
  of the four groups; with an empty group it coloured every point blue.
* `plot.bl_local_result()` colours the target by its own true class. It
  looked the class up in `test_data` by row number, which gave another
  observation's class whenever `bl_select_target()` indexed a different data
  frame. `bl_target` gains a `class_obs` field; a one-row data frame target
  that includes the class column gets its confusion colour too.
* `bl_find_sparse_cf()`: `used_in_sparse` is `FALSE` for a `"fixed"`
  feature, which the sparse counterfactual always reverts. It followed the
  Shapley classification alone and could read `TRUE` for a feature the
  sparse counterfactual did not change.
* `label_offset_dist` and `ticks_n` stay paired with their own variables when
  `label_offset_var` or `ticks_var` names a variable the model does not have,
  and a warning flags a value vector whose length is neither 1 nor the length
  of its variable vector. Affects `plot_biplotEZ()`, `plot.bl_local_result()`
  and `plot.bl_surrogate()` (display only).
* `bl_prepare_data()`, `bl_filter_outliers()`, `bl_build_grid()` and
  `bl_build_result()` no longer change the active graphics device when more
  than one device is open. The internal helpers `.build_hull_polygon()` and
  `.bl_get_plot_bounds()` restore the caller's device after closing their
  temporary `png()`.
* `bl_find_local_cf()`'s per-pair search applies the actionability
  constraints before `train_ranges`, so a segment that fails the more
  selective rule is dropped early.
* `.poly_clip()`, the internal hull filter, is defined once, in
  `R/hull_utils.R`. The copy in `R/boundary.R` was never used (the other
  replaced it at load time) and is removed. Results are unchanged.
* `print.bl_result()` and `summary.bl_result()` label the stored hull
  fraction as the grid hull (`outlie`, default 1), not as the outlier
  filter's fraction, which `bl_prepare_data()` applies earlier.
* `bl_assemble()` no longer fails with `bl_model = NULL`;
  `bl_build_result()` stops with an error where it returned `invisible(NULL)`;
  `print.bl_points()` prints `"(no model)"` when there are no model
  predictions; `.fit_model()` warns on unknown hyperparameters passed in
  `model_params`.
* `DESCRIPTION` lists `grDevices`, `graphics`, `stats` and `utils` under
  `Imports`, matching the `NAMESPACE` imports.
* The help pages and messages of `bl_find_local_cf()` and
  `plot.bl_local_result()` say "component pair", not "eigenvector pair": each
  pair labels a plane through the target, which the Procrustes rotation
  turns, so the selected plane need not coincide with the plane of the two
  components that name it.
* The help of `plot.bl_boundary()` and `bl_robustness()` says that the
  variable importance and robustness measure the in-plane gap to the
  boundary, as the thesis does, while `bl_find_boundary()`'s `dist_x` and
  `dist_mahalanobis` also include the observation's part off the biplot plane
  and are reported only.
* The loan-default vignette is rewritten for package users. It follows the
  development harness on the full data: the same models, the five pruned
  variables, the refit on standardised features with `bl_set_scaling()`, test
  row 13 with the harness's constraints, and the new applicant converted with
  `bl_apply_scaling()`. Each step gives the code, a table of the parameters
  that matter and how to read the output, and it ends with the SHAP
  comparison. It does not call `bl_robustness()`: the pruning reads the
  totals `plot()` returns. The vignette is precomputed:
  `vignettes/precompile.R` knits `Boundary_Logic_loan_default_workflow.Rmd.orig`
  into the shipped `.Rmd` and its figures, so building the package runs no
  vignette code.
* The README's constraint example finds a counterfactual and shows its
  output.
* The technical notes in `documentation/` describe the method and the code:
  the boundary distance and robustness, the biplot rotation and the choice of
  plane, the fixed-feature restore, the sparse counterfactual, the local
  convex hull, the raw-unit display, the nearest-vertex tie-break, and the
  stepwise algorithms (`documentation/algorithms.md`).
* Documentation brought in line with the code across the roxygen help, the
  README, the vignette prose, the implementation summary and the code
  walkthroughs (for example: `bl_build_result(bl_model = NULL)` returns a
  model-free `bl_result`; `bl_filter_outliers()` takes a `bl_data` from
  `bl_wrap_data()`; `dist_z` / `dist_mahalanobis` belong to the winning
  pair only; `bl_predict()` rounds to 3 decimals; the `@return` fields of
  `bl_surrogate()` and `bl_build_projection()` match the objects).

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
