# boundarylogic 0.2.0

## Breaking changes

* The empty-region rule (thesis section 3.10, step 4: a boundary area that
  contains no data is discarded) is judged on the training data, and applied
  in both searches. In `bl_find_boundary()` a contour that is still closed
  after the hull clip is kept only if it encloses training rows; before, the
  rows being searched decided, so a row's counterfactual could depend on which
  other rows were passed in `data` or `tdp`. `bl_find_local_cf()` now applies
  the same rule on each rotated plane when the hull is on. On the loan harness
  the global search keeps 18 contours instead of 12 on the reduced model (346
  of 4,429 test counterfactuals move), and the constrained counterfactual of
  test row 13 moves (an empty closed contour held its nearest boundary point);
  on the Pima example of Paper 1 no global result and no run of either target
  changes. The private helper `.encloses_data()` holds the rule. Decided by
  Adriaan Rowan on 2026-11-06 (`documentation/algorithms.md`, finding F2).
* `bl_find_local_cf()` builds each rotated plane's convex hull from the
  training data alone; it included the target (as Appendix A of the thesis
  does). A target outside the training data no longer enlarges the region its
  counterfactual may come from: the loan harness's new applicant now gets a
  counterfactual inside the training hull (credit score 534.92 instead of
  489.44). `hull = FALSE` is unaffected. The stored `hull_polygon` follows.
  Decided by Adriaan Rowan on 2026-11-06 (finding F4).

* `bl_find_local_cf()` now applies the thesis's plausibility rule (section
  4.3.1): the counterfactual is selected from a local decision boundary within
  the convex hull of the training data. For each component pair, the convex
  hull of the training data and the target is rebuilt on the rotated plane, as
  in Appendix A of the thesis, and only the boundary points inside it are kept,
  before the actionability, `train_ranges` and model filters. The global
  biplot's hull polygon is still not used, since the rotation invalidates it.
  The switch is `set_filters(..., hull = TRUE)`, the default, also when no
  `set_filters()` object is passed; `hull = FALSE` reproduces the earlier
  search exactly. The result gains the fields `hull` and `hull_polygon` (the
  winning plane's hull, stored also with the hull switched off, for the
  plot), and the no-solution message names the hull. Results
  change where the nearest boundary point lay outside the hull: on the Pima
  example (Paper 1), 11 of the 79 test targets move, none loses its
  counterfactual, and the TP target's constrained runs that hold *Age* now
  have none. The loan harness is unchanged. Decided by Adriaan Rowan on
  2026-11-05; see `documentation/local-convex-hull.md`.
* `set_filters()` now holds all three of the thesis's constraints on a
  counterfactual (section 4.3.1, p.116): actionability (the per-variable
  constraints), plausibility (`hull = TRUE`, above) and feasibility
  (`train_ranges = TRUE`, every feature within its training range). The
  training-range rule had no switch before. Both switches are on by default,
  and `set_filters = NULL` in `bl_find_local_cf()` keeps them on, so no
  existing result changes. `set_filters(tgt, hull = FALSE, train_ranges =
  FALSE)` with no per-variable constraint gives the unconstrained
  counterfactual (thesis p.115). With `train_ranges = FALSE` the restored
  `"fixed"` and rounded `"integer"` points are not checked against the ranges
  either. The result gains the logical field `train_ranges`; a search with a
  rule switched off prints `Not applied`, and `print(<bl_filters>)` shows both
  switches. `bl_find_local_cf()`'s `hull` argument (added earlier in this
  version) is removed: pass `set_filters = set_filters(tgt, hull = FALSE)`.
  A feature named `hull` or `train_ranges` cannot be constrained through
  `set_filters()`. Decided by Adriaan Rowan on 2026-11-06.
* `bl_prepare_data()`'s `hull_fraction` default is now `1` (no outlier removal;
  it was `0.9`), as in the thesis, where outlier removal is optional. Pass
  `hull_fraction = 0.9` to remove the outermost 10% of the training rows, as the
  thesis's loan example and the package's loan harness and vignette do. Code
  that relied on the old default keeps every training row now.
  `bl_filter_outliers()` keeps `0.9`, since calling it means filtering is
  wanted.
* The distance now follows the thesis rule by projection (thesis section
  3.9.1): Euclidean distance for a PCA biplot, Mahalanobis distance for a CVA
  biplot. The new default `distance = "auto"` resolves it in
  `bl_find_local_cf()` (the pair comparison and the `"integer"` rounding) and
  in `plot.bl_boundary()` / `bl_robustness()` (the variable-importance
  standardisation). Under PCA the defaults change: the pairs are compared by
  Euclidean distance in X-space, and the importance divides by the total
  standard deviation (`X_sd`), as in the thesis; pass
  `distance = "mahalanobis"` to reproduce earlier PCA results. Under CVA
  nothing changes. Made on explicit instruction (2026-11-05); it resolves the
  deferred PCA selector question (rotation note §4.4).
* `bl_find_local_cf(distance = "euclidean")` now measures in X-space, each
  feature divided by its total standard deviation, as `bl_find_boundary()`'s
  `dist_x` (it was the distance on the rotated 2D plane). Under a
  standardised PCA the two are equal at the search's point, since the rotated
  loadings are orthonormal; the X-space distance also measures the
  `"integer"` candidates correctly, which lie off the plane. Under CVA an
  explicit `"euclidean"` is a different metric from the plane distance, whose
  ranking equals the Mahalanobis one. `bl_find_local_cf()` returns the new
  fields `dist_x` and `all_distances_x`, and its print shows them.
* `bl_find_sparse_cf()`'s `round_to` is deprecated: a non-NULL value stops
  with a message. It rounded in model units (z-scores for a standardised
  model), one value for every feature, without re-checking the constraints,
  and could reverse a move and invalidate the sparse counterfactual. Mark a
  whole-number feature `"integer"` in `set_filters()` instead (see below).
  `round_to = NULL` still runs.
* A `"fixed"` constraint in `set_filters()` is now strict. When the restored
  point (the search's boundary point with each moved fixed feature set back
  to its observed value) no longer crosses the boundary, or leaves
  `train_ranges`, `bl_find_local_cf()` reports no counterfactual:
  `solution_found` is `FALSE`, `fixed_check` is kept with `applied = FALSE`,
  and `blocking_constraint` (also given as a message) says where the search
  moved the fixed feature, why the restored point fails, and that changing
  the constraints on the actionable features may make a counterfactual
  possible. Previously the search's point was kept, so a reported
  counterfactual could move a "fixed" feature by up to 0.5 model units. A
  reported counterfactual now always holds every fixed feature exactly. See
  `documentation/fixed-constraint-restore.md`.
* `bl_find_local_cf()`'s `max_pairs` default is now `21L` (was `10L`), so the
  search tries every eigenvector pair for up to 7 features, as the thesis
  does. Results change for models with 6 or more features: more pairs are
  tried, the selected pair can change, and run time grows with the number of
  pairs. Pass `max_pairs = 10L` to reproduce earlier results.
* `rotate_deg` is removed from `plot.bl_local_result()` and
  `plot.bl_sparse_result()`; passing it stops with an error. `rotate_deg` is
  only for viewing the global biplot (`plot_biplotEZ()`, `plot.bl_surrogate()`
  keep it). The local biplot's orientation is fixed by the optimal rotation
  (`.bl_rotate()`), which is part of the method, and is never turned further.
  This also removes a display bug: with `rotate_deg != 0`, the sparse-CF cross
  was drawn in the unrotated frame.
* `rounding` parameter replaced by `b_margin` in `bl_build_grid()` and
  `bl_build_result()`; `bl_find_local_cf()` reads it from
  `bl_result$b_margin`. `b_margin` sets the decision
  boundary contour band half-width *directly* (default `0.001`, equivalent to
  the old `rounding = 3L`), instead of indirectly via `b_margin = 1/10^rounding`.
  Model predictions are unaffected — they are always floor-rounded to 3 d.p.
* `no_points` renamed to `plot_points` in `plot_biplotEZ()` and
  `plot.bl_local_result()`, with positive semantics (`TRUE` = show points).
  Defaults preserve prior visual behaviour for each function.
  (`plot.bl_surrogate()` has neither argument; it always draws its points.)
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

* `plot.bl_local_result()` (and `plot.bl_sparse_result()` through it) colours
  the prediction grid only inside the convex hull of the training data and
  the target on the plane when the search applied the hull, as the global
  biplot colours only inside its hull: the coloured area is where the search
  could choose boundary points, and boundary lines in the white area lie
  outside it. The new argument `trim_to_hull` (`NULL`, the default, follows
  the search; `TRUE` also shows the hull of a search with the hull switched
  off; `FALSE` colours the whole grid). Display only.
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
* `bl_find_boundary()` gains a `dist_mahalanobis` field: the squared
  Mahalanobis distance of each observation's move to its boundary point,
  `v' W^-1 v` with `v = B_x - x_obs` and `W^-1 = bl_result$metric_inv`. It is
  the same quantity as `bl_find_local_cf()`'s field of that name, so global
  and local distances can be compared directly. `print()` shows the mean of
  its square root. It is `NA`, with a warning, for a `bl_result` without
  `metric_inv`. No search result changes.
* New `bl_apply_scaling(data, bl_result)` converts new data given in raw
  units, such as an external applicant, to the model's units with the
  scaling recorded by `bl_set_scaling()`. `bl_set_scaling()` is
  display-only, so every function that takes new data (`bl_select_target()`,
  `bl_project_points()`, `plot(..., target_point =)`) expects model units.
  `bl_select_target()` now warns when a data-frame target lies outside the
  training ranges, and suggests `bl_apply_scaling()` when the target looks
  like raw units.
* Added Mahalanobis-aware distance measures, replacing per-feature Euclidean
  (divided by total SD) with squared Mahalanobis distance using the
  within-class covariance matrix `W`:
  - `bl_find_local_cf()` gains a `distance = c("mahalanobis", "euclidean")`
    parameter governing the *cross-pair* eigenvector-pair selector. New result
    fields `dist_mahalanobis`, `all_distances_mahalanobis`, `distance`.
  - `bl_robustness()` and `plot.bl_boundary()` gain the same `distance`
    parameter governing the per-feature standardisation denominator
    (`sqrt(diag(W))` instead of total SD `X_sd`). `plot.bl_boundary()` also
    prints a Mahalanobis check, with a note when the per-feature measure
    leaves out at least 25% (in absolute value) of the full squared
    Mahalanobis distance.
  - `bl_build_projection()` now computes and stores the metric matrix and its
    Cholesky-based inverse (`metric`, `metric_inv`, `metric_type`), propagated
    through `bl_assemble()` into `bl_result`. See
    `documentation/boundary-distance-and-robustness.md` for the full derivation.
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
* `bl_find_local_cf()` now restores `"fixed"` features after the search. The
  search still lets a fixed feature move by up to 0.5 model units, because a
  point on the 2D boundary rarely holds a feature at exactly its observed
  value. Each fixed feature the search moved is then set back to its observed
  value; if that point still crosses the boundary (and lies within
  `train_ranges`), it becomes the counterfactual. A new `fixed_check` field
  records the outcome and the search's own point. Results without a moved
  fixed feature are unchanged. See `documentation/fixed-constraint-restore.md`.
* `print()` of a `bl_local_result` now ends with an observed / counterfactual /
  change / constraint table (raw units when a scaling is set), plus a note
  explaining any `"fixed"` feature the search moved. `print()` of a
  `bl_sparse_result` uses the same table layout, with `sparse` and `used`
  columns.
* `bl_shapley()` adds a `"No effect"` class for a Shapley value of 0 (the
  feature's change never moves the prediction). For a class-0 target such a
  feature was previously labelled `"Supports"` and changed in the sparse
  counterfactual; for a class-1 target it was labelled `"Contradicts"`.
  `print()` of a `bl_shapley` uses the same comparison-table layout, with
  `shapley` and `contribute` columns and the Shapley sum.
* `bl_find_sparse_cf()` now selects features greedily: it changes the
  `"Supports"` features one at a time, largest Shapley value first, re-scores
  the model after each, and stops once the boundary is crossed. Features not
  needed stay at their observed values. The new field `sparse_path` records
  each step, and `print()` shows it. Previously every `"Supports"` feature
  was changed. See `documentation/sparse-counterfactual.md`, which also
  records how `round_to` was used and the `"integer"` constraint that
  replaced it.

## Bug fixes and other improvements

* `bl_find_local_cf()` clamps its grid scores to [0, 1], as `bl_build_grid()`
  does, with one warning per call. A custom `predict_fn` returning values
  outside [0, 1] gave missing colours, silently dropped cells, or stopped the
  search with an error (finding F4).
* `plot.bl_local_result()` (and so `plot.bl_sparse_result()`) widens its
  window when the target or the counterfactual lies outside the training
  points, so a new applicant beyond the training data is drawn; other plots
  are unchanged.
* `plot_biplotEZ(grid_col =)` splits the two-colour grid at the model's
  cutoff, not at a fixed 0.5 (finding F6). `bl_build_grid()` no longer
  computes the unused `inside_data`.
* The documentation of `plot.bl_boundary()` and `bl_robustness()` says that
  the variable importance and robustness measure the in-plane gap to the
  boundary, as the thesis does, while `bl_find_boundary()`'s `dist_x` and
  `dist_mahalanobis` also include the observation's part off the biplot plane
  and are reported only (finding F3; no result changes).
* The loan-default vignette is rewritten for package users. It follows the
  development harness on the full data (it drew a 4,000-row subsample): the
  same models, the five pruned variables, the refit on standardised features
  with `bl_set_scaling()`, test row 13 with the harness's constraints, and the
  new applicant converted with `bl_apply_scaling()`. Each step gives the code,
  a table of the parameters that matter and how to read the output, and it
  still ends with the SHAP comparison, now showing feature values in raw
  units. It no longer calls `bl_robustness()`: the pruning reads the totals
  `plot()` returns. The vignette is precomputed: `vignettes/precompile.R`
  knits `Boundary_Logic_loan_default_workflow.Rmd.orig` into the shipped
  `.Rmd` and its figures, so building the package runs no vignette code.
* `.poly_clip()`, the internal hull filter, is defined once, in
  `R/hull_utils.R`. The copy in `R/boundary.R` was never used (the other
  replaced it at load time) and is removed. Results are unchanged.
* Every biplot now uses the thesis colours for the four prediction outcomes:
  TP red, TN blue, FP purple, FN orange. `plot.bl_projection()` (through
  `bl_build_projection()`'s `point_col`) used green, steel blue, orange and
  red, so red meant FN there and TP everywhere else. The colours are held in
  one private constant, `.confusion_cols` (`R/utils.R`).
* `bl_find_local_cf()` warns when the target lies at the centre of the biplot:
  the Procrustes cross-product is then zero, the rotation is undefined, and
  every component pair gives the plane of components (1,2). Results are
  unchanged.
* The help pages and messages of `bl_find_local_cf()` and
  `plot.bl_local_result()` say "component pair", not "eigenvector pair": each
  pair labels a plane through the target, which the Procrustes rotation turns,
  so the selected plane need not coincide with the plane of the two components
  that name it. The `rotate_deg` error names the Procrustes rotation.
* The README's constraint example finds a counterfactual (it found none) and
  shows its output; the loan harness's first biplot is titled "no prior
  default", which matches its data.
* `plot.bl_shapley()` labels every variable `"variable: observed ->
  counterfactual"`, with and without a recorded scaling. Without a scaling it
  showed the observed value and the *change* (e.g. `"Age: 22 -> -0.88"`),
  which read as a negative age. Observed and counterfactual now share the
  decimals of the comparison tables (2, more where the change would print as
  0.00), as do the Shapley values on the bars, which are no longer clipped at
  the plot's edge. The x-axis reads "Impact on prediction (Shapley value)".
  The display rules are written up in
  `documentation/scaling_axis_relabel_note.md`.
* The observed / counterfactual / change tables of `print()` for
  `bl_local_result`, `bl_shapley` and `bl_sparse_result`, the Target and
  Counterfactual blocks, and the sparse path now print in fixed notation with
  2 decimals. They printed with 4 decimals and switched to scientific
  notation when a column mixed large and small values. A feature whose change
  would print as 0.00 gets more decimals (up to 4) in its row, so a small
  move such as a loan-to-income change of -0.0016 stays visible.
  `print.bl_target()` and `print.bl_counterfactual()` gain a `digits`
  argument, which `print.bl_local_result()` uses to keep its blocks
  consistent.
* `plot.bl_boundary()`'s Mahalanobis diagnostic now describes the plotted
  measure. Its "diagonal" used `diag(W^-1)`, which is not what the plot
  divides by (`sqrt(diag(W))`); for correlated features
  `(W^-1)_jj >= 1/W_jj`, so the printed "cross-correlation" share mostly
  measured that mismatch (-42.6% on the loan harness). The renamed
  "Mahalanobis check" compares the squared length of the plotted distances,
  `sum(vec_to_boundary_sd^2)`, with the full `v' W^-1 v` (+17.7% on the same
  run), and returns the numbers as a new `mahalanobis_check` field. See
  `documentation/boundary-distance-and-robustness.md` §4.2.
* `bl_find_boundary()`, `bl_surrogate()` and `bl_find_local_cf()` are now
  deterministic. Their shared helper `.nearest_idx_block()` picked the nearest
  contour vertex with `max.col(-d2)`, whose default `ties.method = "random"`
  treats every vertex within 1e-5 of the *farthest* vertex's squared distance
  as tied and draws one from the session's random-number stream. Results could
  therefore change between repeated calls and depend on what had been called
  before (and on any `set.seed()`). It now takes the exact nearest vertex
  (`ties.method = "first"`). The behaviour was inherited from the PhD code. In
  the loan harness, 1-3% of global counterfactuals could change under the old
  rule. Their 2D distance moved by less than 1e-4, but their feature values
  moved by up to 0.36 training SDs and their model score by up to 0.18. See
  `documentation/nearest-vertex-tie-break.md`.
* Documentation: the help for `bl_find_local_cf()` no longer says the
  `"euclidean"` selector is "not comparable across pairs". Under CVA the two
  `distance` options are the same criterion: the CVA loadings satisfy
  `V'WV = I/(n-g)` for every rotated basis, so the squared Mahalanobis distance
  is `n - g` times the squared in-plane distance for every pair, and both pick
  the same pair. Under PCA they differ, and the thesis's rule is the Euclidean
  one. The rotation and pair-selection material is now one note,
  `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md` (the
  original PhD method, the implementation, its problems, a proposed solution,
  and future work); the Mahalanobis note is renamed
  `documentation/boundary-distance-and-robustness.md` and covers Phase 2 only.
* `bl_build_grid()` now takes its `cutoff` from the model by default
  (`cutoff = bl_model$cutoff`) and stores it as `bl_grid$cutoff`; `bl_assemble()`
  stops if the grid's and the model's cutoffs differ. Previously the grid
  always drew its contours at `0.5` (`bl_build_result()` passed no cutoff),
  while `bl_find_boundary()`, `bl_surrogate()` and `bl_find_local_cf()` used
  `bl_result$cutoff`, so a model built with another cutoff gave a global
  boundary at one threshold and a local one at another, without warning.
  Results at the default `0.5` are unchanged.
* The `"XGB"` model type now accepts models from `xgboost::xgboost()` (class
  `"xgboost"`), which are scored from the data frame because their predict
  method refuses an `xgb.DMatrix`. `bl_wrap_model()` stops unless such a model
  is a `binary:logistic` classifier on levels `"0"`/`"1"`: fitted with a
  numeric response, `xgboost()` gives a regression. `xgb.train()` boosters are
  scored as before.
* `bl_build_projection()` with CVA: if a TP/TN/FP/FN group is empty, a message
  names it and the empty level is dropped, as is any empty level of a
  user-supplied `cva_classes`, so CVA and the metric `W` count only groups
  with rows. `plot.bl_projection()` now uses the confusion colours for any
  subset of the four groups; with an empty group it previously coloured every
  point blue.
* `bl_find_sparse_cf()`: `used_in_sparse` is now `FALSE` for a `"fixed"`
  feature, which the sparse counterfactual always reverts. It previously
  followed the Shapley classification alone and could read `TRUE` for a
  feature the sparse counterfactual did not change.
* `bl_prepare_data()`, `bl_filter_outliers()`, `bl_build_grid()` and
  `bl_build_result()` no longer change the active graphics device when more
  than one device is open. The internal helpers `.build_hull_polygon()` and
  `.bl_get_plot_bounds()` now restore the caller's device after closing their
  temporary `png()`.
* `plot.bl_local_result()` now colours the target by its own true class.
  Previously it looked the class up in `test_data` by row number, which gave
  another observation's class whenever `bl_select_target()` indexed a
  different data frame. `bl_target` gains a `class_obs` field; a one-row data
  frame target that includes the class column now gets its confusion colour
  too.
* `label_offset_dist` and `ticks_n` now stay paired with their own variables
  when `label_offset_var` or `ticks_var` names a variable the model does not
  have. Previously the unknown name was dropped but its value was not, so
  every later value shifted one variable along. A warning now also flags a
  value vector whose length is neither 1 nor the length of its variable
  vector. Affects `plot_biplotEZ()`, `plot.bl_local_result()` and
  `plot.bl_surrogate()` (display only).
* `plot.bl_boundary()` printed `25%%` instead of `25%` in its
  cross-correlation note.
* `print.bl_result()` and `summary.bl_result()` now label the stored hull
  fraction as the grid hull (`outlie`, default 1) rather than as the outlier
  filter's fraction, which `bl_prepare_data()` applies earlier.
* `DESCRIPTION` now lists `grDevices`, `graphics`, `stats` and `utils` under
  `Imports`, matching the `NAMESPACE` imports.
* Documentation brought in line with the code across the roxygen help, the
  README, the vignette prose, the implementation summary and the code
  walkthroughs (for example: `bl_build_result(bl_model = NULL)` returns a
  model-free `bl_result`; `bl_filter_outliers()` takes a `bl_data` from
  `bl_wrap_data()`; `dist_z` / `dist_mahalanobis` belong to the winning
  pair only; the sparse-CF cross is green/yellow; XGB models must come from
  `xgboost::xgb.train()`).
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
  calls on the permutation path (more than `exact_max_vars` features) are no
  longer silently deterministic. The exact path (the default for up to 14
  features) is deterministic and ignores `seed`.
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
