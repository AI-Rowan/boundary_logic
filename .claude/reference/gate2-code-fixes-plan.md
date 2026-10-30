> **IMPLEMENTED 2026-10-30** on branch `paper_nr1_changes`. Archived copy of the approved plan for the four Gate 2 code issues (sparse-CF `rotate_deg`, grid `cutoff`, `xgboost()` models, empty CVA groups); no separate technical note -- the rationale is in this file, NEWS and `2 implementation_summary.txt`.

# Gate 2 code issues: background and suggested fixes

Date created: 2026-10-30
Date last edited: 2026-10-30

## Implementation log (2026-10-30)

- **Decisions:** all four fixed now. Issue 1 removes the argument; Issue 3 is fixed with the
  wrap-time checks. Session 2 was committed first as `7aa1344` (not pushed). These fixes are
  uncommitted.
- **Change from the plan, Issue 4:** the confusion factor keeps alphabetical levels
  (`FN, FP, TN, TP`, exactly what `as.factor()` gave), not `TP, TN, FP, FN`.
  - The first harness run with `TP, TN, FP, FN` changed every Phase 3 number. Target 13's
    constrained CF went from 0.547 to 0.509, and Step 15 went from 0.5554 / 0.542 to
    0.6094 / 0.517.
  - CVA is order-invariant in exact arithmetic. The level order feeds `biplotEZ::CVA()`, though,
    and the resulting change in `V` is enough to move `.bl_rotate()`'s under-determined plane
    (CLAUDE.md §9).
  - A test now pins the alphabetical order. The lesson is recorded in CLAUDE.md §6.
- **Verification:**
  - `devtools::test()`: 278 PASS / 0 FAIL / 4 WARN / 0 SKIP (the 4 baseline WARNs).
  - Harness: the full console log matches a run of the committed `7aa1344` code line for line.
    Only creation timestamps differ (the `HEAD` export was run from the scratchpad).
  - `R CMD check`: Status OK, 0 / 0 / 0.
- **Not as planned:** `.rotation_matrix()` and `.sparse_z()` were not needed once the argument was
  removed.

## Context

The Task C plan (`commit-and-push-the-mutable-hearth.md`, "Issues found while planning") lists four
code issues for the Gate 2 issues table. You decide "fix now" or "defer" for each one. Deferred
issues go into `documentation/known-issues.md` (Task D). This file gives each issue's background,
where it is in the code, whether it affects the Pima paper, and a suggested fix. Each fix comes
with an implementation plan in case you choose "fix now".

At `cutoff = 0.5`, none of the four fixes changes any number in the loan harness or the Pima
example.
- **Fix 1** removes an argument. A call that passes `rotate_deg` to a local plot now stops with an
  error; nothing in the package, harness, vignette or tests makes such a call.
- **Fix 4** changes only what is drawn, plus a message.
- **Fixes 2 and 3** make code paths that currently fail or disagree work correctly.

---

## Issue 1: the sparse-CF cross ignores `rotate_deg`

**Your rule (2026-10-30).** `rotate_deg` is only for viewing the global (optimal) biplot. Whenever
the rotation function `.bl_rotate()` is used, that is in the local biplot for an individual target,
`rotate_deg` must be 0. It is never used to plot a rotated individual point.

**Decision (2026-10-30):** remove `rotate_deg` from `plot.bl_local_result()` and
`plot.bl_sparse_result()`, and stop with an error if it is passed. `plot_biplotEZ()` and
`plot.bl_surrogate()` keep it.

**What `rotate_deg` is.**
- **Display only.** `rotate_deg` turns the finished 2D biplot clockwise by a fixed angle, at plot
  time. The projection, grid, contours, counterfactual search and every stored number are
  untouched.
- **Why turning changes nothing you read.** You read a biplot by dropping a point perpendicularly
  onto a variable's calibrated axis. `.apply_biplot_rotation()` (`R/plot_biplot.R:68`) turns the
  points and the axes by the same orthogonal matrix and recomputes `ax.one.unit`. Axis readings,
  distances and which side of the boundary a point is on are therefore all unchanged.
- **Why it is useful on the global biplot.** The eigenvectors' signs are arbitrary, so the plane's
  "up" carries no meaning. `rotate_deg` lets you choose a readable orientation, for example:
  - stand the decision boundary upright;
  - move crowded axis labels away from a corner, the title or the legend;
  - line a figure up with an earlier run.
- **Limit.** A rotation cannot undo a reflection. If one eigenvector's sign flips, the plot is
  mirrored and no angle turns it back.
- **Why it does not belong on the local biplot.** The local biplot's orientation comes from the
  optimal rotation in `.bl_rotate()`, which is part of the method: it chooses the plane through the
  target, and the counterfactual search runs in that plane. Turning that plane again afterwards is
  not part of the method.

**Background.**
- `plot.bl_sparse_result()` (`R/shapley.R:530`) first draws the local biplot by calling
  `plot.bl_local_result()`, passing `...` through.
- If you pass `rotate_deg = 30`, `plot.bl_local_result()` (`R/local_cf.R:1216-1231`) rotates
  everything it draws by the same 2x2 matrix `R_mat`: the axes, grid, training points, contours,
  target and full-CF cross.
- `plot.bl_sparse_result()` never sees `R_mat`. It projects the sparse CF with `Vr_rot` only
  (`R/shapley.R:550-553`), so its cross is drawn in the unrotated frame.
- The walkthrough already states this (`review_section9_to_15.md:861`).
- **Colours:** the second half of this item ("its colours are green3/yellow2, not the documented
  blue/red") is already settled. The roxygen now documents green3/yellow2. Making the colours
  consistent across plots is Task D item 5, not part of this issue.

**Effect on the paper:** none. Figures 7 and 9 (`plot(bl_local)`, `plot(bl_sparse)`) are drawn
with no extra rotation, and the harness never passes `rotate_deg` to a local plot.

**Fix (no numbers change).**
1. **`plot.bl_local_result()`** (`R/local_cf.R:1097`):
   - remove `rotate_deg = 0` from the signature and its `@param`;
   - remove the `.apply_biplot_rotation()` block (`:1215-1231`). The grid and contours come
     straight from `x$Zgrid` and `x$ct_local`;
   - as the first line of the body, before the `solution_found` check and before any device opens:
     `if ("rotate_deg" %in% ...names()) stop(...)`. `...names()` is available from R 4.1, which
     matches DESCRIPTION's `R (>= 4.1.0)`. The message says that `rotate_deg` is only for viewing
     the global biplot (`plot_biplotEZ()`, `plot(bl_surrogate)`), and that the local biplot's
     orientation is fixed by the optimal rotation.
   - Without this check, `rotate_deg` would fall into `...` and be silently ignored.
2. **`plot.bl_sparse_result()`** needs no rotation code. `rotate_deg` reaches
   `plot.bl_local_result()` through `...` and hits the same error before anything is drawn. The
   `@param ...` text says `rotate_deg` is not accepted. The sparse cross then always sits in the
   same frame as the rest of the plot, which removes the bug.
3. **Unchanged:** `plot_biplotEZ()` and `plot.bl_surrogate()` keep `rotate_deg`, and
   `.apply_biplot_rotation()` (`R/plot_biplot.R:68`) is unchanged. The comment at
   `R/local_cf.R:1203` ("same approach as the `rotate_deg` path in `plot_biplotEZ()`") is about the
   `ax.one.unit` formula and stays.
4. **CLAUDE.md §4:** add one line: "Never add `rotate_deg` back to `plot.bl_local_result()` /
   `plot.bl_sparse_result()`. It is for viewing the global biplot only, and the local biplot's
   orientation comes from `.bl_rotate()`."

**Recommendation:** fix now (decided).

---

## Issue 2: `cutoff` is not passed to `bl_build_grid()`

**Background.**
- The decision threshold lives on the model: `bl_fit_model(cutoff =)` / `bl_wrap_model(cutoff =)`
  sets `bl_mod$cutoff`, and `bl_assemble()` copies it to `bl_result$cutoff` (`R/result.R:161`).
- `bl_build_grid()` has its own `cutoff = 0.5` argument (`R/biplot_grid.R:168`). It extracts the
  boundary contours `ct` at `cutoff +/- b_margin` and the surrogate contours `ct_surrogate`. It
  does not store the value it used.
- `bl_build_result()` calls `bl_build_grid()` without `cutoff` (`R/result.R:348-356`) and has no
  `cutoff` argument of its own. The grid therefore always uses 0.5, whatever the model says.
- **Who reads what:**
  - `bl_find_boundary()` and `bl_surrogate()` use the grid's 0.5 contours, but classify them with
    `bl_result$cutoff` (`R/boundary.R:311,347`, `R/surrogate.R:193`).
  - `bl_find_local_cf()` extracts its own contours at `bl_result$cutoff` (`R/local_cf.R:550`).
  - So with, for example, `cutoff = 0.4`, the global boundary and plot sit at 0.5 while the local
    CF sits at 0.4, and the global search mislabels which side each contour is on. Nothing warns.
- `review_section4_to_6.md:374` states the current behaviour.

**Effect on the paper:** none. Only `cutoff = 0.5` is methodologically valid (CLAUDE.md §4), and
Pima uses 0.5. This is a consistency guard for users who explore other cutoffs.

**Suggested fix (no numbers change at 0.5).**
1. Change `bl_build_grid()`'s default to `cutoff = bl_model$cutoff`. R evaluates the default lazily,
   after `bl_model` has been validated. The effective default is still 0.5, because both model
   constructors default to 0.5, so this does not breach the "never change the cutoff default" rule.
   Direct callers who pass `cutoff` explicitly keep that control.
2. Store `cutoff` on the returned `bl_grid`. This is a new field; no field is removed.
3. In `bl_assemble()`, stop with a clear message if `bl_grid$cutoff` and `bl_model$cutoff` differ.
   Skip the check when `bl_grid$cutoff` is `NULL`, so grids built before this field existed still
   work.
4. `bl_build_result()` then inherits the model's cutoff with no change of its own. Its roxygen gains
   one sentence: "the grid uses `bl_model$cutoff`".

**Recommendation:** fix now. It is small, and it turns silent disagreement into either agreement or
an error.

---

## Issue 3: models from `xgboost::xgboost()` fail in `.pred_function()`

**Background.** Verified this session against the installed xgboost 3.2.1.1.
- The `"XGB"` branch (`R/predict_utils.R:124-128`) expects
  `list(model = <xgb.Booster>, features = ...)`. It predicts with
  `predict(model, newdata = xgb.DMatrix(...))`.
- In xgboost >= 3.0, `xgboost::xgboost()` is the recommended high-level interface. It returns
  class `c("xgboost", "xgb.Booster")`. Its `predict.xgboost()` method refuses DMatrix input:
  > Predictions on 'xgb.DMatrix' objects are not supported with 'xgboost' class. Try 'xgb.train'
  > or 'predict.xgb.Booster'.
- So `bl_wrap_model(list(model = <xgboost() fit>, features = ...), "XGB", train_data = ...)` errors
  inside `bl_wrap_model()`. Without `train_data`, it errors later, in `bl_build_grid()`.
- A second trap (verified): given a numeric 0/1 `y`, `xgboost()` fits a *regression*
  (`objective = "reg:squarederror"`), so its "probabilities" are not probabilities. Given a factor
  `y` it fits `binary:logistic`. `predict(model, newdata = <data.frame>)` then returns P(second
  level), and `attributes(model)$metadata$y_levels` holds `c("0", "1")`.
- Current state: the roxygen (`R/model_fit.R:160-164`), README:234 and
  `review_section4_to_6.md:328` say `xgboost()` models are not supported. The workaround is
  `model_type = "custom"` with a `predict_fn`.

**Effect on the paper:** none. Pima uses a GAM, and the loan harness uses `xgb.train()` and the
custom path.

**Suggested fix (keeps the required `list(model, features)` format, CLAUDE.md §4).**
1. In the `"XGB"` branch of `.pred_function()`:
   - if `inherits(model_use$model, "xgboost")`, call
     `predict(model_use$model, newdata = new_data[, model_use$features, drop = FALSE])` on the
     data frame;
   - otherwise keep the existing DMatrix path for `xgb.train()` boosters.
2. In `bl_wrap_model()`, when `model_type = "XGB"` and the model is of class `"xgboost"`, check once,
   at wrap time:
   - `attr(model$model, "params")$objective` must be `"binary:logistic"`;
   - `attr(model$model, "metadata")$y_levels` must be `c("0", "1")`.

   If either check fails, stop with a message saying to refit with `y = factor(class)`.
3. Update the roxygen note, README:234 and walkthrough :328 to say that both `xgb.train()` and
   `xgboost()` boosters are supported, and that `xgboost()` needs a factor `y`.
4. I rejected calling the internal `xgboost:::predict.xgb.Booster` directly. It works, but it is
   unexported, so `R CMD check` would flag it.

**Caveat:** step 2 reads xgboost's object attributes (`params`, `metadata`). That is how xgboost's
own print method reads them, but it is not a documented API. If you would rather not depend on it,
a lighter variant does step 1 only and documents the factor-`y` requirement. The grid's existing
[0, 1] clamp warning then flags a regression fit.

**Recommendation:** your call; it is independent of the paper. The effort is small (about 30 lines
plus tests). Fixing it is user-facing value, because `xgboost()` is now xgboost's main interface.
Deferring it costs only a known-issues entry.

---

## Issue 4: no guard for an empty outcome group

**Background.**
- With `method = "CVA"` and a model, `bl_build_projection()` labels each training row TP, TN, FP or
  FN (`R/projection.R:231-246`). It uses those labels as the CVA groups. g groups give g - 1
  canonical dimensions, so a 2D plane needs at least 3 non-empty groups.
- `case_when()` creates only the labels that occur. If the model makes, for example, no false
  positives on the training data, the factor has 3 levels.
- **(a) Wrong colours, silently.** `point_col` (`R/projection.R:289-295`) uses the confusion
  palette only when the levels are exactly `{FN, FP, TN, TP}`. Otherwise it falls through to the
  binary rule `cva_classes == "1"`. No confusion label equals "1", so every point is blue.
  - This affects `plot.bl_projection()` only, that is `plot(bl_proj)` called directly.
  - `plot_biplotEZ()`, `plot.bl_boundary()` and `plot.bl_local_result()` colour each observation
    themselves, so they are unaffected.
- **(b) No message.** Nothing says which group is empty. With only 2 groups left (for example, a
  perfect fit on the training data), the second CVA axis is arbitrary. The user sees only the
  generic "CVA works best with at least 3 class levels" warning.
- **(c) Unused factor levels.** If you pass `cva_classes` yourself with an unused level,
  `.compute_metric_inverse()` counts it in g (`R/projection.R:43,51`). The denominator of `W` is
  then n - g - 1 instead of n - g. That rescales every Mahalanobis distance by a constant: the
  ranking is unchanged, but the values shift. Unverified: whether biplotEZ's `CVA()` itself copes
  with an empty level. I will check this before writing the fix.

**Effect on the paper:** none expected. The Pima assertions require all four groups to be non-empty,
and the Task C plan already stops the script if one is empty. This is a package guard.

**Suggested fix.**
1. Build the confusion labels as
   `factor(..., levels = c("FN", "FP", "TN", "TP"))`, the alphabetical order `as.factor()` gave;
   see the implementation log for why the order matters. If any group has zero rows, `message()`
   names it, for example "CVA groups: FP empty (FN 30, TN 200, TP 120); using 3 groups.", then
   call `droplevels()`.
2. Call `droplevels()` on a user-supplied `cva_classes` too, before CVA and before the metric
   (closes (c)).
3. `point_col`: use the confusion palette whenever every level is in `{TP, TN, FP, FN}` (any
   subset), and keep the binary rule for 0/1 labels (closes (a)).
4. Keep the existing `< 3 levels` warning as a warning. Turning it into an error would break the
   model-free CVA path, which uses binary labels, and the tests behind the 4 baseline WARNs.

**Recommendation:** fix now. It is small, it adds a message rather than changing results, and it
closes a silent wrong plot.

---

## Files to read first (at implementation)

- `R/shapley.R` (`plot.bl_sparse_result`), `R/local_cf.R` (`plot.bl_local_result`, rotation
  block), `R/plot_biplot.R` (`.apply_biplot_rotation`), and the existing local-plot tests in
  `tests/testthat/test-local_cf.R` and `test-graphics_device.R`.
- `R/biplot_grid.R` (`bl_build_grid`), `R/result.R` (`bl_assemble`, `bl_build_result`),
  `R/boundary.R:300-350`, `R/surrogate.R:180-200`.
- `R/predict_utils.R`, `R/model_fit.R` (`bl_wrap_model`).
- `R/projection.R` (`bl_build_projection`, `.compute_metric_inverse`, `plot.bl_projection`).
- Callers: `grep` for each changed function first (CLAUDE.md §5).

## Error cases

| case | action |
|---|---|
| `bl_grid$cutoff != bl_model$cutoff` in `bl_assemble()` | stop, naming both values |
| old `bl_grid` with no `cutoff` field | skip the check (back-compatible) |
| `xgboost()` model with a regression objective or levels other than 0/1 | stop in `bl_wrap_model()`, telling the user to refit with `y = factor(class)` |
| `features` not all in `new_data` | the existing subsetting error, unchanged |
| a CVA group is empty | `message()` naming it, then `droplevels()`; still a warning if fewer than 3 groups remain |
| `rotate_deg` passed to `plot(bl_local)` or `plot(bl_sparse)` | stop before any device opens, saying it is for the global biplot only |

## Test cases

- **Issue 1** (`test-local_cf.R`):
  - `plot(bl_local, rotate_deg = 30)` errors with a message naming `rotate_deg`, and
    `dev.cur()` is unchanged (no device opened);
  - `plot(bl_sparse, rotate_deg = 30)` errors the same way;
  - `plot(bl_local)` and `plot(bl_sparse)` still run on `pdf(NULL)`;
  - `plot_biplotEZ(bl_results, rotate_deg = 30)` still runs, so the global path keeps the
    argument (add the test if none exists).
- **Issue 2** (`test-biplot_grid.R`, `test-result.R`):
  - a grid built from a model with `cutoff = 0.4` stores `cutoff = 0.4`, and its `ct` levels are
    0.399 and 0.401;
  - `bl_assemble()` stops when the grid and model cutoffs differ;
  - `bl_build_result()` with a 0.4 model gives a grid at 0.4;
  - the 0.5 path gives byte-identical `ct` to before.
- **Issue 3** (`test-predict_utils.R`), skipped if xgboost is not installed:
  - an `xgboost()` fit with factor `y` wraps, and its predictions match `predict(m, df)`
    floor-rounded;
  - an `xgb.train()` booster still works;
  - an `xgboost()` fit with numeric `y` stops at wrap time.
- **Issue 4** (`test-projection.R`):
  - a model or label set with no FP gives a message naming FP, 3 levels, and confusion colours,
    not all blue;
  - the confusion levels stay alphabetical, as `as.factor()` gives;
  - a user factor with an unused level gives `metric_type = "W_cva"`, and `W` is divided by n - 3;
  - four non-empty groups give no message.

## Documents to update (per fix, as Gate 2 requires)

- **For each fix:**
  - a NEWS entry and roxygen, then `devtools::document()`;
  - the walkthrough passage that now describes the bug, rewritten to present-tense current
    behaviour:
    - `review_section4_to_6.md:374` (issue 2);
    - `review_section4_to_6.md:328` (issue 3);
    - the section 4 `bl_build_projection()` passage (issue 4).
- **Issue 1: the `rotate_deg` note you asked for.** Its substance: "`rotate_deg` is only for
  viewing the global (optimal) biplot. Whenever the rotation function `.bl_rotate()` is used (the
  local biplot for an individual target, `plot(bl_local)` and `plot(bl_sparse)`), the plot is not
  turned further: `rotate_deg` is in effect always 0 there. These functions have no `rotate_deg`
  argument and stop if one is passed. It is never used to plot a rotated individual point."
  - `review_section9_to_15.md`:
    - add the note after the `plot(bl_local)` Parameters table, and delete the table's
      `rotate_deg` row (:534);
    - drop "(and any extra `rotate_deg`)" from the raw-unit relabel paragraph (:544);
    - replace the bug description at :861 with one line: the sparse cross is projected with
      `Vr_rot`, the same frame as the rest of the local plot.
  - `2 implementation_summary.txt`:
    - add the note at §4.4 after the `rotate_deg` sentence (:254);
    - remove `rotate_deg (applied on top of the SVD rotation)` from the list of parameters shared
      with `plot_biplotEZ()` (:425), and add the note there.
  - Also:
    - `review_section4_to_6.md:535`: the `plot_biplotEZ()` `rotate_deg` row gains "global biplot
      only; not available on the local plots";
    - `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md:580`: that bullet
      describes `rotate_deg` as an argument of the local plot. Reword it to say the argument
      exists only on the global plots, for the reason above;
    - the roxygen for `plot_biplotEZ()` and `plot.bl_surrogate()` `rotate_deg` gains "global
      biplot only";
    - NEWS: one entry saying `rotate_deg` was removed from `plot.bl_local_result()` and
      `plot.bl_sparse_result()`, and why;
    - CLAUDE.md §4: the line from Fix 1, step 4.
- **Issue 2:** `2 implementation_summary.txt` gains the new `bl_grid$cutoff` field (an S3 field
  change).
- **Issue 3:** the README:234 sentence.
- **Each deferred issue:** it goes to Task D's `known-issues.md` instead.

## Verification

- `devtools::document()` then `devtools::test()` through a scratchpad `.R` driver (PowerShell +
  Rscript). Expect 250 + new PASS, 0 FAIL, and the same 4 baseline WARNs.
- Run the loan harness end to end from the scratchpad driver. Every number quoted in the
  walkthroughs must be unchanged (all fixes are no-ops at cutoff 0.5, and the harness passes no
  `rotate_deg` to a local plot).
- Run `grep -rn rotate_deg R/ tests/ vignettes/ development/ README.md`. It should show
  `rotate_deg` only on `plot_biplotEZ()`, `plot.bl_surrogate()`, `.apply_biplot_rotation()`, the new
  error check, its tests and the docs.
- `devtools::check(document = FALSE, args = "--no-manual")`: 0 / 0 / 0.
- Stop and report. Commit nothing unless you ask.

## Assumptions

- Session 2's uncommitted work should be committed first (your call), so these fixes land as a
  separate, reviewable commit.
- These fixes go in before Task C's first run, so that run already includes them. If you would
  rather keep the Gate 2 order (run Task C first, then decide), this plan waits.
