> **IMPLEMENTED 2026-10-29** on branch `paper_nr1_changes`. Archived copy of the approved plan; the technical note is `documentation/fixed-constraint-restore.md`.

# Restore "fixed" features after the local search; matching comparison tables in print(bl_local) and print(bl_sparse)

Date created: 2026-10-29
Date last edited: 2026-10-29

(This file earlier held three finished 2026-10-29 tasks: the walkthrough Parameters tables, the
surrogate move plus the pruning section, and the Step 11 script table. The file now holds this
new task.)

## Context

1. **Comparison table in the print method.** The harness builds an observed / counterfactual /
   change / constraint table after `print(bl_local)` (script Step 11, lines 539–560). The user
   wants `print(bl_local)` to produce it, shown last, after the existing summary, Target,
   Counterfactual and per-pair blocks.
2. **Why the `"fixed"` `loan_amnt` moved (12,000 → 12,398.5).** `.apply_actionability()`
   implements `"fixed"` as `|B_x - x_obs| <= 0.5` in model units (`R/local_cf.R:89`). That is a
   tolerance, not an equality, for a geometric reason: the counterfactual must be a vertex of a
   boundary contour in the 2D rotated biplot, and back-projecting it sets every feature from two
   coordinates. A vertex almost never sits on one feature's observed value exactly, so an
   equality would usually leave no candidate. For z-scored `bl_results_v2`, 0.5 is 0.5 training
   SD: ±3,001 for `loan_amnt`.
3. **User's rule.** Keep the search as it is. Afterwards, set each `"fixed"` feature back to its
   observed value and re-score. If that point still crosses the boundary, it becomes the selected
   counterfactual. Otherwise keep the search's point, and print a note explaining the move.
4. **`print(bl_sparse)`** gets the same table layout and headings as `print(bl_local)`, keeping
   the full CF column.
5. **`used_in_sparse` fix.** The flag is `TRUE` for a `"fixed"` feature even after the sparse CF
   puts it back (a Gate 2 known issue). It is fixed now in the stored object, so the object and
   the print agree.

Item 3 changes which counterfactual `bl_find_local_cf()` returns whenever a `"fixed"` feature
moved, so it gets the CLAUDE.md §2 triplet: technical note, plan archive and implementation
summary.

User's decisions (2026-10-29): the restore rule as in item 3; a note below the table; the
sparse table keeps the full CF column; fix the stored `used_in_sparse` now.

## Changes

### 1. `R/local_cf.R`: the restore step

- **One tolerance constant:** a private `.fixed_tol <- 0.5`, used by `.apply_actionability()` in
  place of the literal. The search is unchanged.
- **New private helper `.restore_fixed(best_result, x_obs, var_names, set_filters, bl_result,
  bl_target, sv)`**, called in `bl_find_local_cf()` after the pair loop and only when a solution
  exists:
  1. `moved` = the `"fixed"` features with `|B_x - x_obs| > 1e-8` (model units). If none, return
     `NULL` and leave everything unchanged.
  2. `B_x_r` = `B_x` with the `moved` features set to `x_obs`. `pred_r` =
     `.pred_function(model, model_type, B_x_r)`.
  3. `applied` is `TRUE` when the restored point crosses the boundary
     (`as.integer(pred_r >= cutoff) != bl_target$pred_class`) and passes the search's
     `train_ranges` check (`get_filter_logical_vector()`). Other constraints cannot break,
     because only the fixed features change.
  4. If applied: `B_x <- B_x_r`, `B_pred <- pred_r`, and `B_z <- ((B_x_r - X_center) / sv) %*%
     Vr_rot`. That is the same projection `bl_select_target()` uses; a point on the plane
     re-projects to itself.
- **New field `bl_local$fixed_check`** (`NULL` when no solution, or when no fixed feature moved):
  `list(vars, applied, pred_restored, search_point = list(B_x, B_z, B_pred))`. `dist_z`,
  `dist_mahalanobis` and the per-pair tables stay the search's values: they describe the pair
  selection. `bl_counterfactual` is built from the selected point, after restoring.
- **Verbose message** when applied: "Fixed feature(s) loan_amnt restored to observed value(s);
  the restored point still crosses the boundary (pred 0.547) and is used as the counterfactual."

### 2. One shared table builder for both prints (in `R/local_cf.R`)

- **`.constraint_labels(constraints, var_names)`:** returns each feature's constraint as given,
  a range as `"[lo, hi]"` (`%.4g`, as in `print.bl_filters()`), otherwise `"none"`.
- **`.cf_compare_table(x_obs, cols, var_names, scaling, constraints)`.** Rows are the features in
  `var_names` order. Columns:
  - `observed`, then each model-unit vector in `cols` (for example
    `list(counterfactual = B_x)`, or adding `sparse = x_sparse`), converted with
    `.scale_to_raw(..., "level")`;
  - `change` = last column minus observed, converted with `.scale_to_raw(..., "delta")`;
  - `constraint`.

  All numbers are rounded to 4 d.p. It returns the table plus `units` (`"raw"` / `"model"`).
- **Fixed-move notes**, built from `x$fixed_check` (§4).

### 3. The print methods

**`print.bl_local_result()`.** The existing blocks are unchanged and stay in order. Add the
summary line `Fixed restored : loan_amnt` when a restore happened. When a solution exists, finish
with:
```
  Observed vs counterfactual (raw units):
                      observed counterfactual   change constraint
  loan_amnt           12000.00     12000.0000   0.0000      fixed
  ...
  Note: ...                                   (only when a fixed feature moved)
```

**`print.bl_sparse_result()`** (`R/shapley.R`) uses the same builder with
`cols = list(counterfactual = B_x, sparse = x_sparse)`, adds a `used` column (the stored
`used_in_sparse`, matched by feature name), and ends with a footer line:
```
-- bl_sparse_result --
  Solution valid : TRUE
  Observed pred  : 0.1110  (class 0)
  Full CF pred   : 0.5470
  Sparse CF pred : 0.5470

  Observed vs sparse counterfactual (raw units):
                      observed counterfactual     sparse   change constraint  used
  loan_amnt           12000.00     12000.0000 12000.0000   0.0000      fixed FALSE
  loan_int_rate          12.69        14.1460    14.1460   1.4560   increase  TRUE
  loan_percent_income     0.05         0.0484     0.0484  -0.0016       none  TRUE
  credit_score          634.00       634.2206   634.2206   0.2206   increase  TRUE
  (change = sparse - observed)
```
The summary labels match the `plot(bl_sparse)` console summary. Rows now follow `var_names`
order, as in `print(bl_local)`, instead of the Shapley order.

### 4. The notes (one per fixed feature the search moved)

**Restore applied** (the expected harness case):
```
  Note: loan_amnt is "fixed", but the boundary search moved it to 12398.5032
        (+398.5032). The search lets a "fixed" feature move by up to 0.5 model
        units (+/- 3001.1775 here), because a point on the 2D decision boundary
        rarely holds a feature at exactly its observed value. Restoring
        loan_amnt to its observed value still crosses the boundary (prediction
        0.5470), so the restored point is the selected counterfactual.
```

**Restore not applied:** the same first three lines, then "Restoring loan_amnt to its observed
value gives prediction 0.4880, which does not cross the boundary, so the search's point is
kept." When the `train_ranges` check is what failed, that reason is given instead. Without a
scaling, the values are in model units and the bracket reads `(+/- 0.5 here)`.

### 5. `used_in_sparse` fix (`R/shapley.R`, `bl_find_sparse_cf()`)

`df_out$used_in_sparse <- df_out$Contribute == "Supports" & !(varnames %in% fixed_vars)`, where
`fixed_vars` are the `"fixed"` features in `bl_local$set_filters`. Update the roxygen field
description to match: `TRUE` when the sparse CF takes the feature's counterfactual value.

### 6. Roxygen

- `bl_find_local_cf()`: a **Fixed constraints** section (the tolerance, why, the restore rule),
  the `fixed_check` field, and the `B_x` / `B_z` / `B_pred` meaning (the selected
  counterfactual, after restoring).
- `set_filters()` `"fixed"` item: the tolerance, why, the restore step.
- `print.bl_local_result()` and `print.bl_sparse_result()`: the block order and the table.
- `bl_find_sparse_cf()`: `used_in_sparse`.

ASCII only in R source.

### 7. Tests (`tests/testthat/test-local_cf.R`)

**`.restore_fixed()` is tested deterministically.** Each test builds a small `bl_result`-like
list with a `"custom"` model whose `predict_fn` depends on chosen features:
- the prediction ignores the fixed feature: the restore is applied, `B_x[fixed] == x_obs[fixed]`,
  `B_pred` is re-scored, and `search_point` keeps the original;
- the prediction depends on the fixed feature: not applied, and the result is unchanged;
- no fixed feature moved: returns `NULL`;
- the restored value is outside `train_ranges`: not applied.

**The table builder and notes, on synthetic inputs:**
- the levels and the delta use the right conversion, and the rows follow `var_names` order;
- the constraint labels are right, including `"[lo, hi]"`;
- applied and not-applied notes have the right text, and there is no note without a move;
- the header says "(model units)" without scaling.

**Integration (iris pipeline):**
- `print(bl_local)` ends with the table, after "Per-pair distances";
- `print(bl_sparse)` shows "Observed vs sparse counterfactual" with the six headings;
- with a `"fixed"` filter, `used_in_sparse` is `FALSE` for the fixed feature;
- the no-solution case prints the Reason and no table;
- `fixed_check` is `NULL` without filters.

Integration tests that need a solution use `skip_if_not(solution_found)`, as the existing tests
do.

### 8. Script, documents, changelog

- **Script:** delete the Step 11 `cf_tab` block; add a one-line comment that `print(bl_local)`
  ends with the table and explains any restored `"fixed"` feature.
- **Technical note** `documentation/fixed-constraint-restore.md`:
  - why `"fixed"` is a tolerance in the search;
  - the restore rule and why it is post-selection;
  - limits: the restored point is off the 2D boundary and not necessarily the nearest valid
    point; the pair choice is unchanged; 0.5 model units means different things for raw and
    standardised features.
- **Plan archive:** after implementation, copy this plan to `.claude/reference/` with an
  IMPLEMENTED banner.
- **`2 implementation_summary.txt`:** the `bl_local_result` fields (`fixed_check`, `B_*` meaning),
  the line-363 `"fixed"` row, the print descriptions (line 401 and the sparse print), the
  raw-unit list at line 451, and `used_in_sparse`.
- **`1 Foundation intro documents.txt` line 224, README line 206, vignette line 453:** the
  `"fixed"` description gains the restore step.
- **CLAUDE.md:**
  - §5 raw-unit list: add `print.bl_local_result()`;
  - §6: one key-fact line on the restore step.
- **`NEWS.md`:** one entry covering all three changes.
- **`progress.md`:** remove the `used_in_sparse` item from the Gate 2 issues list (fixed).
- **`review_section9_to_15.md`:**
  - a new "Stage H — Restore fixed features" in the Step 11 algorithm;
  - the Step 11 console sample, with the table and note, from a new harness run (replacing the
    script-snippet subsection);
  - Steps 13 and 14 (Shapley table, the new sparse print and its text, without the
    `used_in_sparse` caveat) from the same run;
  - the object flow and key point 5.
- **`devtools::document()`** for the man pages.

## Files to read first

- `R/local_cf.R`: 77–102, 223–289, 292–716 (`bl_find_local_cf()`), 735, 1060–1167.
- `R/shapley.R`: `bl_find_sparse_cf()` (350–450), `print.bl_sparse_result()` (535–573), and the
  `plot.bl_sparse_result()` summary (about 520).
- `R/data_prepare.R` (`.scale_to_raw()`), `R/feasibility_utils.R:35`, `R/predict_utils.R:76`.
- `tests/testthat/test-local_cf.R`, `test-set_scaling.R` (print tests at 196–197 and 234–240).

## Error cases

- **No solution:** no restore, `fixed_check` is `NULL`, and print shows the Reason and no table.
- **No `"fixed"` constraint, or no move:** `fixed_check` is `NULL` and the results are identical
  to before.
- **The restored point does not cross the boundary, or is outside `train_ranges`:** the search
  point is kept, and the note says why.
- **Older objects** without `fixed_check` / `set_filters`: the prints treat missing fields as
  `NULL`.
- **Incomplete scaling:** `.scale_to_raw()` warns and returns model units; the header follows.
- **Existing tests:** `test-set_scaling.R:197` expects "raw units" in `print(sparse)`; the new
  header keeps "(raw units)".

## Verification

1. `devtools::document()` then `devtools::test()`. Baseline 176 PASS / 0 FAIL / 4 WARN plus the
   new tests, all passing. Report the counts and stop (no commit).
2. Harness run. Expected:
   - the Step 11 counterfactual becomes `loan_amnt` 12,000, `loan_int_rate` 14.146,
     `loan_percent_income` 0.0484, `credit_score` 634.2206, prediction 0.547;
   - the `print(bl_local)` table and the "restored" note print last;
   - `print(bl_sparse)` shows the new table with `loan_amnt` `used = FALSE`;
   - Step 15 (no filters) is unchanged: 0.5554 / 0.542.
3. Update the walkthrough numbers from that log; grep for any left-over 12398.5 and
   "Variable summary" in Steps 11–14.

## Implementation notes (2026-10-29)

Implemented as planned. Details settled during implementation:

- **Notes text:** built with `strwrap()` (width 78) rather than fixed line breaks. Numbers use
  `formatC(..., drop0trailing = TRUE)`, so the tolerance reads `+/- 3001.1777` (0.5 × 6002.355)
  and a model-unit tolerance reads `+/- 0.5`. With several moved fixed features, the tolerance
  reads `+/- a for x, +/- b for y`.
- **`fixed_check` gained `in_range`,** so the note can say why a restore was not applied.
- **Tests:** 48 new expectations in `test-local_cf.R`.
  - `.restore_fixed()`: four deterministic tests with a custom-model fixture.
  - Table builder: two tests.
  - Notes: one test covering three outcomes plus `NULL`.
  - Iris integration: three tests.
- **Extra edits beyond the plan:** CLAUDE.md §9 (Phase 3 test coverage); the
  `2 implementation_summary.txt` §5.x `"fixed"` rationale paragraph (line 544).

**Verification**
- `devtools::document()` + `devtools::test()`: **224 PASS / 0 FAIL / 4 WARN / 0 SKIP** (was 176).
- `R CMD check`: 0 / 0 / 0.
- Harness run (log `harness_log4.txt` in the session scratchpad):
  - Step 11 restores `loan_amnt` to 12,000 at 0.547, and the table and note print last;
  - Step 14's sparse table shows `loan_amnt` `used = FALSE`;
  - Step 15 is unchanged (0.5554 / 0.542).
