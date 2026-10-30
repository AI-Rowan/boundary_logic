> **IMPLEMENTED 2026-10-29** on branch `paper_nr1_changes`. Archived copy of the approved plan; the technical note is `documentation/sparse-counterfactual.md`.

# Sparse CF: "No effect" class, greedy feature selection, and a round_to note

Date created: 2026-10-29
Date last edited: 2026-10-29

(This file earlier held four finished 2026-10-29 tasks. The last one, the fixed-feature restore,
is archived at `.claude/reference/fixed-constraint-restore-plan.md`. The file now holds this new
task.)

## Context

The user ran `bl_find_sparse_cf(bl_shapley_values, round_to = NULL)` (script Step 14) and asked
why `credit_score` still changes (634 → 634.22) when its Shapley value is 0, and why
`round_to = 1` moves it further.

**1. Why a zero-effect feature changes.**
- `bl_shapley()` labels a class-0 target's feature `"Supports"` when `shapley_cause >= 0`
  (`R/shapley.R:218`). The rule is inherited from the PhD script
  (`origin/original_PhD_code:scripts/4.4 Local Interpretation Shapley.R:102`).
- `bl_find_sparse_cf()` then changes every "Supports" feature.
- The value is exactly 0 because the move crosses no XGBoost split, and predictions are
  floor-rounded to 3 d.p.

**2. Why `round_to = 1` moves it further.** `round_to` rounds the stored values, which are
z-scores for `bl_results_v2`. `credit_score`'s z of −0.221 rounds to 0, the training mean,
643.38.

**User's decisions (2026-10-29):**
1. Add a `"No effect"` class for zero Shapley values.
2. Keep `round_to` as it is. Write a technical note on it for future investigation.
3. Remove `round_to = NULL` from the script's example calls (`NULL` is the default).
4. **Greedy sparse CF:** change the "Supports" features one at a time, largest Shapley value
   first. Re-score the model after each, and stop at the first point that crosses the boundary.
   Features not needed stay at their observed values.
5. Store the path and print it under the sparse table.
6. `print(bl_shapley_values)` takes the same layout as `print(bl_local)` and `print(bl_sparse)`.

## Changes

### 1. `bl_shapley()` (`R/shapley.R`)

- **The rule:** first `case_when()` branch: `abs(df$shapley_cause) < 1e-8 ~ "No effect"` (zero
  up to floating-point noise). The Supports/Contradicts branches follow unchanged.
- **`plot.bl_shapley()`:** add `"No effect" = "grey45"` to both manual scales. It is distinct
  from Contradicts' `grey70`, and the white value label stays readable.
- **Roxygen:** the three classes, and that zero was "Supports" in the PhD rule.

### 2. `bl_find_sparse_cf()` (`R/shapley.R`): greedy selection

1. `fixed_vars` = the `"fixed"` features in `bl_local$set_filters`.
2. **Candidates:** "Supports" features not in `fixed_vars`, ordered by `abs(shapley_cause)`,
   largest first (`order()` is stable, so ties keep the Shapley-table order). "No effect",
   "Contradicts", "Unknown" and fixed features never change.
3. **Loop:** start from `x_sparse = x_obs`. For each candidate:
   - set its counterfactual value, rounded by `round_to` exactly as now (model units, nearest
     multiple);
   - score with `.pred_function()`;
   - append a row to `sparse_path` (`step`, `variable`, `shapley`, `pred`, `crosses`);
   - stop when the predicted class differs from the target's.
4. **Result:**
   - `pred_sparse` is the last score; with no candidates it is the score at `x_obs`;
   - `solution_valid` comes from the last score;
   - if no step crosses, every candidate has been changed (the previous behaviour) and
     `solution_valid` is `FALSE`;
   - `used_in_sparse` = the feature appears in `sparse_path`;
   - `x_sparse_val` is unchanged in meaning.
5. **New field `sparse_path`** (data frame, zero rows when there are no candidates).
6. The old post-hoc revert loop for fixed features goes: they are never changed now. The fixed
   features are still read from `set_filters`.
7. **Roxygen:**
   - the greedy rule, and that it is not guaranteed minimal (a smaller subset outside the
     Shapley order may exist);
   - `used_in_sparse` and `sparse_path`;
   - `round_to` acts on model units (z-scores for a standardised model), see the technical note.

At most `p` extra model calls.

### 3. `print.bl_sparse_result()` (`R/shapley.R`)

After the table and the `(change = sparse - observed)` line:
```
  Path (Supports features, largest Shapley value first):
 step            variable shapley   pred crosses
    1       loan_int_rate  0.2915 0.xxxx   FALSE
    2 loan_percent_income  0.1445 0.5470    TRUE
```
With no candidates, it prints: `  Path: no "Supports" feature to change.` The shapley and pred
columns are rounded to 4 d.p.

### 3b. `print.bl_shapley()` (`R/shapley.R`): the same layout as `print(bl_local)`

The summary labels follow `print(bl_sparse)`, and the table uses the shared
`.cf_compare_table(x_obs, list(counterfactual = B_x), ...)` plus two columns: `shapley`
(`shapley_cause`, 4 d.p., prediction units, not converted) and `contribute`.
```
-- bl_shapley --
  Row ID         : 13
  Observed pred  : 0.1110  (class 0)
  Full CF pred   : 0.5470
  Shapley sum    : 0.4360

  Observed vs counterfactual, with Shapley values (raw units):
                    observed counterfactual  change constraint shapley contribute
loan_amnt           12000.00     12000.0000  0.0000      fixed  0.0000  No effect
loan_int_rate          12.69        14.1460  1.4560   increase  0.2915   Supports
loan_percent_income     0.05         0.0484 -0.0016       none  0.1445   Supports
credit_score          634.00       634.2206  0.2206   increase  0.0000  No effect
```
- **Rows** follow `var_names`, like the other two tables; the plot keeps its Shapley order.
- **`Shapley sum`** lets the reader check that the values add up to Full CF pred minus
  Observed pred.
- **Constraints** come from `bl_local_result$set_filters`.
- **Header:** "(model units)" without a scaling. The existing tests (`test-set_scaling.R:196`,
  `206`) check "raw units" present or absent, which still holds.
- **Roxygen:** `print.bl_shapley()` describes the layout.

### 4. Tests (`tests/testthat/test-local_cf.R`)

These use synthetic `bl_local_result` objects with custom models (the restore-fixture pattern),
run through the real `bl_shapley()` and `bl_find_sparse_cf()`:
- **"No effect":** a model `plogis(4 * a)` ignores `b`. `b` is "No effect", stays at observed,
  and is not used.
- **Greedy stops early:** a model `plogis(4 * a + c)`, observed `(-0.2, -0.3)`, end `(0.3, 0.2)`.
  `a` alone crosses (0.711), so the path has 1 step, `c` stays at observed, and `used` is
  `FALSE` for `c`.
- **Greedy needs both:** a model `plogis(a + c)`, observed `(-0.4, -0.4)`, end `(0.35, 0.25)`.
  `a` alone gives 0.487 (no cross); both give 0.645. The path has 2 steps, and the last one
  crosses.
- **Print:** `print(sparse)` contains "Path (Supports features".
- **`print(shp)`:**
  - it has the header "Observed vs counterfactual, with Shapley values" and the six column
    names;
  - its rows follow `var_names`;
  - the "Shapley sum" equals `sum(shapley_cause)` to 4 d.p.
- The `plot()` of the Shapley object builds without error (`skip_if_not(ggplot2)`).
- The existing iris test (`used_in_sparse` `FALSE` for a fixed feature) still holds.

### 5. Script (`development/03_loan_status_Boundary_Logic.R`)

Remove `round_to = NULL` from the `bl_find_sparse_cf()` calls in Steps 14 and 15.

### 6. Technical note `documentation/sparse-counterfactual.md`

- **§1 "No effect" class (implemented):** the rule, why an exact 0 arises, the departure from
  the PhD `>= 0` rule, and the effects.
- **§2 Greedy selection (implemented):**
  - the rule, and why the model is re-scored rather than Shapley values summed (a Shapley value
    is an average over coalitions, not the effect when added first);
  - it is not guaranteed minimal;
  - the harness path as the example.
- **§3 `round_to` (open, for future investigation):**
  - current behaviour;
  - a real scratch run with `round_to = 1`: values, `pred_sparse`, `solution_valid`;
  - the issues: z-score units; one value for all scales; nearest rounding can cancel or reverse
    a move; constraints are not re-checked;
  - the options: raw units, a per-feature `round_to`, rounding away from observed, a constraint
    re-check.

### 7. Documents

- **`review_section9_to_15.md`:**
  - Step 13: the Contribute table gains "No effect"; the Shapley print and its sentence come from
    the new run;
  - Step 14: the Call block without `round_to`, the Parameters row (default; links the note), the
    rewritten algorithm (greedy loop), `sparse_path` in the returned-object table, and the print
    sample, path and text from the new run;
  - Step 15 code and the object flow;
  - Key point 4 (greedy by Shapley order with model re-scoring).
- **`2 implementation_summary.txt` §4.9:** the classification block; the `bl_find_sparse_cf()`
  steps; `sparse_path`; both prints (`print.bl_shapley()` and `print.bl_sparse_result()`); the
  `round_to` pointer.
- **CLAUDE.md §5:** the shared-builder sentence covers all three tables (`print.bl_local_result()`,
  `print.bl_shapley()`, `print.bl_sparse_result()`).
- **`1 Foundation intro documents.txt` (line ~266) and the vignette (line ~530):** the sparse CF
  description.
- **`documentation/fixed-constraint-restore.md` §5:** after a restore, a fixed feature reads
  "No effect", and the sparse CF follows the greedy path.
- **`NEWS.md`:** one entry for the two behaviour changes.
- **CLAUDE.md:**
  - §6: a one-line key fact on the greedy sparse CF;
  - §9: a Deferred item for `round_to`, and the test-coverage line.

## Files to read first

- `R/shapley.R` 140–600.
- `tests/testthat/test-local_cf.R`.
- `review_section9_to_15.md` Steps 13–15.

## Error cases

- **No candidates** (all features are "No effect", "Contradicts" or fixed): `x_sparse = x_obs`,
  `sparse_path` has 0 rows, `solution_valid` is `FALSE`, and the print says so.
- **No step crosses:** every candidate is changed, and `solution_valid` is `FALSE` (as before).
- **Rounding moves a value back to observed:** the step still counts and is re-scored, so the
  path shows the effect.
- **Older `bl_sparse_result` objects without `sparse_path`:** the print skips the path block.

## Verification

1. `devtools::document()` + `devtools::test()` (baseline 224 PASS / 0 FAIL / 4 WARN), then
   `R CMD check` 0/0/0. Report and stop (no commit).
2. Harness:
   - Step 13: `loan_amnt` and `credit_score` read "No effect";
   - Step 14: the path starts with `loan_int_rate`, and `credit_score` stays at 634. Read whether
     `loan_int_rate` alone crosses from the log.
3. A scratch run of `bl_find_sparse_cf(bl_shapley_values, round_to = 1)` on the harness objects
   supplies the note's §3 example.

## Implementation notes (2026-10-29)

Implemented as planned.

**Results**
- Tests: **250 PASS / 0 FAIL / 4 WARN / 0 SKIP** (was 224).
- `R CMD check`: 0 / 0 / 0.
- Harness (`harness_log5.txt` in the session scratchpad):
  - Step 13: `loan_amnt` and `credit_score` read "No effect";
  - Step 14 greedy path: `loan_int_rate` reaches 0.360 (no cross), then `loan_percent_income`
    reaches 0.547 (crosses); `credit_score` stays at 634;
  - `round_to = 1` scratch run: `loan_int_rate` 14.6504, `loan_percent_income` 0.0666 (direction
    reversed), path 0.384 → 0.483, `solution_valid = FALSE`.

**Display detail.** At the harness log's 80-column width, the `print(bl_shapley)` table is one
character too wide, so R prints `contribute` as a second block. A wider console shows it on one
line. The walkthrough shows the logged output and says so.

**Extra edits.** The README function table and the vignette's Step 15 prose now mention
"No effect".
