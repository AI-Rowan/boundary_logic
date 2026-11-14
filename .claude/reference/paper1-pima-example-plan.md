> **ARCHIVED 2026-11-07.** The Task C plan and revision log of the Paper 1 Pima example
> (`scripts/paper1-pima/`, untracked): the only written specification of those scripts. Implemented;
> the open items (the paper sentences for issues 5, 9, 11, 12 and 13, and the port to the paper repo,
> T4) are in `progress.md`, "Outstanding items".

# Task C: the Paper 1 Pima example, on the current code

Date created: 2026-10-31
Date last edited: 2026-11-04

The authoritative Task C plan. Reviewed in plan mode on 2026-11-02 (as a plan-mode file (removed 2026-11-07))
and approved with the revisions below.

This plan details Task C of `paper1-task-spec.md`, as updated for the code at `49c9067`
(sessions 3-5). It supersedes the Task C section of `paper1-working-plan.md` where the
two differ. Branch: `paper_nr1_changes`. No package code changes before Gate 2.

**Revised 2026-11-02 after your review:**
- **The CVA distance identity is no longer recorded or asserted.** The paper explains why Z-space and
  X-space distances agree, and gives the mathematical reason.
- **A comment on the `"fixed"` restore is added:** the restored counterfactual lies off the rotated
  plane. This is accepted because the point is tested to be a counterfactual. The comment includes
  how much the restore changes the distance and the prediction.
- **Run 6 is added:** `Glucose = "decrease", BMI = "decrease"`, for both targets, with a Table 1
  column and figure panel 07c.

## Context

Paper 1 (JCGS) takes every figure and number from this package. Task C writes the paper's Pima
example as standalone scripts, runs them against the working tree, and stops at Gate 2 so you can
decide, issue by issue, what to fix in the package. The scripts are later copied to the paper's repo
and re-run against the installed package, so they must run unchanged there apart from the one line
that loads the package.

### What changed since the prompt was written

| prompt says | now | consequence for Task C |
|---|---|---|
| `max_pairs = 21L` must be passed | `21L` is the default (`49c9067`) | scripts rely on the default; run 5 still passes `max_pairs = 1L` |
| results depend on call history; quote first calls only; fresh `bl_results` per run | the searches are deterministic (`f0e3dd4`) | **your decision 2026-10-31:** one shared `bl_results`; `tests.R` proves identity three ways (below) |
| targets nearest the median of "a" global distance (`dist_x` planned) | Paper 1 compares by Mahalanobis only | **your decision 2026-10-31:** rank by the global counterfactual's Mahalanobis distance (below) |
| Task E measures the pair selector (Mahalanobis vs Euclidean) | under CVA the two are identical, d²_M = (n−g)·dist_z² | **your decision 2026-11-02:** not measured, recorded or asserted. The paper explains that Z-space and X-space distances are consistent and gives the derivation (rotation note §3.4); Task E cites it. |
| `"fixed"` may move ±0.5 | the search allows ±0.5, then `.restore_fixed()` puts it back when the point still crosses (`fixed_check`) | Table 1 and `numbers.rds` record `fixed_check$applied` per run; **your decision 2026-11-02:** a comment on the restored point lying off the rotated plane (see "The `"fixed"` restore and the rotated plane") |
| `bl_robustness()` for global importance | CLAUDE.md §8: redundant after `plot(bl_bnd)` | **CLAUDE.md wins:** take `sum_of_distance`, `robustness` and `mahalanobis_check` from `plot(bl_bnd)`'s return value; no `bl_robustness()` call |
| Shapley 0 counts as "Supports"; sparse CF = all Supports | "No effect" class; greedy sparse path (`sparse_path`) | Table 1 and `numbers.rds` record the sparse path |
| print output | fixed 2-dp tables | console logs quote cleanly |

## Files to read first

- `.claude/reference/paper1-task-spec.md`, Task C (the fixed set-up, runs, figures, assertions).
- `documentation/select-target-class/select_target_class_example.R`: the Pima set-up already runs
  there (392 rows, 313 train / 79 test, GAM via `model_type = "GAM"`). Reuse its set-up block.
- `R/boundary.R` `bl_find_boundary()` (fields `B_x`, `x_obs`, `pred_obs`, `class_obs`, `dist_x`);
  `R/boundary_plot.R` `plot.bl_boundary()` (return list, `mahalanobis_check`).
- `R/local_cf.R`: `bl_select_target()`, `set_filters()`, `bl_find_local_cf()` (return fields,
  `fixed_check`, `all_distances*`), `plot.bl_local_result()`.
- `R/shapley.R`: `bl_shapley()` (stops when `solution_found = FALSE`), `bl_find_sparse_cf()`
  (`sparse_path`, fixed handling), `plot.bl_shapley()` (returns a ggplot), `plot.bl_sparse_result()`.
- `R/plot_biplot.R` `plot_biplotEZ()` (`target_point`, `target_label`, `new_title`).
- `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md` §3.4: the CVA identity, which
  the paper explains (no test in Task C).
- `documentation/fixed-constraint-restore.md`: why the restored point leaves the plane, for the
  restore comment.

## Set-up (fixed; unchanged from the prompt)

- `inst/extdata/pima_diabetes.csv`; drop rows with a zero in Insulin, SkinThickness, Glucose or BMI;
  assert 392 rows.
- `bl_prepare_data(data, class_col = "Outcome", feature_cols = <7>, train_fraction = 0.8,
  seed = 121L, hull_fraction = 1)`; p = 7 (BloodPressure dropped via `feature_cols`).
- `mgcv::gam(class ~ s(Pregnancies) + s(Glucose) + s(SkinThickness) + s(Insulin) + s(BMI) +
  s(DiabetesPedigreeFunction) + s(Age), family = binomial)` on `train_data`, raw units; wrapped with
  `bl_wrap_model(model_type = "GAM")`; cutoff 0.5.
- `bl_build_result(bl_dat, bl_mod, method = "CVA")`: CVA on TP/TN/FP/FN, K = 3.
- Local search: package defaults (`max_pairs = 21L`, `distance = "mahalanobis"`).

## Files (all in `scripts/paper1-pima/`, gitignored, local only)

| file | does |
|---|---|
| `00_setup.R` | `library(mgcv)`; **the one repo-specific line** `devtools::load_all("../..", quiet = TRUE)` (paper repo: `library(boundarylogic)`); `out_dir <- Sys.getenv("PAPER1_OUT", "output")`, created if missing; `meta` list (package version, git SHA + dirty flag or `"installed"`, R and mgcv versions, time); helpers below. Guarded so sourcing it twice is harmless. |
| `01_data-model.R` | filter (assert 392), `bl_prepare_data()`, GAM, `bl_wrap_model()`; asserts the wrapped predictions equal `floor(predict.gam(type = "response") * 1000) / 1000` on `test_data` (else switch to `"custom"` + `predict_fn` and record which); saves `output/data-model.rds` |
| `02_global.R` | `bl_results` (built once); `bl_find_boundary()` on `train_data` and `test_data`; `plot(bl_bnd_train, type = "boxplot")` return value for importance; saves `output/global.rds` |
| `03_local.R` | target choice; `run_local()`; runs 1-6 for each target; builds Table 1 and `numbers.rds`; saves `output/local.rds` |
| `fig-03.R` ... `fig-09.R` | one figure each; each panel written as its own PDF (6 x 6 in) and PNG (300 dpi) to `output/`, e.g. `fig-05a.pdf`, `fig-05b.pdf`, so the paper lays out sub-figures itself. Each file starts with `source("00_setup.R")` (guarded) and `readRDS()`s what it needs, so Quarto's `#| file:` can include it alone. |
| `run_all.R` | sources 01, 02, 03, then fig-03 ... fig-09, in order; writes `output/run_log.txt` via `sink()` |
| `tests.R` | the assertions (below); writes `output/tests.rds`; stops with a non-zero exit if any failed |
| `issues.csv` | the issues table (symptom, where, effect on the paper, proposed fix, decision blank); filled after the first run |
| `review.qmd` → `review.html` | self-contained page: meta (SHA), figures, Table 1 with the restore comment, one pass/fail line per assertion, the issues table. Reads `output/` only; computes nothing. |

### Helpers in `00_setup.R`

- `fresh_result()`: rebuilds `bl_results` from `output/data-model.rds` (`bl_build_result(bl_dat,
  bl_mod, method = "CVA")`). Used only by `tests.R` check 8(b).
- `run_local(bl_results, target_row, filters = NULL, max_pairs = 21L)` returns a list:
  1. `tgt <- bl_select_target(bl_results, target = target_row)` (an index into `test_data`, the default);
  2. `flt <- if (length(filters)) do.call(set_filters, c(list(tgt), filters))`;
  3. `bl_local <- bl_find_local_cf(bl_results, tgt, set_filters = flt, max_pairs = max_pairs)`;
  4. if `bl_local$solution_found`: `bl_shapley_values <- bl_shapley(bl_local)` and
     `bl_sparse <- bl_find_sparse_cf(bl_shapley_values)`; otherwise both `NULL` (`bl_shapley()`
     stops on no solution).
- `extract_numbers(run)`: the quotable numbers of one run (below), so `03_local.R` and `tests.R`
  extract the same fields.
- `mahal(x_from, x_to, metric_inv)`: per-row `sqrt(rowSums((v %*% metric_inv) * v))` with
  `v = as.matrix(x_to) - as.matrix(x_from)`. A distance from package fields, not a re-derived plane
  or projection. It has two uses:
  - the target rule, with `bl_bnd_test$x_obs`, `bl_bnd_test$B_x` and
    `bl_results$metric_inv`;
  - the restore comment, the target to the search point (`fixed_check$search_point$B_x`) against the
    target to the restored point (`B_x`).
- `save_fig(name, expr, gg = FALSE)`: opens `pdf()`, evaluates the plot (wrapped in `print()` for a
  ggplot), closes **its own** device, then repeats for `png()`.

## Logic steps

### `02_global.R`

1. `bl_results <- bl_build_result(bl_dat, bl_mod, method = "CVA")`; assert all four outcome groups
   are non-empty (`table()` of the confusion group on `train_data`).
2. `bl_bnd_train <- bl_find_boundary(bl_results, data = bl_results$train_data)` and
   `bl_bnd_test <- bl_find_boundary(bl_results, data = bl_results$test_data)` (CLAUDE.md §6: pass
   `data =` explicitly).
3. `imp <- plot(bl_bnd_train, type = "boxplot")` inside `save_fig("fig-04b", ...)` if the method
   draws itself, else store `imp$plot` for `fig-04.R` (checked at implementation). Keep
   `sum_of_distance`, `robustness`, `mahalanobis_check`.
4. Record: confusion counts (train, test), counterfactuals found / NA per data set, how many global
   `B_x` fall outside `train_ranges` (reported, not asserted: `bl_find_boundary()` keeps them by
   design), median `B_pred`.

### `03_local.R`

1. **Targets.** From `bl_bnd_test`: `pred_class = pred_obs >= cutoff`; TP = class 1 & pred 1,
   TN = class 0 & pred 0; `d_M <- mahal(bl_bnd_test$x_obs, bl_bnd_test$B_x, bl_results$metric_inv)`;
   drop `NA`; in each group take
   `which.min(abs(d_M - median(d_M)))` (first on ties, so deterministic). Record the row, `x_obs`,
   `pred_prob`, `d_M` and the group median.
2. **Runs** for each target (12 in all):

   | run | filters | `max_pairs` |
   |---|---|---|
   | 1 unconstrained | none | 21 |
   | 2 cannot decrease | `Pregnancies = "increase", Age = "increase"` | 21 |
   | 3 both fixed | `Pregnancies = "fixed", Age = "fixed"` | 21 |
   | 4 mixed (the paper's main constrained example) | `Pregnancies = "fixed", Age = "increase"` | 21 |
   | 5 pinned pair (1,2) | none | 1 |
   | 6 cannot increase (the `"decrease"` example) | `Glucose = "decrease", BMI = "decrease"` | 21 |

   Run 6 is numbered after run 5 so that the prompt's run numbers (figure 6(b) uses run 5) stay
   valid. For the TN target the counterfactual moves towards diabetic, against the direction run 6
   allows, so it may find no solution. That is a result to report, not a reason to change the run.
3. **`numbers.rds`**, one list:
   - `meta` (excluded from the reproducibility comparison);
   - `data`: 768 / 392 / n_train / n_test / p, class counts;
   - `model`: wrap path used, test confusion counts;
   - `cva`: the outcome groups and their sizes;
   - `global`: per data set the counts above, `sum_of_distance`, `robustness`, `mahalanobis_check`;
   - `targets`: TP and TN as in step 1;
   - `runs[[target]][[run]]`: `solution_found`, `blocking_constraint`, `best_pair`,
     `sqrt(dist_mahalanobis)`, `dist_z`, `all_distances_mahalanobis` (21 values), `B_x`, `B_pred`,
     `fixed_check` (vars, applied, pred_restored, the search point's `B_x` and `B_pred`),
     `restore_impact` (below), Shapley table (`varnames`, `shapley_cause`, `Contribute`),
     `x_sparse`, `pred_sparse`, `solution_valid`, `sparse_path`.
4. **Table 1** (`output/table1.csv`, plus a formatted version in `review.html`). Per target:
   - **Rows:** the 7 features, then `B_pred`, √d²_M of the search, √d²_M of the restored point (blank
     unless a restore was applied), `best_pair`, `solution_found`, and whether the fixed restore was
     applied.
   - **Columns:** observed, runs 1-4 and 6, and the sparse counterfactual of run 1.
   - **Units and rounding:** raw clinical units (no scaling is set, so model units are raw units);
     stored unrounded, shown to 2 dp.
   - **Under the table:** the restore comment below.

### The `"fixed"` restore and the rotated plane

Runs 3 and 4 hold Pregnancies (and in run 3 also Age) fixed. The search accepts a fixed feature
within ±0.5. `.restore_fixed()` then sets it back to its observed value.

**What happens when `fixed_check$applied` is `TRUE`:**
- **`B_x` is the restored point, and it lies off the rotated plane.** Changing one feature moves the
  point out of the two-dimensional plane the search worked in.
- **The local biplot shows its projection.** `B_z` is the projection onto the plane, and it need not
  sit on the drawn boundary contour. Figure 7(b) (run 4) shows this.
- **The distances are the search's.** `best_pair` and √d²_M come from the on-plane search point,
  not the restored point.

**Why this is accepted:** the restored point is tested to be a counterfactual.
- The model re-scores it on the other side of the cutoff (assertion 3).
- It lies within `train_ranges` (assertion 7).
- It holds the fixed features exactly (assertion 5).

So it is a valid counterfactual, though not necessarily the nearest valid one.

**The impact is quantified per restored run in `restore_impact`:**
- how far each fixed feature moved back (search value against observed);
- √d²_M from the target to the search point and to the restored point;
- `B_pred` at the search point and at the restored point.

`review.html` prints this comment and the numbers under Table 1 and beside figure 7(b). Sources:
`documentation/fixed-constraint-restore.md` and `bl_find_local_cf()`'s "Fixed constraints"
roxygen section.

### Figures (plot the package's own objects only)

| file | panels | call |
|---|---|---|
| `fig-03.R` | 03 | `plot(bl_results, new_title = ...)`: four outcome groups, the decision surface, the boundary contours |
| `fig-04.R` | 04a, 04b | (a) `plot_biplotEZ(bl_results)`, training points only (**no arrows**: see issue 1); (b) `plot(bl_bnd_train, type = "boxplot")` |
| `fig-05.R` | 05a, 05b | (a) `plot_biplotEZ(bl_results, target_point = <TP x_obs>, target_label = "TP")`; (b) `plot(<TP run 1 bl_local>)` |
| `fig-06.R` | 06a, 06b | (a) base-R `dotchart()` of `sqrt(all_distances_mahalanobis)` over the 21 pairs, TP and TN, shortest marked; (b) `plot()` of TP run 5 beside run 1, **only if** run 1's `best_pair` is not (1,2), otherwise a one-line note in `review.html` |
| `fig-07.R` | 07a, 07b, 07c | `plot()` of TP run 1 (unconstrained), TP run 4 (mixed; shows the restored point's projection) and TP run 6 (`"decrease"`), each on its own plane (never overlaid) |
| `fig-08.R` | 08a, 08b | `print(plot(bl_shapley_values))` for TP and TN, run 1 |
| `fig-09.R` | 09a, 09b | `plot(bl_sparse)` for TP and TN, run 1 |

**Assumption for approval:** figures 8 and 9 and the sparse column of Table 1 use run 1
(unconstrained). Run 4's Shapley values and sparse counterfactual are still in `numbers.rds`.

## Assertions (`tests.R`)

Each assertion is recorded by a `check(name, expr)` helper, so every one gets a line in
`review.html`. The script stops with a non-zero exit at the end if any failed. (The prompt says
"stops on the first failure"; recording all of them first is needed for the review page.)

1. 392 rows, p = 7, all four outcome groups non-empty in `train_data`.
2. The GAM wrapper's predictions match `predict.gam()` floored to 3 dp.
3. For every run with a solution: `B_pred` is on the other side of 0.5 from the target's
   `pred_prob`.
4. Shapley efficiency: `sum(shapley_cause) == B_pred - pred_prob` within `1e-8`.
5. Constraints: each `"fixed"` feature is within ±0.5 of observed (exactly observed when
   `fixed_check$applied`); each `"increase"` feature is ≥ observed; each `"decrease"` feature
   (run 6) is ≤ observed.
6. In each sparse counterfactual, every `"fixed"` feature equals its observed value.
7. Every local `B_x` lies within `train_ranges`. For the sparse counterfactual, every changed feature
   equals `B_x`'s value and every unchanged one equals observed. A test-set target may itself lie
   outside `train_ranges`, so the sparse point is not range-checked as a whole.
8. **Reproducibility:**
   - (a) repeating `run_local()` for every run in the same session gives `identical()` numbers;
   - (b) a `fresh_result()` gives identical numbers for TP and TN, run 1;
   - (c) `run_all.R` run twice through `system2(file.path(R.home("bin"), "Rscript"), ...)`, each in
     a fresh process with its own `PAPER1_OUT` (set with `Sys.setenv()` before the call, so the child
     inherits it), gives `identical()` `numbers.rds` apart from `meta`.

## Issues table, pre-seeded

| # | symptom | where | proposed fix |
|---|---|---|---|
| 1 | Figure 4(a) cannot show boundary arrows non-interactively; arrows exist only in `bl_pick_point()` | `R/pick_point.R` | a `type = "biplot"` option in `plot.bl_boundary()`, or a non-interactive `rows =` mode in `bl_pick_point()`; never `boundary =` in `plot_biplotEZ()` (CLAUDE.md §4) |
| 2 | `bl_find_boundary()` returns no Mahalanobis distance, so the target rule computes it in the script | `R/boundary.R` | a `dist_mahalanobis` field (display-neutral; does not change any search) |

Earlier planning issues are already fixed and are not repeated: `used_in_sparse`, the sparse cross
and `rotate_deg`, the grid cutoff, empty CVA groups, the `seed` docs and the stale implementation
summary. The two items you parked on 2026-10-29 go straight to Task D.

**Watch for in the first run** (log as issues if they appear): non-integer counterfactuals for
Pregnancies and Age (the `round_to` question, CLAUDE.md §9); the target-range warning for a test-set
target; `solution_found = FALSE` in a constrained run (most likely the TN target's run 6); a restored
point whose √d²_M is much larger than the search point's.

## Error cases

| case | action |
|---|---|
| the zero filter does not leave 392 rows | stop; do not change the set-up |
| an outcome group is empty | stop with a message; do not change the set-up |
| `s(Pregnancies)` fails | lower its `k` and flag it as a set-up change for approval (not expected: 17 distinct values) |
| the GAM path's predictions mismatch | switch to `"custom"` + `predict_fn` and record it |
| a constrained run returns `solution_found = FALSE` | report it in Table 1 with `blocking_constraint`; skip its Shapley and sparse steps; do not re-pick the target |
| a group has no non-`NA` global counterfactual | stop with a message (the target rule cannot be applied) |
| check 8 (a), (b) or (c) fails | record it as an issue; never re-run to hide it |
| no run applies a fixed restore | say so in the restore comment; no `restore_impact` numbers to report |
| `plot.bl_boundary()` draws on its own and returns invisibly | call it inside the open device and do not `print()` it again |
| a figure needs something the package cannot draw | do not hand-roll it; log it as an issue |
| Quarto fails to render | check `rmarkdown` and report; the review falls back to the PNGs and `table1.csv` |

## Verification (running it)

- R from the PowerShell tool: `& "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" run_all.R`, then
  `tests.R`, with the working directory `scripts/paper1-pima/` (both in the background; `tests.R`
  runs the pipeline twice more).
- Render: `& "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe" render review.qmd`.
- Nothing is written outside `scripts/paper1-pima/` (no `Rplots.pdf`: every figure opens its own
  device; `run_all.R` opens `pdf(NULL)` first as a guard).

## Gate 2 (stop)

Report the path of `review.html`, the assertion results, Table 1 in brief and the issues table. For
each issue you mark "fix now": a short plan, the fix, a test, a NEWS entry, the walkthrough and
harness updates, `devtools::test()`, then re-run Task C. Each "defer" goes to Task D.

## After implementation

- **Before any script is written:** copy this revised plan over the plan-mode working copy, which
  still holds the 2026-10-31 version, so the two match.
- **The Claude Doc copy** (https://claude.ai/artifact/A3XvG3gydFu8DUm2umqaLc) is half-filled and
  now out of date. I leave it alone unless you ask me to finish or delete it.
- Update the plan-mode working copy with any changes made during implementation, and add a pointer from
  the Task C section of `paper1-working-plan.md` to it.
- No package code changes, so `devtools::test()` is not re-run before Gate 2 (348 / 0 / 4 stands).
- Update `progress.md` at session close.

## Implementation log (2026-11-02, first run, Gate 2 reached)

- **Built as planned.** All the files above are in `scripts/paper1-pima/`. The GAM wrapper path
  is `"GAM"` (max |diff| against `predict.gam()` = 0).
- **Changes from the plan:**
  - `run_local()` passes `max_pairs` only when it is set, so runs 1-4 and 6 use the package default
    (21). The default is recorded as `pairs_tried` = 21.
  - The wrapper check runs in `02_global.R` through the exported `bl_predict()`, not in
    `01_data-model.R`, which has no `bl_result` yet.
  - `bl_find_boundary()`'s global importance comes from `plot(bl_bnd_train)`, which prints the
    ggplot itself. `fig-04.R` re-prints the stored `importance$plot`.
  - `00_setup.R` opens a `pdf(NULL)` device in non-interactive runs, so no `Rplots.pdf` is ever
    written.
  - Figure 6(a): the pair labels were drawn in the TP colour; the dot chart now draws its frame
    and labels first, then both series.
  - The restore comment on the review page also covers a restore that was not applied (the
    TP target's runs 3 and 4).
- **Targets:** TP = test row 68 (19 in the group), TN = test row 73 (41 in the group).
- **Results:** all 12 runs found a counterfactual. `tests.R`: 10 of 10 checks pass, including
  8a-8c (identical numbers on repeat, after a fresh rebuild, and across two fresh processes).
- **Issues for Gate 2:** 10, in `scripts/paper1-pima/issues.csv`, shown on `review.html`.

## Revision (2026-11-04, after your review of the Gate 2 re-run)

Plan: `.claude/reference/known-issues-taskD-plan.md`. Scripts only; no package changes.

- **Pruning.** After the 7-feature global search, DiabetesPedigreeFunction, SkinThickness and BMI
  (the three smallest per-variable totals: 48.68, 44.11, 40.05) are dropped. The new
  `02b_prune-refit.R` checks that they are the three weakest, re-splits the same 392 rows with
  the same seed (same 313 / 79 rows, checked) and refits the GAM on Pregnancies, Glucose, Insulin
  and Age. **Phase 3 runs on this 4-feature model**: targets, all runs, Table 1, figures 5-9.
- **Run 6** is `Glucose = "decrease", Insulin = "decrease"` (BMI is pruned).
- **Figures:** 3 = the 7-feature biplot; 4(a) = the importance boxplot (was 4(b)), with a
  "global search" subtitle; 4(b) = figure 3 redrawn on the 4-feature model. The old 4(a)
  (training points without arrows) and every mention of the arrow issue are removed.
- **Ticks:** every biplot passes `ticks_v = 4L` (`biplot_ticks`); figure 7(c) raises Age,
  Glucose and Insulin with `ticks_var`/`ticks_n`. Every axis shows at least 2 ticks.
- **Shapley subtitles** give the target's prediction and the counterfactual's.
- **Figure 7** takes the target whose runs 4 and 6 both bind (TN, else TP, else TN's drawable
  panels). On the 4-feature model neither qualifies, so it shows TN runs 1 and 6 (issue 4).
- **`review.qmd`** follows the workflow: 7-feature data, global search, figures 3 and 4(a); the
  pruning and figure 4(b); then targets, local runs, Table 1 (four feature rows) and figures 5-9;
  assertions; open issues only.
- **`issues.csv`** has a `status` column; issue 1 is removed; open: 4, 5, 9, 11, 12, 13.
- **`run_all.R`** removes the run's old `fig-*` files first, so an undrawn panel cannot show a
  stale image. **`tests.R`** gains check 1b (pruning) and checks the wrapper for both models.
- **Results (4-feature model):** TP = test row 23 (pred 0.963), TN = test row 25 (pred 0.374).
  TP run 3 and TN runs 3 and 4 have no counterfactual.

## Revision (2026-11-04, your decisions on `pruned-model-issues-review.md`)

Plan: `.claude/reference/known-issues-taskD-plan.md` (archived as
`.claude/reference/known-issues-taskD-plan.md`).

- **Issue 4 (a):** figure 7 is the TP target's run 4 (its run 1 is figure 5(b)), the TN target's
  run 1 and the TN target's run 6. Figure 7(a) shows the same point as figure 6(b): TP runs 4 and 5
  give the same counterfactual on pair (1,2).
- **Issue 5 (b):** run 2 is out of Table 1 (it equals run 1 for both targets); it stays in
  `numbers.rds` and the local-runs table.
- **Issue 9 (b):** `02b_prune-refit.R` records the 4-feature model's Mahalanobis check (-35.2%,
  against -27.4% for 7 features); the review page shows both.
- **Issues 11, 12, 13 (a):** paper text (Task E). The package options of 12 and 13 are recorded in
  `documentation/known-issues.md` (Task D).
- `tests.R`: 13 of 13 PASS.

## Revision (2026-11-04, the duplicate figure and run 5)

Plan: `.claude/reference/known-issues-taskD-plan.md` (this version).

- **Figure 6(b) is removed.** Figure 6 is the per-pair dot chart alone (`fig-06`), whose (1,2) row is
  run 5's search distance. Figure 7(a)'s caption says it is also run 5's counterfactual: the TP
  target's run 4 picks pair (1,2), and setting Pregnancies 7.42 back to 7 gives the point that run 5
  reaches by rounding it to 7.
- **Table 1 has a run 5 column** ("pair (1,2) only"), with a sentence: run 5 is run 1 without the
  pair search, and shows what the search gains (√d²_M TP 1.50 against 1.92, TN 0.14 against 0.17).
  Run 5 had been specified only to draw figure 6(b), so it was never a Table 1 column.
- `issues.csv`: issue 14, closed. `documentation/known-issues.md` §1.2 no longer cites figure 6(b).
