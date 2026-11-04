# Paper 1 Gate 2 fixes: the package and the Pima example

Date created: 2026-11-04
Date last edited: 2026-11-04

> **IMPLEMENTED 2026-11-04** on `paper_nr1_changes` (uncommitted at archive time). The
> implementation log at the end records the results and the changes from the plan.


Branch: `paper_nr1_changes`. The authoritative plan, approved 2026-11-04. It is copied from the
plan-mode working copy `.claude/plans/graceful-watching-shannon.md`.

## Context

Task C ran the Paper 1 Pima example and stopped at Gate 2 with 10 issues. You marked your decisions
in `scripts/paper1-pima/gate2-issues-review.md`. `.claude/reference/gate2-code-fixes-plan.md` is
unchanged since 2026-10-30, so this plan works from the review file and from your answers of
2026-11-04. The aim is to fix what the paper needs and re-run Task C, so that the numbers in the
paper come from the fixed package.

| # | Your decision | Work |
|---|---|---|
| 1 | Defer the arrow layer: the CLAUDE.md §4 rule stands | Figure 4(a) stays without arrows, and figure 4(b) carries the distances. Add a `bl_pick_point()` example for the Phase 3 test point to the harness and `review_section9_to_15.md`. The arrow layer goes to Task D. |
| 2 | Fix now, (a) | `dist_mahalanobis` field in `bl_find_boundary()` |
| 3 | Fix now, (b) only | Strict `"fixed"`: a failed restore means no solution, with a message. Update `fixed-constraint-restore.md`, and record (c) and (d) as future developments. |
| 4 | Fix now, (a) | Figure 7(c) becomes the TN target's run 6. Your answer of 2026-11-04 makes figure 7 the TN target's runs 1, 4 and 6. |
| 5 | Paper text, (a) | none |
| 6 | Fix now, (c), **via `set_filters()`** | A new `"integer"` constraint, e.g. `credit_score = c("increase", "integer")`. It rounds the local counterfactual, so the Shapley values and the sparse counterfactual inherit it. It replaces `round_to`. |
| 7 | Fix now, (a) | `"observed -> counterfactual"` everywhere, with the decimals of the comparison tables. The display rules are written up in `documentation/scaling_axis_relabel_note.md`. |
| 8 | Fix now, (a) | Tune the label offsets in the Pima figure scripts |
| 9 | Paper text, (a) | none |
| 10 | Paper text, (a) | none |

## Assumptions for you to confirm at approval

1. **The rounding rule for `"integer"`** (your answer pointed to `set_filters()`, so the rule itself
   is mine):
   - For each integer feature the search moved off a whole number, try the whole number just below
     and just above, in raw units. Every combination of these is a candidate.
   - Keep the candidates that still cross the cutoff, meet every `set_filters()` constraint and lie
     within `train_ranges`.
   - Take the one nearest the observed point, by the search's own selector (Mahalanobis, or
     Euclidean in the rotated plane); on a tie, the first.
   - If none is valid, report no solution, as for a failed fixed restore.
2. **Restore first, then round.** A failed fixed restore ends the search before rounding.
3. **`dist_mahalanobis` in `bl_find_boundary()` is squared**, `v' W^-1 v`, as the field of the same
   name in `bl_find_local_cf()` is; `print()` shows the mean of its square root. It uses `W`, not a
   rotated `W`; the note explains why (Fix 1, Docs).
4. **A numeric range with `"integer"`** is given as `list(c(min, max), "integer")`, since
   `c(600, 700, "integer")` would coerce to character.
5. **`round_to` is kept as a deprecated argument that stops when used.** A non-NULL value stops with
   a message pointing to `set_filters(..., "integer")`, so `round_to = NULL` in old code still runs.
   `sparse-counterfactual.md` §3 shows, with the harness example, how `round_to` was used and what
   it gave, and the `"integer"` call that replaces it, with its result (Fix 3b, Docs).
6. **Pima:** only Pregnancies is `"integer"`, in every run where it is not `"fixed"`. Age and the
   other features stay continuous. Run 1 becomes "Pregnancies integer only".
7. **Table 1** keeps the TP run 6 and TN run 2 columns and adds a row, `binds`, showing whether the
   counterfactual differs from run 1's.
8. **The `bl_pick_point()` example** goes in harness Step 9, for `tdp = 13`, wrapped in
   `if (interactive())` so that sourcing the script skips it.

## Files to read first

- `R/local_cf.R`:
  - `.fixed_tol` (L73), `.apply_actionability()` (L84), `set_filters()` (L318);
  - the post-search block of `bl_find_local_cf()` (L744-835), `.restore_fixed()` (L892);
  - `.constraint_labels()` (L956), `.cf_compare_table()`, `.fmt_digits()`, `.fmt_fixed()`
    (L984-1057), `.fixed_notes()` (L1070);
  - `print.bl_local_result()` (L1448), `print.bl_filters()` (L1580).
- `R/shapley.R`:
  - the labels in `bl_shapley()` (L234) and in `plot.bl_shapley()` (L276-321);
  - the fixed detection and `round_to` in `bl_find_sparse_cf()` (L416-495).
- `R/boundary.R`: the distances in `bl_find_boundary()` (L443-449), the return value (L452), and
  `print.bl_boundary()` (L479).
- `R/data_prepare.R`: `.scale_to_raw()` (L381) and `bl_apply_scaling()` (L436).
- `tests/testthat/test-local_cf.R`: `rf_fixture()` (L78) and the restore tests (L108-250).
  `test-boundary.R`, `test-set_scaling.R` (L177-233).
- The notes: `documentation/fixed-constraint-restore.md`, `documentation/sparse-counterfactual.md`
  §3, `documentation/scaling_axis_relabel_note.md`.

## Fix 1: `dist_mahalanobis` in `bl_find_boundary()` (issue 2)

**Design.**
- After `dist_x` (`R/boundary.R` L449), add:
  `v <- B_x_mat - x_obs_mat; dist_mahalanobis <- rowSums((v %*% metric_inv) * v)`.
  It is the squared Mahalanobis distance in model units, the same as `bl_find_local_cf()`
  (L687-692).
- `metric_inv` is read from `bl_result`. If it is NULL, `dist_mahalanobis` is `NA` and a warning is
  given (CLAUDE.md §5 convention). `.empty_return()` gets an NA vector.
- Add the field to the return list after `dist_x`, and to the `@return` roxygen.
- `print.bl_boundary()` gains `Dist Mahal (mean sqrt): %.4f`.
- No search result changes.

**Tests (`test-boundary.R`):**
- the field equals a hand computation of `v' W^-1 v` on the iris fixture;
- an NA row stays NA;
- with `metric_inv` set to NULL, it is NA with a warning;
- the print shows the line.

**Docs:**
- `2 implementation_summary.txt` §4.5;
- the field table and a print sample from a real run in `review_section7_to_8.md`;
- `documentation/boundary-distance-and-robustness.md`:
  - the new field, in the Phase 2 distances;
  - a subsection, **why the distance uses `W` and not a rotated `W`**:
    - The distance is computed on X-space vectors, `v = B_x - x_obs`. `W` is the within-class
      scatter of the features, a property of the data, not of the plane.
    - A rotation `A` (orthogonal, `A A' = I`) acts only on the Z coordinates: the local search's
      plane `V A`, for each pair. It decides *which* boundary point is found, not how an X-space
      move is measured. `bl_find_boundary()` uses the global plane and has no rotation at all.
    - Measured in Z coordinates instead, the metric would be `A' V' W V A`. Under CVA, `V'WV` is a
      multiple of the identity, so `A` cancels: this is the identity proved in the rotation note
      §3.4, cited there and not repeated. Under PCA it is not, so a Z-space distance depends on the
      rotation. Computing in X-space avoids that.
- NEWS "New features".

## Fix 2: Shapley labels and decimals (issue 7)

**Design.**
- One private label builder, `.shapley_labels(var_names, observed, cf, scaling)`, in `R/shapley.R`:
  - it converts both values to display units with `.scale_to_raw(kind = "level")`;
  - it sets each row's decimals with `.fmt_digits(cf - observed)`, the comparison tables' rule;
  - it formats both numbers with `.fmt_fixed()`;
  - it returns, for example, `"Glucose: 179.00 -> 170.63"`.
- `bl_shapley()` stores `varnames_p` from it in model units. `plot.bl_shapley()` rebuilds it with
  `scaling` (identical when there is no scaling). Both paths then read "observed -> counterfactual".
- The y title is "Variable (observed -> counterfactual)", with ", raw units" when scaled.
- The x label becomes "Impact on prediction (Shapley value)".
- The bar labels use `.fmt_fixed(s, .fmt_digits(s))`.
- Widen the value axis (`ggplot2::expansion()`, with a single-line `@importFrom`) so that labels are
  not clipped.
- Grep `R/` for any other "->" label that shows a change.

**Tests:**
- the default-path label for one feature is exact, for example
  `paste0(v, ": ", .fmt_fixed(obs, d), " -> ", .fmt_fixed(cf, d))`;
- the raw-unit path label;
- a change below 0.005 widens to 4 dp;
- no scientific notation;
- the titles;
- update `test-set_scaling.R` L177-207 (the old "change required" expectation).

**Docs:**
- `scaling_axis_relabel_note.md` gets a "Display conventions" section and becomes the one place for
  these rules: raw-unit levels and deltas, the 2-dp rule widened up to 4, the Shapley labels and the
  bar values.
- These point to it instead of restating the rules:
  - CLAUDE.md §5 (and the §9 decimal-place item is marked resolved);
  - `2 implementation_summary.txt` (L409, L460-462);
  - Step 13 of `review_section9_to_15.md` (text and a sample from a real run);
  - the vignette prose (L519-520, L769);
  - `1 Foundation intro documents.txt` L264 (the label format only).
- NEWS "Bug fixes and other improvements".

## Fix 3a: strict `"fixed"` (issue 3, option b)

**Design** (in `bl_find_local_cf()`, after the `.restore_fixed()` call at L785; `.restore_fixed()`
itself is unchanged, so its unit tests stay valid):
- When `rf$fixed_check$applied` is FALSE:
  - set `solution_found <- FALSE` and `best_result <- NULL`, and keep `fixed_check`;
  - build `blocking` from a new `.fixed_blocking(fixed_check, x_obs, var_names, scaling, cutoff)`.
    It reuses `.fixed_notes()`' display-unit formatting and reads like this:
    > No counterfactual holds the fixed feature(s) Pregnancies, Age at their observed values. The
    > nearest boundary point moved them to 0.32, 21.61; with them restored the model predicts 0.561,
    > which does not cross the cutoff (0.5) [or: which lies outside the training ranges]. Changing
    > the constraints on the actionable features in set_filters() may make a counterfactual that
    > holds Pregnancies, Age possible.
  - `message()` it as the existing no-solution path does (L762).
- The return value uses `null_pair_fields`, and keeps `all_distances*`, `fixed_check` and
  `blocking_constraint`.
- `print.bl_local_result()` wraps the Reason line with `strwrap()`.
- The "kept" branches of `.fixed_notes()` and the verbose message (L795-801) become unreachable and
  are removed. The applied note stays.
- `.fixed_tol` stays 0.5: the search still needs it to find a vertex.

**Tests (`test-local_cf.R`):**
- A new `bl_find_local_cf()` fixture where the restore must fail: two-class iris and a custom
  `predict_fn` that depends on the fixed feature only.
  - It gives `solution_found` FALSE, `fixed_check$applied` FALSE, a blocking message matching the
    feature and "set_filters()", and `B_x` NULL;
  - `print()` shows Reason and no table;
  - `bl_shapley()` stops.
- Update `.fixed_notes` (L198): the kept-note expectations go.
- L250 must no longer skip silently. It asserts either a restored point or a no-solution result
  whose message names the fixed feature.

**Docs:**
- `fixed-constraint-restore.md` is rewritten:
  - §1 what changed;
  - §2 the tolerance (kept);
  - §3 why restore after the search (kept);
  - §4 the strict rule, why a "fixed" that moves is not a constraint, and the message;
  - §5 what the adjusted point is and is not (it lies off the plane and is accepted because it is
    tested);
  - §6 Shapley and the sparse counterfactual;
  - §7 future developments: (c) a follow-up 1-D search holding the fixed features exactly,
    (d) searching in the free features only, and retrying the restore on the next-best pair;
  - §8 reproduction.
- The roxygen for the `set_filters()` "fixed" item, `@section Fixed constraints`, and the
  `fixed_check` and `blocking_constraint` fields.
- CLAUDE.md §6, the fixed bullet.
- `2 implementation_summary.txt` (L407-409).
- `review_section9_to_15.md` L541: the not-restored case is now no solution.
- The prose in README L206, the vignette L453-455 and L533, and Foundation L224-226.
- The stale line references in `biplot_rotation_and_pair_selection.md` L624-634.
- NEWS "Breaking changes".
- The harness numbers do not change: its `loan_amnt` restore is applied.

## Fix 3b: the `"integer"` constraint, and `round_to` removed (issue 6)

**Constraint parser.** One private `.parse_constraint(v, nm)` returns
`list(type = "decrease"|"increase"|"fixed"|"range"|"none", range, integer)`. It accepts:
- one direction;
- `"integer"`;
- a character vector of one direction plus `"integer"`, in either order;
- `c(min, max)`;
- `list(c(min, max), "integer")`.

Anything else stops with a message listing the forms. Every reader of `set_filters$constraints`
switches to it, and none tests `v == "fixed"` on a vector any longer:
- the validation in `set_filters()`;
- `.apply_actionability()`, which ignores `integer`;
- `.restore_fixed()` (L896);
- `.constraint_labels()`, e.g. `"increase, integer"`;
- `print.bl_filters()`;
- the fixed detection in `bl_find_sparse_cf()` (L436-439).

Two small helpers sit on top of it: `.fixed_vars()` and `.integer_vars()`.

**Rounding**, a new `.round_integers(best_result, x_obs, var_names, set_filters, bl_result,
bl_target, sv, distance, metric_inv)`, called after a successful restore (or on the search point when
no restore was needed):
1. Find the integer features that are not fixed. Skip a feature whose value already equals observed,
   or is a whole number in raw units, within 1e-8.
2. Convert to raw with `.scale_to_raw()`, and take `floor` and `ceiling`. Convert back as
   `(raw - center) / scale`, the vector form of `bl_apply_scaling()`. With no scaling, model units
   are raw units.
3. Build the `expand.grid` of candidates; 2^k rows, and stop with an error if k > 12. Score them in
   one batch through the existing scoring call.
4. Keep the rows that pass three checks:
   - they cross (`pred class != bl_target$pred_class`, as `.restore_fixed()` does);
   - `.apply_actionability()` passes;
   - `get_filter_logical_vector(train_ranges)` passes.
5. Pick the nearest row with the search's selector:
   - `"mahalanobis"`: `v' W^-1 v`;
   - `"euclidean"`: the distance to `Z_target` after projecting with `Vr_rot`. `B_z` is projected
     as in `.restore_fixed()` (L933-942); factor that projection into one helper.
   - Ties go to the first row.
6. Return the new `best_result` (with `B_x`, `B_z` and `B_pred`) and an `integer_check`: `vars`,
   `applied`, `n_candidates`, `n_valid`, `search_point`.

   If no row is valid, the search reports no solution with this message:
   > No whole-number values of credit_score keep the counterfactual across the boundary within the
   > constraints. Changing the constraints on the other features in set_filters(), or removing
   > "integer", may make one possible.

`dist_*` and the per-pair tables stay the search's, as for the restore.

**Downstream:**
- `print.bl_local_result()` gains an `Integer rounded :` line and a note in display units (the
  boundary values, the whole values chosen and the prediction).
- The Shapley values and the sparse counterfactual use the rounded `B_x`. The sparse counterfactual
  takes features from `B_x` or `x_obs`, so it holds whole numbers for whole-number observed values.
- `bl_find_sparse_cf(round_to = )`: a non-NULL value stops (assumption 5).
- Remove `round_to = NULL` from README L124 and the vignette (L540, L564).

**Tests:**
- `set_filters()` accepts each form and rejects bad ones (two directions, an unknown string, a range
  of length 3);
- `.parse_constraint()`;
- `.apply_actionability()` with `c("increase", "integer")`;
- `.round_integers()` on an `rf_fixture()`-style fixture: the chosen candidate is the nearest valid
  one, a direction-violating candidate is dropped, no valid candidate gives NULL or no solution, and
  an already-whole value is untouched;
- `test-set_scaling.R`: rounding happens in raw units (raw `B_x` whole);
- `bl_find_sparse_cf(round_to = 1)` stops;
- `"fixed"` with `"integer"` is not rounded.

**Docs:**
- `fixed-constraint-restore.md`: a new section on integer rounding, the second post-search
  adjustment, and a title that covers both adjustments (the filename stays);
- `sparse-counterfactual.md` §3, rewritten as "`round_to`, replaced by the `"integer"` constraint":
  - how `round_to` was used, with the harness example as it stands: the call
    `bl_find_sparse_cf(bl_shapley_values, round_to = 1)` and its table (z-score rounding reverses
    `loan_percent_income`, `solution_valid = FALSE`);
  - what using it does now (it stops, with the message quoted);
  - the replacing call `set_filters(..., credit_score = c("increase", "integer"))` and its result
    from the Fix 3b harness run;
  - how the integer constraint answers each of the five issues (units, one value for all features,
    direction, constraints, validity);
  - a pointer to the restore note's integer section;
- the roxygen for the `set_filters()` "integer" item, `@section Integer features` and the
  `integer_check` field;
- CLAUDE.md §6 gains a bullet, and the `round_to` item in §9 is marked resolved;
- `2 implementation_summary.txt` (L443-479);
- Foundation L267-275 (`round_to = 0.5` "recommended for clinical data" becomes `"integer"`);
- NEWS: "Breaking changes" (`round_to`) and "New features" (`"integer"`).
- **Harness Step 10**, your filter:
  `flt <- set_filters(tgt, loan_amnt = "fixed", loan_int_rate = "increase", credit_score = c("increase", "integer"))`,
  with its comments.
- **Walkthrough Steps 10-14** re-quoted from a real run.

## Fix 4: the `bl_pick_point()` example (issue 1)

- **Harness Step 9**: after the target biplot is drawn, add
  `if (interactive()) bl_pick_point(bl_results_v2, data = bl_results_v2$test_data, bl_boundary = bl_bnd_v2)`.
  The comment says to click the target, row 13, and press Esc. The script has no
  `interactive()` guard yet, so this is the first.
- **`review_section9_to_15.md` Step 9**:
  - a subsection with the call and what a click prints and draws;
  - row 13's global counterfactual, printed by the harness driver from `bl_bnd_v2` (`B_x`,
    `B_pred`, `dist_z`, √`dist_mahalanobis`);
  - a Parameters table;
  - the `rotate_deg = 0` caveat.
- Task D records the deferred non-interactive arrow layer.

## The Pima scripts (`scripts/paper1-pima/`)

- **`00_setup.R`**:
  - `runs_spec` with Pregnancies `"integer"` (assumption 6):

    | run | filters |
    |---|---|
    | 1 | `Pregnancies = "integer"` |
    | 2 | `Pregnancies = c("increase", "integer"), Age = "increase"` |
    | 3 | `Pregnancies = "fixed", Age = "fixed"` |
    | 4 | `Pregnancies = "fixed", Age = "increase"` |
    | 5 | `Pregnancies = "integer"`, `max_pairs = 1L` |
    | 6 | `Glucose = "decrease", BMI = "decrease", Pregnancies = "integer"` |

  - `extract_numbers()` records `integer_check`, an `integer_impact` (search against rounded: √d²_M
    and `B_pred`), `blocking_constraint` and `binds`.
- **`03_local.R`**:
  - targets ranked by `sqrt(bl_bnd_test$dist_mahalanobis)`;
  - Table 1 gains the `binds` and `integer_rounded` rows;
  - TP runs 3 and 4 show "no solution" with the reason.
- **`fig-07.R`**: TN runs 1, 4 and 6.
- **Figure scripts**: `label_offset_var`, `label_offset_dist` and `label_dir` tuned per figure (03,
  04a, 05, 06b, 07, 09) by viewing the PNGs.
- **`tests.R`**:
  - check 2b: `dist_mahalanobis` equals the script's `mahal()^2`;
  - check 5: every solution holds each `"fixed"` feature exactly and meets each direction;
  - new check 9: Pregnancies is a whole number in every `B_x` and sparse counterfactual;
  - checks 8a-c are unchanged.
- **`review.qmd`**:
  - the restore comment is rewritten for strict fixed and integer rounding (both adjusted points lie
    off the plane and are accepted because they are tested);
  - figure 4(a)'s caption: "figure 4(b) carries the distances";
  - figure 7's captions.
- **`issues.csv`**: the decision column is filled in from the review file.
- **Paper text**, recorded for Task E: issue 1 (c), issue 3 (TP runs 3 and 4 have no solution, and
  why), issues 4 and 5 (non-binding constraints), issue 6, issue 9 (one sentence on the check) and
  issue 10.

## Order of work

Fix 1, then 2, 3a, 3b and 4, then the Pima scripts. Each fix follows the same steps:
1. code;
2. `devtools::document()`;
3. tests;
4. `devtools::test()`;
5. NEWS and docs;
6. a harness run (scratchpad driver, `pdf()` opened first);
7. the walkthrough re-quoted.

After that come the Pima `run_all.R` and `tests.R`, and the `review.html` render. Then
`devtools::check(document = FALSE, args = "--no-manual", error_on = "never")`.

## Error cases

| case | action |
|---|---|
| `metric_inv` NULL in `bl_find_boundary()` | `dist_mahalanobis` NA, with a warning |
| a malformed constraint (e.g. two directions, `c(1, 2, 3)`, an unknown word) | `set_filters()` stops, listing the valid forms, `call. = FALSE` |
| the fixed restore fails (does not cross, or leaves `train_ranges`) | no solution; the message names the features and the reason, and suggests changing the other constraints |
| no valid whole-number candidate | no solution, with the integer message |
| more than 12 integer features to round at once | stop, with a message |
| scaling missing for an integer feature | `.scale_to_raw()` warns; round in model units |
| an observed value that is not whole, for an unmoved integer feature | left as observed (not forced to move) |
| `round_to` non-NULL | stop, pointing to `set_filters(..., "integer")` |
| a Pima panel with no solution | `plot()` prints "nothing to plot"; log it, and do not re-pick the run |
| the harness walkthrough numbers change | re-quote from the log; never edit numbers by hand |

## Verification

- `devtools::test()` after each fix. The baseline is 348 PASS / 0 FAIL / 4 SKIP; report the new
  counts.
- **Harness:** a scratchpad driver (`setwd` to the root, `pdf()` before `source()`) after Fixes 1,
  2, 3a and 3b. Fix 3a must leave Steps 10-14 identical; Fix 3b changes them, and the walkthrough is
  re-quoted.
- **Pima:**
  1. `& "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" run_all.R`, then `tests.R`, in
     `scripts/paper1-pima/`;
  2. all checks PASS, including 8c;
  3. `review.html` re-rendered with the RStudio Quarto.
- `devtools::check()`: 0 errors, 0 warnings, 0 notes.
- `git status` against the session start: only the intended files have changed, and there is no
  `Rplots.pdf`.

## After implementation

- Copy the plan to `.claude/reference/paper1-gate2-fixes-plan.md` with an IMPLEMENTED banner, and
  add it to the CLAUDE.md §2 plan archive. (The other two parts of the triplet are the notes and the
  implementation summary above.)
- Update `progress.md` and stop at Gate 2 again, with the test counts, the new Table 1 and the
  checks. No commit or push without your instruction.
- Next: Task D (`documentation/known-issues.md`, including issue 1) and Task E.

## Implementation log (2026-11-04)

All fixes implemented in the planned order. Tests 434 PASS / 0 FAIL / 0 SKIP / 4 WARN (the four
existing 2-group CVA warnings; 348 PASS before). `devtools::check()`: 0 errors, 0 warnings, 0 notes.
Pima: `run_all.R` clean, `tests.R` 12 of 12 checks PASS (8a, 8b and 8c identical), `review.html`
re-rendered.

**Harness.** Fixes 1, 2 and 3a leave the harness log identical apart from timestamps (strict
`"fixed"` does not bite: the `loan_amnt` restore is applied). Fix 3b (`credit_score =
c("increase", "integer")`) rounds 634.2206 to 635: both 634 and 635 cross at 0.547, and the
Mahalanobis selector prefers 635 (√d²_M 0.5914 against 0.5949), because `W` correlates the credit
score with the features that move. The Shapley values and the sparse CF are unchanged. Row 13's
global counterfactual (Step 9c): √`dist_mahalanobis` 0.9106, `B_pred` 0.538.

**Pima.** Targets unchanged (TP row 68, TN row 73). The TP target's runs 3 and 4 have no
counterfactual (restored predictions 0.561 and 0.513). Pregnancies is whole in every counterfactual
(TP 0.11 → 1, which moves √d²_M from 0.43 to 0.62 because 0 does not cross; TN 2.28 → 2, 4.27 → 4).
Figure 7 is the TN target's runs 1, 4 and 6.

**Changes from the plan.**
- `print.bl_boundary()` prints `Dist M (mean)  : <x> (sqrt of dist_mahalanobis)`, aligned with the
  other lines, instead of `Dist Mahal (mean sqrt)`.
- The "kept" branches of the fixed note were not deleted but rewritten: the private `.fixed_text()`
  gives one text, used wrapped as the note when the restore is applied and unwrapped as
  `blocking_constraint` when it fails. It says "its"/"their" by the number of features.
- `.integer_text()` and `.integer_notes()` are the counterparts for the rounding; `integer_check`
  also records `pred_rounded`. The projection of a restored or rounded point is one helper,
  `.project_local()`.
- Pima Table 1: `sqrt_d2M_restored` became `sqrt_d2M_cf` (the reported counterfactual's distance,
  after the restore and the rounding), and a `no_solution` row was added. The new checks are 2b
  (`dist_mahalanobis`) and 9 (whole Pregnancies).
- Figure label offsets: the global plane staggers the left labels (Glucose 1.2, Pregnancies 2.4,
  Age 3.6) and moves SkinThickness inside the top border (-1.2); each local panel has its own.
- README and the vignette got the `"integer"` row and lost `round_to = NULL`; the vignette's own
  `set_filters()` example (`credit_score = "increase"`) is unchanged, and the pkgdown site is not
  rebuilt.
