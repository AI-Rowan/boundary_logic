> **IMPLEMENTED 2026-10-31** on branch `paper_nr1_changes`. Archived copy of the approved plan.

# Plan: Step 16 raw units, `max_pairs = 21`, the Phase 2 Mahalanobis check, 2-dp tables, PCA note

Date created: 2026-10-31
Date last edited: 2026-10-31

> **IMPLEMENTED 2026-10-31.** Uncommitted; see Implementation notes at the end.

## Context

The session-4 notes are approved. The user asked for five follow-ups (2026-10-31):

1. **The Step 16 external applicant is in raw units.** `development/03_loan_status_Boundary_Logic.R`
   builds `new_applicant` in raw loan units (credit score 420, loan 12000) and passes it to
   `bl_results_v2`, whose model works in z-scores.
   - `bl_set_scaling()` is display-only, and `bl_select_target()` takes a data-frame target as model
     units. So the search explains a point hundreds of SDs out (Z distance about 115.76).
   - The vignette is not affected: it does not standardise, so its raw applicant is correct.
2. **The PCA default departs from the thesis** (rotation note §4.4). Record the decision: Paper 1
   compares by Mahalanobis distance only (CVA, where both selectors are identical, note §3.4). The
   PCA default is addressed only once the rotation problem (note §4.1, §5) is fully addressed.
3. **Raise `max_pairs`.** The user chose **21** after being told that 15 never reaches component 7:
   - 21 covers every pair for p ≤ 7, so Paper 1's Pima set-up (p = 7) searches all 21 pairs, as the
     thesis does, including its winner (1,7);
   - from p = 8 the cap still leaves pairs out.
4. **The Phase 2 diagnostic does not describe the plot** (Phase 2 note §4.2).
   - **The plot** standardises each feature by √W_jj.
   - **The diagnostic's "diagonal"** uses (W⁻¹)_jj, which is ≥ 1/W_jj.
   - **So** on the harness it reports a −42.6% "cross-correlation" that partly measures this
     mismatch, not the plot.
5. **Scientific notation in the comparison tables.** The observed / counterfactual / change tables
   (`.cf_compare_table()`) are numeric data frames rounded to 4 dp. `print.data.frame()` switches to
   scientific notation when one column mixes large and small values.
   - **The user's choices:**
     - 2 dp in fixed notation;
     - a row is widened only where a change would otherwise print as 0.00;
     - the same format for the Target/Counterfactual blocks, the Shapley column and the sparse-path
       table.

**Order of work.** Tasks 1–5 run in sequence. Each one updates its code, tests, harness and docs
before the next starts (CLAUDE.md §2). Walkthrough numbers come only from real harness runs.

## Files to read first

- `R/local_cf.R`:
  - `.generate_pairs()` (12–32);
  - `bl_select_target()` (116–226);
  - `bl_find_local_cf()` roxygen and signature (370–444);
  - `.cf_compare_table()` (914–948);
  - `print.bl_local_result()` (≈1320–1383);
  - `print.bl_target()` (1386–1405);
  - `print.bl_counterfactual()` (1408–1430).
- `R/shapley.R`: `print.bl_shapley()` (339–363) and `print.bl_sparse_result()` (608–643).
- `R/boundary_plot.R`: `plot.bl_boundary()` (94–246).
- `R/data_prepare.R`: `bl_set_scaling()` (301–321) and `.scale_to_raw()` (381–394).
- `R/feasibility_utils.R`: `get_filter_logical_vector()`.
- `tests/testthat/`:
  - `test-local_cf.R` (133–155: the table tests);
  - `test-set_scaling.R` (153–249: `build_shap_pipeline()`; 231 asserts a 4-dp string).
- `documentation/boundary-distance-and-robustness.md` §4.2.
- `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md` §3.2, §3.6, §4.4, §4.7 and
  §6.2.
- The walkthroughs:
  - `development/reference/review_section7_to_8.md` (378–385, return value);
  - `development/reference/review_section9_to_15.md` (131–188, 266, 396, 484–493, 738–747, 833–841,
    911–919, 961).

## Task 1: record the PCA-default decision (documentation only)

- **Rotation note §4.4:** add a "Decision (2026-10-31)" paragraph.
  - Paper 1 compares counterfactuals by Mahalanobis distance only, under CVA, where it equals the
    in-plane distance (§3.4).
  - The PCA default (`W_binary`) stays as it is until the rotation problem (§4.1) is fully
    addressed (§5); §5.4's choice of metric is taken then.
  - Add the same ordering to §6.2.
- **CLAUDE.md:**
  - **§6:** extend the PCA sentence of the CVA-identity bullet with "deferred until the rotation is
    fully addressed; Paper 1 uses Mahalanobis comparisons only".
  - **§9:** add a one-line "PCA pair-selector default" bullet.
- **`progress.md`:** move it from "decide" to "deferred".

## Task 2: `max_pairs` default 10 → 21

**Code** (`R/local_cf.R`):
- `bl_find_local_cf(..., max_pairs = 21L, ...)`.
- Roxygen `@param max_pairs`: "Default `21L`, every pair for up to 7 features (choose(7, 2) = 21).
  With more features, the pairs with the highest components are left out."
- `devtools::document()` regenerates `man/bl_find_local_cf.Rd`.

**Tests** (`test-local_cf.R`):
- the default is `21L`;
- `.generate_pairs(7L, 21L)` has 21 rows, row 16 is (1,7) and the last row is (6,7);
- `.generate_pairs(4L, 21L)` has 6 rows.

**Docs:**
- `1 Foundation intro documents.txt` §10.7 ("default 10");
- `2 implementation_summary.txt` (352, 367: "a maximum of 10 combinations");
- rotation note §3.2 (cap), §3.6 table row, §4.7 and §6.2 item 7:
  - **§4.7:** the default now covers p ≤ 7, including Pima's (1,7);
  - **§6.2 item 7:** a default tied to choose(p,2) remains future work for p ≥ 8;
- walkthrough `review_section9_to_15.md` 266 (Parameters row) and 961;
- the `17_1` comment at 282;
- `NEWS.md` (a behaviour change: results change for p ≥ 6, and run time grows about linearly with
  the number of pairs);
- `progress.md`.

**Local plan files** (gitignored):
- `boundary-logic-fix-prompt.md`:
  - 235: "`max_pairs = 21L` (the default)";
  - 280: "each of the 21 pairs";
  - a one-line note at Task E item 3 that the default changed on 2026-10-31.
- grep `commit-and-push-the-mutable-hearth.md` for "10 pairs" or "ten pairs".

**Effect:**
- the loan harness (p = 4, 6 pairs), iris tests and `17_*` scripts are unchanged;
- the README Pima example (p = 6) now searches 15 pairs instead of 10.

## Task 3: make the Phase 2 Mahalanobis check describe the plot

**Design.**
- **The plotted measure** is v_j/√W_jj per feature. Its squared length,
  Σ_j v_j²/W_jj = `rowSums(vec_to_boundary_sd^2)`, is the diagonal approximation the plot implies.
- **The check compares that with the full** v′W⁻¹v. The difference is exactly what the per-feature
  plot leaves out by ignoring the correlations in W.
- **The sign:**
  - negative: correlated features are double-counted, so the plot overstates;
  - positive: the plot understates.
- **Units:** the check stays in squared units, because d²_M is a quadratic form. The plot's
  per-variable totals T_j are L1 sums of the same standardised components; the docs say so.

**Code** (`R/boundary_plot.R`, the block at 215–237):
- `diag_per_row <- rowSums(vec_to_boundary_sd^2)` replaces the `diag(metric_inv)` line.
- **The print:**
  ```
  Mahalanobis check (squared distances, summed across obs):
    plotted (per feature) = %.2f | full d_M^2 = %.2f | difference = %.2f (%+.1f%% of full)
  ```
- The 25% note is kept, reworded: "the per-feature plot leaves out …".
- **The return value** gains `mahalanobis_check = list(plotted, full, difference, pct)`, or `NULL`
  for `"euclidean"` or when `metric_inv` is missing. Add it to `@return`. No field is removed.
- Update the roxygen "Distance computation" wording ("cross-correlation diagnostic" → "Mahalanobis
  check").

**Tests** (new `tests/testthat/test-boundary_plot.R`; iris GLM, CVA, `bl_find_boundary()`; plot under
`pdf(NULL)` with `capture.output()`):
- `mahalanobis_check$plotted == sum(res$vec_to_boundary_sd^2)`: the check describes the plot;
- `full` equals a direct `v W^-1 v'` sum;
- with `metric <- diag(diag(W))` and `metric_inv <- diag(1/diag(W))` substituted on a copy, the
  difference is about 0;
- `distance = "euclidean"` gives `NULL`.

**Docs:**
- **Phase 2 note §4.2:** the open issue is resolved; describe the check, its sign and units, and
  replace the illustrative block with the harness output. Update the §6 code line numbers.
- `2 implementation_summary.txt` §4.6 (314);
- `review_section7_to_8.md` 378–385, its return-value table, and Step 7b if it is quoted;
- **CLAUDE.md:**
  - §2: the Phase 2 note entry ("its open mismatch with the plotted measure" → resolved);
  - §6: the `plot.bl_boundary()` bullet;
- `NEWS.md`.
- Harness run → quote the new block.

## Task 4: fixed 2-dp tables, widened only where a change would hide

**Helpers** (`R/local_cf.R`, beside `.cf_compare_table()`, `@noRd`):
- **`.fmt_digits(x, digits = 2L, max_digits = 4L)`** gives the decimals per element.
  - It is 2, unless `round(x, 2) == 0` and `round(x, 4) != 0`; then it is
    `min(4, 1 - floor(log10(abs(x))))`.
  - Examples: −0.0016 → 4; 0.004 → 4; 1.456 → 2; 1e-9 → 2 (prints 0.00); NA → 2.
- **`.fmt_fixed(x, d)`** gives fixed-notation strings: `sprintf("%.*f", d, round(x, d))`.
  - Values that round to 0 are set to +0 first, so it never prints "-0.00".
  - NA becomes "NA".
  - It never uses scientific notation.

**Comparison table:**
- **`.cf_compare_table()`** drops `round(tab, 4L)` and also returns `digits`, one per feature row:
  the maximum of `.fmt_digits()` over that row's differences (each level column − observed, and
  `change`), in display units.
  - In the harness, `loan_percent_income` gets 4 (0.0500 / 0.0484 / −0.0016); the other rows get 2.
- **New `.format_cf_table(cmp)`** formats observed, the level columns and change with
  `.fmt_fixed(col, cmp$digits)`. `constraint` is unchanged.
- **The three prints** call it before adding their extra columns:
  - **`print.bl_shapley()`:** shapley = `.fmt_fixed(s, .fmt_digits(s))`; "No effect" values print
    as 0.00.
  - **`print.bl_sparse_result()`:** the path table's `shapley` follows the same rule. `pred` uses
    `.fmt_fixed(pred, .fmt_digits(pred - cutoff))`, so a prediction within 0.005 of the cutoff
    shows the side it is on (0.4996, not 0.50).

**Target/CF blocks:**
- `print.bl_target()` and `print.bl_counterfactual()` gain `digits = NULL`.
  - By default: `.fmt_digits(values)`.
  - The values print as a named character vector with `quote = FALSE`.
- `print.bl_local_result()` builds `cmp` first, then passes `cmp$digits` to both blocks, so the
  three blocks agree. Roxygen `@param digits`.
- The header lines (`%.4f` predictions and distances) and the per-pair distance table are unchanged:
  they were not in scope, and they never print in scientific notation.

**Tests:**
- **`test-local_cf.R`:**
  - `.fmt_digits()` on the examples above;
  - `.fmt_fixed(-0.001, 2)` is "0.00";
  - `.fmt_fixed(36000000.123, 2)` is "36000000.12";
  - `.cf_compare_table()`'s `digits` is 4 for a −0.0016 row and 2 otherwise;
  - `.format_cf_table()` on mixed magnitudes contains no `[0-9]e[+-][0-9]`.
- **`test-set_scaling.R:231`:** expect the 2-dp string, e.g.
  `formatC(round(raw_cf[1L], 2L), format = "f", digits = 2L)`.

**Docs:**
- harness run → refresh every quoted block in `review_section9_to_15.md`: print(tgt), print(bl_local),
  print(bl_shapley_values), print(bl_sparse), the Step 15 sparse table and the path tables;
- CLAUDE.md §5 (the raw-unit values bullet: one sentence on the format and the shared helpers);
- `2 implementation_summary.txt` (the print-method description);
- `NEWS.md`.

## Task 5: raw-unit external targets (`bl_apply_scaling()` and a target-range warning)

**`bl_apply_scaling(data, bl_result)`** (exported; `R/data_prepare.R`, after `.scale_to_raw()`):
- **What it does:** converts new raw-unit data to the model's units using the recorded transform:
  each feature column becomes `(raw - center) / scale`. Other columns (e.g. `class`) are untouched.
  It is the input-side counterpart of `.scale_to_raw(kind = "level")`.
- **What it accepts:** a `bl_result`, `bl_data` or `bl_filter_result`.
- **Errors** (`call. = FALSE`):
  - `data` is not a data frame (`stop_if_not_data_frame()`);
  - the wrong class for `bl_result`;
  - no scaling is recorded ("nothing to convert: the model already uses the data's own units");
  - feature columns are missing;
  - a feature column is non-numeric.
- Full roxygen, with an `@examples` round trip.

**Warning in `bl_select_target()`** (data-frame targets only; index targets come from the package's
own data). It comes from a new private helper, `.warn_target_range(x_obs, bl_result)`:
- **When it fires:** it lists the features outside `bl_result$train_ranges`, and does nothing when
  `train_ranges` is not a list or all values are in range.
- **The raw-units hint:** if a scaling is recorded and converting the target leaves *fewer*
  features out of range, it adds: "This looks like raw units: bl_set_scaling() is display-only,
  so convert the target first with bl_apply_scaling(target, bl_result)."
- **Why a warning, not an error:** an out-of-range applicant is legitimate. The search then finds
  counterfactuals only within the training range.

**Harness Step 16:** after the `intersect()`, add
`new_applicant <- bl_apply_scaling(new_applicant, bl_results_v2)`, with a two-line comment.
- Update the section header comment: "give the values in raw units".
- The applicant's values are not changed. If credit score 420 lies below the (hull-filtered)
  training range, the warning fires without the hint; the walkthrough then shows it.

**Tests** (`test-set_scaling.R`, `make_std_data()`):
- **`bl_apply_scaling()`:**
  - the round trip `bl_apply_scaling(raw, bl)` equals the standardised features;
  - `class` is untouched;
  - column order does not matter;
  - it works on `bl_data` and on `bl_result`;
  - each error case.
- **The warning:**
  - a raw-unit row on the scaled pipeline warns, matching "bl_apply_scaling";
  - the converted row gives no warning;
  - an index target gives no warning;
  - on an unscaled pipeline, an out-of-range data frame warns without "bl_apply_scaling".

**Docs:**
- `review_section9_to_15.md`:
  - Step 16 (158–165): the conversion call, plus a Parameters table for `bl_apply_scaling()`;
  - its console output from the harness run;
  - re-check line 396's ratio for the Step 16 search;
- `review_section4_to_6.md`: one line where `bl_set_scaling()` is described, pointing to
  `bl_apply_scaling()` for new raw-unit data;
- `2 implementation_summary.txt`: inventory row, data-prep section, `bl_select_target()` warning;
- **CLAUDE.md:**
  - §3 file map: `data_prepare.R` gains `bl_apply_scaling()`;
  - §5 raw-unit bullets: one line on input conversion and the warning;
- `NEWS.md`;
- `progress.md`: the Step 16 bug fixed.

## Error cases (summary)

- **`bl_apply_scaling()`:** a non-data-frame, a wrong object class, no scaling, missing or
  non-numeric features → errors.
- **`.warn_target_range()`:** no `train_ranges` → skip; NA values → skipped for that feature.
- **`.fmt_fixed()`:** NA → "NA"; −0 → "0.00".
- **`plot.bl_boundary()`:** a missing `metric_inv`, or `"euclidean"` → no check, and
  `mahalanobis_check = NULL`. The existing `metric` fallback warning is unchanged.

## Verification

1. **After each task:** `devtools::document()` then `devtools::test()` through a scratchpad `.R`
   driver (CLAUDE.md §5), expecting 0 FAIL (baseline 314 PASS / 4 WARN, plus the new tests).
2. **Harness end to end** through a scratchpad driver after Tasks 3, 4 and 5, taking the walkthrough
   numbers from those logs. Check that:
   - Steps 11 and 15 keep their pairs and distances (`max_pairs` does not touch p = 4);
   - no comparison table, target block or path table contains `e+` or `e-` (grep the log);
   - Step 16's in-plane distance is no longer about 115.76;
   - the driver prints the Step 16 ratio `all_distances_mahalanobis / all_distances^2` (expect 15992).
3. **`17_3`** end to end (the default `max_pairs` changed): all assertions still pass.
4. **The README Pima example** (p = 6, now 15 pairs) runs.
5. **`devtools::check(document = FALSE, args = "--no-manual", error_on = "never")`**, expecting
   0 / 0 / 0. It covers the new export and the print-method arguments.
6. **No non-ASCII characters** added to `R/` (grep the changed files).
7. **Stop and report.** No commit or push; session 3's Parts B and C and session 4 are still
   uncommitted.

## Implementation notes (2026-10-31)

**Done as planned:** Tasks 1-5 and verification steps 1-7.

**Additions and deviations:**
- `.print_values()` (private) prints the Target/Counterfactual blocks.
- The walkthrough's Step 16 passage gains:
  - the converted applicant;
  - the warning, quoted from the log;
  - the Step 16 result.
- One stale sentence in `review_section9_to_15.md` ("the reported distance can differ between
  repeated calls") was replaced: the search has been deterministic since 2026-10-30.
- The local plan file `commit-and-push-the-mutable-hearth.md` also had two `max_pairs = 10L`
  mentions; both are updated.

**Measured (harness):**
- **The Mahalanobis check:**
  - Step 7: plotted 9675.25, full 11754.70, difference +17.7% (the old diagnostic: -42.6%);
  - Step 7b: +15.9%.
- **Step 16 after conversion:**
  - Z 0.0135 (was 115.76);
  - predicted 0.296, class 0;
  - pair (1,2) wins, at sqrt(d2_M) 1.7026;
  - the counterfactual's credit score is 489.44, at the training minimum (489).
- **The applicant warns:** credit score 420 is below the training range (489-754).
- **Steps 11 and 15** keep pair (3,4). The ratio is 15992 in all three searches.

**Results:**
- **Tests:** 348 PASS / 0 FAIL / 4 WARN / 0 SKIP.
- **`17_3`:** every assertion passes.
- **README Pima:** runs, 15 pairs.
- **`R CMD check`:** 0 / 0 / 0.
- **Non-ASCII:** none added to `R/`.
