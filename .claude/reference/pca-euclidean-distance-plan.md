# The thesis distance rule under PCA: Euclidean in X-space, chosen by method

Date created: 2026-11-05
Date last edited: 2026-11-05

> **IMPLEMENTED 2026-11-05 (session 8, continued).** Approved plan, verbatim below the implementation notes.
>
> **Implementation notes:**
> - **Code:**
>   - `.resolve_distance()` (`R/utils.R`);
>   - `.cf_sq_distance()` (`R/local_cf.R`);
>   - `bl_find_local_cf()`: `"auto"`, `dist_x`, `all_distances_x`, selection by `dist_x` under Euclidean; the print
>     and plot show them;
>   - `.round_integers()`: X-space Euclidean, with a new `X_sd` argument;
>   - `plot.bl_boundary()` / `bl_robustness()`: `"auto"`.
>
>   The pair loop keeps its existing Mahalanobis expression, so the CVA results are bit-identical.
> - **Tests:** 465 PASS / 0 FAIL / 4 WARN (was 434). New `test-distance_by_method.R`, plus an X_sd rounding test in
>   `test-local_cf.R`.
> - **Paper 1 (CVA):**
>   - `numbers.rds` is identical to before (meta excluded);
>   - `tests.R` 13/13;
>   - `04_method-change-effects.R` passes all 6 self-checks (new: 6, the PCA default chooses the smallest `dist_x`).
> - **The loan harness (CVA):** numbers identical to the `HEAD` baseline. The only output changes are the
>   "Euclidean (X)" line and the `dist_x` column; the walkthrough samples are refreshed from that log.
> - **Rotation scripts `17_1`-`17_3`:**
>   - pinned to `distance = "mahalanobis"`;
>   - their console output is identical to the `HEAD` baseline, apart from the new print line, `17_3`'s timings and
>     its §11 identity line;
>   - `17_3` §11 now compares the Mahalanobis choice with the smallest `all_distances`, instead of an explicit
>     `"euclidean"` run.
> - **The PCA measurements on Pima** (4-feature, run 1):
>   - TP: (1,3) under the thesis rule, against no whole-number counterfactual under the old Mahalanobis default;
>   - TN: (1,4) against (3,4);
>   - the rules choose different pairs for 57 of 79 test rows;
>   - `dist_x` - `dist_z` < 1e-15.
> - **README:** its PCA example quotes no output, so it is unchanged.


(The previous plan in this file, Task E, is archived as `.claude/reference/methodology-changes-taskE-plan.md`.)

## Context

The thesis compares the candidate counterfactuals of the different component pairs by "Euclidean
distance for PCA biplots and Mahalanobis distance for CVA biplots" (§3.9.1, p.86), with the change
vector $\mathbf{d} = (\mathbf{z}_{cf} - \mathbf{z}_0)\mathbf{V}_2'$ measured in the full feature space.
Its Phase 2 variable importance divides by the overall SD (§4.2.2, p.103).

The package uses Mahalanobis distance by default under both PCA and CVA, in two places:
- **Phase 3:** the pair comparison in `bl_find_local_cf()`. Under PCA the metric is `W_binary`.
- **Phase 2:** the importance standardisation in `plot.bl_boundary()` and `bl_robustness()`, which
  divides by the within-class SD.

Under PCA this departs from the thesis. It is recorded as deferred in CLAUDE.md §9 and the rotation
note §4.4, and as a departure in the methodology note §3.5. The user has now asked for the thesis rule
under PCA, and for every note to be updated.

**Decisions (user, 2026-11-05):**
1. **Euclidean distance in X-space is standardised by the total SD:**
   $d_x = \lVert(\mathbf{B}_x - \mathbf{x}_{obs}) / \mathbf{X}_{sd}\rVert$. This is how
   `bl_find_boundary()` already computes `dist_x` (`R/boundary.R:454-456`).
   - For the default (standardised) PCA it equals the thesis's $\mathbf{d}$, and it equals the plane
     distance at the search point, because the rotated loadings are orthonormal.
   - It still measures correctly once `"fixed"` or `"integer"` moves the point off the plane.
2. **The new default is `distance = "auto"`:** Euclidean for PCA, Mahalanobis for CVA.
   `"mahalanobis"` and `"euclidean"` stay available when passed explicitly.
3. **It applies to both phases:** the Phase 3 pair comparison (and the integer-rounding choice, which
   uses the same distance), and the Phase 2 importance (`plot.bl_boundary()`, `bl_robustness()`).
   Under PCA the importance then divides by the total SD, as the thesis does.

**Effect:**
- **Paper 1 (CVA) is unchanged:** `"auto"` resolves to `"mahalanobis"`. The Pima scripts and the
  loan harness must reproduce their numbers.
- **PCA results change.**
- **An explicit `"euclidean"` changes meaning.** It used to be the 2D plane distance; it becomes the
  X-space distance. Under standardised PCA the two are identical at the search point. Under CVA it is
  now a third metric, not the plane distance.

## Files to read first

- `R/local_cf.R`:
  - `bl_find_local_cf()` (the signature at :590, the pair loop at :760-830, the roxygen sections
    "Cross-pair selection criterion" and `@return`);
  - `.round_integers()` (its selector at about :1180, `sel <- if (distance == "mahalanobis" ...)`);
  - `print.bl_local_result()` (the per-pair table at about :1835);
  - `plot.bl_local_result()` (the summary block at about :1750).
- `R/boundary_plot.R`: `plot.bl_boundary()` (:102, `denom` at :138, the Mahalanobis check at :222-249)
  and `bl_robustness()` (:282).
- `R/boundary.R:450-462`: the `dist_x` pattern to reuse.
- `R/utils.R`: home for the new resolver.
- **Tests:** `tests/testthat/test-local_cf.R`, `test-boundary_plot.R`, `test-boundary.R`,
  `test-set_scaling.R` (their fixtures' `method`).
- **Docs:**
  - `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md` §3.4, §4.4;
  - `documentation/boundary-distance-and-robustness.md`;
  - `documentation/methodology-changes/methodology-changes.qmd` §3.1, §3.5, Tables 1, 2 and 5;
  - CLAUDE.md §5, §6 and §9;
  - `2 implementation_summary.txt`;
  - the walkthroughs `development/reference/review_section7_to_8.md` and
    `review_section9_to_15.md`.

## Design

### 1. Shared resolver, `R/utils.R`

```r
.resolve_distance <- function(distance, bl_result) {
  distance <- match.arg(distance, c("auto", "mahalanobis", "euclidean"))
  if (distance != "auto") return(distance)
  if (identical(bl_result$method, "PCA")) "euclidean" else "mahalanobis"
}
```

- CVA, and any older `bl_result` without `method`, resolves to `"mahalanobis"`, the previous default.
- The existing fallbacks stay in each caller: a missing `metric_inv` or `metric` falls back to
  `"euclidean"` with a warning.

### 2. Phase 3, `bl_find_local_cf()` (`R/local_cf.R`)

- **Signature:** `distance = c("auto", "mahalanobis", "euclidean")`, then
  `distance <- .resolve_distance(distance, bl_result)`.
- **One distance helper** for the pair loop and the rounding:

  ```r
  .cf_sq_distance <- function(v, distance, metric_inv, X_sd)
  ```

  Here `v` is an n x p matrix of `B_x - x_obs`. It returns squared distances: `rowSums((v %*%
  metric_inv) * v)` for Mahalanobis, and `rowSums(sweep(v, 2L, X_sd, "/")^2)` for Euclidean.
- **The pair loop:**
  - compute `dist_x` (the square root of the Euclidean value) for every pair, beside `dist_z` and
    `dist_mah`, and store `all_distances_x`;
  - the selector is `dist_mah` for Mahalanobis and `dist_x^2` for Euclidean (it was `dist_z`);
  - the verbose message shows the active distance.
- **The result** gains `dist_x` and `all_distances_x`. `distance` stores the resolved metric. The
  no-solution path sets `dist_x = NA`. `dist_*` stay the search point's, as now.
- **`.round_integers()`:** in its Euclidean branch, rank the candidates by
  `.cf_sq_distance(cand - x_obs, "euclidean", ...)` in X-space, instead of their projected Z distance.
  It receives `X_sd` from `bl_result`, not `sv`, because `sv` is 1 under CVA.
- **Print and plot:**
  - `print.bl_local_result()`'s per-pair table gains a `dist_x` column;
  - the plot summary adds "Euclidean (X)", and its "Selector" line shows the resolved metric (plus
    "(auto)" when it was resolved by method).
- **Roxygen:**
  - `@param distance`;
  - the "Cross-pair selection criterion" section: the `"auto"` rule with the thesis citation; the
    definition of `"euclidean"`; equality with the plane distance under standardised PCA; the CVA
    identity $d^2_M = (n-g)\,d_z^2$ kept;
  - the "Integer features" selector sentence;
  - `@return`: `dist_x`, `all_distances_x`, `distance`.

### 3. Phase 2, `plot.bl_boundary()` and `bl_robustness()` (`R/boundary_plot.R`)

- `distance = c("auto", "mahalanobis", "euclidean")`, resolved with `.resolve_distance()`.
  - Under PCA the default divides by `X_sd` (the total SD, thesis p.103), and no Mahalanobis check is
    printed.
  - Under CVA nothing changes.
- Roxygen: `@param distance` and the `@return` `distance` and `mahalanobis_check` entries.

### 4. Unchanged

- `bl_find_boundary()`: its selection is the plane distance, and it already reports `dist_x` and
  `dist_mahalanobis`.
- `bl_shapley()` and `bl_find_sparse_cf()`.
- The `"fixed"` restore.
- `bl_surrogate()`.

## Tests (`tests/testthat/`)

**Existing tests that are affected:**
- **`test-local_cf.R`:** most fixtures are PCA (`lcf_bl_results`, `:5-13`), and none asserts a pair or
  distance value, so they should pass unchanged.
  - The synthetic `run_round()` (`:301-305`) calls `.round_integers()` with `distance` and
    `metric_inv`; add the new `X_sd` argument.
  - Its Euclidean branch test (`:310-327`, comment `:323`) is re-checked against X-space distances.
- **`test-boundary.R:62-78`:** the RNG-independence test on a PCA fixture still holds, since it compares
  two runs. Its pair moves from (3,4) to (1,4) (rotation note §3.4), which is no test failure.
- **`test-set_scaling.R:154-165,220-234`** (PCA): no pair values are asserted, so they should pass.
- **Use the `formals()` pattern of `test-local_cf.R:61-70`** to test the new default, as for
  `max_pairs`.

**New tests:**

1. **`.resolve_distance()`:**
   - `"auto"` gives `"euclidean"` under PCA and `"mahalanobis"` under CVA;
   - `"mahalanobis"` when `method` is missing;
   - explicit values pass through;
   - an invalid value errors.
2. **`.cf_sq_distance()`:** matches hand-computed Mahalanobis and total-SD Euclidean values on a small
   matrix.
3. **PCA fixture (iris two-class, standardised), default call:**
   - `distance == "euclidean"`;
   - `all_distances_x` equals `all_distances` within 1e-8 for every pair (no filters; the
     orthonormality identity);
   - `best_pair` is the argmin of `all_distances_x`;
   - `dist_x` equals `sqrt(sum(((B_x - x_obs) / X_sd)^2))`.
4. **PCA with `distance = "mahalanobis"`:** `best_pair` is the argmin of `all_distances_mahalanobis`.
5. **CVA fixture:** the default stores `"mahalanobis"` and gives the same `best_pair`, `B_x` and
   `B_pred` as an explicit `"mahalanobis"` (Paper 1 unchanged).
6. **Integer rounding under `"euclidean"`:** with a PCA fixture and an `"integer"` filter, the reported
   point minimises `.cf_sq_distance()` among the valid whole-number candidates.
7. **Phase 2:**
   - under PCA, `plot.bl_boundary(type = "boxplot")` and `bl_robustness()` default to `"euclidean"`,
     with `mahalanobis_check` NULL;
   - under CVA, `"mahalanobis"` with the check;
   - `test-boundary_plot.R`'s existing check test passes `distance = "mahalanobis"` explicitly if its
     fixture is PCA.

## Scripts and measurements

- **Paper 1 Pima (CVA):** run `run_all.R`, then `tests.R` (13/13). `numbers.rds` must be identical
  apart from the version fields.
- **`04_method-change-effects.R`:**
  - **M5:** keep the CVA identity (the per-pair ratio). Drop the explicit `"euclidean"` CVA
    comparison, which now compares a different metric.
  - **Add M5-PCA:** a PCA `bl_results` on the 4-feature model.
    - Run 1 for both targets under `"auto"` (Euclidean) and `"mahalanobis"`: `best_pair`, distance,
      `B_x`.
    - An unconstrained sweep over the test rows that counts the rows where the argmin of
      `all_distances_x` and the argmin of `all_distances_mahalanobis` differ. One run per row, using the
      stored per-pair values.
  - **Add M1-PCA:** the PCA Phase 2 totals under the new default (total SD) against the within-class SD.
- **Loan harness** (`development/03_loan_status_Boundary_Logic.R`, CVA): re-run. The numbers are
  unchanged; refresh the console samples that show the new `dist_x` column.
- **Rotation scripts `documentation/biplot-rotation/17_1`, `17_2`, `17_3` (PCA, no `distance`
  passed):** their measurements were made under the Mahalanobis default.
  - `17_2`'s row choices (`:163-172`) depend on it, and so do `17_3`'s expected values
    (`:123-126`, 1.147879 / 1.312277).
  - Pin `distance = "mahalanobis"` explicitly in every `bl_find_local_cf()` call of the three scripts.
    That covers `17_1:522-523`, `17_2:264-265,301` and `17_3:123-124,242-244,394-395,422,575,593,604`;
    `17_3:576` already passes `"euclidean"`.
  - Update their comments: "the default" becomes "pinned to the selector the note measured".
  - Re-run all three and confirm every number in the note is unchanged.
  - The note says in one line (§3.4 or Appendix B) that the scripts pin the old selector.
- **README.md:** the PCA example (`:88-118`) now selects by Euclidean distance. Re-run it and update any
  quoted output; the CVA example (`:195`) is unchanged.
- **Vignette:** the loan vignette is CVA and its text stays correct. The pkgdown `docs/` is already
  stale and is not rebuilt here.
- **Policy:** CLAUDE.md §4 forbids changing default conventions "without explicit instruction". This
  change is explicitly instructed (2026-11-05), so record that in NEWS.

## Notes to update ("all the notes")

Line numbers come from the 2026-11-05 inventory. If a line has moved, find the passage by its text.

- **Roxygen, then `man/`** through `devtools::document()`.
  - Includes the `.compute_metric_inverse` and `bl_build_projection` text that ties `W_binary` to the
    local selector (`R/projection.R:13-16,134-137,150-156`).
  - Includes `R/boundary.R:147-149` ("as the field of the same name in `bl_find_local_cf()`"), which
    stays true.
- **`1 Foundation intro documents.txt:239`:** the `bl_find_local_cf()` signature default.
- **`documentation/fixed-constraint-restore.md:159-161`:** the integer-rounding selector ("Euclidean
  distance in the rotated plane") becomes X-space.
- **CLAUDE.md:**
  - `:159`: the rounding selector;
  - `:162`: the selector bullet;
  - `:200`: §9, the "PCA pair-selector default" bullet → RESOLVED 2026-11-05;
  - `:204`: the "still untested" list, since selection by argmin is now tested.
- **`documentation/known-issues.md:184`:** the PCA selector entry → resolved.
- **The rotation note:** §2.5 end (`:283`), §3.4 (`:470-594`, the PCA subsection and table `:576-594`),
  §3.5 (`:598-642`), §3.6 (`:660-665`), §4.4 (`:742-762` → resolved), §5.4 (`:1025-1035`) and §6.2
  item 4 (`:1095-1098`).
- **`review_section9_to_15.md`:** `:465-486` (the selector text), `:518-533`, `:548-550`, `:1104` and
  `:1114-1117`.
- **`NEWS.md`:**
  - the new `"auto"` default (thesis rule by method, both phases);
  - `"euclidean"` is now the total-SD-standardised X-space distance (it was the plane distance; they
    are identical under standardised PCA at the search point);
  - new fields `dist_x` and `all_distances_x`;
  - PCA results change, CVA does not.
- **`2 implementation_summary.txt`:**
  - the `bl_local_result` S3 fields;
  - the cross-pair selector description;
  - the Phase 2 standardisation by method;
  - the integer-rounding selector;
  - §4.6.1 where relevant.
- **CLAUDE.md:**
  - §5, the "Distance-measure functions" convention: `"auto"` default by method;
  - §6, the selector bullet: PCA now follows the thesis, and the CVA identity stays;
  - §6, the Phase 2 bullet: by method;
  - §9: mark the "PCA pair-selector default" bullet resolved, with the date.
- **`documentation/biplot-rotation/biplot_rotation_and_pair_selection.md`:** §3.4 (the implementation)
  and §4.4 (problem → resolved 2026-11-05), plus the folder's README if it lists §4.4.
- **`documentation/boundary-distance-and-robustness.md`:** the Phase 2 default by method; the pair
  selector section.
- **`documentation/known-issues.md`:** move the PCA selector entry to the resolved table, if it is
  listed.
- **`documentation/methodology-changes/methodology-changes.qmd`, then re-render:**
  - the Summary: two departures from the text remain;
  - Table 1: a 2026-11-05 row;
  - Table 2: "Comparison across pairs" becomes restored; "Variable importance" changed under CVA
    only;
  - §3.1: PCA now uses the total SD;
  - §3.5: PCA restored, with the M5-PCA numbers;
  - Table 5: drop "Euclidean pair comparison under PCA".
  - In the `claims` chunk, replace the M5 `same_pair` / `same_counterfactual` assertions (CVA explicit
    `"euclidean"`, which now means a different metric) with assertions on the M5-PCA rows. Under CVA
    the identity rests on the per-pair ratio alone.
- **The walkthroughs:** `review_section9_to_15.md` (the `bl_find_local_cf()` Parameters table and the
  per-pair sample) and `review_section7_to_8.md` (the `plot.bl_boundary()` Parameters table), from a
  real harness run.
- `progress.md`; this plan archived to `.claude/reference/` with an IMPLEMENTED banner.

## Error cases

| Case | Behaviour |
|---|---|
| An older `bl_result` without `method` | `"auto"` resolves to `"mahalanobis"` (the previous default) |
| `metric_inv` / `metric` missing with Mahalanobis | The existing fallback to `"euclidean"`, with a warning |
| An invalid `distance` | A `match.arg()` error |
| An explicit `"euclidean"` under CVA | Allowed; documented as not the thesis's CVA rule |
| `X_sd` is 0 for a constant feature | As in `bl_find_boundary()`: not guarded, because PCA cannot be built on a constant column |
| A Pima or harness number changes under CVA | Stop and investigate: `"auto"` must reproduce the old CVA results exactly |

## Verification

1. `devtools::document()`, then `devtools::test()`: all pass; report the counts. `R CMD check`:
   0/0/0.
2. Harness end to end: the CVA numbers are unchanged; refresh the walkthrough samples from its log.
3. Pima: `run_all.R`, `tests.R` 13/13, `04_method-change-effects.R` self-checks; re-render the note.
4. Rotation scripts `17_1`-`17_3`, with `distance = "mahalanobis"` pinned: re-run them through the
   scratchpad driver (CLAUDE.md §5). Their console numbers equal the note's.
5. README PCA example: re-run it; the quoted output matches.
6. Stop and report. No commit until you say so.
