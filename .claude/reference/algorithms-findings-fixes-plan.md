> **IMPLEMENTED (2026-11-06).** The approved plan, with the implementation notes recorded at the top.
> M-d (switching B2 to B1's gap) was measured and reverted at the step 5 gate, on the user's decision.

# Plan: address the findings F1–F8 of `documentation/algorithms.md`

Date created: 2026-11-06
Date last edited: 2026-11-06

> **Housekeeping (done 2026-11-06):** written in plan mode into
> `before-updating-the-paper-wiggly-wave.md` (the algorithms-document plan); step 0 moved it here and
> restored that file.

## Implementation notes (2026-11-06)

**The step 5 gate was hit, and M-d was reverted.** With B2 switched to B1's gap, Pima's three weakest
variables became SkinThickness, BMI and Insulin, and `02b_prune-refit.R` stopped.
- **What the measurement showed:** B1's gap is mostly the reset of the observation's off-plane part:
  - 80% of the summed squared gap on the loan 9-feature model;
  - 59% on the loan 4-feature model;
  - 78% on Pima's 7-feature model, and 63% on its 4-feature model.

  The ranking therefore measured the reconstruction error: `person_gender` went from weakest (47) to
  strongest (4,418).
- **The user's decision:** B1 and B2 differ by design. B1's distances are reported only, and pick the
  Pima targets. Robustness is B2's (`bl_robustness()`).
- **What was done:** `R/boundary_plot.R` was restored from `HEAD`. The roxygen of `plot.bl_boundary()`
  and `bl_robustness()` now explains the difference. The test asserts that `bl_robustness()` equals the
  plot's totals and uses the in-plane gap.

**Done:**
- **Code:**
  - M-a and M-b, through `.encloses_data()` in `R/hull_utils.R`;
  - M-c;
  - the local clamp;
  - the local plot window (`.local_plot_limits()`);
  - the `grid_col` cutoff split;
  - `inside_data` removed;
  - the F1 comment.
- **Tests:** 561 PASS / 0 FAIL / 4 WARN, before the final check.
- **Regenerated:**
  - the harness log;
  - Pima `run_all.R` (log and Table 1 identical), `tests.R` (13 of 13), `04` (F3 rows added; every other
    value unchanged), `05` and the new `06_findings-effects.R`;
  - the methodology note, rendered;
  - Pima `review.qmd`;
  - the vignette, precompiled (25 figures; Step 16 now shows the applicant).
- **Documents:**
  - `algorithms.md` §7, "Findings and decisions";
  - the boundary note (§2 rule and F3; line numbers);
  - the local-hull note;
  - the rotation note (line numbers);
  - both walkthroughs (the ·S step; review_section7_to_8 Step 7b numbers; review_section9_to_15 Steps 11, 13,
    14, 15 and 16 numbers);
  - the implementation summary, CLAUDE.md, NEWS and README;
  - the Pima `fig-07.R` and `review.qmd` wording.

**Measured effects:**
- **M-a** (B1 rule on the training data):
  - the loan 4-feature model keeps 18 contours instead of 12, and 346 of 4,429 counterfactuals move
    (Step 7b robustness 6813.80 → 6790.20);
  - the loan 9-feature model and Pima are unchanged.
- **M-b** (C4 empty-region rule):
  - harness row 13's constrained counterfactual moves from √d²_M 0.5838 to 0.6934 (prediction
    0.547 → 0.693; verified to lie on a closed level-0.510 contour with no training data);
  - its sparse counterfactual now changes only the interest rate.
- **M-c** (training-only hull): the new applicant's counterfactual moves from √d²_M 1.70 to 2.82
  (credit score 489.44 → 534.92).
- **Pima:** no run of either target changes. 4 of the 79 test rows' run-1 counterfactuals move
  (rows 16, 46, 57 and 63; median distance +4.5%).

## Context

`documentation/algorithms.md` §7 lists eight code behaviours that no other document recorded. The user
has decided each one (2026-11-06). Four decisions change what the method computes (M-a to M-d below).
These move Phase 2 and Phase 3 results, so the loan harness, the vignette, the Pima example (Paper 1)
and the methodology note must be re-run. The rest are notes, display fixes and documentation.

**What the code check found** (the facts behind the decisions):
- **F2:** B1's closed-contour rule tests the *searched rows* (`Z_obs`, from `data`/`tdp`), not the
  training data, although its warning says "enclose training observations". The PhD code did the same
  (`1.4 Biplot Boundary Search_optimized.R:136`, `X$Z.st[tdp, ]`). C4 has no such rule: an empty closed
  contour on a rotated plane is kept.
- **F3:** B2's gap is $(\mathbf z_i-\mathbf b_i)\mathbf V_r^{-1}\mathbf S = \hat{\mathbf x}_i-\mathbf x^{cf}_i$, and
  B1's is $\mathbf x^{cf}_i-\mathbf x_i$. They differ by $\mathbf x_i-\hat{\mathbf x}_i$, the part of the
  observation off the biplot plane. That is zero only when the observation lies on the plane (for example
  when $p=2$), so the two are not always equal. The size of the difference is measured in step 4.
- **F4:** a target outside the training cloud does not stop the search: it still compares distances to the
  boundary points inside the hull. But the plane's hull includes the target while its grid covers only the
  training cloud plus 10%, and the plot window, set by biplotEZ from the training points (`plot.biplot`
  `xlim`/`ylim`), leaves the target out.
- **F6:** the PhD code (`1.1 biplot_plane_optimized.R:166-232`) used `inside_data` to drop the training
  points outside a (peeled) hull from the plotted points. The package moved point plotting to
  `bl_project_points(filter_to_polygon =)`, so the variable became dead code.

## Decisions (user, 2026-11-06)

| Finding | Decision | Kind |
|---|---|---|
| F1 B1 zero test | Keep unchanged; note it in the code and in `algorithms.md` | Note |
| F2 closed contours | **M-a:** B1 keeps a closed contour only if it holds *training* rows. **M-b:** C4 drops closed contours with no training rows on that plane, only when `hull = TRUE` | Method |
| F3 gap | Measured, then **M-d reverted** (user, at the step 5 gate). B1 and B2 differ by design: B1's gap also resets the observation's off-plane part (80%, 59% and 78% of the summed squared gap on the loan 9-feature, loan 4-feature and Pima models). B1's `dist_x`/`dist_mahalanobis` are reported only (and pick the Pima targets); the importance and robustness come from B2 (`bl_robustness()`, `plot()`'s totals). The roxygen says so; the notes record it | Docs |
| F4 clamp | Clamp the local grid's scores to [0, 1], as A3 does | Code |
| F4 target | **M-c:** each plane's hull is built from the **training data only** (the target is no longer included). The local plot window widens to show the target and the counterfactual | Method + display |
| F5 Gini | Accept; note it | Note |
| F6 plot split | `plot_biplotEZ()`'s two-colour grid splits at `cutoff`, not at 0.5 | Display |
| F6 `inside_data` | Remove it; record its history | Code |
| F7 Shapley re-scoring | Note it; no change | Note |
| F8 documentation | Recommendations (below), applied | Docs |

## Files to read first

- `R/boundary.R` (B1, `.is_closed()`), `R/boundary_plot.R` (B2: `plot.bl_boundary()`, `bl_robustness()`,
  the Mahalanobis check), `R/local_cf.R` (the pair loop: grid, hull, contour filters, `col_value`; `set_filters()`
  roxygen; `plot.bl_local_result()`), `R/hull_utils.R`, `R/biplot_grid.R` (the clamp, `inside_data`),
  `R/plot_biplot.R` (the `grid_col` split), `R/shapley.R` (`plot.bl_sparse_result()`).
- `documentation/local-convex-hull.md`, `documentation/boundary-distance-and-robustness.md` §2–4,
  `documentation/methodology-changes/methodology-changes.qmd` (the hull and importance sections) and
  `scripts/paper1-pima/04_method-change-effects.R` (the M6 block is the pattern for new measurements).
- `documentation/algorithms.md` (B1, B2, C4, A3 and §7), the walkthroughs `review_section7_to_8.md` and
  `review_section9_to_15.md`, and the tests `test-boundary.R`, `test-boundary_plot.R`, `test-local_cf.R`.

## Design

**Shared helper (new, `R/hull_utils.R`):** `.encloses_data(Mi, Z_ref)`. It returns `TRUE` when the polyline
`Mi` is not closed (`.is_closed()`) or when at least one row of `Z_ref` lies in it
(`sp::point.in.polygon(...) > 0`). B1 and C4 both use it, so the rule is defined once.

| ID | Where | Change |
|---|---|---|
| M-a | `bl_find_boundary()` | Project `bl_result$train_data` as the search projects `data` (same centre, `sv`, `Vr`), and test closed contours against these training rows, not `Z_obs`. A row's result no longer depends on which other rows are searched (F2) |
| M-b | `bl_find_local_cf()`, pair loop | With `hull = TRUE`, after `.poly_clip()`: `if (!.encloses_data(Mi, Z_train_rot)) next` |
| M-c | `bl_find_local_cf()` | `.build_hull_polygon(Z_train_rot, ...)`, without the target, for the search and for the stored `hull_polygon` (also when the hull is off) |
| M-d | `plot.bl_boundary()`, `bl_robustness()` | `vec_to_boundary = x_obs - B_x`, in model units, divided per variable by $\sigma_j$ ($\sqrt{W_{jj}}$ or `X_sd`, as now); the step through $\mathbf V_r^{-1}$ and $\cdot\mathbf S$ goes. The Mahalanobis check's "full" value then equals $\sum_i$ `dist_mahalanobis` from B1, and under `"euclidean"` the per-variable squares sum to `dist_x`² |
| F4 clamp | `bl_find_local_cf()` | Before `contourLines()`: `grid_prob <- pmin(pmax(grid_prob, 0), 1)`, with one warning per call counting the clamped cells over all pairs (A3's wording) |
| F4 window | `plot.bl_local_result()` | When the target or the counterfactual lies outside biplotEZ's default window (`range(Z * 1.2)`), pass `xlim`/`ylim` covering them to both `plot()` calls (base and overlay). Otherwise nothing is passed, so plots inside the window stay byte-identical. `plot.bl_sparse_result()` inherits it |
| F6 split | `plot_biplotEZ()` | `gr$grid_prob >= bl_result$cutoff`, and the roxygen of `grid_col` updated to match |
| F6 dead code | `bl_build_grid()` | Remove `Z_df` and `inside_data` |
| F1 note | `bl_find_boundary()` | A comment at the zero test: the known stop when every contour lies on one side, kept by decision (`algorithms.md` F1) |

**Not changed:** `.empty_return()`, the Gini (F5), the Shapley re-scoring (F7), and the global grid's
range. With `hull = FALSE`, C4 still reproduces the search before 2026-11-05 exactly (M-b and M-c apply
only with the hull on), so the rotation scripts `17_1`–`17_3`, which set `hull = FALSE`, are unaffected.

**F8 recommendations (applied):**
1. **The two walkthroughs:** write the back-projection in general form,
   $\hat{\mathbf x}=\mathbf z\mathbf V_r^{-1}\mathbf S+\bar{\mathbf x}$, adding that $\mathbf S=\mathbf I$ under CVA.
   Both examples are CVA, so no number changes.
2. **Stale line numbers:** in the rotation note §3.3 and the boundary note §3.2, replace them with
   function names (with a section name where useful), as `algorithms.md` does.
3. **CLAUDE.md §5:** one line: "In notes and walkthroughs, cite functions, not line numbers; line numbers go
   stale."

## Steps

0. **Housekeeping:** restore the algorithms plan to this file, word for word, and save this plan as
   `.claude/plans/algorithms-findings-fixes.md`.
1. **Baseline:** export `HEAD` (`git archive HEAD DESCRIPTION NAMESPACE LICENSE R inst development .Rbuildignore`,
   per CLAUDE.md) into `tempdir()`; `R/` has not changed since `HEAD`. Run the harness against it and keep
   the log.
2. **Code changes,** as in the Design table. Then `devtools::document()`.
3. **Tests.** Add:
   - `.encloses_data()`: an open line is kept; a closed square is kept with a point inside and dropped
     without one;
   - B1: each row's result equals its row in a search over all rows (`tdp` independence, M-a);
   - C4: the stored `hull_polygon` is the hull of `Z_train_rot` alone: a far-out external target leaves it
     unchanged and its counterfactual lies inside it (M-c);
   - the local clamp: a custom `predict_fn` returning values in [−0.005, 1.2] gives a warning and no
     error, and `length(col_value) == nrow(Zgrid)`;
   - B2: under `"euclidean"`, `rowSums(vec_to_boundary_sd^2)` equals `dist_x^2`; under `"mahalanobis"`, the
     check's "full" equals `sum(dist_mahalanobis)` (M-d);
   - the local plot window: the limits cover an outside target (a small internal helper,
     `.local_plot_limits()`, tested directly);
   - `plot_biplotEZ(grid_col =)` runs for a model wrapped with `cutoff = 0.6`.

   Update the tests that assume the old hull or B2 values. Then run `devtools::test()`.
4. **Measure,** with a scratchpad driver that `load_all()`s the baseline root and the changed root in turn,
   on the loan harness (9 and 4 features) and Pima (7 and 4 features):
   - **F3:** per row, the off-plane part $\lVert(\mathbf x_i-\hat{\mathbf x}_i)/\sigma\rVert$ against the
     B1 gap $\lVert(\mathbf x_i-\mathbf x^{cf}_i)/\sigma\rVert$ (median and maximum share); the per-variable
     totals and their ranking, old and new; the Mahalanobis check, old and new, against
     $\sum$ `dist_mahalanobis`;
   - **M-a:** closed contours kept under the old and the new rule, and the rows whose global counterfactual
     changes;
   - **M-b and M-c:** the Phase 3 counterfactuals that change: harness row 13 and the new applicant, Pima's
     79 test targets and its Table 1 runs.
5. **Gate, only if a pruning set changes.** If M-d changes the harness's four weakest variables or Pima's
   three weakest (`02b_prune-refit.R` stops on that), stop and report the new rankings. Choosing the
   variables to prune is the analyst's decision, and it changes Paper 1's model. Otherwise continue.
6. **Regenerate:**
   - run the harness and update the walkthrough numbers (`review_section7_to_8.md`,
     `review_section9_to_15.md`);
   - Pima: `run_all.R`, `tests.R`, `review.qmd`;
   - extend `04_method-change-effects.R` with blocks for M-a to M-d (old variant computed inline where
     possible, as K5 does; otherwise from step 4's two-root run, recorded by a new
     `06_findings-effects.R` that exports `HEAD` itself); re-run `04` and `05`; add the four method changes
     and the code improvements to `methodology-changes.qmd` and render it;
   - re-run the vignette precompile: restore `plot(bl_sparse_ext)` (now visible); update the hull wording
     and the Step 7, Step 16 and pruning text and claims.
7. **Documents:**
   - **`algorithms.md`:** B1 step 3, B2's formula, C4's hull and contour steps, A3; §7 becomes "Findings and
     decisions", with the F3 measurements;
   - **technical notes:** `boundary-distance-and-robustness.md` (§2 the closed-contour rule; §3–4 the gap
     and the check, with the reason for the change); `local-convex-hull.md` (training-only hull, the empty
     contour rule, the plot window);
   - **summaries and rules:** `2 implementation_summary.txt`; CLAUDE.md §4 and §6 (the hull bullets) and
     §5 (the F8 line);
   - **NEWS:** under 0.2.0, "Breaking changes" for M-a to M-d and "Bug fixes and other improvements" for
     the rest;
   - **wording sweep:** the roxygen and print text that says "hull of the training data and the target",
     and the harness Step 10 comment;
   - **F8:** the walkthrough $\cdot\mathbf S$ fix and the stale line numbers;
   - **archive:** this plan in `.claude/reference/` with an IMPLEMENTED banner.
8. **Checks:**
   - `devtools::document()` and `devtools::test()`;
   - `devtools::check()`: 0 errors, 0 warnings, 0 notes;
   - rebuild pkgdown (the CLAUDE.md procedure);
   - check every changed figure by eye (Step 16's local and sparse plots now show the applicant).
9. **Stop and report:**
   - the test and check counts;
   - the measured effects;
   - the changed files.

   Do not commit. Remind the user about `progress.md` and `/revise-claude-md`.

## Error cases and risks

- **The pruning ranking changes (M-d):** step 5 gate; nothing downstream is regenerated until the user
  chooses.
- **Training-only hull (M-c) leaves a target with no counterfactual on some plane, or on every plane:**
  this is the rule working (the counterfactual must lie in the training hull). Report the count; the
  no-solution message already names the hull.
- **No closed contour encloses training rows in B1:** the existing fallback (keep all contours, with a
  warning) still applies, and its wording is now accurate.
- **Plot limits and aspect:** biplotEZ's plot keeps an equal aspect, so `xlim`/`ylim` may be widened
  further. Check the Step 16 figures by eye.
- **Phase 3 values in Paper 1 move (M-b, M-c):** the methodology note and `numbers.rds` are regenerated in
  step 6. The Paper 1 plan itself is updated later, as already planned.
- **Long runs:** harness about 1 min, `run_all.R` 1 min, `tests.R` several minutes, `04` about 10 min,
  `05` about 15 min, the precompile 1 min, check about 5 min, pkgdown a few minutes.

## Verification

- `devtools::test()` passes, with the new tests; `devtools::check()` gives 0/0/0.
- **M-a:** B1 results are independent of `tdp`.
- **M-d:** B2 totals decompose B1's distances (the identities in step 3).
- **Unchanged where they should be:** with `hull = FALSE`, the local search equals the baseline's (diff
  the harness run with `set_filters(tgt, hull = FALSE)` against the baseline root), and a plot whose target
  lies inside the window is byte-identical.
- **Re-run without breaking:** the methodology note renders with its `claims` chunk passing, and the
  vignette precompile passes its `*-claims` chunks.
- `git status` shows only the files listed in steps 2–7, plus the earlier uncommitted files of this session.
