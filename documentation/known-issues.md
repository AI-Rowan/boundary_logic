# Known Issues: Open Defects and Limitations

**Author:** Adriaan Rowan
**Date:** 2026-11-04
**Package:** `boundarylogic`
**Branch:** `paper_nr1_changes`
**Location:** `documentation/known-issues.md`
**Companion documents (paths relative to this file):**
- `../CLAUDE.md` §9 -- deferred and future work (the methodological items this note does not repeat)
- `fixed-constraint-restore.md` §8 -- the future developments for `"fixed"` and `"integer"`
- `biplot-rotation/biplot_rotation_and_pair_selection.md` §4 and §6 -- the rotation's problems and further work
- `../.claude/reference/known-issues-taskD-plan.md` -- the approved plan that wrote this note

---

This note records the open defects and limitations found while preparing Paper 1 (Task D of the
Paper 1 fixes). None is fixed for the paper. Each entry gives the symptom, where it is in the code,
its effect on results, and a suggested fix. Line numbers are those of 2026-11-04.

The Pima numbers come from the Paper 1 example (`scripts/paper1-pima/`, local), run on the GAM
refitted on the four features kept after the global search: Pregnancies, Glucose, Insulin and Age.

## 1. Open issues

| # | Issue | Kind | Moves results if fixed? |
|---|---|---|---|
| 1 | The 2D-boundary approximation gap has not been measured | experiment | no (it measures) |
| 2 | Only the winning pair's plane is stored | package | no |
| 3 | `plot.bl_projection()` uses its own confusion colours | display | no |
| 4 | The `"fixed"` tolerance is hard-coded | package | yes, for `"fixed"` runs |
| 5 | No follow-up search when the `"fixed"` restore fails | method | yes, for `"fixed"` runs |
| 6 | The local search does not work in the free features only | method | yes, for `"fixed"` runs |
| 7 | Integer rounding cannot let the other features compensate | method | yes, for `"integer"` runs |
| 8 | The README's constraint example finds no counterfactual | documentation | no |
| 9 | The harness's CVA biplot title says "prior defaulters" | documentation | no |

### 1.1 The 2D-boundary approximation gap has not been measured

- **Symptom.** The counterfactual searches find the boundary as a contour of the model's
  predictions on an `m x m` grid in a 2D plane (`grDevices::contourLines()`), not in the
  p-dimensional feature space. How far the 2D boundary point lies from the nearest point of the
  p-dimensional boundary has never been measured.
- **Where.** `R/biplot_grid.R:277` (the global grid's contour, used by `bl_find_boundary()`);
  `R/local_cf.R:698` (the local grid's contour, per eigenvector pair, in `bl_find_local_cf()`).
- **Effect on results.** Validity is not affected: every local counterfactual is re-scored by the
  model and must cross the cutoff. The distance is: a counterfactual is the nearest boundary point
  *within the searched planes*, which may be further from the target than the nearest point of the
  p-dimensional boundary. Paper 1's distances are therefore upper bounds of unknown tightness.
- **Suggested fix.** An experiment, not a code change: for a sample of targets, compare the
  package's counterfactual distance with the distance from a p-dimensional search (e.g. a
  constrained optimiser minimising the Mahalanobis distance subject to the model crossing the
  cutoff), and report the ratio by model and data set.

### 1.2 Only the winning pair's plane is stored

- **Symptom.** `bl_find_local_cf()` searches every eigenvector pair but keeps the plane, the grid
  and the boundary point of the winning pair only. The other pairs keep only their distances.
- **Where.** `R/local_cf.R:645` (`all_distances` and `all_distances_mahalanobis`, the only
  per-pair record) and `R/local_cf.R:819` (`best_result`, which overwrites `Vr_rot`, `B_x` and
  the grid when a nearer pair is found).
- **Effect on results.** None on the counterfactual. Two consequences:
  - another pair's plane cannot be plotted without a second search pinned to it (Pima Figure 6(b)
    needs run 5, `max_pairs = 1L`);
  - the `"fixed"` restore and the `"integer"` rounding are tried on the winning pair's point only,
    so a target can lose a counterfactual that another pair's point would have given
    (`fixed-constraint-restore.md` §8, third bullet).
- **Suggested fix.** Keep each pair's plane (`Vr_rot`, `tVr_rot`) and boundary point (`B_x`, `B_z`,
  `B_pred`, distances) in a list on the result; keep the grids out, and rebuild one on demand for
  plotting. An argument (e.g. `keep_pairs = FALSE`) would keep the default result small.

### 1.3 `plot.bl_projection()` uses its own confusion colours

- **Symptom.** `plot.bl_projection()` colours the outcome groups TP green, TN steel blue, FP orange
  and FN red. Every other biplot uses TP red, TN blue, FP purple and FN orange.
- **Where.** `R/projection.R:317` (`col_map` in `bl_build_projection()`, used by
  `plot.bl_projection()`), against `R/boundary_plot.R:163` and the local plot in
  `R/local_cf.R:1629`.
- **Effect on results.** Display only. Red means FN in one plot and TP in the others, so a reader
  comparing a projection plot with any other biplot can misread the groups.
- **Suggested fix.** One shared colour map (a private helper) used by every plot method.

### 1.4 The `"fixed"` tolerance is hard-coded

- **Symptom.** The local search accepts a `"fixed"` feature within ±0.5 of its observed value, in
  model units, for every feature. There is no argument to change it.
- **Where.** `R/local_cf.R:73` (`.fixed_tol <- 0.5`).
- **Effect on results.** Its meaning depends on each feature's scale: 0.5 of a training SD for a
  standardised feature (±3,001 for `loan_amnt` in the loan harness), but 0.5 of a pregnancy or of a
  year in the raw-unit Pima model. Too tight a tolerance leaves the search without candidates: in
  the Pima example the TP target with Pregnancies and Age fixed (run 3) has no boundary point
  within it on any of the 6 pairs. Too loose a tolerance makes the restore fail more often (§1.5).
- **Suggested fix.** A `fixed_tol` argument to `bl_find_local_cf()`, per feature, in raw units or
  as a fraction of the within-group SD (`sqrt(diag(W))`). The restore keeps the counterfactual
  exact whatever the tolerance.

### 1.5 No follow-up search when the `"fixed"` restore fails

- **Symptom.** When the search's point, with its fixed features set back to their observed values,
  no longer crosses the boundary, `bl_find_local_cf()` reports no counterfactual (strict
  `"fixed"`), even when one that holds the fixed features exists.
- **Where.** `R/local_cf.R:1018` (`.restore_fixed()`); `fixed-constraint-restore.md` §4 and §8
  (first bullet).
- **Effect on results.** In the Pima example, the TN target's runs 3 (Pregnancies and Age fixed) and
  4 (Pregnancies fixed, Age cannot decrease) find boundary points, but the restored points score
  0.475 and 0.487, so neither run has a counterfactual. On the 7-feature model, before the pruning,
  the TP target's runs 3 and 4 failed the same way (0.561 and 0.513).
- **Suggested fix.** From the restored point, extend the move in the free features along the
  search's own direction until the model crosses the cutoff, within `train_ranges` and the other
  constraints. This is a 1-D, deterministic search; its point lies off the plane, as a restored
  point already does.

### 1.6 The local search does not work in the free features only

- **Symptom.** The rotation and the pair search use every feature, fixed ones included, and the
  fixed features are held only by the tolerance and the restore afterwards.
- **Where.** `R/local_cf.R:45` (`.bl_rotate()`) and the pair loop at `R/local_cf.R:650`;
  `fixed-constraint-restore.md` §8 (second bullet).
- **Effect on results.** With half the features fixed, few points on a 2D plane hold them. In the
  Pima example, with two of four features fixed, only the TP target's run 4 has a counterfactual
  (TP run 3 has no boundary point; TN runs 3 and 4 fail the restore).
- **Suggested fix.** Hold the fixed features at their observed values *before* the rotation, so
  that the local biplot and its boundary are built in the subspace of the actionable features;
  every boundary point then holds the fixed features exactly and no restore is needed. This changes
  the projection the local search uses, so it is a change of method and needs the same approval as
  the other changes since the thesis.

### 1.7 Integer rounding cannot let the other features compensate

- **Symptom.** `.round_integers()` tries only the whole values just below and above each
  `"integer"` feature's counterfactual value, with the other features as the search left them.
  When the nearer whole value does not cross the boundary, the other one is taken, even when a
  small change in the free features would let the nearer one cross.
- **Where.** `R/local_cf.R:1110` (`.round_integers()`); `fixed-constraint-restore.md` §5.
- **Effect on results.** In the Pima example, the TN target's run 1 search point has Pregnancies
  10.21 (√d²_M 0.14, prediction 0.500). Pregnancies 10 does not cross (0.487), so 11 is taken: √d²_M
  rises to 0.42 and the prediction to 0.548. Run 5 behaves the same way (√d²_M 0.18 to 0.44). The
  counterfactual distance is then mostly the cost of a whole pregnancy.
- **Suggested fix.** When a whole value does not cross, hold it and run the follow-up search of
  §1.5 along the free features; keep the result if it is nearer than the other whole value.

### 1.8 The README's constraint example finds no counterfactual

- **Symptom.** The README's "Actionability constraints" example takes the first test patient the
  Pima GAM predicts as diabetic and sets Glucose to `"decrease"`, Age and Pregnancies to `"fixed"`.
  `bl_find_local_cf()` returns `solution_found = FALSE`: every segment on all 15 pairs is
  eliminated by the constraints (re-run 2026-11-04).
- **Where.** `README.md:184` (the example) and `README.md:137` (its set-up: the GAM and
  `bl_results_pima`).
- **Effect on results.** None on the package. A reader who runs the README's demonstration of
  `set_filters()` gets only "No counterfactual satisfying all constraints was found".
- **Suggested fix.** Choose constraints or a patient for which a counterfactual exists (checked by
  running the example), or keep the example and say in the README that it finds none and why.

### 1.9 The harness's CVA biplot title says "prior defaulters"

- **Symptom.** The first-pass result in the loan harness is titled "Loan default (prior defaulters)
  -- XGB, CVA biplot", but Step 3 keeps only the applicants **without** a prior default on file.
- **Where.** `development/03_loan_status_Boundary_Logic.R:251` (the title);
  `development/03_loan_status_Boundary_Logic.R:142` (the filter).
- **Effect on results.** None; the title misdescribes the population on the plot.
- **Suggested fix.** "Loan default (no prior default) -- XGB, CVA biplot", and the same wording
  where the walkthrough quotes it (`development/reference/review_section4_to_6.md:358`).

## 2. Resolved since the Task D list was written

The Paper 1 task list (2026-10-29) and its planning named these for this note; they were fixed
before it was written.

| Issue | How it was fixed | Commit |
|---|---|---|
| `best_pair` not reproducible across identical calls | The boundary searches take the exact nearest contour vertex (`max.col(..., ties.method = "first")`), not a random one among near-ties; see `nearest-vertex-tie-break.md`. The rotation plane's under-determination remains, but it is deterministic: CLAUDE.md §9 and the rotation note §4.1. | `f0e3dd4` |
| Results depend on earlier calls in the session | The same fix: no search draws from the session's random numbers. | `f0e3dd4` |
| The cutoff is not passed to `bl_build_grid()` | The grid takes the model's cutoff (`bl_model$cutoff`). | `f4f3647` |
| `rotate_deg` is not applied to the sparse counterfactual's cross | `rotate_deg` was removed from the local plots, whose orientation the optimal rotation fixes. | `f4f3647` |
| `xgboost::xgboost()` models fail in prediction | `.pred_function()` accepts them. | `f4f3647` |
| An empty outcome group breaks four-group CVA | `bl_build_projection()` drops empty groups, with a message. | `f4f3647` |
| `used_in_sparse` is `TRUE` for a fixed feature after the restore | The flag is `FALSE` for fixed features, matching what the sparse counterfactual does. | `7aa1344` |
| Stale documentation: `dist_z`, the sparse colours and the Shapley `seed` in roxygen; the README; the "stratified" split; the walkthrough's `set_filters(person_age = "fixed")` examples after `person_age` was pruned | Corrected in the documentation sweep and the walkthrough review, and re-checked in the documentation audit of 2026-11-04. | `f2a076b`, `7aa1344` |

## 3. Other future work

The methodological future work is in `../CLAUDE.md` §9 and is not repeated here: the minimal
Givens rotation with a plane-family search (rotation note §5 and §6), the PCA pair-selector default,
the full-Mahalanobis Shapley decomposition in Phase 2, per-variable standardisation methods, and
LDA multi-class biplots.
