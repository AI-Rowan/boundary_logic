# Technical Note: Restoring "fixed" Features after the Local Search

**Author:** Adriaan Rowan
**Date:** 2026-10-29
**Package:** `boundarylogic`
**Branch:** `paper_nr1_changes`
**Location:** `documentation/fixed-constraint-restore.md`
**Companion documents (paths relative to this file):**
- `../2 implementation_summary.txt` Section 4.8 (fixed-feature restore) and 4.9 -- code mechanics
- `../development/reference/review_section9_to_15.md` Step 11, Stage H -- the harness walkthrough
- `../.claude/reference/fixed-constraint-restore-plan.md` -- the approved plan that was executed

---

## 1. What changed

A `"fixed"` constraint in `set_filters()` marks a feature the applicant cannot or
should not change. Before this change, `bl_find_local_cf()` treated it only as a
tolerance: a boundary point qualified when the fixed feature lay within ±0.5
(model units) of its observed value. The reported counterfactual could therefore
move a "fixed" feature. In the loan harness, `loan_amnt` moved from 12,000 to
12,398.5.

`bl_find_local_cf()` now adds a post-selection step. After the best eigenvector
pair is chosen, each fixed feature the search moved is set back to its observed
value and the point is re-scored. If that restored point still crosses the
decision boundary and lies within `train_ranges`, it becomes the counterfactual.
Otherwise the search's point is kept. `print()` explains the outcome either way.

| Situation | Counterfactual returned |
|---|---|
| No `"fixed"` constraint, or the search did not move the fixed feature | The search's boundary point (unchanged behaviour) |
| Restored point crosses the boundary and is within `train_ranges` | The restored point: fixed features at their observed values |
| Restored point does not cross, or leaves `train_ranges` | The search's boundary point, with a note giving the reason |

## 2. Why the search uses a tolerance

The local search finds counterfactuals on a 2D boundary contour in the rotated
biplot. Each contour vertex `B_z` back-projects to a full feature vector,
`B_x = B_z %*% tVr_rot (* X_sd) + X_center`, in which every feature is a linear
function of the two plane coordinates. The condition "feature j equals its
observed value" is a single line in that plane. The contour is a polyline traced
on an `m x m` grid, so it crosses that line at isolated points that almost never
coincide with a vertex. An equality filter would therefore remove almost every
candidate.

The tolerance (`.fixed_tol = 0.5` model units, `R/local_cf.R`) keeps the vertices
near that line, so the search still has candidates. The implementation summary
gives the same reason in terms of back-projection precision.

The tolerance is absolute in model units, so its meaning depends on the
feature's scale. For a standardised feature it is 0.5 training SD: ±3,001 for
`loan_amnt` in the harness (training SD 6,002). For an unstandardised feature it
is 0.5 in the feature's own unit. This note does not change the tolerance.

## 3. Why restore after the search, and not inside it

Holding a fixed feature exactly at its observed value moves the point off the 2D
plane of the chosen pair, so it cannot be searched for on that plane's contour.
The restore therefore tests one specific off-plane point, the search's boundary
point with its fixed features put back, and keeps it only if the model still
predicts the opposite class. This:

- leaves the pair search, its distances and its per-pair tables exactly as
  before (they describe how the pair was chosen);
- changes only the features marked fixed, so the `"increase"`, `"decrease"` and
  range constraints, which the search point already met, still hold;
- applies the same `train_ranges` feasibility check as the search.

## 4. What the restored point is, and is not

- **It satisfies the fixed constraint exactly.** The fixed features equal the
  target's observed values.
- **It crosses the boundary, but is not on it.** Its prediction can lie further
  past the cutoff than the search point's. In the harness both score 0.547,
  because the model is insensitive to the loan-amount change there.
- **It need not be the nearest valid counterfactual.** The pair was chosen for
  the search's point. A different pair could give a restored point closer to the
  target, but the selection step does not compare restored points.
- **Its biplot position is a projection.** `B_z` becomes
  `((B_x - X_center) / sv) %*% Vr_rot`, the same projection `bl_select_target()`
  uses for the target. A point on the plane re-projects to itself; the restored
  point lands at its projection, slightly off the contour.

## 5. Interaction with Shapley values and the sparse counterfactual

`bl_shapley()` attributes the move from the observed point to `B_x`. After a
restore, the fixed features do not change, so their Shapley values are 0, and
`bl_shapley()` labels them "No effect". `bl_find_sparse_cf()` never changes
fixed features. It changes the "Supports" features one at a time, largest
Shapley value first, until the boundary is crossed (see
`sparse-counterfactual.md`). `used_in_sparse` is `FALSE` for fixed features,
so the stored flag matches what the sparse counterfactual does.

## 6. Reproducing the harness result

`development/03_loan_status_Boundary_Logic.R`, Steps 10–14, with
`loan_amnt = "fixed"`. `print(bl_local)` reports `Fixed restored : loan_amnt`,
ends with the observed vs counterfactual table, and explains the move in a note.
The counterfactual keeps `loan_amnt` at 12,000, raises `loan_int_rate` to 14.146
and `credit_score` to 634.22, and scores 0.547.
