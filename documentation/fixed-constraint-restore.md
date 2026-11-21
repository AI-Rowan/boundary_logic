# Technical Note: Restoring "fixed" Features and Rounding "integer" Features after the Local Search

Date created: 2026-10-29
Date last edited: 2026-11-16

**Author:** Adriaan Rowan
**Package:** `boundarylogic`
**Location:** `documentation/fixed-constraint-restore.md`
**Companion documents (paths relative to this file):**
- `../2 implementation_summary.txt` Section 4.8 (fixed-feature restore) and 4.9 -- code mechanics
- `../development/reference/review_section9_to_15.md` Step 11, Stage H -- the harness walkthrough
- `methodology-changes/methodology-changes.md` §3.7 and §3.8 -- the history of both rules

---

## 1. What the rule does

A `"fixed"` constraint in `set_filters()` marks a feature the applicant cannot or
should not change. `bl_find_local_cf()` handles it in two steps:

1. **The search** accepts a boundary point when the fixed feature lies within
   ±0.5 (model units) of its observed value (§2).
2. **The restore**, after the best component pair is chosen, sets each fixed
   feature the search moved back to its observed value and re-scores the point.

| Situation | Result |
|---|---|
| No `"fixed"` constraint, or the search did not move the fixed feature | The search's boundary point |
| The restored point crosses the boundary and is within `train_ranges` | The restored point: fixed features at their observed values |
| The restored point does not cross, or leaves `train_ranges` | **No solution** (`solution_found = FALSE`), with the reason |

The `train_ranges` checks in this table, and in the rounding (§5), apply while
the feasibility switch of `set_filters()` is on (`train_ranges = TRUE`, the
default; also with no filter object). With `train_ranges = FALSE` the search,
the restore and the rounding all skip the training ranges.

The last row is the strict rule (§4): a failed restore means no solution, so a
reported counterfactual always holds every fixed feature exactly. Without it, the
search's point would be reported, and a counterfactual could move a "fixed"
feature by up to 0.5 model units.

A second post-search adjustment follows the restore: features marked
`"integer"` are rounded to whole values and the point is checked again (§5).
Both adjustments move the counterfactual off the rotated plane, and both are
accepted for the same reason: the adjusted point is re-scored and re-checked
(§6).

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
is 0.5 in the feature's own unit.

## 3. Why restore after the search, and not inside it

Holding a fixed feature exactly at its observed value moves the point off the 2D
plane of the chosen pair, so it cannot be searched for on that plane's contour.
The restore therefore tests one specific off-plane point, the search's boundary
point with its fixed features put back. This:

- leaves the pair search, its distances and its per-pair tables unchanged
  (they describe how the pair was chosen);
- changes only the features marked fixed, so the `"increase"`, `"decrease"` and
  range constraints, which the search point already met, still hold;
- applies the same `train_ranges` feasibility check as the search (none when
  `set_filters(..., train_ranges = FALSE)` switches it off).

## 4. The strict rule: a failed restore means no solution

**The rule.** When the restored point does not cross the boundary, or leaves
`train_ranges`, `bl_find_local_cf()` reports no counterfactual:

- `solution_found` is `FALSE`, and the counterfactual fields (`B_x`, `B_z`,
  `B_pred`, `best_pair`, `bl_counterfactual`) are empty, as for a search that
  found no boundary point;
- `fixed_check` is kept, with `applied = FALSE` and the search's point in
  `search_point`, so the failure can be inspected;
- `all_distances`, `all_distances_x` and `all_distances_mahalanobis` keep the
  per-pair search distances;
- `bl_shapley()`, and so `bl_find_sparse_cf()`, stop with "No solution found".

**Why.** A "fixed" feature is one the applicant cannot change. A counterfactual
that changes it, even by a fraction of the tolerance, is not one the applicant
can act on, and a table that shows it contradicts the constraint the user set.
Keeping the search's point instead would always return an answer, at the cost of
a counterfactual that breaks the constraint. Reporting no solution is honest
about what the search found: on the chosen plane, every boundary point near the
target needs the fixed feature to move.

**The message.** The reason is given as `blocking_constraint`, shown on the
`Reason` line of `print()`, and as a message from `bl_find_local_cf()`. It says:

1. where the search moved each fixed feature, in display units;
2. why the search allows that (the tolerance);
3. why the restored point fails: its prediction does not cross the boundary, or
   it leaves the training feature ranges;
4. that changing the constraints on the actionable features in `set_filters()`
   may make a counterfactual that holds the fixed features possible. Relaxing a
   direction or a range on another feature changes which boundary points the
   search accepts, and so which point is restored.

**An example.** In the Paper 1 Pima example (`scripts/paper1-pima/`; Phase 3
uses the GAM on four features, whose model units are raw units, so the tolerance
is half a pregnancy or half a year), the true-negative target (test row 25) with
Pregnancies and Age both fixed (run 3) gives:

> Pregnancies is "fixed", but the boundary search moved it to 10.0697 (+0.0697).
> Age is "fixed", but the boundary search moved it to 39.4768 (+0.4768). [...]
> Restoring Pregnancies, Age to their observed values gives prediction 0.4750,
> which does not cross the boundary, so no counterfactual that holds
> Pregnancies, Age fixed was found. Changing the constraints on the actionable
> features in set_filters() may make a counterfactual that holds Pregnancies,
> Age fixed possible.

The search's point scores 0.501, just past the cutoff; putting Pregnancies and
Age back lowers the prediction to 0.475. With Pregnancies alone fixed and Age
`"increase"` (run 4) the restore fails too: the search moves Pregnancies to
10.2082, and restoring it lowers the prediction from 0.500 to 0.487. The
true-positive target (test row 23) has no counterfactual in runs 3 and 4 for a
different reason: on none of the 6 pairs does it have a boundary point inside
the hull that meets the constraints, so there is no point to restore. In run 7
(Pregnancies fixed, Glucose `"decrease"`) the true negative's restore succeeds
(0.501 to 0.512) and the true positive's fails (0.499 to 0.526). The
false-positive target (test row 26, prediction 0.688) shows the restore
succeeding in all three runs: its search points score 0.499, and the restored
points 0.487 (run 3), 0.479 (run 4) and 0.486 (run 7), which still cross, so
each restored point is the counterfactual
(`scripts/paper1-pima/output/numbers.rds`). A follow-up search that holds the
fixed features is a future development (`future-improvements-review.md`
§1.3.1, "A follow-up search from the restored point").

## 5. Rounding "integer" features

**The constraint.** `set_filters()` accepts `"integer"` for a feature that takes
whole values, such as a count or a score. It combines with the other forms:
`"integer"` alone, `c("increase", "integer")` (either order) for a direction or
`"fixed"`, and `list(c(min, max), "integer")` for a range, because
`c(600, 700, "integer")` would coerce the range to character. Without
`"integer"` a feature is continuous. Whole means whole in raw units: when
`bl_result$scaling` is set, the value is converted to raw units before rounding
and back afterwards; otherwise model units are raw units.

**Why after the search.** The same argument as for `"fixed"` (§2): the search
finds boundary points on a 2D contour, and a contour vertex almost never has a
whole-number value for a given feature. The search therefore treats the feature
as continuous, and the counterfactual it finds is then put onto whole values.

**The rule** (private `.round_integers()`, after the restore):

1. **Which features.** Every `"integer"` feature that is not `"fixed"` (a fixed
   feature already holds its observed value), whose counterfactual value moved
   from the observed value and is not already a whole number (both within
   1e-8).
2. **The candidates.** For each such feature, the whole numbers just below and
   just above its value. Every combination is a candidate point, so k features
   give 2^k candidates; more than 12 at once stops with an error. A whole value
   equal to the observed raw value takes the observed model value exactly, so
   that an `"increase"` or `"decrease"` check is not upset by rounding error.
3. **The checks.** Each candidate is re-scored and kept only when its predicted
   class differs from the target's, it meets every `set_filters()` constraint
   (directions, ranges, `"fixed"`) and it lies within `train_ranges` (when that
   switch is on).
4. **The choice.** The valid candidate nearest the target in X-space by the
   search's own distance (`.cf_sq_distance()`: squared Mahalanobis, or for
   `distance = "euclidean"`, the default under PCA, the Euclidean distance
   in X-space with each feature divided by its total SD); on a tie, the
   first.
5. **No valid candidate.** No counterfactual is reported, as for a failed
   restore (§4). `blocking_constraint` gives the boundary values, says that none
   of the combinations crosses the boundary within the constraints, and that
   changing the other constraints, or removing `"integer"`, may make one
   possible.

`integer_check` records the features rounded, whether a candidate was used, the
number of candidates and of valid ones, the chosen point's prediction and the
point before rounding. `print()` adds an `Integer rounded :` line and a note.

**Why both neighbours, and the search's selector.** Rounding to the nearest
whole value can cancel a move (back onto the observed value) or, when the
observed value is not whole, reverse it, and either can undo the crossing. Trying
both neighbours and re-checking keeps every valid option. Choosing among them by
the selector that chose the pair keeps the counterfactual nearest in the same
sense as the search. In the Mahalanobis metric that need not be the whole value
nearest the observed one. In the loan harness (Step 11) `credit_score` reaches
634.0454. Both 634 (the observed value) and 635 cross, at 0.693, and 635 is
chosen because its distance from the target is √d²_M 0.6996 against 0.7031 for
634: `W` correlates the credit score with the features that move.

**The sparse counterfactual does not round.** It takes each changed feature's
value from the rounded point (§7). `bl_find_sparse_cf()`'s `round_to` argument
is deprecated: any value other than `NULL` stops with an error that points to
the `"integer"` constraint (`sparse-counterfactual.md` §3).

## 6. What the restored or rounded point is, and is not

- **It satisfies the fixed constraint exactly.** The fixed features equal the
  target's observed values.
- **It crosses the boundary, but is not on it.** Its prediction can lie further
  past the cutoff than the search point's. In the harness (Step 11) the search
  point, the restored point and the rounded point all score 0.693, because the
  model is insensitive to the loan-amount and credit-score changes there.
- **It need not be the nearest valid counterfactual.** The pair was chosen for
  the search's point. A different pair could give a restored point closer to the
  target, but the selection step does not compare restored points.
- **Its biplot position is a projection.** `B_z` becomes
  `((B_x - X_center) / sv) %*% Vr_rot`, the same projection `bl_select_target()`
  uses for the target. A point on the plane re-projects to itself; the restored
  point lands at its projection, slightly off the contour.
- **It is accepted because it is tested.** It lies off the rotated plane the
  search worked in, so the plane does not certify it. The model does: the
  point is re-scored and must lie on the other side of the cutoff, and within
  `train_ranges`. A rounded point must also meet every `set_filters()`
  constraint.
- **The same holds for a rounded point.** Its `"integer"` features are whole,
  it lies off the plane, `B_z` is its projection, and `dist_z`, `dist_x`,
  `dist_mahalanobis` and the per-pair tables still describe the search's point.

## 7. Interaction with Shapley values and the sparse counterfactual

`bl_shapley()` attributes the move from the observed point to `B_x`. After a
restore, the fixed features do not change, so their Shapley values are 0, and
`bl_shapley()` labels them "No effect". `bl_find_sparse_cf()` never changes
fixed features. It changes the "Supports" features one at a time, largest
Shapley value first, until the boundary is crossed (see
`sparse-counterfactual.md`). `used_in_sparse` is `FALSE` for fixed features,
so the stored flag matches what the sparse counterfactual does. With the strict
rule a failed restore yields no Shapley values or sparse counterfactual, so
neither can start from a point that moved a fixed feature.

The Shapley values attribute the move to the rounded point, and the sparse
counterfactual takes each changed feature's value from it, so an `"integer"`
feature is whole there too (when its observed value is whole).

## 8. Future developments

The future developments of the restore and the rounding are in `future-improvements-review.md`,
with the reason each waits: a follow-up search from the restored point (§1.3.1), a search in the
free features only (§1.3.2), retrying the restore on the next-best pair (§1.3.3) and other
rounding step sizes (§1.4.2).

## 9. Reproducing the harness result

`development/03_loan_status_Boundary_Logic.R`, Steps 10–14, with
`loan_amnt = "fixed"`, `loan_int_rate = "increase"` and
`credit_score = c("increase", "integer")`. `print(bl_local)` reports
`Fixed restored : loan_amnt` and `Integer rounded : credit_score`, ends with the
observed vs counterfactual table, and explains both steps in notes. The search's
point (pair (3,4), √d²_M 0.6934) moves `loan_amnt` to 12,523.4 and
`credit_score` to 634.0454. The counterfactual keeps `loan_amnt` at 12,000,
raises `loan_int_rate` from 12.69 to 14.44, lowers `loan_percent_income` from
0.0500 to 0.0457, raises `credit_score` to 635, and scores 0.693.
