# The convex hull on each rotated plane

Date created: 2026-11-05
Date last edited: 2026-11-15

This note explains why `bl_find_local_cf()` keeps only the boundary points inside the convex hull of
the training data on each rotated plane (`set_filters(hull = TRUE)`, the default), how it does so,
and what the rule changes on the Pima example of Paper 1 and on the loan harness. Section 6 explains
why this switch and the one for the training ranges sit in `set_filters()`, and why a search without
a filter object keeps both rules on. `documentation/algorithms.md` measures the rule's two parts on
their own: closed contours that enclose no training data (Finding F2) and a hull built from the
training data alone (Finding F4). The history of the rule is in the methodology note
(`methodology-changes/methodology-changes.md`).

Thesis pages are the printed pages of the final thesis (printed page = PDF page − 16).

## 1. The thesis rule

- **The text.** "The counterfactual must be selected from a local decision boundary that is within
  the convex hull of the training data" (§4.3.1, p.117). The hull ensures that the counterfactual is
  plausible: it lies in a region the training data covers (p.116). The rule is dropped only when an
  unconstrained solution is acceptable (p.115).
- **Appendix A7** rebuilds the hull on each rotated plane from the target and the training data
  (pp.207, 209).
- **Empty regions.** "Any local decision boundary area that does not contain data points is
  discarded" (§3.10 step 4, p.93).

## 2. Why the hull is rebuilt on each plane

The global biplot's hull polygon (from `bl_build_grid()`) has coordinates on the unrotated plane.
The local search turns the plane for every component pair (`.bl_rotate()`), so that polygon does not
describe the training data on a rotated plane, and the local search never uses it. Instead it
projects the training data onto each rotated plane and builds the hull there, as Appendix A7 does.

## 3. What the code does

For each component pair, after the rotation:

1. **The hull.** The training data are projected onto the rotated plane (`Z_train_rot`), and their
   convex hull is built with `.build_hull_polygon(..., outlie = 1)` (`R/hull_utils.R`). The fraction 1
   gives the full hull, as the global grid's default (`bl_build_grid(outlie = 1)`). The target is not
   part of it, so a target outside the training data does not enlarge the region its counterfactual
   may come from (`algorithms.md` F4).
2. **The filter.** Each boundary contour on the plane keeps only its vertices inside the hull
   (`.poly_clip()`, the vertex filter the global search uses). A contour still closed after this clip
   is kept only if it encloses training data on the plane (`.encloses_data()`, the global search's
   empty-region rule; `algorithms.md` F2). The actionability constraints (`set_filters()`), the
   training ranges and the model's re-score then run on the remaining vertices.
3. **The selection:** the nearest valid vertex per pair, and the pair with the smallest distance by
   the thesis rule.

The result stores `hull` (whether the rule applied) and `hull_polygon` (the winning plane's hull
vertices; with the rule off it is built after the search, for the plot). When no pair has a valid
point, the message names the hull and says that `set_filters(..., hull = FALSE)` drops the rule.

**The plot.** `plot.bl_local_result()` colours the prediction grid only inside the winning plane's
hull when the search applied it (`trim_to_hull`), as the global biplot colours only inside its hull
(`bl_build_grid(calc_hull = TRUE)`): the coloured area is where the search could choose boundary
points, and contour lines in the white area lie outside the hull. `trim_to_hull = TRUE` shows the
hull of a search with the rule off; the Pima Figure 7(a) uses it to show the TP target's run 4,
which has no counterfactual inside the hull, searched without it: its counterfactual is the point the
hull removes.

**Design choices:**

| Choice | What the code does | Reason |
|---|---|---|
| The points the hull is built from | The training data (Appendix A7 also includes the target) | The valid region is set by the training data alone, as the text says ("the convex hull of the training data", p.117) and as in the global search. A target outside the training data does not enlarge the region its counterfactual may come from |
| A switch | `set_filters(hull = TRUE)` by default; `hull = FALSE` gives the search without the rule (§6) | The thesis rule by default, with the unconstrained search one call away |
| The rotation scripts (`17_1` to `17_3`) | Pin `set_filters(tgt, hull = FALSE)` | Their measurements concern the rotation and the pair search, so the rotation note's numbers do not depend on the hull |
| An enclosure rule for closed contours | With the hull on, a closed contour is kept only if it encloses training data on the plane | The same empty-region rule as the global search (§3.10 step 4, p.93), tested on the training data in both searches |

**Off the plane.** The hull applies to the boundary points on the plane. The restored point of a
`"fixed"` feature and the rounded point of an `"integer"` feature lie off the plane, and are checked
against the constraints and the training ranges only. The sparse counterfactual
(`bl_find_sparse_cf()`) also lies off the plane; it is re-scored by the model and checked against no
rule.

**A target inside every training range can lie outside the hull.** This follows from how the hull
is built, and is not an error:

- The hull is built from the training rows projected onto the rotated plane. The projection is
  linear, so the hull on a plane is the projection of the convex hull of the training data in
  X-space.
- The training ranges check each feature on its own. Together they form a box, and the training data
  fill only part of it, so a target inside the box can lie outside the data's convex hull.
- The rotation places the target in the plane, so its position on the plane is exact, while the
  training rows are projections. On a plane that separates the target from the data, the target lies
  outside the plane's hull; on another plane it can lie inside. A target inside the convex hull in
  X-space lies inside the hull on every plane.
- The biplot axes do not show this. A value read off an axis is that of the projected point, so the
  axes show the training ranges only approximately.

The hull applies to the counterfactual, not the target, so the search still returns a counterfactual
inside the hull. On the Pima example, the made-up patient C (*Pregnancies* 4, *Glucose* 95,
*Insulin* 400, *Age* 35; `rule_patients` in `scripts/paper1-pima/00_setup.R`) is inside every
training range but outside the training data's convex hull in four dimensions: no training row has
*Insulin* of 350 or more with *Glucose* of 110 or less. On pair (3,4), which the search chooses, the
patient lies outside the plane's hull; on pair (1,3) it lies inside. Its counterfactual lies inside
the hull on both. Patient D (*Pregnancies* 4, *Glucose* 130, *Insulin* 30, *Age* 35) lies inside the
four-dimensional hull, and so inside the hull on every plane.

## 4. The effect

**The Pima example of Paper 1** (four-feature GAM, CVA; the methodology note's measurement M6, from
`scripts/paper1-pima/04_method-change-effects.R`, which runs every search with and without the hull):

- Over the 79 test rows (run 1), every row has a counterfactual with and without the hull. Without
  it, 11 rows would get a different one; with it, their search distance is 4.5% longer at the median.
- The TN target's counterfactual is the same with and without the hull in every run.
- The TP target's run 1 lies on pair (3,4) at a search distance of 1.84. Without the hull it would lie
  on pair (2,3) at 1.50, raising *Age* from 51 to 59.6 and *Pregnancies* from 7 to 10, a point
  outside the training data on its plane. Inside the hull every TP counterfactual lowers *Age*, to
  about 41, so run 2 (*Pregnancies* and *Age* unable to decrease) and run 4 (*Pregnancies* fixed,
  *Age* unable to decrease) have no counterfactual; without the hull each would have one (run 4's is
  Figure 7(a)). Run 3 (both fixed) has none either way. Run 5 (pair
  (1,2) only) has a boundary point inside the hull, but no whole-number *Pregnancies* value still
  crosses the boundary.

The TP target shows the trade-off between the constraints: the plausible counterfactuals need a
younger patient, which is not actionable, so no counterfactual meets both rules.

**The loan harness** (`development/03_loan_status_Boundary_Logic.R`): every number its walkthrough
quotes comes from the search with the hull on. The Step 16 applicant, outside the training data,
is the case where the rule matters (§6).

## 5. Reproduction

- `set_filters(tgt, hull = FALSE)` gives the search without the rule; the rotation scripts use it.
- The measurement M6 in `04_method-change-effects.R` records both settings for every run of both
  targets, and the sweep over the 79 test rows.

## 6. The plausibility and feasibility switches in `set_filters()`

**Where the rules run.** For each component pair, every opposing-class boundary segment on the
rotated plane passes four filters in this order (`R/local_cf.R`, `bl_find_local_cf()`):

| # | Filter | Thesis term | Switch |
|---|---|---|---|
| 1 | Inside the hull of the training data on this plane (`.poly_clip()`, 2D), and a closed contour only if it encloses training data (`.encloses_data()`) | plausibility (pp.116-117); empty regions (p.93) | `set_filters(hull = TRUE)` |
| 2 | The per-variable constraints (`.apply_actionability()`, X-space) | actionability (p.116) | the constraints in `set_filters()` |
| 3 | Every feature within its training range (`get_filter_logical_vector()`, X-space) | feasibility (pp.116-117) | `set_filters(train_ranges = TRUE)` |
| 4 | The model scores the point on the other side of the cutoff | not a constraint | always on |

After the search, the restored `"fixed"` point and the rounded `"integer"` point lie off the plane, so
they are checked against the training ranges (when that switch is on) but not the hull.

**Why one place.** The thesis asks a counterfactual to be actionable, feasible and plausible (defined
in §3.9.2, pp.88–90; applied on p.116), and drops the hull only when an unconstrained solution is
acceptable (p.115). Holding all three in `set_filters()` follows that framing, gives the training
ranges a switch of their own (otherwise the only way to drop them would be to set
`bl_result$train_ranges` to `NULL`), and makes the unconstrained counterfactual one call:
`set_filters(tgt, hull = FALSE, train_ranges = FALSE)` with no per-variable constraint.
`bl_find_local_cf()` has no `hull` argument.

**Why `set_filters = NULL` keeps both rules on.** The switches live in the `bl_filters` object, and a
`NULL` has none, so the code must give `NULL` a meaning. The private `.filter_switches()` reads it as
both on, for four reasons:

1. **The thesis's counterfactual by default.** Reading `NULL` as both off would change every search
   without filters wherever the nearest boundary point lies outside the plane's hull or outside a
   training range: the loan harness's runs without filters (Steps 15 and 16, quoted in the
   walkthrough), the vignette's and README's equivalents, the Pima runs without filters and the
   approximation-gap experiment (`05_full-space-gap.R`, whose ratio is measured with the hull on and
   whose check that every search point lies inside its plane's hull would fail).
2. **The size of that change.** Without the hull, 11 of the 79 Pima test rows would get a different
   counterfactual (§4); the TP target's would lie on pair (2,3) at 1.50, raising *Age* to 59.6,
   outside the training data. Dropping the ranges as well can move more.
3. **Consistency.** `set_filters(tgt)` with no constraint fills in both switches as `TRUE`. If `NULL`
   meant both off, "no filter object" and "an empty filter object" would give different
   counterfactuals for the same target.
4. **The thesis.** Its counterfactual is the constrained one (p.117); the unconstrained one is the
   stated exception (p.115), so it is the explicit opt-out.

A `bl_filters` object without the switch fields (built by hand) is read the same way.

**What the switches change on the loan harness and the Pima example.** The loan rows come from a
driver run after the harness, with the searches without actionability constraints. The Pima rows
come from the made-up external patients of `scripts/paper1-pima/03_local.R`, on the 4-feature GAM
with run 1's filter (`Pregnancies = "integer"`). Distances are the search's √d²_M.

| Target | Default (both on) | `hull = FALSE` | `train_ranges = FALSE` | Both off |
|---|---|---|---|---|
| Loan test row 13, Step 15 (no actionability constraints) | pair (3,4), √d²_M 0.5554 | the same | the same | the same |
| Loan Step 16 external applicant | pair (1,2), √d²_M 2.8241; credit score 534.92 | pair (1,2), 1.7026; credit score 489.44 | the same as the default | pair (1,4), 0.4750; credit score 421.26 |
| Pima patient A (Pregnancies 17, Glucose 150, Insulin 200, Age 45) | pair (2,4), √d²_M 3.1973; Pregnancies 9 | pair (1,4), 2.5145; Pregnancies 13 | the same as the default | the same as `hull = FALSE` |
| Pima patient B (Pregnancies 4, Glucose 180, Insulin 900, Age 25) | pair (2,4), √d²_M 2.4535; Insulin 641.45 | pair (2,3), 0.7658; Insulin 819.01 | the same as the default | the same as `hull = FALSE` |

The external applicant's credit score (420) lies below the training range (489 to 754), and the
applicant lies outside the training data on the chosen plane. With the hull on, the counterfactual
must lie inside the training data's hull, which already keeps the credit score within its range, so
dropping the ranges changes nothing. With the hull off and the ranges on, the nearest counterfactual
sits at the range's lower end; with both off it leaves the credit score near 420, outside the
training data, at a sixth of the default's distance. (If the hull also included the applicant, as
Appendix A7 builds it, the default would equal the `hull = FALSE` column.)

The two Pima patients lie outside the training data on every plane chosen. Patient A's *Pregnancies*
(17) is above its range (0 to 15), and patient B's *Insulin* (900) is above its range (15 to 846).
- **The ranges alone do not keep the point within the data.** With the hull off, both
  counterfactuals keep every feature within its training range (*Pregnancies* 13, *Insulin* 819.01)
  but leave the training data. The ranges hold each feature separately, not the point.
- **Dropping the ranges with the hull on changes nothing,** as for the loan applicant.
- **If the hull also included the patient** (as Appendix A7 builds it):
  - Patient B's default would equal the `hull = FALSE` column, as the loan applicant's would: that
    hull would reach out to the point.
  - Patient A's default would not change, because that hull would not reach a nearer boundary
    point.

The review page of the Pima example (`review.html`) draws both patients with the hull on and off
(`fig-ext.R`).

**The training ranges on a target inside the data.** The made-up Pima patient D (*Pregnancies* 4,
*Glucose* 130, *Insulin* 30, *Age* 35), inside every range and inside the training data, shows what
the ranges add to the hull. With both rules on, the counterfactual lies on pair (1,4) at √d²_M 0.8407,
with *Insulin* 111.72. With `train_ranges = FALSE` it would lie on pair (3,4) at 0.5686, with
*Insulin* −15.41: inside the plane's hull, but a value no patient can have. Over the 79 Pima test
rows (run 1's filter), the ranges change 6 counterfactuals (`scripts/paper1-pima/07_rule-effects.R`).
