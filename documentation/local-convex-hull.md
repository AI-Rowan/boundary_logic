# The convex hull on each rotated plane

Date created: 2026-11-05
Date last edited: 2026-11-06

This note explains why `bl_find_local_cf()` keeps only the boundary points inside the convex hull of
the training data on each rotated plane (`set_filters(hull = TRUE)`, the default since 2026-11-05),
how it does so, and what the rule changes on the Pima example of Paper 1 and on the loan harness.
Section 6 explains why the switch, and the one for the training ranges, sit in `set_filters()`
(since 2026-11-06), and why a search without a filter object keeps both rules on. Section 7 records
the two changes of 2026-11-06 (`documentation/algorithms.md`, findings F2 and F4): the hull is built
from the training data alone, and a closed contour must enclose training data. The decision is
recorded in `documentation/methodology-changes/future-improvements-review.md` (item B10, round 4)
and in the methodology note's section on the convex hull in the local search.

Thesis pages are the printed pages of the final thesis (printed page = PDF page − 16).

## 1. The thesis rule

- **The text.** "The counterfactual must be selected from a local decision boundary that is within
  the convex hull of the training data" (§4.3.1, p.117). The hull ensures that the counterfactual is
  plausible: it lies in a region the training data covers (p.116). The rule is dropped only when an
  unconstrained solution is acceptable (p.115).
- **Appendix A** rebuilds the hull on each rotated plane from the target and the training data
  (pp.207, 209).
- **The post-thesis branch** `original_PhD_code` builds it from the test data, which contains the
  target (`scripts/3.2 Biplot search - Test data.R`, through `biplot_plane(polygon = NA,
  no_polygon = FALSE)`), and keeps a closed contour only if it encloses the target.

## 2. The package before 2026-11-05

`bl_find_local_cf()` applied no hull. The reason recorded in the code was that the rotation
invalidates the hull polygon of the global biplot (from `bl_build_grid()`). That is correct for the
global polygon, whose coordinates belong to the unrotated plane, but it does not cover a hull rebuilt
on each rotated plane, as Appendix A does. The methodology note listed this as the largest departure
from the thesis text (its section "No convex hull in the local search").

## 3. What the package does now

For each component pair, after the rotation:

1. **The hull.** The training data are projected onto the rotated plane (`Z_train_rot`), and their
   convex hull is built with `.build_hull_polygon(..., outlie = 1)` (`R/hull_utils.R`). The fraction 1
   gives the full hull, as the global grid's default (`bl_build_grid(outlie = 1)`). The target is not
   part of it (since 2026-11-06, §7).
2. **The filter.** Each boundary contour on the plane keeps only its vertices inside the hull
   (`.poly_clip()`, the vertex filter the global search uses). A contour still closed after this clip
   is kept only if it encloses training data on the plane (`.encloses_data()`, the global search's
   empty-region rule; since 2026-11-06, §7). The actionability constraints (`set_filters()`), the
   training ranges and the model's re-score then run on the remaining vertices, as before.
3. **The selection** is unchanged: the nearest valid vertex per pair, and the pair with the smallest
   distance by the thesis rule.

The result stores `hull` (whether the rule applied) and `hull_polygon` (the winning plane's hull
vertices; since 2026-11-06 also when the rule is off, built after the search for the plot). When no
pair has a valid point, the message names the hull and says that `set_filters(..., hull = FALSE)`
drops the rule.

**The plot.** `plot.bl_local_result()` colours the prediction grid only inside the winning plane's
hull when the search applied it (`trim_to_hull`, 2026-11-06), as the global biplot colours only inside
its hull (`bl_build_grid(calc_hull = TRUE)`): the coloured area is where the search could choose
boundary points, and contour lines in the white area lie outside the hull. `trim_to_hull = TRUE`
shows the hull of a search with the rule off; the Pima Figure 7(a) uses it to show the TP target's
run 4, which has no counterfactual inside the hull, searched without it: its counterfactual is the
point the hull removes.

**Design choices** (Adriaan Rowan, 2026-11-05; the first and last rows revised on 2026-11-06, §7):

| Choice | Decision | Reason |
|---|---|---|
| The points the hull is built from | The training data (since 2026-11-06; before, the training data and the target, as Appendix A) | The valid region is set by the training data alone, as the text says ("the convex hull of the training data", p.117) and as in the global search. A target outside the training data no longer enlarges the region its counterfactual may come from. |
| A switch | `hull = TRUE` by default; `hull = FALSE` gives the earlier search exactly. Since 2026-11-06 the switch is in `set_filters()` (§6) | The thesis rule by default, with the earlier results reproducible (the methodology note's measurement M6) |
| The rotation scripts (`17_1` to `17_3`) | Pin `set_filters(tgt, hull = FALSE)` | Their measurements concern the rotation and the pair search, so the rotation note's numbers stay valid |
| An enclosure rule for closed contours | Since 2026-11-06, a closed contour is kept only if it encloses training data on the plane (with the hull on); before, none | The same empty-region rule as the global search ("any local decision boundary area that does not contain data points is discarded", §3.10 step 4, p.93), tested on the training data in both searches |

**Off the plane.** The hull applies to the boundary points on the plane. The restored point of a
`"fixed"` feature, the rounded point of an `"integer"` feature and the sparse counterfactual lie off
the plane, and are checked against the constraints and the training ranges only, as before.

## 4. The effect

**The Pima example of Paper 1** (four-variable GAM, CVA; the methodology note's M6, from
`scripts/paper1-pima/04_method-change-effects.R`):

- Over the 79 test rows (run 1), every row keeps a counterfactual: 68 are unchanged and 11 move, with
  a median search distance 4.5% longer among the moved rows.
- The TN target is unchanged in every run.
- The TP target moves from pair (2,3) to pair (3,4), and its search distance from 1.50 to 1.84. The
  counterfactual without the hull raised *Age* from 51 to 59.6 and *Pregnancies* from 7 to 10, a point
  outside the training data on its plane. Inside the hull, every TP counterfactual lowers *Age*, to
  about 41, so run 4 (*Pregnancies* fixed, *Age* unable to decrease) now has no counterfactual; run 3
  (both fixed) had none before either. Run 5 (pair (1,2) only) has a boundary point inside the hull,
  but no whole-number *Pregnancies* value still crosses the boundary.

The TP target shows the trade-off between the constraints: the plausible counterfactuals need a
younger patient, which is not actionable, so no counterfactual meets both rules.

**The loan harness** (`development/03_loan_status_Boundary_Logic.R`) is unchanged: every number its
walkthrough quotes is reproduced with the hull on.

## 5. Reproduction

- `hull = FALSE` reproduces the Pima results of commit `12b38e4` exactly (checked 2026-11-05 for
  every run of both targets).
- The measurement M6 in `04_method-change-effects.R` records both settings for every run and target,
  and the sweep over the 79 test rows.

## 6. The plausibility and feasibility switches in `set_filters()` (2026-11-06)

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
acceptable (p.115). Holding all three in
`set_filters()` follows that framing, gives the training ranges a switch they lacked (before, the only
way off was to set `bl_result$train_ranges` to `NULL`), and makes the unconstrained counterfactual one
call: `set_filters(tgt, hull = FALSE, train_ranges = FALSE)` with no per-variable constraint. The
argument `bl_find_local_cf(hull = )` of 2026-11-05 was removed (version 0.2.0 was not yet tagged).

**Why `set_filters = NULL` keeps both rules on.** The switches live in the `bl_filters` object, and a
`NULL` has none, so the code must give `NULL` a meaning. The private `.filter_switches()` reads it as
both on, for four reasons:

1. **No result changes.** Before the move, a call without filters already applied both rules (the
   hull by its own default, the ranges always). Reading `NULL` as both off would silently change every
   such call wherever the nearest boundary point lies outside the plane's hull or outside a training
   range: the loan harness's runs without filters (Steps 15 and 16, quoted in the walkthrough), the
   vignette's and README's equivalents, the Pima runs without filters and the approximation-gap
   experiment (`05_full-space-gap.R`, whose ratio is measured with the hull on and whose check that
   every search point lies inside its plane's hull would fail).
2. **The size of that change.** The hull alone moves 11 of the 79 Pima test targets (§4); the TP
   target would go from pair (3,4) at 1.84 to pair (2,3) at 1.50, raising *Age* to 59.6, outside the
   training data. Dropping the ranges as well can move more.
3. **Consistency.** `set_filters(tgt)` with no constraint fills in both switches as `TRUE`. If `NULL`
   meant both off, "no filter object" and "an empty filter object" would give different
   counterfactuals for the same target.
4. **The thesis.** Its counterfactual is the constrained one (p.117); the unconstrained one is the
   stated exception (p.115), so it is the explicit opt-out.

A `bl_filters` object built by hand or before the switches existed (no fields) is read the same way.

**Checks (2026-11-06).** With the defaults, the loan harness, the rotation scripts `17_1` to `17_3`
and the Pima pipeline give the same results as before the move.

**What the switches change on the loan harness** (a driver run after the harness, re-run on
2026-11-06 after the changes of §7; the searches without actionability constraints):

| Target | Default (both on) | `hull = FALSE` | `train_ranges = FALSE` | Both off |
|---|---|---|---|---|
| Step 10 test row | pair (3,4), √d²_M 0.5554 | the same | the same | the same |
| Step 16 external applicant | pair (1,2), √d²_M 2.8241; credit score 534.92 | pair (1,2), 1.7026; credit score 489.44 | the same as the default | pair (1,4), 0.4750; credit score 421.26 |

The external applicant's credit score (420) lies below the training range (489 to 754), and the
applicant lies outside the training data on the chosen plane. With the hull on, the counterfactual must lie
inside the training data's hull, which already keeps the credit score within its range, so dropping
the ranges changes nothing. With the hull off and the ranges on, the nearest counterfactual sits at the
range's lower end; with both off it leaves the credit score near 420, outside the training data, at a
sixth of the default's distance. (Before §7 the hull included the applicant, and the default equalled
the `hull = FALSE` column.)
