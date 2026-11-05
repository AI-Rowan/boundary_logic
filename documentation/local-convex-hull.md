# The convex hull on each rotated plane

Date created: 2026-11-05
Date last edited: 2026-11-05

This note explains why `bl_find_local_cf()` keeps only the boundary points inside the convex hull of
the training data on each rotated plane (`hull = TRUE`, the default since 2026-11-05), how it does so,
and what the rule changes on the Pima example of Paper 1 and on the loan harness. The decision is
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

1. **The hull.** The training data and the target are projected onto the rotated plane
   (`Z_train_rot`, `Z_target`), and their convex hull is built with `.build_hull_polygon(...,
   outlie = 1)` (`R/hull_utils.R`). The fraction 1 gives the full hull, as the global grid's default
   (`bl_build_grid(outlie = 1)`).
2. **The filter.** Each boundary contour on the plane keeps only its vertices inside the hull
   (`.poly_clip()`, the vertex filter the global search uses). The actionability constraints
   (`set_filters()`), the training ranges and the model's re-score then run on the remaining vertices,
   as before.
3. **The selection** is unchanged: the nearest valid vertex per pair, and the pair with the smallest
   distance by the thesis rule.

The result stores `hull` (the argument) and `hull_polygon` (the winning plane's hull vertices). When no
pair has a valid point, the message names the hull and says that `hull = FALSE` drops the rule.

**Design choices** (Adriaan Rowan, 2026-11-05):

| Choice | Decision | Reason |
|---|---|---|
| The points the hull is built from | The training data and the target | As Appendix A. The target always lies inside its own hull, so a target at the edge of the training cloud can still reach a counterfactual near it. |
| A switch | `hull = TRUE` by default; `hull = FALSE` gives the earlier search exactly | The thesis rule by default, with the earlier results reproducible (the methodology note's measurement M6) |
| The rotation scripts (`17_1` to `17_3`) | Pin `hull = FALSE` | Their measurements concern the rotation and the pair search, so the rotation note's numbers stay valid |
| An enclosure rule for closed contours | None | The text requires only the hull for the local search (p.117) |

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
