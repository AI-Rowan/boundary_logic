# Methodology changes since the thesis

*The boundarylogic package against the PhD method: the changes made and the future work*

Adriaan Rowan · 6 November 2026

**The record.** This note records the changes from the thesis to `boundarylogic` 0.2.0 (tag `v0.2.0`), as reviewed by the supervisor. Its values were measured at commit `e41f22f`, whose R code equals `v0.2.0`. They come from `method-change-effects.csv` and `findings-effects.csv` in this folder and from the Pima example's `output/full_space_gap.rds`. The note is plain Markdown, edited by hand; it is the one place for the package's history, with `NEWS.md`.

## Contents

- [Summary](#summary)
- [1 How to read this note](#1-how-to-read-this-note)
- [2 Changes made: history](#2-changes-made-history)
- [3 Changes made: core methodology](#3-changes-made-core-methodology)
  - [3.1 Variable importance standardised by the within-class standard deviation (CVA)](#31-variable-importance-standardised-by-the-within-class-standard-deviation-cva)
  - [3.2 The metric W as a covariance matrix](#32-the-metric-w-as-a-covariance-matrix)
  - [3.3 Global counterfactuals outside the training ranges](#33-global-counterfactuals-outside-the-training-ranges)
  - [3.4 A cap on the number of pairs of components](#34-a-cap-on-the-number-of-pairs-of-components)
  - [3.5 The distance that compares the pairs (restored under PCA)](#35-the-distance-that-compares-the-pairs-restored-under-pca)
  - [3.6 The convex hull in the local search (restored)](#36-the-convex-hull-in-the-local-search-restored)
  - [3.7 Fixed variables: a search tolerance and a strict restore](#37-fixed-variables-a-search-tolerance-and-a-strict-restore)
  - [3.8 Whole-number variables: the integer constraint](#38-whole-number-variables-the-integer-constraint)
  - [3.9 A “No effect” class for zero Shapley values](#39-a-no-effect-class-for-zero-shapley-values)
  - [3.10 The sparse counterfactual: a greedy rule with a stopping point](#310-the-sparse-counterfactual-a-greedy-rule-with-a-stopping-point)
  - [3.11 Outlier removal before the projection](#311-outlier-removal-before-the-projection)
  - [3.12 Three properties the thesis does not state](#312-three-properties-the-thesis-does-not-state)
- [4 Changes made: code improvements](#4-changes-made-code-improvements)
- [5 Future work](#5-future-work)
- [Reproduction](#reproduction)

# Summary

This note records the changes made to the method of the thesis (Rowan,
2026) in the `boundarylogic` package since the thesis was submitted
(Sections 2 to 4), and the future work that remains (Section 5). It
separates the **core methodology**, which determines what the method
computes, from the **code improvements** used to implement it, and
measures the effect of each change on the Pima diabetes example of Paper
1 wherever the package can switch the change back.

Eleven core changes are recorded in Section 3. One of them departs from
the thesis text itself, not only from its code: **global counterfactuals
outside the training ranges are kept and reported** rather than
discarded (Section 3.3), whereas the thesis filters them (p.93). The
package keeps them so that the distance to the nearest boundary stays
correct, at the cost of reporting some infeasible counterfactuals.

Two other departures were removed on 5 November 2026:

- the comparison of the pairs of components by Mahalanobis rather than
  Euclidean distance under PCA: the package now applies the thesis’s
  rule by projection (p.86), in the pair comparison and in the variable
  importance (Sections 3.1 and 3.5). Paper 1, which uses CVA, is
  unaffected;
- **the local search without a convex hull**: the package again keeps
  only the boundary points inside the convex hull of the training data
  on each rotated plane, as the thesis requires (p.117; Section 3.6).
  This changes the Phase 3 results of the true-positive target.

On 6 November 2026 two rules were aligned with the text and between the
two searches: the hull on each rotated plane is built from the training
data alone, and a closed boundary contour is kept only if it encloses
training data, in the global search (where the rows being searched
decided before) and in the local search (which had no such rule). The
Paper 1 results are unchanged (Sections 3.3 and 3.6).

The remaining core changes extend the thesis method (the `"fixed"` and
`"integer"` constraints, the greedy sparse counterfactual, the “No
effect” class) or change the scale or reach of a measure (the
within-class standard deviation in the variable importance, the
covariance form of $\mathbf{W}$, the cap on the number of pairs).
Section 3.12 adds three properties of the thesis method that the thesis
does not state, all of which the paper relies on, including the measured
gap between the distances on the biplot planes and those in the full
space. The code improvements are summarised in Section 4, and Section 5
lists the future work, each item of which Paper 1 names as an area of
future improvement.

# 1 How to read this note

**The three baselines.** The *thesis version* of the method is defined
by the thesis text (Rowan, 2026, Stellenbosch University, submitted 27
February 2026). Two versions of the thesis code exist and are reported
where they differ from the text or from each other:

- **Appendix A** of the thesis (pp.158–218), which prints the code used
  for the thesis;
- **the branch `original_PhD_code`** of the package repository,
  committed on 9 March 2026, after submission, as the starting code for
  the package. Its `*_optimized.R` scripts are not always the same as
  Appendix A.

Where the text, Appendix A and the branch disagree, the text defines the
thesis version. Page numbers are the printed pages of the final thesis,
as in the package’s other notes.

**The Pima set-up** is that of Paper 1. The 392 rows without zeros in
*Insulin*, *SkinThickness*, *Glucose* or *BMI* are split 80:20 (seed
121) into 313 training and 79 testing rows. A generalised additive model
(GAM) with a smooth term per feature is fitted to the training data, and
the CVA biplot is built on the four prediction outcomes (TP, TN, FP,
FN). Phases 1 and 2 use the seven features; the global variable
importance then drops *DiabetesPedigreeFunction*, *SkinThickness* and
*BMI*, and Phase 3 uses the GAM refitted on *Pregnancies*, *Glucose*,
*Insulin* and *Age*. The local examples use a true positive (test row
23) and a true negative (test row 25), each the test row whose global
counterfactual distance lies nearest its group’s median.

**Measurements.** Every measured number in this note is read from
`method-change-effects.csv`, written by the script
`04_method-change-effects.R` at package commit e41f22f, except the
approximation gap of Section 3.12, which `05_full-space-gap.R` measures,
and the effects of the changes of 6 November 2026, which
`06_findings-effects.R` measures against the package before them
(`findings-effects.csv`). A *toggle* measurement switches one change
back through a package argument and reruns the search; a *stored*
measurement reads a diagnostic the package keeps (for example the search
point before a constraint is applied). Changes that cannot be switched
back are described and marked “not measured”. Distances are quoted as
$\sqrt{d^2_M}$, the Mahalanobis distance under the package’s
$\mathbf{W}$ (Section 3.2).

**Proposed paper sentences.** Where a change affects how Paper 1 reports
a Pima result, the section ends with the sentence the paper proposes to
use.

# 2 Changes made: history

The package was started on 15 January 2026, while the thesis was being
finalised, and was built from the thesis code over March 2026. Table 1
lists the dated changes that bear on the method. The category column
says whether a change belongs to the core methodology (Section 3), is a
code improvement (Section 4) or is a bug fix.

| Date | Source | What changed | Category |
|:---|:---|:---|:---|
| 27 Feb 2026 | Thesis text and Appendix A | The baseline |  |
| 9 Mar 2026 | Branch `original_PhD_code` | Optimised thesis scripts committed as the package’s starting code; the boundary search takes a random tie-break (`max.col()`) where Appendix A uses `which.min` | Code |
| 10 Mar 2026 | `b8d4965` | Phase 1 port; outlier removal before the projection, default fraction 0.9 | Core (3.11) |
| 13 Mar 2026 | `e42fa03` | Phase 2 port; training ranges applied before the nearest boundary point is chosen | Core (3.3) |
| 16 Mar 2026 | `4058730` | Phase 3 port: pair search capped at 10 pairs; `.bl_rotate()` returned the wrong columns for pairs other than (1,2); no hull in the local search | Core (3.4, 3.6); fix |
| 16 Mar 2026 | `f030a64` | Global grid hull fraction set to 1 (thesis code 0.9) | Code |
| 17 Mar 2026 | `bb96efe` | The `"fixed"` constraint, with a 0.5 tolerance | Core (3.7) |
| 20 Mar 2026 | `fe47d75`, `dffbeb8` | Predictions floored to three decimals (thesis code two) | Code |
| 24 Mar 2026 | `ea6a12f` | Training ranges checked after the nearest point is chosen; points outside are kept and reported (prompted by a PCA biplot, Section 3.3) | Core (3.3) |
| 23 May 2026 | `271b71e` | $\mathbf{W}$ divided by $(n-g)$; within-class standard deviation in the variable importance; Mahalanobis pair selection by default; contour band width 0.001; `.bl_rotate()` fixed | Core (3.1, 3.2, 3.5); code; fix |
| 25 May 2026 | `d509201` | Outlier removal moved into `bl_prepare_data()` | Code |
| 5 Jun 2026 | `1719c73` | Grid probabilities clamped to \[0, 1\] | Code |
| 29–30 Oct 2026 | `f2a076b`, `7aa1344`, `f4f3647` | Paper 1 Tasks A and B; the `"fixed"` restore; the “No effect” class; the greedy sparse counterfactual; Gate 2 code fixes | Core (3.7, 3.9, 3.10); fix |
| 31 Oct 2026 | `f0e3dd4` | The nearest boundary vertex taken deterministically | Code |
| 31 Oct 2026 | `49c9067` | Pair cap raised to 21; the Phase 2 Mahalanobis check | Core (3.4); code |
| 4 Nov 2026 | `67a1424` | Strict `"fixed"`; the `"integer"` constraint; Mahalanobis distance returned by the global search | Core (3.7, 3.8); code |
| 5 Nov 2026 | `12b38e4` | The thesis’s distance rule by projection (`distance = "auto"`): Euclidean under PCA, Mahalanobis under CVA, in the pair comparison and the variable importance | Core (3.1, 3.5) |
| 5 Nov 2026 | `c13fe44` | The convex hull of the training data and the target rebuilt on each rotated plane in the local search (`hull = TRUE`), as the thesis requires; outlier removal off by default (`hull_fraction = 1`); one colour scheme for every biplot; a warning for a target at the centre of the biplot | Core (3.6, 3.11); code |
| 6 Nov 2026 | `9bc8e1c` | The hull and training-range switches held in `set_filters()` with the actionability constraints (`hull = TRUE`, `train_ranges = TRUE`, both on by default); the training ranges gain a switch; no result changes | Code |
| 6 Nov 2026 | after `ef014ca` | The findings of the package’s algorithms document: the empty-region rule tested on the training data in the global search and added to the local search; the hull on each rotated plane built from the training data alone; the local grid clamped to \[0, 1\]; the local plot widened to show a target outside the training data | Core (3.3, 3.6); code |

Table 1. Dated changes to the method since the thesis.

# 3 Changes made: core methodology

Table 2 follows the method through its three phases and marks each step
as unchanged or changed. The steps marked as changed are discussed in
Sections 3.1 to 3.11; the code-level differences behind the unchanged
steps are listed in Section 4.

| Step | Thesis (text) | Package | Status |
|:---|:---|:---|:---|
| **Phase 1** |  |  |  |
| Outlier removal | Optional; none for Pima (p.99), a 90% hull for the loans (p.127) | Optional, off by default: `bl_prepare_data(hull_fraction =)`, on a stand-in projection | Changed (3.11) |
| Projection | PCA or CVA; CVA on the four prediction outcomes (p.92) | `bl_build_projection()` | Unchanged |
| Metric $\mathbf{W}$ | Within-groups SSP (§3.4, p.48) | Within-groups covariance, SSP$/(n-g)$ | Changed (3.2) |
| Decision boundary | Contours at $b \pm \varepsilon$ on a prediction grid (§3.6.1, p.62; $\varepsilon = 0.001$, p.78) | `bl_build_grid(b_margin = 0.001)` | Unchanged |
| Plausibility on the plane | Grid within the convex hull of the training data (p.77) | Hull polygon of the training data | Unchanged |
| **Phase 2** |  |  |  |
| Global counterfactual | Nearest boundary point on the plane, by Euclidean distance (p.75, p.101) | `bl_find_boundary()` | Unchanged |
| Feasibility of the global counterfactual | Training bounds, convex hull, no empty regions (§3.10 step 4, p.93) | Hull and empty-region rules (empty regions judged by the training data); training bounds reported, not applied | Changed (3.3) |
| Re-scoring by the model | Appendix A (pp.183–184) | The same rule | Unchanged |
| Variable importance $T_j$ | Standardised by the overall standard deviation (§4.2.2, p.103) | Within-class standard deviation under CVA; overall standard deviation under PCA | Changed under CVA (3.1) |
| Robustness $R = \sum_j T_j$ | §4.2.3, p.105 | The same sum | Follows 3.1 |
| **Phase 3** |  |  |  |
| Procrustes rotation through the target | First two components of the rotated matrix (§3.5, pp.52–54; p.115) | `.bl_rotate()` | Unchanged (see 3.12) |
| Pairs of components | All $p(p-1)/2$ pairs (p.117) | At most `max_pairs = 21` | Changed (3.4) |
| Nearest point on one plane | Euclidean distance on the plane (p.75) | The same | Unchanged |
| Comparison across pairs | Mahalanobis (CVA) or Euclidean (PCA) in the full space (p.86, p.117) | The same, by `distance = "auto"` (since 5 Nov 2026) | Restored (3.5) |
| Plausibility of the local counterfactual | Within the convex hull of the training data (p.117) | The convex hull of the training data, rebuilt on each rotated plane, with the empty-region rule (`set_filters(hull = TRUE)`) | Restored (3.6) |
| Actionability | Boundary points filtered by the constraints (p.88) | `set_filters()`, with a `"fixed"` restore | Changed (3.7) |
| Whole-number variables | Treated as continuous (p.87, p.119) | `"integer"` constraint | Changed (3.8) |
| Shapley values | All orderings of the change (§4.3.2, p.121) | Exact computation, the same values | Unchanged |
| Supports and Contradicts | By the sign of the Shapley value (pp.119–120) | Adds “No effect” | Changed (3.9) |
| Sparse counterfactual | Supporting variables changed, contradicting kept (p.120, p.122) | Greedy, with a stopping rule | Changed (3.10) |

Table 2. The thesis method step by step against the package.

## 3.1 Variable importance standardised by the within-class standard deviation (CVA)

**Thesis.** The global variable importance decomposes each data point’s
change vector to its counterfactual along the variable axes. Each
component “is standardised by dividing it by the overall standard
deviation of that variable” (§4.2.2, p.103), and the importance of
variable $j$ is the sum of the absolute standardised components over the
data points, $T_j$. The robustness of the model is $R = \sum_j T_j$
(§4.2.3, p.105). Appendix A and the branch use the total standard
deviation of the training data.

**Package.** Under CVA, `plot.bl_boundary()` and `bl_robustness()`
divide each component by the within-class standard deviation,
$\sqrt{W_{jj}}$, by default; `distance = "euclidean"` restores the
thesis version. Under PCA the default has been the thesis’s overall
standard deviation since 5 November 2026 (`distance = "auto"`); before,
the within-class standard deviation was the default for both
projections.

**Reason.** The total standard deviation includes the spread between the
groups. A variable that separates the groups well therefore has a large
total standard deviation, which shrinks its standardised distance. This
under-weights the variables that matter most for the classification. The
within-class standard deviation measures each component in the same
units as the CVA metric.

**Effect on Pima.** The ranking is unchanged on both models, and the
thesis rule would have pruned the same three variables (*BMI,
DiabetesPedigreeFunction, SkinThickness*). The totals are larger under
the package’s rule, with the robustness rising from 528.82 to 602.37 on
the seven-feature model.

| Variable | Total SD | Rank (total SD) | Within-class SD | Rank (within-class SD) |
|:---|---:|---:|---:|---:|
| Glucose | 139.6 | 1 | 181.2 | 1 |
| Age | 107.2 | 2 | 123.9 | 2 |
| Insulin | 78.4 | 3 | 83.6 | 3 |
| Pregnancies | 75.8 | 4 | 80.8 | 4 |
| DiabetesPedigreeFunction | 47.5 | 5 | 48.7 | 5 |
| SkinThickness | 41.8 | 6 | 44.1 | 6 |
| BMI | 38.5 | 7 | 40.0 | 7 |

Table 3. Variable importance totals on the seven-feature model (training
rows).

**The Mahalanobis check (Pima issue 9).** The per-variable measure
ignores the correlations in $\mathbf{W}$. Summed over the training rows,
the squared per-variable measure exceeds the full squared Mahalanobis
distance by 27.4% on the seven-feature model and by 35.2% on the
four-feature model. *Proposed paper sentence:* “The per-variable
standardisation ignores the correlations between the variables; on the
Pima data the summed squared per-variable distances exceed the full
squared Mahalanobis distance by 27.4% (seven variables) and 35.2% (four
variables), so the importance ranks the variables but does not decompose
the distance.” A full decomposition is listed in Section 5.

**The importance against the global distances.** The change vector of
the thesis lies on the plane,
$\boldsymbol{\delta}_i = (\mathbf{z}_{cf}-\mathbf{z}_i)\mathbf{M}_2^{-}$,
and the package keeps it. The distances `bl_find_boundary()` reports
(`dist_x`, `dist_mahalanobis`) measure instead from the observation
itself to its global counterfactual. Since the counterfactual is a point
of the plane (Section 3.3), that gap also resets the part of the
observation that lies off the plane, its other $p-2$ canonical
coordinates, to the mean. That part has nothing to do with the boundary
and is most of the gap: on the training rows it is 78.1% of the summed
squared global gap on the seven-feature model and 62.6% on the
four-feature model (within-class standard deviations). The importance
and the robustness therefore measure the in-plane change, as the thesis
does, and the global distances are reported only; Paper 1 uses the
global $\sqrt{d^2_M}$ only to choose its two targets.

## 3.2 The metric W as a covariance matrix

**Thesis.** CVA partitions the total sums of squares and products (SSP)
into the between-groups and within-groups SSP, and $\mathbf{W}$ is the
within-groups SSP (§3.4, p.48; p.118), computed on the four prediction
outcomes (p.92). The text also refers to “the common within-group
covariance matrix” (p.84). Appendix A and the branch take biplotEZ’s
`Wmat`, an SSP with no divisor, invert it with `solve()`, and use it
only to compare the pairs of components.

**Package.** `bl_build_projection()` computes $\mathbf{W}$ as the pooled
within-groups covariance, the SSP divided by $(n-g)$, inverts it with a
Cholesky decomposition and stores both. Under CVA it is used by default
in the pair search (Section 3.5) and the variable importance (Section
3.1), and in the Mahalanobis distances the package reports. Under PCA,
$\mathbf{W}$ is computed on the two classes and is used only when
`distance = "mahalanobis"` is passed.

**Reason.** As a covariance, $\mathbf{W}$ measures distances in
within-class standard deviations, so $\sqrt{d^2_M}$ can be read as a
number of standard deviations, and $\sqrt{W_{jj}}$ is a standard
deviation that Section 3.1 can use. This also follows the covariance
reading on p.84.

**Effect on Pima.** The change is a fixed rescaling. With $n = 313$
training rows and $g = 4$ groups, every squared distance is
$n - g = 309$ times the thesis value, and every distance is
$\sqrt{309} = 17.58$ times larger. No ranking or selection changes. The
distances quoted in Paper 1 are on the package’s scale.

## 3.3 Global counterfactuals outside the training ranges

**Thesis.** The local decision boundaries “are filtered to only contain
feasible counterfactuals within the upper and lower bounds for each
feature of the training data”, only boundaries within the convex hull of
the training data are considered, and “any local decision boundary area
that does not contain data points is discarded” (§3.10 step 4, p.93).
Page 101 is less strict, keeping regions “that contain training data or
fall within the convex hull”. The thesis code does not apply the ranges
in the global search: the branch’s driver scripts pass no training
ranges to it, and instead remove the test rows that lie outside the
training ranges before the search.

**Package.** `bl_find_boundary()` keeps the hull and the empty-region
rules: contours are clipped to the training hull, and a closed contour
is kept only if it encloses training rows. Until 6 November 2026 the
rows being searched decided, as in the branch, so one row’s
counterfactual could depend on which other rows were searched with it;
the training data now decide, so it does not. On Pima no contour or
counterfactual changes (6 contours on the four-feature model before and
after); on the package’s loan example the rule keeps 18 contours instead
of 12 on the four-feature model and moves 346 of its 4,429 global
counterfactuals. The package first followed the thesis text and filtered
the boundary points by the training ranges before choosing the nearest
one (13 March 2026, `e42fa03`). Since 24 March 2026 (`ea6a12f`) the
ranges are checked after the nearest point is chosen: points outside
them are kept, and only their number is reported in a message, with no
flag on the individual rows. The test rows are also left in place, where
the branch removed them.

**Reason.** The change was prompted by a PCA biplot of the iris data. A
global counterfactual is a point on the biplot plane, mapped back to the
$p$ features through the two plotted components only:

$$\hat{\mathbf{x}}_{cf} = \bar{\mathbf{x}} + \mathbf{z}_{cf}\mathbf{M}_2^{-},$$

where $\mathbf{z}_{cf}$ is the boundary point on the plane,
$\mathbf{M}_2^{-}: 2 \times p$ holds the rows of the inverse projection
matrix for the two plotted components, and $\bar{\mathbf{x}}$ is the
vector of training means (on the standardised scale for a standardised
PCA biplot). In other words, every component that is not plotted is set
to its mean. The individual feature values of this reconstruction can
therefore fall outside the training ranges, even for a boundary point
well inside the convex hull on the plane. Filtering by the ranges before
the selection then discarded the nearest boundary point and selected a
farther one. This changed the distance that the variable importance and
the robustness aggregate, which are meant to measure the distance to the
nearest boundary. The package therefore keeps the nearest point and
reports its feasibility.

The commit records that CVA was unaffected. The same reconstruction
applies under CVA, however, and the Pima results below show that
counterfactuals on a CVA biplot can also fall outside the training
ranges.

**The trade-off.** The two rules trade correct distances against
feasibility. Keeping the nearest point gives the correct distance to the
approximated boundary, which the global measures need, but may report a
counterfactual with feature values outside the training data. The thesis
rule reports only feasible counterfactuals, but measures the distance to
a farther point, and a data point may be left without one. A middle
course keeps the nearest point for the global measures and flags each
infeasible counterfactual by row, which the package does not yet do; it
and the thesis rule as an option are listed in Section 5.

**Effect on Pima.** Not switchable, so the effect on the importance is
not measured. On the seven-feature model no global counterfactual lies
outside the training ranges (0 of 313 training rows and 0 of 79 test
rows), so the pruning is unaffected. On the four-feature model, also a
CVA biplot, 10 of 313 training rows and 5 of 79 test rows have a
counterfactual outside them; the thesis rule would have replaced these
with farther points. The global counterfactuals of both Paper 1 targets
lie within the ranges, so their own global distances are the same under
either rule; the group medians used to choose them could still move,
since the distances of other rows would change. The branch would also
have removed 6 test rows on the seven-feature model before the search.
The local counterfactuals of Paper 1 are unaffected, since the local
search applies the training ranges (Section 3.6).

## 3.4 A cap on the number of pairs of components

**Thesis.** “the projection matrix is recalculated for all possible
two-dimensional biplots formed from pairs of components”, that is
$p!/(2!(p-2)!)$ pairs (§4.3.1, p.117). The pair search matters in
practice: the seven-variable Pima counterfactual of the thesis “was
created from the biplot generated by components numbers 1 and 7” (p.118,
Figure 4.11), and the true-positive loan example used the second and
fourth components (§5.3.1, p.137).

**Package.** `bl_find_local_cf()` searches the pairs in the order (1,2),
(1,3), (2,3), (1,4), …, up to `max_pairs`. The default was 10 pairs from
March to October 2026 and is 21 pairs since 31 October 2026. Twenty-one
pairs cover every pair for $p \le 7$; for eight or more variables the
pairs with the highest components are left out.

**Reason.** Each pair requires a rotation, a prediction grid and a
boundary search, so the run time grows with $p^2$. The earlier cap of 10
would have missed the thesis’s own Pima counterfactual: pair (1,7) is
number 16 in the package’s order. This is the reason the cap was raised.

**Effect on Pima.** The four-feature model has 6 pairs, so a cap of 10
or 21 searches all of them, and the results are identical. The gain of
the pair search is measured instead against pair (1,2) alone
(`max_pairs = 1`). Over all pairs, the search distance falls from 1.93
to 1.84 for the true positive and from 0.17 to 0.14 for the true
negative. After the `"integer"` rounding of Section 3.8, the true
positive has no counterfactual on pair (1,2) alone (no whole number of
pregnancies near its search point crosses the boundary), against a
reported distance of 1.82 over all pairs; for the true negative the
reported distances are 0.44 against 0.42. The winning pairs are (3,4)
and (2,4). With more variables the number of pairs grows quickly (21 for
seven variables, 28 for eight), so a larger model may need a higher cap,
or no cap, to find the nearest counterfactual. Removing the default cap
would return the package to the thesis version.

## 3.5 The distance that compares the pairs (restored under PCA)

**Thesis.** The counterfactual is “the shortest change vector measured
in the full dimensional space using Euclidean distance for PCA biplots
and Mahalanobis distance for CVA biplots” (§3.9.1, p.86), and the CVA
rule is $\min \, \mathbf{d}'\mathbf{W}^{-1}\mathbf{d}$ (p.117). The text
also notes that Euclidean distances on a CVA biplot correspond to
Mahalanobis distances in the feature space (p.91). Appendix A and the
branch use the Mahalanobis rule only, and only under CVA.

**Package.** Since 5 November 2026 `bl_find_local_cf()` follows the
thesis rule by default (`distance = "auto"`): Mahalanobis distance under
CVA, and under PCA the Euclidean distance in the feature space with each
feature divided by its total standard deviation. Before, Mahalanobis
distance was the default for both projections, and
`distance = "euclidean"` compared the pairs by the distance on the
rotated plane.

**Reason.** The default now follows the thesis rule for each projection.
Under CVA the package’s Mahalanobis rule agrees with the distance on the
plane: its $\mathbf{W}$ makes the rotated loadings satisfy
$\mathbf{V}'\mathbf{W}\mathbf{V} = \mathbf{I}/(n-g)$, so for a
counterfactual on the plane:

$$d^2_M = (n-g)\, d^2_z,$$

where $d^2_M$ is the squared Mahalanobis distance in the feature space
and $d_z$ is the distance on the rotated plane. In other words, the
Mahalanobis rule ranks the pairs as the plane distance does. This
identity is the package’s derivation, not the thesis’s: with the SSP
form of $\mathbf{W}$ in the thesis the factor is 1. Under a standardised
PCA the rotated loadings are orthonormal, so the Euclidean distance in
the feature space equals the plane distance at the search point;
measuring it in the feature space also covers the whole-number
candidates of Section 3.8, which lie off the plane.

**Effect on Pima.** Paper 1 uses CVA, so its results are unchanged: the
ratio $d^2_M / d^2_z$ equals $n - g = 309$ for every pair and both
targets, to a relative error below 6e-15, and the Mahalanobis choice is
the pair with the smallest plane distance for both targets, with and
without the `"integer"` constraint (pairs (3,4) and (2,4)).

On a PCA biplot of the same four-feature model, outside Paper 1’s
set-up, the two rules differ. For the true positive the thesis rule
selects pair (3,4) (distance 2.05), while the earlier Mahalanobis
default finds no whole-number counterfactual; for the true negative the
rules select (1,4) and (3,4). Over the 79 test rows the two rules choose
different pairs for 51. On the seven-feature PCA biplot the two
standardisations of Section 3.1 rank the variables differently but
identify the same three weakest.

## 3.6 The convex hull in the local search (restored)

**Thesis.** “The counterfactual must be selected from a local decision
boundary that is within the convex hull of the training data” (§4.3.1,
p.117), which ensures plausibility (p.116). The hull is dropped only
when an unconstrained solution is acceptable (p.115). Appendix A
rebuilds the hull on each rotated plane from the target and the training
data (pp.207, 209); the branch builds it from the target and the test
data, and keeps a closed contour only if it encloses the target.

**Package.** From March to November 2026 `bl_find_local_cf()` applied no
hull: the candidate boundary points were filtered by the actionability
constraints and the training ranges only. The recorded reason, that the
rotation invalidates the hull of the global biplot, is correct for that
hull, but it does not cover a hull rebuilt on each rotated plane. Since
5 November 2026 a convex hull is rebuilt on each rotated plane and only
the boundary points inside it are kept, before the other filters
(`set_filters(hull = TRUE)`, the default; `hull = FALSE` gives the
earlier search). From 5 to 6 November the hull was built from the
training data and the target, as in Appendix A; since 6 November 2026 it
is built from the training data alone, as the text says, and a contour
that is still closed after the clip is kept only if it encloses training
data on the plane, the empty-region rule of the global search (Section
3.3; the thesis states it for the decision boundaries in general, §3.10
step 4, p.93). Since 6 November 2026 `set_filters()` holds the switches
for all three of the thesis’s rules (actionability, feasibility and
plausibility, p.116): `train_ranges = FALSE` drops the training ranges,
and both switches off give the unconstrained counterfactual (p.115). A
search without a filter object keeps both rules.

**Reason.** The text requires the hull for plausibility, and Appendix A
builds it on each rotated plane; the earlier reason did not cover such a
hull. The hull is the training data’s alone, so the region a
counterfactual may come from does not grow with a target that lies
outside the training data, and the empty-region rule is the same in both
searches, judged on the same data.

**Effect on Pima.** The true negative is unchanged in every run. The
true positive’s search moves from pair (2,3) to (3,4), and its search
distance from 1.50 to 1.84: without the hull its counterfactual raised
*Age* from 51 to 59.58, a point outside the training data on its plane,
while inside the hull it lowers *Age* to 41.20. Every true-positive
counterfactual inside the hull lowers *Age*, so the runs that keep *Age*
from falling (runs 2, 3 and 4) have no counterfactual; run 4 had one
without the hull. Over the 79 test rows (run 1), every row keeps a
counterfactual: 68 are unchanged and 11 move, with a median search
distance 4.5% longer among the moved rows. See the package’s
`documentation/local-convex-hull.md`.

*The changes of 6 November* (measured by `06_findings-effects.R` against
the package before them) leave every run of both targets unchanged. Over
the 79 test rows (run 1), every row keeps a counterfactual and 4 move
(test rows 16, 46, 57, 63), with a median distance 4.5% longer. On the
package’s loan example the empty-region rule moves the harness target’s
constrained counterfactual (an empty closed contour held its nearest
boundary point), and the training-only hull moves the counterfactual of
an applicant who lies outside the training data.

**A constraint that does not bind (Pima issue 5).** With *Pregnancies*
and *Age* not allowed to decrease (run 2), the true negative’s
counterfactual equals the unconstrained one, because its unconstrained
counterfactual already increases these variables. The true positive’s
unconstrained counterfactual lowers both, so the constraint leaves it
without a counterfactual. *Proposed paper sentence:* “A directional
constraint changes the counterfactual only when it opposes the
unconstrained move; for the true negative the unconstrained
counterfactual already increases *Pregnancies* and *Age*, so the
constraint does not bind.”

**Plausibility against actionability (Pima issue 11).** Inside the
convex hull, the true positive (prediction 0.963) reaches the boundary
by lowering every variable: *Pregnancies* from 7 to 6, *Glucose* from
181 to 148.91, *Insulin* from 192 to 106.10 and *Age* from 51 to 41.20.
Without the hull it raised *Age* to 59.58, where the GAM’s non-monotone
smooth of *Age* (effective degrees of freedom 7.2) predicts a lower risk
but only 19 of the 313 training patients are older than 50. *Proposed
paper sentence:* “Without the plausibility rule, the nearest
counterfactual for the true positive raised *Age* from 51 to 59.6,
outside the training data; inside the convex hull every counterfactual
lowers *Age*, which is not actionable, so the constrained runs report
none.”

## 3.7 Fixed variables: a search tolerance and a strict restore

**Thesis.** Two approaches are described (§3.9.2, p.88): building the
biplot from the actionable variables only, with the others fixed at the
target’s values, or building it on all the data and filtering the
boundary points by the constraints. The examples use the second
approach, and a variable that must not change is held in the sparse
counterfactual: “the biplot cannot identify counterfactual solutions
that are precisely whole numbers. This is required to ensure both
Pregnancies and Age remain unchanged” (p.119), and “The number of
Pregnancies can remain unchanged” (p.122).

**Package.** `set_filters(tgt, x = "fixed")` lets the search accept a
boundary point that moves the variable by up to 0.5 model units, since a
point on a 2D boundary rarely holds a variable exactly. After the
search, `.restore_fixed()` sets each fixed variable back to its observed
value and re-scores the point. If the restored point no longer crosses
the boundary, no counterfactual is reported (strict since 4 November
2026), and the message says that changing the constraints on the other
variables may make one possible.

**Reason.** A reported counterfactual must hold every fixed variable
exactly, which the thesis achieved only in the sparse counterfactual.

**Effect on Pima (issue 12).** With *Pregnancies* and *Age* fixed (run
3), and with *Pregnancies* fixed and *Age* not allowed to decrease (run
4), the true positive has no boundary point within the tolerance and
inside the hull on any of the 6 pairs (Section 3.6: every true-positive
counterfactual inside the hull lowers *Age*). For the true negative the
restored points score 0.475 (run 3) and 0.487 (run 4), which do not
cross, so neither run reports a counterfactual. With *Pregnancies* fixed
and *Glucose* not allowed to increase (run 7), the true-negative search
moves *Pregnancies* from 10 to 9.83; restoring it changes the prediction
from 0.501 to 0.512, which still crosses, so the counterfactual is kept.
The true positive’s run 7 restore fails (0.526). *Proposed paper
sentence:* “When variables are held fixed, the two-dimensional search
rarely finds a boundary point that also holds them; five of the six runs
with a fixed variable report no counterfactual.” A follow-up search that
holds the fixed variables exactly is listed in Section 5.

## 3.8 Whole-number variables: the integer constraint

**Thesis.** “Integer data is included as a continuous variable data
type” (p.87), and the biplot “cannot identify counterfactual solutions
that are precisely whole numbers” (p.119).

**Package.** `set_filters(tgt, x = "integer")` rounds the counterfactual
after the search. Each moved integer variable is tried at the whole
values just below and above (in raw units); the combinations that still
cross the boundary and meet every constraint and the training ranges are
kept, and the nearest by the search’s distance is taken. If none is
valid, no counterfactual is reported. This replaces an earlier rounding
of the sparse counterfactual (`round_to`), which rounded in model units
without a re-check.

**Reason.** A count such as *Pregnancies* can only change in whole
steps, so a fractional counterfactual is not actionable.

**Effect on Pima (issue 13).** For the true negative, the search point
has *Pregnancies* 10.21. The value 10 does not cross the boundary, so 11
is taken: the distance rises from 0.14 to 0.42 and the prediction from
0.500 to 0.548. For the true positive the rounding moves *Pregnancies*
from 6.41 to 6 and the distance falls from 1.84 to 1.82. *Proposed paper
sentence:* “A count variable can only change in whole steps, so the
nearest valid counterfactual may lie well past the boundary; for the
true negative, the whole pregnancy raises the distance from 0.14 to
0.42.” Letting the other variables compensate after the rounding is
listed in Section 5.

## 3.9 A “No effect” class for zero Shapley values

**Thesis.** A variable *supports* the counterfactual when its Shapley
value moves the prediction towards the boundary and *contradicts* it
otherwise (pp.119–120); variables with no impact can be excluded from
the explanation (p.92). Appendix A counts a value of exactly zero as
Supports for a class-0 target and as Contradicts for a class-1 target
(p.216).

**Package.** `bl_shapley()` labels a value with $|\phi_j| < 10^{-8}$ “No
effect”, and `bl_find_sparse_cf()` keeps such a variable at its observed
value.

**Reason.** A variable whose change never moves the prediction should
not be changed in the sparse counterfactual. This follows the thesis
text, which the code’s labelling did not.

**Effect on Pima.** “No effect” occurs only where *Pregnancies* ends at
its observed value: the true negative’s run 7, after the restore. The
true negative is a class-0 target, so Appendix A would have labelled it
Supports, and the sparse counterfactual is the same under both rules.

## 3.10 The sparse counterfactual: a greedy rule with a stopping point

**Thesis.** “only changing the variables with combined Shapley values
equal to the change needed” reduces the distance and the number of
variables (p.120). The Pima example changes *BMI*, *Glucose*, *Insulin*
and *Age* and keeps *Pregnancies* (p.122). No stopping rule is given; in
Appendix A the variables to change are chosen by hand (p.218).

**Package.** `bl_find_sparse_cf()` takes the Supports variables in
decreasing order of $|\phi_j|$, changes one at a time to its
counterfactual value, re-scores the point, and stops at the first point
that crosses the boundary. The path is stored. Until 30 October 2026 the
package changed all Supports variables at once.

**Reason.** The Shapley values add up to the full change only at the
full counterfactual, so a subset must be re-scored. This rule is
reproducible and gives the shortest prefix of the ranked Supports
variables that changes the class.

**Effect on Pima.** For the true positive, the greedy rule changes 4
variables (*Glucose, Insulin, Age, Pregnancies*, prediction 0.481),
against 4 when every Supports variable changes (prediction 0.481). For
the true negative it changes 3 (prediction 0.526) against 4 (prediction
0.548). Both versions cross the boundary for both targets. For the true
positive every Supports variable is needed, so the two rules agree; for
the true negative the greedy rule is sparser by one variable.

## 3.11 Outlier removal before the projection

**Thesis.** Outlier removal is optional. “No outlier data is removed” in
the Pima example (p.99); in the loan example “all data points lying
outside the 90% convex hull of the training data were removed” and “the
convex hull fraction was adjusted” (p.127). Appendix A computes the hull
on the training biplot itself.

**Package.** `bl_prepare_data(hull_fraction = )` removes the training
rows outside a convex hull that encloses the given fraction of them,
computing the hull on a two-component PCA of the standardised training
features, before the model is fitted. The test data is never filtered.
The default was 0.9 from March to November 2026; since 5 November 2026
it is 1, which removes nothing, so outlier removal is optional, as in
the thesis. Each application sets its own value: the loan harness and
vignette use 0.9, as the thesis’s loan example does, and the Pima
example of Paper 1 uses 1, as the thesis’s Pima example does.

**Reason.** The filter runs as part of the data preparation, before the
model and the final projection exist, so it uses a stand-in projection.
A value below 1 is worth choosing when extreme training rows would
stretch the biplot and its prediction grid; the default keeps the data
as given.

**Effect on Pima.** None: Paper 1 sets `hull_fraction = 1`, which is now
also the default. Not measured.

## 3.12 Three properties the thesis does not state

These properties are unchanged from the thesis, but Paper 1 relies on
them, so they are recorded here.

**Components beyond the canonical dimensions.** With $g$ groups, CVA has
$K = \min(g - 1, p)$ non-zero eigenvalues (§3.4, p.49), so $K = 3$ for
the four prediction outcomes. The remaining $p - K$ components span the
null space of the between-groups matrix. When that space has two or more
dimensions, its basis is not unique and is set by the numerical routine;
this holds for the seven-feature model ($p - K = 4$), where the thesis’s
pair (1,7) uses such a component. On the four-feature model of Paper 1,
$p - K = 1$, so the fourth component is fixed up to its sign. Both
counterfactuals use it (pairs (3,4) and (2,4)).

**The rotation is not unique.** The Procrustes rotation aligns the plane
with a single target point, so the cross-product matrix has rank 1. This
fixes the direction of the target, but the second direction of the plane
is left to a numerical tie-break. The thesis does not discuss this
(§3.5, pp.52–54); it is inherited from the thesis code and documented in
the package’s rotation note (§4.1). The results are deterministic for
identical inputs, but a perturbation of the loadings at machine
precision can move the plane. On Pima, ordering the CVA groups as the
thesis code did changes the loadings by no more than 5.0e-17, yet the
true positive loses its counterfactual: its search point moves, and no
whole number of pregnancies near the new point crosses the boundary
(against pair (3,4) and distance 1.82 with the package’s order). A
deterministic replacement is proposed in Section 5.

**The reported distance is an upper bound (the approximation gap).** The
counterfactual is the nearest boundary point within the planes that are
searched, so its distance is at least the distance to the nearest point
of the decision boundary in the full space of the four variables. A
constrained search in the full space measures the gap: for each of the
79 test rows, the smallest Mahalanobis distance at which the GAM’s
probability equals 0.5 within the training ranges, found by a
gradient-based solver from the package’s point and from 50 random
directions. The package’s search distance is a median 1.19 times the
full-space distance (quartiles 1.09 and 1.31, range 1.00 to 2.75), and
within 5% of it for 15% of the rows. For the two examples the ratios are
1.52 (true positive) and 1.05 (true negative). 59 of the full-space
points lie inside the convex hull of the training data in four
dimensions, and their median ratio is 1.18, so the gap is not explained
by implausible full-space points. *Proposed paper sentence:* “The biplot
search finds a boundary point a median 19% further from the target than
the nearest boundary point in the full four-dimensional space, so the
distances it reports are upper bounds (quartiles 9% and 31% further).”

# 4 Changes made: code improvements

These changes alter how the method is computed, its default parameters
or its numerical behaviour, not what it computes. Table 4 summarises
them with their measured effect on Pima. Run 1 is the paper’s
unconstrained run (*Pregnancies* `"integer"`).

| Item | Thesis code → package | Why | Effect on Pima |
|:---|:---|:---|:---|
| Contour band and prediction rounding | $\varepsilon = 0.01$ and predictions floored to 2 decimals → $\varepsilon = 0.001$ (as in the text, p.78) and 3 decimals | Separate the band from the rounding; a narrower band | With $\varepsilon = 0.01$: true-positive distance 1.828 against 1.819, prediction 0.472 against 0.481; true-negative distance 0.418 against 0.419; the target rule picks test row 1 for the true positive (row 23 at 0.001). The 2-decimal floor cannot be switched back. |
| Global grid hull fraction | 0.9 → 1 | Every training point inside the hull | Run 1 unchanged for both targets; the target rule picks test row 69 for the true negative (row 25 at 1) |
| Order of the CVA groups | TN, FN, FP, TP (colour names) → FN, FP, TN, TP, fixed | Reproducibility | True positive: no counterfactual with the thesis order (no whole-number point), against pair (3,4) and distance 1.82; true negative: 0.39 against 0.42 (Section 3.12) |
| Inverse of $\mathbf{W}$ | `solve()` → Cholesky, with a ridge fallback | Numerical stability; detects a singular $\mathbf{W}$ | The ridge fallback is not needed on either model |
| Grid probabilities | No clamp → clamped to \[0, 1\], with a warning, on the global grid (June 2026) and on each local grid (6 Nov 2026) | Out-of-range model outputs gave missing colours, and on a local grid could stop the search | No warning raised |
| The local plot’s window | Set by the training data → widened to show a target or counterfactual outside it (6 Nov 2026) | A target outside the training data was not drawn | Display only |
| The two-colour grid of the global plot | Split at 0.5 → at the model’s cutoff (6 Nov 2026) | Consistency with the cutoff | Display only; the same at the cutoff of 0.5 |
| Nearest boundary vertex | Appendix A `which.min`; branch random tie-break → the first nearest vertex | Determinism, as the text states (p.92) | Not measured; the scale on the loan data is in the package’s tie-break note (§5) |
| Data behind the local plane | Test data and plot bounds → training data, range plus 10% | The training data defines the projection | Not measured |
| Confusion colours | The thesis scheme in every plot except `plot.bl_projection()` → one shared scheme (TP red, TN blue, FP purple, FN orange) | One colour means one outcome in every biplot | Display only |
| A target at the centre of the biplot | No check → a warning (the rotation is undefined there) | The local plane is arbitrary at the centre | No result changes |
| Hull for the global test search | Test-data hull (branch) → training hull | As the text states (p.93) | Not measured |
| Constraint syntax and pair order | Free expressions; pairs in `combn` order → named forms (`"increase"`, `"decrease"`, ranges); pairs by the second component | Usability | Only ties are affected |
| New diagnostics | None → Mahalanobis distance in the global search; the Phase 2 Mahalanobis check | Reporting | No search result changes |
| Minor items | Grid size; the source of the cutoff; the fallback when no boundary is found; the order of the filters; the Shapley permutation seed ($p > 14$ only); centring on the filtered data; a pair without a boundary point is skipped; PCA allowed in the local search; empty CVA groups dropped |  | No effect on the Pima results |

Table 4. Code improvements and their effect on Pima.

**Bug fixes.** One line each:

- `.bl_rotate()` returned the wrong two columns for pairs other than
  (1,2), placing the target at the origin; fixed on 23 May 2026,
  restoring the thesis behaviour.
- Task A: two internal helpers switched the user’s active graphics
  device.
- Task B: the local plot coloured the target by another observation’s
  class when the target came from a data frame other than the test data.
- Gate 2:
  - fixed variables were marked as used in the sparse counterfactual;
  - the local plots accepted a rotation that misplaced the sparse
    counterfactual;
  - the grid cutoff was not checked against the model’s;
  - models from `xgboost::xgboost()` were rejected;
  - empty CVA groups were coloured wrongly;
  - raw-unit data could be passed to functions that expect model units
    (`bl_apply_scaling()` and a warning were added);
  - the Shapley labels read “observed -\> change”.

**Future work.** The open items that bear on the method are listed in
Section 5.

# 5 Future work

None of these is implemented. Each is left as future work, and Paper 1
names it as an area of future improvement. Table 5 summarises them; each
is described below. The complete list for later projects, with the code
each item touches, is in `../future-improvements-review.md`.

| Change | Why | Expected effect on results | Effort | Decision |
|:---|:---|:---|:---|:---|
| Minimal rotation with a search over planes | Makes the plane unique (Section 3.12) | Moves every Phase 3 result | High | Future work, after Paper 1 |
| A check on other platforms | The plane depends on the linear algebra routine (Section 3.12) | Tests whether the results reproduce | Low | Future work; Paper 1 states its platform |
| Follow-up searches for `"fixed"` and `"integer"`, and a search in the actionable variables only | Fewer runs without a counterfactual; smaller rounding cost | Changes constrained runs only | Medium to high | Future work |
| Feasibility as a constraint inside the search | Finds the nearest feasible point directly | Changes constrained runs | Medium | Future work |
| A configurable `"fixed"` tolerance | The 0.5 tolerance is hard-coded | None at the default | Low | Future work |
| A per-row feasibility flag, and the training-range filter as an option, in the global search | Shows which counterfactuals are infeasible; offers the text’s rule (Section 3.3) | None at the default | Low | Future work |
| Storing every pair’s plane | Shows the search across planes; allows a failed restore to be retried on the next-best pair | None | Low to medium | Future work |
| Full Mahalanobis decomposition in Phase 2 | Removes the gap of Section 3.1 | Changes the importance values | Medium | Future work |
| Aggregated Shapley values | Global importance from local attributions (thesis §6.2) | New output | High | Future work |
| CVA biplot for a singular $\mathbf{W}$ ($p > n$) | Wide data | New capability | Medium to high | Future work |
| Multi-class models | One boundary for each pair of classes | New capability | High | Future work, a key area |
| A different transformation per variable | Skewed or mixed-scale variables (display only) | None | Medium | Future work |

Table 5. Possible future changes.

**Minimal rotation with a search over planes.** The package’s rotation
note (§5) proposes a closed-form minimal (Givens) rotation that carries
the plane onto the target, followed by a reproducible search over the
family of planes that contain the target. This removes the tie-break of
Section 3.12 and broadens the search beyond fixed pairs, in the spirit
of the grand tour proposed in the thesis (§6.3, p.144). Since every
Phase 3 number would move, the change should follow Paper 1.

**A check on other platforms.** The Procrustes rotation leaves the
plane’s second direction to the linear algebra routine (Section 3.12),
so another operating system or linear algebra library may select a
different plane, and hence a different counterfactual. The package has
been run on one platform (R 4.6.0 on Windows); Paper 1 states this, and
its supplement records the session information. Rerunning the Pima
example on other platforms would test whether its results reproduce; the
minimal rotation above would remove the dependence.

**Follow-up searches for `"fixed"` and `"integer"`.** When the restore
of a fixed variable fails, a second search could move the free variables
along the search direction, with the fixed variables held exactly, until
the model crosses the boundary; building the local biplot in the
actionable variables only is the larger alternative. After an integer
rounding, the other variables could compensate, reducing the distance
cost of Section 3.8, and a step size other than 1 could be allowed.

**Feasibility as a constraint inside the search.** The constraints and
the training ranges filter the vertices of the boundary contour, so the
nearest feasible point may lie between two vertices. Treating
feasibility as a constraint inside the search would find it directly.

**A configurable `"fixed"` tolerance and the global training ranges.**
The 0.5 tolerance suits variables measured in whole units, but not
variables on a small scale. For the global search, a flag on each
counterfactual outside the training ranges would show which rows are
affected, since the package reports only their number. An option to
filter by the ranges before the selection would offer the thesis version
of Section 3.3 without changing the default.

**Storing every pair’s plane.** Only the winning plane is kept, so the
other candidates cannot be drawn. Storing all planes would allow a
figure of the search across planes, and the restore of a fixed variable
to be retried on the next-best pair, without changing any other result.

**Full Mahalanobis decomposition in Phase 2.** A Shapley decomposition
of the full Mahalanobis distance would attribute the correlated part of
the distance that the per-variable measure ignores. Its cost grows with
$2^p$ per data point.

**Aggregated Shapley values.** The thesis proposes a global importance
measure that aggregates the Shapley contributions over a representative
set of data points (§6.2, p.143). This needs a rotation and a Shapley
computation for each data point.

**A CVA biplot for a singular $\mathbf{W}$.** When $p > n$, $\mathbf{W}$
is singular. The generalised singular value decomposition of Ganey &
Gardner-Lubbe (2026, *Statistical Methods & Applications* 35, 825–847; R
package `wideRhino`) gives a CVA biplot in this case.

**Multi-class models.** The method explains binary classifiers. A
multi-class model can be explained one class at a time, by coding the
class of interest as 1 and the other classes as 0, so that the
counterfactual moves the observation out of that class rather than
towards a specific other class. The CVA biplot already separates more
than two groups (here the four outcome groups), so extending the
boundaries and the search to one boundary for each pair of classes, for
example through linear discriminant analysis, is a key area of future
research.

**A different transformation per variable.** The package records one
affine scaling for display. Allowing a different transformation for each
variable (z-score, log, none) would help skewed variables. A non-affine
transformation needs calibrated non-linear axes, such as the
spline-based biplots of Groenen, le Roux & Gardner-Lubbe (2015,
*Advances in Data Analysis and Classification* 9, 219–238).

# Reproduction

The measurements were made at package commit e41f22f (`boundarylogic`
0.2.0, R version 4.6.0 (2026-04-24 ucrt), mgcv 1.9.4, run 2026-11-15
07:05:40). The scripts of the Pima example are kept outside the package
repository. From that folder:

1.  `Rscript run_all.R` builds the data, the models and the Paper 1
    results;
2.  `Rscript 04_method-change-effects.R` reruns the searches with each
    change switched back and writes `method-change-effects.csv` to this
    folder; it stops if the baseline does not reproduce the Paper 1
    results;
3.  `Rscript 05_full-space-gap.R` measures the approximation gap of
    Section 3.12;
4.  `Rscript 06_findings-effects.R` compares the package before the
    changes of 6 November 2026 (commit ef014ca, exported from git) with
    the current one and writes `findings-effects.csv` to this folder.

The values in this note are those of the two CSVs and of `full_space_gap.rds` at the commit
above; after a re-run, update the values here by hand.
