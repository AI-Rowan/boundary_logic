# Future developments

Date created: 2026-11-15

This note lists the developments of `boundarylogic` that later projects can take up: changes to the
method, to the code and its interface, to the tests, and to the documentation and release. None is
implemented, and none is part of Paper 1. Each needs an explicit instruction (CLAUDE.md §9).

Each item gives what it is and why, the code it touches, its expected effect on results, a rough
effort, and its source. Thesis pages are the printed pages of the final thesis.

The sources are:
- the methodology note's "Future work" (`methodology-changes/methodology-changes.md`, Table 5);
- the review form of the open items (`methodology-changes/future-improvements-review.md`, items A
  and B);
- the future sections of the technical notes: the rotation note (`biplot-rotation/`, §5 and §6),
  `fixed-constraint-restore.md` §8 and `boundary-distance-and-robustness.md` §7.

---

## 1. Method

These change what the method computes, or add a new result.

### 1.1 A unique local plane: minimal rotation with a search over planes

- **What and why.** `.bl_rotate()` aligns the target with the chosen component pair through a
  Procrustes rotation whose cross-product matrix has rank 1. The minimiser is fixed only along the
  target's direction, and the plane's second direction is the linear-algebra routine's tie-break.
  - **On the loan data:** bases that differ by 4.65e-15 give planes 65.91° apart, and only 21% of a
    selected (2,3) plane's second direction lies in the (2,3) eigenplane.
  - **The results** reproduce for identical inputs on one platform; across platforms this is
    untested (§3.3).
  - **The proposal:** a closed-form minimal (Givens) rotation that carries the plane onto the target,
    followed by a reproducible search over the family of planes that contain it. This makes the
    plane unique and broadens the search beyond fixed pairs, in the spirit of the thesis's grand
    tour (§6.3, p.144).
- **Frame convention.** The same work fixes the orientation of the axes within the local plane,
  which is now also the routine's choice. This changes the figures' orientation, not their numbers.
- **Code:** `.bl_rotate()`, the pair loop of `bl_find_local_cf()`, `plot.bl_local_result()`.
- **Effect:** moves every Phase 3 result (counterfactuals, Shapley values, sparse counterfactuals).
- **Effort:** high. The rotation note's acceptance tests 4 (monotonicity, like for like), 5 (variance
  retained) and 7 (cost for $p\ge4$) are incomplete (§3.4).
- **Source:** rotation note §4.1, §4.3, §5, §6.1; methodology note Table 5; review items B1, B4.

### 1.2 Feasibility as a constraint inside the search

- **What and why.** The constraints and the training ranges filter the vertices of the boundary
  contour, so the nearest feasible point can lie between two vertices. On the loan data the
  training ranges remove 18–56% of the candidates per pair. Treating feasibility as a constraint
  inside the search would find the nearest feasible point directly.
- **Code:** the filter cascade of `bl_find_local_cf()`.
- **Effect:** changes the constrained runs.
- **Effort:** medium.
- **Source:** rotation note §4.5, §5.4; methodology note Table 5; review item B6.

### 1.3 Counterfactuals that hold `"fixed"` features exactly

A failed restore of a `"fixed"` feature means no counterfactual (`fixed-constraint-restore.md` §4),
although one may exist. Three remedies:

| Remedy | What it does | Effect | Effort | Source |
|---|---|---|---|---|
| A follow-up search from the restored point | Holds the fixed features at their observed values and extends the move in the free features along the search's direction until the model crosses the cutoff, within the other constraints and the training ranges. A 1-D, deterministic search; its point lies off the plane | Fixed runs (the Pima TN target's runs 3 and 4) | Medium | review item A5 |
| Searching in the actionable features only | Builds the local biplot in the subspace of the free features, so every boundary point holds the fixed features exactly and no restore is needed | Every fixed run; a different projection for the local search | High | review item A6 |
| Retrying the restore on the next-best pair | Tries the restore on the other pairs' points, nearest first, and takes the first that still crosses | Fixed runs | Low, once §2.1 exists | review item B12 |

**Code:** `.restore_fixed()`, `bl_find_local_cf()`. **Source:** `fixed-constraint-restore.md` §8.

### 1.4 Integer rounding

| Development | What it does | Effect | Effort | Source |
|---|---|---|---|---|
| Free features compensate after rounding | After an `"integer"` feature is rounded, the other features move to reduce the distance cost the rounding adds | Integer runs (the Pima TN target's runs 1 and 5) | Medium | review item A7 |
| Other step sizes | A step size per feature in `set_filters()`, for example 0.5, or 100 for a loan amount; only whole numbers are supported | None at the default | Low | review item B13; `fixed-constraint-restore.md` §8; `sparse-counterfactual.md` §3 |

**Code:** `.round_integers()`, `.parse_constraint()`, `set_filters()`.

### 1.5 Off-plane counterfactuals

- **What and why.** The restored, rounded and sparse counterfactuals lie off the plane, while
  `dist_*` stay the search's. The plot could mark these points as off-plane, or the distances could
  be recomputed for the reported point.
- **Code:** `bl_find_local_cf()`, `plot.bl_local_result()`, `plot.bl_sparse_result()`.
- **Effect:** the printed distances only.
- **Effort:** low to medium.
- **Source:** rotation note §4.6; review item B7.

### 1.6 Full Mahalanobis decomposition in Phase 2

- **What and why.** `plot.bl_boundary()` and `bl_robustness()` standardise each feature by its own
  within-class SD (the diagonal of **W**), so the plotted measure ignores the correlations in **W**.
  On Pima the gap to the full Mahalanobis distance is 27.4% on seven features and 35.2% on four.
  A Shapley decomposition of the full distance attributes the correlated part too.
- **Code:** `plot.bl_boundary()`, `bl_robustness()`.
- **Effect:** changes the importance values.
- **Effort:** medium. The cost grows as $O(n\cdot 2^p)$.
- **Source:** `boundary-distance-and-robustness.md` §7; `2 implementation_summary.txt` §4.6.1;
  methodology note Table 5; review item B14.

### 1.7 New results

| Development | What and why | Code | Effort | Source |
|---|---|---|---|---|
| Aggregated Shapley values | Global importance built from the local attributions over a representative set of data points (thesis §6.2, p.143). Needs a rotation and a Shapley computation per point | new function after `bl_shapley()` | High | methodology note Table 5; review item B15 |
| Calibrated probabilities from the boundary distance | Platt scaling on the length of the change vector (Euclidean or Mahalanobis) gives a model-agnostic posterior class probability (thesis §6.4, pp.145–146) | new function after `bl_find_boundary()` | Medium | `boundary-distance-and-robustness.md` §7 |
| A Cholesky-whitened diagnostic plot | Shows the whitened decomposition ($u_k^2$) across the training set, without claiming a per-feature attribution; useful when the Mahalanobis check's percentage is large | `plot.bl_boundary()` | Low | `boundary-distance-and-robustness.md` §7 |

### 1.8 New capabilities

| Development | What and why | Code | Effort | Source |
|---|---|---|---|---|
| Multi-class models | The method explains binary classifiers: a multi-class model is explained one class at a time (the class of interest as 1, the others as 0), so the counterfactual leaves that class rather than moving towards a specific other class. One boundary per pair of classes, for example through linear discriminant analysis (`MASS::lda()` in `bl_wrap_model()`, or a `bl_fit_lda()`), is a key area of future research. The CVA biplot already separates more than two groups | `bl_wrap_model()`, `.pred_function()`, the grid, both searches | High | methodology note Table 5; review item B17; `2 implementation_summary.txt` §4.1 |
| A CVA biplot for a singular **W** ($p>n$) | **W** is singular for wide data. The generalised SVD of Ganey & Gardner-Lubbe (2026, *Statistical Methods & Applications* 35, 825–847; R package `wideRhino`) gives a CVA biplot in this case | `bl_build_projection()` | Medium to high | methodology note Table 5; review item B16 |
| A different transformation per feature | `bl_set_scaling()` records one affine method for all features (each with its own centre and scale). A transformation per feature (z-score, min-max, log or Box-Cox, none) would help skewed or mixed-scale features. A non-affine inverse needs calibrated non-linear axes, such as the spline-based biplots of Groenen, le Roux & Gardner-Lubbe (2015, *Advances in Data Analysis and Classification* 9, 219–238; biplotEZ's spline axes), or a stated "affine only" limit. Display only | `bl_set_scaling()`, `.bl_rescale_biplot_axes()`, `.scale_to_raw()`, `bl_apply_scaling()` | Medium | methodology note Table 5; review item B18 |

---

## 2. Code and interface

These change no result at the defaults.

### 2.1 Keeping every pair's plane and boundary point

- **What and why.** `bl_find_local_cf()` keeps only the winning plane. Storing every pair's plane
  and boundary point allows a figure of the search across planes, and §1.3's retry on the next-best
  pair.
- **Code:** the result of `bl_find_local_cf()`.
- **Effect:** none.
- **Effort:** low to medium.
- **Source:** methodology note Table 5; review item A2.

### 2.2 Smaller interface changes

| Development | What and why | Code | Effort | Source |
|---|---|---|---|---|
| A configurable `"fixed"` tolerance | `.fixed_tol` is 0.5 model units, which suits features measured in whole units but not features on a small scale | `bl_find_local_cf()`, `.apply_actionability()` | Low | methodology note Table 5; review item A4 |
| A per-row feasibility flag in the global search | `bl_find_boundary()` keeps counterfactuals outside the training ranges and reports only their number (10 of 313 training rows and 5 of 79 test rows on the 4-feature Pima model). A flag per row would show which; an option to filter by the ranges before the selection would offer the thesis's rule (§3.10 step 4, p.93) | `bl_find_boundary()` | Low | methodology note Table 5; review item B11 |
| Pair control | A `pairs =` argument to force any pair, and a `max_pairs` default tied to $\binom{p}{2}$ so the default covers every pair for any $p$ (21 covers every pair up to $p=7$) | `bl_find_local_cf()` | Low | rotation note §4.7, §6.2; review item B8 |
| The pair label in output | Report how far the selected plane lies from the eigenplane its pair names (the angle $\theta$, or the overlap) | `print.bl_local_result()` | Low | rotation note §4.2, §6.2 |

---

## 3. Tests

The tests cover the code paths and the display rules; these would pin computed values.

| Development | What it pins | Note | Source |
|---|---|---|---|
| 3.1 Values of computed results | The `bl_find_local_cf()` boundary point; the Shapley values on a real model; the values of `bl_robustness()` and `bl_surrogate()`; `bl_build_result()` | The local plane depends on the linear-algebra routine (§1.1), so values need a tolerance and a stated platform, or §1.1 first | review item B21 |
| 3.2 The CVA identity | $d^2_M=(n-g)\,d_Z^2$ for every pair, including two-group CVA (rotation note §3.4) | `17_3` §11 asserts it; `tests/testthat/` does not | rotation note §6.2; review item B5 |
| 3.3 Cross-platform reproducibility | Re-run the Pima pipeline (`run_all.R`, `tests.R`) and `17_3` on Linux and macOS and on other BLAS builds, for example through GitHub Actions, and compare `numbers.rds` | Because of §1.1, small numerical differences matter: listing the CVA groups in another order changes **V** by 5.0e-17, and the Pima TP target's run 1 then has no counterfactual (methodology note §3.12, measurement K3) | rotation note §4.1, §6.1 test 2, §6.2; review item B2 |
| 3.4 The proposal's acceptance tests | Tests 4, 5 and 7 of the rotation note §6.1 | Needed before §1.1 can replace the current rotation | rotation note §6.1 |

---

## 4. Documentation and release

| Development | What it gives | Source |
|---|---|---|
| A vignette on the mathematical method | Derives the biplot projection and the counterfactual geometry from first principles | review item B22 |
| CRAN preparation | Checks, examples within time limits, and the submission | review item B23 |
