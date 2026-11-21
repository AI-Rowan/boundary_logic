# Future improvements

Date created: 2026-11-05
Date last edited: 2026-11-16

This file is the only place for the future developments of `boundarylogic`: changes to the method,
to the code and its interface, to the tests, and to the documentation and release. Each item says
why it is future work and not done now. None is implemented, and none is part of the package that
Paper 1 uses; Paper 1 names some of them as areas of future improvement (each item's **Decision**
says which). Each needs an explicit instruction (CLAUDE.md §9). CLAUDE.md and the technical notes
point here and do not list future developments themselves.

**Part 1** is the list. Each item gives:
- **What and why:** what the development is, and what it would give.
- **Decision:** the decision on it, from the replies of 2026-11-05 (the Appendix) or a later record.
- **Reason:** why it is future work and not done now, as the decision records it. "Reason: not
  recorded" means that the record gives no reason.
- **Code**, **Effect** (the expected effect on results), **Effort** (a rough estimate) and
  **Source** (the reply, A1–A9 or B1–B24, or the note section it comes from).

**The Appendix** keeps the replies of 2026-11-05: the form in which every open item was decided (fix
now, future improvement or no action), with its later rounds. It is the dated record of those
decisions.

Paths are relative to `documentation/`. "The rotation note" is
`biplot-rotation/biplot_rotation_and_pair_selection.md`, and "the methodology note" is
`methodology-changes/methodology-changes.md`. Thesis pages are the printed pages of the final thesis.

---

# Part 1: The developments

## 1. Method

These change what the method computes, or add a new result.

### 1.1 A unique local plane: minimal rotation with a search over planes

- **What and why.** `.bl_rotate()` aligns the target with the chosen component pair through a
  Procrustes rotation. Its target matrix, `Y <- rbind(-x, 0, x)`, has rank 1, so the Procrustes
  cross-product has rank 1 too. The minimiser is then fixed only along the target's direction, and
  the plane's second direction is LAPACK's tie-break for a degenerate subspace.
  - **On the loan data** (GAM, PCA, $p=4$): bases that differ by 4.65e-15 give planes 65.91° apart,
    and only 21% of a selected (2,3) plane's second direction lies in the (2,3) eigenplane. The
    target lies in every such plane, so the check a reader would make, whether the target lands
    correctly, passes either way.
  - **Reproducibility.** The plane is deterministic for identical inputs on one platform. Across
    platforms it is untested (§3.3). The thesis method has the same property: it is inherited, not a
    porting error.
  - **The proposal.** A closed-form minimal (Givens) rotation carries the named eigenplane onto the
    target. A reproducible search over the family of planes that contain the target follows. This
    makes the plane unique and broadens the search beyond fixed pairs, in the spirit of the thesis's
    grand tour (§6.3, p.144). The rotation note §5 has the derivation and the measurements (its
    A.9, A.10). It stays there.
  - **Adoption path** (rotation note §5.5, §5.6). The proposal lands as new options with the
    defaults unchanged: `rotation = c("procrustes", "minimal")` and
    `plane_search = c("pairs", "family")`. The defaults switch only after the acceptance tests
    (§3.4) pass across a full training set. The output then reports the objective over the searched
    family with its runner-up minima, and keeps $\theta$ and the nearest named pair.
  - **Frame convention** (B4). The same work fixes the orientation of the axes within the local
    plane, which is now also LAPACK's choice (rotation note §4.3, §5.2). This changes the figures'
    orientation, not their numbers.
  - **What it does not fix**. The search still finds the best plane within the
    biplot framework, not the nearest boundary in the full space. The Shapley values and the sparse
    counterfactual still inherit the plane that is drawn. Restored, rounded and sparse
    counterfactuals still leave the plane (§1.5).
- **Decision:** future improvement, after Paper 1, named in the paper (B1; B4 goes with it). Do not
  change `.bl_rotate()` or implement the proposal without an explicit instruction.
- **Reason:** it moves every Phase 3 number, so it follows Paper 1 (B1: "Future (after Paper 1)";
  methodology note §5). Because it moves every Phase 3 counterfactual, it is not covered by the rule
  of CLAUDE.md §4 that the maths may be replaced only when the outputs are proven equivalent. Its
  acceptance tests 4, 5 and 7 are also incomplete (§3.4).
- **Code:** `.bl_rotate()`, the pair loop of `bl_find_local_cf()`, `plot.bl_local_result()`.
- **Effect:** moves every Phase 3 result (counterfactuals, Shapley values, sparse counterfactuals).
- **Effort:** high.
- **Source:** B1, B4; CLAUDE.md §9; rotation note §4.1, §4.3, §5; methodology note §3.12,
  §5 (Table 5).

### 1.2 Feasibility as a constraint inside the search

- **What and why.** The constraints and the training ranges filter the vertices of the boundary
  contour after it is extracted. The reported point is the nearest among the survivors, so the
  nearest feasible point can lie between two vertices. On the loan data the training ranges remove
  18–56% of the candidates per pair, and a replication of the pair search without that filter
  agrees with the package on only 76 of 221 rows. Restricting the grid and the contour extraction to
  the feasible region, with the actionability constraints as box or half-space constraints, would
  find the nearest feasible point directly.
- **Decision:** future improvement, named in the paper (B6).
- **Reason:** not recorded. The form notes that it would move Paper 1's numbers.
- **Code:** the filter cascade of `bl_find_local_cf()` (`.apply_actionability()`, the training-range
  filter).
- **Effect:** changes the constrained runs.
- **Effort:** medium.
- **Source:** B6; rotation note §4.5, §5.4; methodology note §5 (Table 5).

### 1.3 Counterfactuals that hold `"fixed"` features exactly

A failed restore of a `"fixed"` feature means no counterfactual (`fixed-constraint-restore.md` §4),
although one may exist. On Pima, the TN target's runs 3 and 4 find boundary points, but the
restored points score 0.475 and 0.487, so neither run reports a counterfactual (methodology note
§3.7). Three remedies follow.

#### 1.3.1 A follow-up search from the restored point

- **What and why.** Hold the fixed features at their observed values, and extend the move in the
  free features along the search's own direction until the model crosses the cutoff, within the
  training ranges and the other constraints. It is a 1-D search, so it is cheap and deterministic.
  Its point lies off the plane, as a restored point already does, and its distance is at least the
  restored point's.
- **Decision:** future improvement, named in the paper (A5, option (b)). Paper 1 reports the empty
  runs as a result of the strict rule (Pima issue 12).
- **Reason:** it would move Paper 1's numbers (the TN target's runs 3 and 4 may gain
  counterfactuals). It is a new method change, which needs the supervisor's approval, and option (b)
  keeps Paper 1 on the method the supervisor was reviewing (A5).
- **Code:** `.restore_fixed()`, `bl_find_local_cf()`.
- **Effect:** the fixed runs (the Pima TN target's runs 3 and 4).
- **Effort:** medium.
- **Source:** A5; `fixed-constraint-restore.md` §8; methodology note §3.7, §5 (Table 5); Pima
  issue 12.

#### 1.3.2 Searching in the actionable features only

- **What and why.** Hold the fixed features at their observed values before the rotation, so that
  the local biplot and its boundary are built in the subspace of the free features. Every boundary
  point then holds the fixed features exactly, and no restore is needed. This is the first of the
  thesis's two approaches to constraints (§3.9.2, p.88). It changes the projection the local search
  uses: a different one for each set of constraints.
- **Decision:** future improvement, named in the paper (A6, option (b)).
- **Reason:** it is the cleanest solution, but too large a change to make during Paper 1. It is a
  change of method, which needs the supervisor's approval, and it would move every fixed run (A6).
- **Code:** `.bl_rotate()`, the pair loop of `bl_find_local_cf()`.
- **Effect:** every fixed run; a different projection for the local search.
- **Effort:** high.
- **Source:** A6; `fixed-constraint-restore.md` §8; methodology note §3.7, §5 (Table 5).

#### 1.3.3 Retrying the restore on the next-best pair

- **What and why.** The restore is tried only on the winning pair's point. Trying it on the other
  pairs' points, nearest first, and taking the first that still crosses would rescue some targets.
- **Decision:** future improvement, named in the paper (B12).
- **Reason:** it needs each pair's boundary point to be kept (§2.1), which the search does not do and
  which is itself future work (B12: "Future (needs A2)").
- **Code:** `.restore_fixed()`, `bl_find_local_cf()`.
- **Effect:** the fixed runs.
- **Effort:** low, once §2.1 exists.
- **Source:** B12; `fixed-constraint-restore.md` §8.

### 1.4 Integer rounding

#### 1.4.1 Free features compensate after rounding

- **What and why.** `.round_integers()` tries only the whole values just below and just above each
  moved `"integer"` feature, with the other features where the search put them. When the nearer
  whole value does not cross, the further one is taken, even when a small change in the free
  features would let the nearer one cross. Holding the nearer value and running §1.3.1's follow-up
  search along the free features would reduce the distance cost of the rounding. On Pima, the TN
  target's search point has *Pregnancies* 10.21 ($\sqrt{d^2_M}$ 0.14); 10 does not cross, so 11 is
  taken, and $\sqrt{d^2_M}$ rises to 0.42.
- **Decision:** future improvement, named in the paper (A7, option (b)).
- **Reason:** it needs §1.3.1's follow-up search, which is itself future work: the two belong
  together (A7). It would also move Paper 1's numbers (the TN target's runs 1 and 5).
- **Code:** `.round_integers()`, `bl_find_local_cf()`.
- **Effect:** the integer runs (the Pima TN target's runs 1 and 5).
- **Effort:** medium, or low once §1.3.1 exists.
- **Source:** A7; `fixed-constraint-restore.md` §5; methodology note §3.8, §5 (Table 5); Pima
  issue 13.

#### 1.4.2 Other step sizes

- **What and why.** A step size per feature in `set_filters()`, for example 0.5, or 100 for a loan
  amount. Only whole numbers are supported.
- **Decision:** future improvement, named in the paper (B13).
- **Reason:** not recorded.
- **Code:** `.round_integers()`, `.parse_constraint()`, `set_filters()`.
- **Effect:** none at the default.
- **Effort:** low.
- **Source:** B13; `fixed-constraint-restore.md` §8; `sparse-counterfactual.md` §3; methodology note
  §5.

### 1.5 Off-plane counterfactuals

- **What and why.** The restored, rounded and sparse counterfactuals lie off the plane. The plot
  draws them at their projections, which need not lie on the drawn boundary. For a restored or
  rounded counterfactual, `dist_z`, `dist_x` and `dist_mahalanobis` describe the search's point, not
  the reported one. The plot could mark these points as off-plane, or the distances could be
  recomputed for the reported point.
- **Decision:** future improvement (B7). Paper 1 states which distances Table 1 quotes (the search
  point's and the counterfactual's $\sqrt{d^2_M}$), as the method stands.
- **Reason:** not recorded. The form notes that only the printed distances would change.
- **Code:** `bl_find_local_cf()`, `plot.bl_local_result()`, `plot.bl_sparse_result()`.
- **Effect:** the printed distances only.
- **Effort:** low to medium.
- **Source:** B7; rotation note §4.6.

### 1.6 Full Mahalanobis decomposition in Phase 2

- **What and why.** `plot.bl_boundary()` and `bl_robustness()` standardise each feature by its own
  within-class SD (the diagonal of **W**), so the plotted measure ignores the correlations in **W**.
  On Pima the gap to the full Mahalanobis distance is 27.4% on seven features and 35.2% on four. A
  Shapley decomposition of the full distance would attribute the correlated part too: it is the
  theoretically correct per-feature decomposition.
  - **A closed form.** With the value function $d^2_M=v'\mathbf{W}^{-1}v$ (a feature outside the
    coalition keeps its observed value), the exact Shapley values are
    $\phi_j=v_j\,(\mathbf{W}^{-1}v)_j$. Each cross term $a_{ij}v_iv_j$ ($a_{ij}$ an entry of
    $\mathbf{W}^{-1}$) is split equally between $i$ and $j$. The values sum to $v'\mathbf{W}^{-1}v$,
    need no model call, and cost $O(p^2)$ per row.
  - **A caveat.** A $\phi_j$ can be negative, when its cross terms are negative, so it cannot always
    be plotted as a distance.
- **Decision:** future improvement, named in the paper (B14).
- **Reason:** the recorded reason is the cost: $O(n\cdot 2^p)$, growing exponentially with the number
  of features (B14's paper sentence; `boundary-distance-and-robustness.md` §7). That argument does
  not hold for this value function: the closed form costs $O(p^2)$ per row, with no model call. The
  decision is unchanged until the user decides on it; the negative $\phi_j$ is a possible reason not
  to plot the decomposition as a per-feature distance.
- **Code:** `plot.bl_boundary()`, `bl_robustness()`.
- **Effect:** changes the importance values.
- **Effort:** recorded as medium; low with the closed form.
- **Source:** B14; `boundary-distance-and-robustness.md` §4, §7; `2 implementation_summary.txt`
  §4.6.1; methodology note §3.1, §5 (Table 5); Pima issue 9.

### 1.7 New results

#### 1.7.1 Aggregated Shapley values

- **What and why.** A global importance built from the local attributions over a representative set
  of data points (thesis §6.2, p.143). It needs a rotation and a Shapley computation for each point.
- **Decision:** future improvement, named in the paper (B15).
- **Reason:** not recorded.
- **Code:** a new function after `bl_shapley()`.
- **Effect:** a new output; no current result changes.
- **Effort:** high.
- **Source:** B15; methodology note §5 (Table 5).

#### 1.7.2 Calibrated probabilities from the boundary distance

- **What and why.** Platt scaling on the length of the change vector (Euclidean or Mahalanobis) would
  give a model-agnostic posterior class probability from the distance to the boundary (thesis §6.4,
  pp.145–146).
- **Decision:** future work in the boundary note; no reply covers it.
- **Reason:** not recorded.
- **Code:** a new function after `bl_find_boundary()`.
- **Effect:** a new output; no current result changes.
- **Effort:** medium.
- **Source:** `boundary-distance-and-robustness.md` §7.

#### 1.7.3 A Cholesky-whitened diagnostic plot

- **What and why.** A plot of the whitened decomposition ($u_k^2$) across the training set, without
  a claim of per-feature attribution. It is useful when the Mahalanobis check's percentage is large.
- **Decision:** future work in the boundary note; no reply covers it.
- **Reason:** it has not been requested (`boundary-distance-and-robustness.md` §7).
- **Code:** `plot.bl_boundary()`.
- **Effect:** a new display; no result changes.
- **Effort:** low.
- **Source:** `boundary-distance-and-robustness.md` §4.2, §7.

### 1.8 New capabilities

#### 1.8.1 Multi-class models

- **What and why.** The method explains binary classifiers. A multi-class model is explained one
  class at a time (the class of interest as 1, the others as 0), so the counterfactual leaves that
  class rather than moving towards a specific other class. One boundary per pair of classes, for
  example through linear discriminant analysis (`MASS::lda()` in `bl_wrap_model()`, or a
  `bl_fit_lda()`), would extend the method. The CVA biplot already separates more than two groups.
- **Decision:** future improvement, named in the paper as a key area of future research (B17).
- **Reason:** a multi-class model is "currently resolved by only targeting one of the multi classes
  as 1, the rest as 0, so counterfactual is out of the target class, not towards a specific class"
  (the reply on B17).
- **Code:** `bl_wrap_model()`, `.pred_function()`, the grid, both searches.
- **Effect:** a new capability; binary results do not change.
- **Effort:** high.
- **Source:** B17; methodology note §5 (Table 5); `2 implementation_summary.txt` §4.1.

#### 1.8.2 A CVA biplot for a singular W ($p>n$)

- **What and why.** **W** is singular for wide data. The generalised SVD of Ganey & Gardner-Lubbe
  (2026, *Statistical Methods & Applications* 35, 825–847; R package `wideRhino`) gives a CVA biplot
  in this case.
- **Decision:** future improvement, named in the paper (B16).
- **Reason:** not recorded.
- **Code:** `bl_build_projection()`.
- **Effect:** a new capability; no current result changes.
- **Effort:** medium to high.
- **Source:** B16; methodology note §5 (Table 5).

#### 1.8.3 A different transformation per feature

- **What and why.** `bl_set_scaling()` records one affine method for all features (each with its own
  centre and scale). A transformation per feature (z-score, min-max, log or Box-Cox, none) would help
  skewed or mixed-scale features. A non-affine inverse needs calibrated non-linear axes, such as the
  spline-based biplots of Groenen, le Roux & Gardner-Lubbe (2015, *Advances in Data Analysis and
  Classification* 9, 219–238; biplotEZ's spline axes, `PCOaxes = "splines"`), or a stated "affine
  only" limit. Display only.
- **Decision:** future improvement, named in the paper (B18). The reply names spline axes as the
  future area to incorporate.
- **Reason:** not recorded.
- **Code:** `bl_set_scaling()`, `.bl_rescale_biplot_axes()`, `.scale_to_raw()`, `bl_apply_scaling()`.
- **Effect:** none (display only).
- **Effort:** medium.
- **Source:** B18; methodology note §5 (Table 5); `scaling_axis_relabel_note.md`.

---

## 2. Code and interface

These change no result at the defaults.

### 2.1 Keeping every pair's plane and boundary point

- **What and why.** `bl_find_local_cf()` searches every pair, but keeps the plane, grid and boundary
  point of the winning pair only; the other pairs keep only their distances. Storing each pair's
  `Vr_rot`, `tVr_rot`, `B_x`, `B_z`, `B_pred` and distances (for example when `keep_pairs = TRUE`),
  with a pair's grid rebuilt on demand for plotting, would allow a figure of the search across
  planes, the alternative counterfactuals on their own planes, and §1.3.3's retry.
- **Decision:** future improvement, named in the paper (A2, option (b), sentence (ii)).
- **Reason:** nothing in Paper 1 needs the other planes (A2). Its figure 6 already shows the search as
  a dot chart of the per-pair distances.
- **Code:** the result of `bl_find_local_cf()`.
- **Effect:** none.
- **Effort:** low to medium.
- **Source:** A2; methodology note §5 (Table 5).

### 2.2 Smaller interface changes

#### 2.2.1 A configurable `"fixed"` tolerance

- **What and why.** `.fixed_tol` is 0.5 model units for every feature, and no argument changes it. In
  the raw-unit Pima model that is half a pregnancy or half a year; on a standardised model it is half
  a training SD (±3,001 for `loan_amnt` on the loan harness). So it suits features measured in whole
  units, but not features on a small scale. Too tight leaves no candidates; too loose makes the
  restore fail more often. A `fixed_tol` argument (one number, or a named vector per
  feature; possibly in raw units or as a fraction of the within-class SD, `sqrt(diag(W))`), with a
  default equal to 0.5 model units, would suit features measured on different scales.
- **Decision:** future improvement, named in the paper (A4, option (c)), against the recommendation
  to fix it now.
- **Reason:** not recorded.
- **Code:** `.fixed_tol`, `.apply_actionability()`, `bl_find_local_cf()`.
- **Effect:** none at the default.
- **Effort:** low.
- **Source:** A4; methodology note §3.7, §5 (Table 5).

#### 2.2.2 A per-row feasibility flag in the global search

- **What and why.** `bl_find_boundary()` keeps counterfactuals outside the training ranges and
  reports only their number (10 of 313 training rows and 5 of 79 test rows on the 4-feature Pima
  model). A flag per row would show which. An option to filter by the ranges before the selection
  would offer the thesis's rule (§3.10 step 4, p.93).
- **Decision:** future improvement, on its own, not in the paper (B11).
- **Reason:** on the seven-feature model that Paper 1 uses for Phase 2, no global counterfactual falls
  outside the training ranges (0 of 313 training and 0 of 79 test rows), so the paper needs no
  sentence on it (B11).
- **Code:** `bl_find_boundary()`.
- **Effect:** none at the default.
- **Effort:** low.
- **Source:** B11; methodology note §3.3, §5 (Table 5).

#### 2.2.3 Pair control

- **What and why.** `max_pairs = 1L` forces pair (1,2), but no other pair can be forced. A
  `pairs =` argument would force any pair. A `max_pairs` default tied to $\binom{p}{2}$ would make
  the default cover every pair for any $p$: the default of 21 covers every pair up to $p=7$, and from
  $p=8$ it leaves out every pair that involves the highest components.
- **Decision:** future improvement, not in the paper (B8).
- **Reason:** not recorded. The form notes that it changes nothing for $p\le7$, where the default
  already searches every pair.
- **Code:** `bl_find_local_cf()`, `.generate_pairs()`.
- **Effect:** none at the default for $p\le7$.
- **Effort:** low.
- **Source:** B8; rotation note §4.7.

#### 2.2.4 The pair label in output

- **What and why.** Report how far the selected plane lies from the eigenplane its pair names: the
  angle $\theta$, or the overlap. The angle does not predict the overlap: among 74 cases with
  $\theta<5°$ the overlap spans 0% to 100% (rotation note §4.2).
- **Decision:** further work (this item). B3 settled the wording only: the help
  pages and Paper 1 say "component pair", and that the selected plane need not coincide with the
  plane of the two components that name it.
- **Reason:** not recorded.
- **Code:** `print.bl_local_result()`, and `bl_find_local_cf()` to compute it.
- **Effect:** none (output only).
- **Effort:** low.
- **Source:** rotation note §4.2; B3.

#### 2.2.5 A switch for the empty-region rule

- **What and why.** With the hull on, the local search drops a closed boundary line that encloses no
  training data on the plane (`.encloses_data()`; thesis §3.10 step 4, p.93), and `set_filters()` has
  no switch to turn that part off on its own. A switch (for example
  `set_filters(empty_region = TRUE)`) would let a user show what the rule changes, as the Pima
  example's patient C does by replacing the helper inside the script
  (`without_empty_region_rule()` in `scripts/paper1-pima/00_setup.R`).
- **Decision:** a possible future development, listed with Pima issue 17 (2026-11-15).
- **Reason:** the current method is complete without one (Pima issue 17), and a switch in
  `set_filters()` would be a package change after tag v0.2.0, which Paper 1 pins (the review plan of
  2026-11-15, assumption 6).
- **Code:** `set_filters()`, `.filter_switches()`, `bl_find_local_cf()`.
- **Effect:** none at the default.
- **Effort:** low.
- **Source:** Pima issue 17 (2026-11-15); `local-convex-hull.md`; `algorithms.md` F2.

#### 2.2.6 `new_title` for `plot.bl_projection()`

- **What and why.** Every other biplot plot method accepts `new_title = NA` for a custom title
  (CLAUDE.md §5). `plot.bl_projection()` does not: its `...` is ignored. Adding it would make the
  plot methods uniform.
- **Decision:** future development (2026-11-15).
- **Reason:** adding it after tag v0.2.0 would change the API that Paper 1 pins.
- **Code:** `plot.bl_projection()`.
- **Effect:** none.
- **Effort:** low.
- **Source:** the review plan of 2026-11-15 (`.claude/plans/before-the-next-steps-temporal-storm.md`,
  WP4 and assumption 2); CLAUDE.md §5 (the `new_title` rule).

### 2.3 Behaviours of the code (`algorithms.md` §7)

On 2026-11-06 these three behaviours were recorded as notes, with no code change
(`.claude/reference/algorithms-findings-fixes-plan.md`). Each describes something that could be
developed.

#### 2.3.1 The global search's zero test (F1)

- **What and why.** `bl_find_boundary()` combines the two opposite-side choices by adding matrices
  with zeros for the side not used, then sets any row whose `B_z` is exactly $(0,0)$ to NA. A row
  with no opposite-side contour is already NA, so the test never marks a missing row; it can only
  turn a real boundary point at the biplot's origin into NA. When every surviving contour lies on one
  side, the rows predicted on that side get NA, the test is NA for them, and `if (any(zero_rows))`
  stops with "missing value where TRUE/FALSE needed". The empty result (`.empty_return()`) returns
  `B_z` as zeros, not NA. Recording the rows without an opposite-side contour explicitly, as NA, would
  let the search return a result in that case.
- **Decision:** keep unchanged; a comment at the test in `bl_find_boundary()` and `algorithms.md` F1
  record the stop (2026-11-06).
- **Reason:** not recorded.
- **Code:** `bl_find_boundary()` (the zero test, `.empty_return()`).
- **Effect:** none where the search runs now; the case that stops would return NA rows instead.
- **Effort:** low.
- **Source:** `algorithms.md` §7 F1; `.claude/reference/algorithms-findings-fixes-plan.md`.

#### 2.3.2 A Gini that does not depend on the row order (F5)

- **What and why.** `calc_gini()` sorts by score and does not average tied scores. With scores
  floored to 3 decimals, ties are common, so the Gini can change with the row order. Averaging over
  tied scores would make it independent of the order.
- **Decision:** accept, and note it (2026-11-06).
- **Reason:** not recorded.
- **Code:** `calc_gini()` (called by `bl_fit_model()` and `bl_wrap_model()`).
- **Effect:** the reported Gini only. No algorithm uses it: it is stored in `bl_model` and
  `bl_result`, and printed.
- **Effort:** low.
- **Source:** `algorithms.md` §7 F5; `.claude/reference/algorithms-findings-fixes-plan.md`.

#### 2.3.3 Caching coalition scores in the exact Shapley method (F7)

- **What and why.** The exact method scores each coalition $T$ again for every $S\cup\{j\}=T$, with no
  caching: for $p=7$, about 575 scorings for 128 coalitions. Scoring each coalition once would save
  time.
- **Decision:** note it; no change (2026-11-06).
- **Reason:** not recorded.
- **Code:** `.shapley_exact_one()`, `.build_rows_from_masks()`.
- **Effect:** none on the values; time only.
- **Effort:** low.
- **Source:** `algorithms.md` §7 F7; `.claude/reference/algorithms-findings-fixes-plan.md`.

### 2.4 biplotEZ display add-ons

- **What and why.** `plot_biplotEZ()` uses only a thin slice of biplotEZ's rendering API: `axes()` on
  a blank canvas (`samples(opacity = 0)`), with the package's own grid, points, contours and target
  marker drawn on top. biplotEZ offers more, none of it wired into `plot_biplotEZ()`,
  `plot.bl_local_result()` or `plot.bl_surrogate()`. The candidates, with the walkthrough's notes:

  | biplotEZ function | What it would add | Note |
  |---|---|---|
  | `means()` | Class or group mean markers and labels | Needs grouping by outcome group (TP/TN/FP/FN) to fit the package's colours |
  | `alpha.bags()` | Alpha-bag contours: the region holding alpha% of a class's points | A density-based summary region per class |
  | `ellipses()` | Kappa concentration ellipses per class | A parametric alternative to alpha-bags |
  | `density1D()` / `density2D()` | Kernel density curves or heatmaps over the biplot | An alternative or complement to the probability-coloured grid |
  | `legend.type()` | An in-plot legend for samples, means, bags, ellipses and regions | The package prints a text summary to the console instead |
  | `newsamples()` / `newaxes()` | Supplementary points or axes, styled apart from the primary biplot | Could set external data apart from the training and test overlays |
  | `interpolate()` | Projects a new dataset onto a fitted biplot object | A different mechanism from `bl_project_points()`, which uses `bl_result$V` and `X_center` directly |
  | `classification()` / `prediction()` | biplotEZ's own classification regions and prediction overlays | An alternative to the package's prediction grid and contours |
  | `fit.measures()` | Quality-of-fit diagnostics (overall quality, axis, sample and class predictivity) | Not visual, but could help decide which additions are worth making |

  The walkthrough rates three groups as not useful: `regress()` (an alternative construction, not an
  addition to the PCA/CVA biplot); `rotate()`, `reflect()` and `translate_axes()` (the package has
  its own rotation, `.apply_biplot_rotation()`, which keeps the grid, contours and target in step);
  and `CA()`, `CATPCA()`, `PCO()`, `AoD()`, `CLPs()`, `CLRs()` (for categorical or distance data, not
  continuous features).
- **Decision:** candidates for future additions in the Phase 1 walkthrough; no reply covers them.
- **Reason:** each would be new, scoped feature work. `plot_biplotEZ()` uses a thin slice of
  biplotEZ on purpose: the package needs full control over the outcome-group colours, the
  prediction-probability surface, the boundary contours and the synchronised rotation, and none of
  these map cleanly onto biplotEZ's group-driven sample styling.
- **Code:** `plot_biplotEZ()`, `plot.bl_local_result()`, `plot.bl_surrogate()`.
- **Effect:** display only; no result changes.
- **Effort:** not estimated.
- **Source:** Phase 1 walkthrough (`development/reference/review_section4_to_6.md`, "biplotEZ visual
  features").

---

## 3. Tests

The tests cover the code paths and the display rules. These would pin computed values, or check
them on other platforms.

### 3.1 Values of computed results

- **What and why.** Tests of the values of the `bl_find_local_cf()` boundary point, the Shapley
  values on a real model, `bl_robustness()`, `bl_surrogate()` and `bl_build_result()`. The local
  plane depends on the linear-algebra routine (§1.1), so values need a tolerance and a stated
  platform, or §1.1 first. The design agreed before the tests were deferred (B21, option (a)):
  - property tests everywhere: the counterfactual crosses the cutoff, meets every constraint and the
    training ranges, its distances recompute from its values, and the chosen pair minimises the
    distance;
  - exact values for everything that does not depend on the plane;
  - exact values of the local search on Windows only (`skip_on_os()`), with a comment that points to
    §1.1;
  - the data: iris as two classes, and one test on the bundled Pima file.
- **Decision:** future improvement, not in the paper; it goes with CRAN preparation (§4.2) (B21,
  Round 4).
- **Reason:** the user's instruction (Round 4): "I don't want the tests. It must only work on my
  computer. There can be a note in the paper that this may be a problem." The tests are needed only
  if the package is used beyond that computer.
- **Code:** `tests/testthat/`.
- **Effect:** none.
- **Effort:** not estimated.
- **Source:** B21.

### 3.2 The CVA identity

- **What and why.** A test that $d^2_M=(n-g)\,d_Z^2$ for every pair, including two-group CVA
  (rotation note §3.4). The identity follows from biplotEZ's CVA scaling,
  $\mathbf{V}'\mathbf{W}\mathbf{V}=\mathbf{I}/(n-g)$; the per-pair rotation cancels, so it holds
  for the current rotation and for the minimal one alike. `17_3` §11 asserts it; `tests/testthat/`
  does not.
- **Decision:** future improvement, not in the paper, with §3.1 (B5, Round 4).
- **Reason:** as §3.1: no tests now, since the package only has to work on the user's computer
  (Round 4).
- **Code:** `tests/testthat/`.
- **Effect:** none.
- **Effort:** low (about ten lines).
- **Source:** B5; rotation note §3.4.

### 3.3 Cross-platform reproducibility

- **What and why.** Re-run the Pima pipeline (`run_all.R`, `tests.R`), the loan harness and `17_3` on
  Linux and macOS, and on other BLAS/LAPACK builds (for example OpenBLAS), for example through GitHub
  Actions, and compare `numbers.rds`. Because of §1.1, small numerical differences matter: ordering
  the CVA groups as the thesis code did changes **V** by no more than 5.0e-17, and the Pima TP target's
  run 1 then has no counterfactual (no whole number of pregnancies near the new search point crosses),
  against pair (3,4) at 1.82 with the package's order (methodology note §3.12 and Table 4,
  measurement K3).
- **Decision:** future improvement, named in the paper, with a mitigation (B2, option (a)): the paper
  states the platform (R 4.6.0 on Windows) and that the local plane is under-determined, and the
  supplement includes `sessionInfo()`.
- **Reason:** the user's instruction (Round 4): the package only has to work on the user's computer,
  with a note in the paper that this may be a problem. A check was also not practical (B2, Round 2):
  no Linux distribution is installed under WSL on this machine, and `scripts/` is gitignored, so
  GitHub Actions cannot run the Pima scripts until they are in a repository.
- **Code:** none (a check, not a change).
- **Effect:** none; it tests whether the results reproduce.
- **Effort:** low: about half a session for one GitHub Actions run, once the scripts are in the
  paper repository (B2, option (b)).
- **Source:** B2; rotation note §4.1; methodology note §3.12, §5 (Table 5).

### 3.4 The acceptance tests of §1.1's proposal

- **What and why.** Seven acceptance tests decide whether §1.1's minimal rotation can replace the
  current one. `17_3_minimal-rotation-verification.R` runs them all; tests 1, 2, 3 and 6 are asserted,
  so the script fails if they stop holding. Tests 4, 5 and 7 are incomplete:

  | # | Test | Runs in | Status |
  |---|---|---|---|
  | 1 | **Exactness:** $x_t \in \operatorname{span}(\mathbf{V}_r)$ to about 1e-12 for every candidate | `17_3` §2 | Asserted; residuals 1.1e-15 to 4.2e-15, all six pairs |
  | 2 | **Determinism:** identical across at least 5 calls and across BLAS/LAPACK builds | `17_3` §3 | Asserted bit-identical across 5 calls; across builds not yet run (§3.3) |
  | 3 | **Minimality / label integrity:** principal angles $(0,\theta)$ exactly | `17_3` §2 | Asserted; the first angle 0.00e+00 and the second $\theta$, all six pairs |
  | 4 | **Monotonicity:** the family optimum is at most the current distance, per target | `17_3` §7, §8 | Asserted and holding at both targets, but the family figure comes from the stage-1 ray scan and the package figure from its own 200×200 grid. Needs stage 2 re-run on the winning plane for a like-for-like number, and a sweep over all rows to show what the label sample costs |
  | 5 | **Variance retained:** adequacy, axis predictivity, RSS at $\mathbf{J}=\operatorname{diag}(1,1,0,\ldots)$, minimal against current | — | Not run. Reuse the diagnostic block of `17_1` §13 / `17_2` §11 |
  | 6 | **Thesis reproduction:** the 28.14° at $p=2$, $r=1$ | `17_3` §10 | Asserted; 28.1411°, matching Procrustes to 6.75e-14 |
  | 7 | **Cost** | `17_3` §9 | Measured at $p=3$ and reported, not asserted (wall clock depends on the machine); $p\ge4$ not yet timed |

  `17_3` also asserts two properties of the current package: within-platform determinism (§5, §8b;
  rotation note A.2) and the CVA identity (§11; §3.2).
- **Decision:** future work, with §1.1 (B1: "The acceptance tests 4, 5 and 7 (§6.1) are
  incomplete").
- **Reason:** they are needed only before §1.1 replaces the current rotation, and §1.1 follows
  Paper 1 (B1).
- **Code:** `biplot-rotation/17_3_minimal-rotation-verification.R`.
- **Effect:** none.
- **Effort:** not estimated.
- **Source:** B1; rotation note §5.

---

## 4. Documentation and release

### 4.1 A vignette on the mathematical method

- **What and why.** A vignette that derives the biplot projection and the counterfactual geometry
  from first principles.
- **Decision:** future improvement, not in the paper (B22).
- **Reason:** not recorded.
- **Code:** `vignettes/`.
- **Effect:** none.
- **Effort:** not estimated.
- **Source:** B22.

### 4.2 CRAN preparation

- **What and why.** The checks, examples within the time limits, and the submission, with §3.1's
  tests.
- **Decision:** future improvement, not in the paper (B23).
- **Reason:** not recorded at the decision (Round 2). Round 4 records that the package only has to
  work on the user's computer, and puts the tests (§3.1, §3.2) with CRAN preparation.
- **Code:** the whole package.
- **Effect:** none.
- **Effort:** not estimated.
- **Source:** B23.

### 4.3 DESCRIPTION and `.Rbuildignore` at the next version

- **What and why.**
  - DESCRIPTION has no `URL` or `BugReports` field. Both could point to the GitHub repository
    (`AI-Rowan/boundary_logic`, and its issues page), so that the help index, pkgdown and a CRAN page
    link to the source and to the place to report problems. Its `Date` (2026-11-05) is also stale.
  - `.Rbuildignore` has a redundant line (`^scripts/pima\ diabetes\.csv`, already covered by
    `^scripts`) and lacks one; the review plan does not name the missing line
    (`outputs_dev_oos/`, a gitignored folder at the root, has none).
- **Decision:** not edited now; they go to the next version (the review plan of 2026-11-15,
  assumption 4).
- **Reason:** DESCRIPTION is not edited after tag v0.2.0, which Paper 1 pins: "its date and the
  missing `URL` go to the next version" (assumption 4).
- **Code:** DESCRIPTION, `.Rbuildignore`.
- **Effect:** none.
- **Effort:** low.
- **Source:** the review plan of 2026-11-15 (`.claude/plans/before-the-next-steps-temporal-storm.md`,
  WP4, assumption 4 and Appendix A).

### 4.4 A help page for the bundled Pima data

- **What and why.** The README's Pima example reads `inst/extdata/pima_diabetes.csv` (768 rows; the
  class column `Outcome`), but `R/datasets.R` documents only `loan_data`. A help page would give the
  data's source, its columns, and the zeros that mark missing values (the README keeps complete rows
  only).
- **Decision:** future development (2026-11-15).
- **Reason:** not recorded.
- **Code:** `R/datasets.R`, `man/`.
- **Effect:** none.
- **Effort:** low.
- **Source:** the review plan of 2026-11-15 (WP4); README ("Bringing your own model").

### 4.5 The plot display arguments

- **What and why.** No note explains the display arguments of the plot methods; only their roxygen
  help and the walkthroughs do:
  - `ticks_v`, `ticks_var`, `ticks_n` (tick marks per axis) and `label_offset_var`,
    `label_offset_dist` (labels moved outward), in `plot_biplotEZ()`, `plot.bl_local_result()` and
    `plot.bl_surrogate()`, with their pairing rule (each distance or count pairs with its variable by
    position; `.make_label_line_vec()`, `.make_ticks_vec()`);
  - `grid_col` (a two-colour grid split at the model's cutoff) and `confusion_cols` (outcome-group
    colours or predicted-class colours), in `plot_biplotEZ()`;
  - `show_arrows`, in `plot.bl_local_result()`, `plot.bl_sparse_result()` and `bl_pick_point()`;
  - `no_grid` and `plot_points`, in `plot_biplotEZ()` and `plot.bl_local_result()` (where
    `plot_points` has different defaults).

  A section of a note, for example next to the "Display conventions" of
  `scaling_axis_relabel_note.md`, would explain them, with the advice for figures (CLAUDE.md §5:
  `ticks_v = 4L`, and `ticks_var`/`ticks_n` for a short axis).
- **Decision:** future development (2026-11-15).
- **Reason:** not recorded.
- **Code:** none (documentation); the roxygen is in `R/plot_biplot.R`, `R/local_cf.R`,
  `R/surrogate.R`, `R/shapley.R` and `R/pick_point.R`.
- **Effect:** none.
- **Effort:** low.
- **Source:** the review plan of 2026-11-15 (WP4).

---

# Appendix: The replies of 2026-11-05

**Status (2026-11-16).** This is the reply form "Future Improvements Review: Replies on the Open
Items", written and answered on 2026-11-05. On 2026-11-16 it was moved here from
`documentation/methodology-changes/future-improvements-review.md`, and its code line numbers were
replaced by function and section names; its wording is otherwise unchanged. It is the dated record
of the decisions; Part 1 holds the current list. Later work superseded these statements:
- **The local hull (Round 4's summary and B10).** Round 4 records the hull on each rotated plane as
  built from the training data and the target, with no enclosure rule for closed contours, and
  switched by `hull =` in `bl_find_local_cf()`. Since 2026-11-06 it is built from the training data
  alone. A contour that is still closed after the clip is kept only if it encloses training data on
  the plane (the empty-region rule, `.encloses_data()`). The hull's switch is `set_filters(hull =)`,
  not an argument of `bl_find_local_cf()` (methodology note §3.6; `documentation/local-convex-hull.md`).
- **B2's measurement K3.** Its TP numbers (pair (2,3) to (1,3), distance 1.42 to 1.70) predate the
  hull. With the hull, ordering the CVA groups as the thesis code did leaves the TP target with no
  counterfactual, against pair (3,4) at 1.82 with the package's order (methodology note §3.12 and
  Table 4).
- **A5 and A6.** "Only TP run 4 has a counterfactual" predates the hull: with the hull, TP run 4 has
  none either. Of the six runs with a fixed variable for the TP and TN targets, only the TN target's
  run 7 (added later) has one (methodology note §3.7).

**Author of the replies:** Adriaan Rowan
**Prepared by:** Claude, from the package's notes as of commit `12b38e4` (branch `paper_nr1_changes`)
**Source of Part A:** the former `documentation/known-issues.md` (removed on 2026-11-05 once every
item was decided; in git history at `c13fe44`).
**Companion documents** (paths are relative to the package root):
- `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md`, §4 and §6: the rotation's
  problems and further work;
- `documentation/fixed-constraint-restore.md`, §8: the follow-up searches for `"fixed"` and
  `"integer"`;
- `documentation/methodology-changes/methodology-changes.md`: §5 (future work); the interpretation
  points are in §3.6 and §3.12;
- Part 1 of this file: the future developments, formerly CLAUDE.md §9's deferred work;
- `.claude/reference/future-improvements-hull-plan.md`: the plan behind this form.

The future items of these notes are now in Part 1.

---

## Purpose

This form collects your decision on every open item, so that none is left as an "issue". Each item
ends in one of three ways:

- **Fix now:** implemented before Paper 1. Each item says whether this changes Paper 1's numbers.
- **Future improvement:** the package is unchanged, and the paper describes the item as an area of
  future improvement.
- **No action:** the item is dropped, and is not in the paper.

On 2026-11-04 you decided that none of the known issues would be implemented for Paper 1. This form
lets you revisit that item by item.

Once every item is resolved and the fixes are implemented, the methodology note
(`methodology-changes.md`) becomes the single record of the design decisions from the thesis to the
package. Its section on future work is then rebuilt from the items you mark "Future improvement",
with the paper sentences you agree here.

## How to reply

1. For each item in Part A, tick one decision box (change `[ ]` to `[x]`).
2. Write your recommended solution in the quote block under **Your solution**. A few words are
   enough if you agree with an option as written.
3. Where a draft paper sentence is given, tick "OK as drafted" or write your edit.
4. For Part B, write `Now`, `Future` or `Drop` in the **Your reply** column. Longer comments go in
   the block below the table.
5. Save the file and tell me.

I will then read your replies and add follow-up questions under each item's **Follow-up** heading.
The rounds continue until every item's **Status** reads "Resolved".

**Notes on the numbers.**
- All Pima numbers come from the Paper 1 example: `scripts/paper1-pima/output/table1.csv`, and
  `documentation/methodology-changes/method-change-effects.csv`.
- Phase 3 runs on the GAM refitted on four features (*Pregnancies*, *Glucose*, *Insulin*, *Age*).
- √d²_M is the Mahalanobis distance from the target.
- Code is cited by function (the form's line numbers of 2026-11-05 were replaced on 2026-11-16).

**The runs in Table 1:**

| Run | Constraints |
|---|---|
| 1 | unconstrained (*Pregnancies* `"integer"` only) |
| 3 | *Pregnancies* and *Age* fixed |
| 4 | *Pregnancies* fixed, *Age* cannot decrease |
| 5 | run 1 searched on pair (1,2) only |
| 6 | *Glucose* and *Insulin* cannot increase |

Run 2 is not in Table 1.

---

## Summary

| # | Improvement area | Fixing it moves Paper 1's numbers? | Effort | My recommendation | Your decision |
|---|---|---|---|---|---|
| A1 | Measuring how far the 2D boundary lies from the full boundary | No (it adds a result) | Medium | Fix now (an experiment in the scripts) | Fix now (a) |
| A2 | Keeping every pair's plane and boundary point | No | Low to medium | Future improvement | Future (b) |
| A3 | One colour scheme for every biplot | No | Low | Fix now | Fix now (a) |
| A4 | A configurable tolerance for fixed variables | No, at the default | Low | Fix now | Future (c) |
| A5 | A follow-up search when a fixed variable cannot be restored | Yes: TN runs 3 and 4 | Medium | Future improvement | Future (b) |
| A6 | Searching in the actionable variables only | Yes: every fixed run | High | Future improvement | Future (b) |
| A7 | Letting the free variables compensate after integer rounding | Yes: TN runs 1 and 5 | Medium | Future improvement | Future (b) |
| A8 | A README constraint example that finds a counterfactual | No | Low | Fix now, before the release tag | Fix now (a) |
| A9 | The loan harness's biplot title | No | Low | Fix now | Fix now (a) |
| B1–B24 | The other open items (Part B) | see Part B | | see Part B | see Part B |

---

## Round 4 (2026-11-05): the final decisions

Every item now has a final decision. Round 3's answers were given in the chat (AskUserQuestion)
and are recorded under each item's follow-up. Round 4 changed two things:
- **No tests now.** The package only has to work on your computer, so B21 (tests of computed
  values) and B5 (the CVA identity test) are future work, in the design document only, with B23.
  B2's paper sentence notes the risk (revised in Part B's Round 4 note).
- **B10 is done now.** You decided the methodology note's §3.6 yourself: the local search keeps
  only boundary points inside the convex hull of the training data and the target on each rotated
  plane, as the thesis requires (p.117). This changes Paper 1's Phase 3 numbers. The supervisor has
  not seen the methodology note yet; the version they will receive presents §3.6 as restored.

The implementation follows the plan `.claude/reference/future-improvements-hull-plan.md` (approved
2026-11-05).

| Step | Items | What it does | Moves Paper 1's numbers? |
|---|---|---|---|
| 1. The hull on each rotated plane | B10 | `hull = TRUE` in `bl_find_local_cf()`; `hull = FALSE` gives the earlier results | Yes: Phase 3 |
| 2. Small package fixes | A3, A8, A9, B3, B9 | One colour map; the README constraint example (run 4, with its output); the harness title; the pair wording; a warning for a target at the centroid | No |
| 3. Re-runs | | The loan harness, the Pima pipeline, the methodology effects (a new M6 for the hull) and the notes | |
| 4. The approximation-gap experiment | A1 | `scripts/paper1-pima/05_full-space-gap.R`, local search only, with the hull checks | No; it adds a result |
| 5. Release housekeeping | B24 | Version 0.2.0, the final thesis page citations, the pkgdown site; merge and tag on your instruction | No |

**Paper only, no package change:** the agreed sentences, grouped as in G1, and B2's note (the
platform statement in the paper, and `sessionInfo()` in the supplement).

**Future improvement, in the paper:** A2, A4, A5, A6, A7, B1, B2, B6, B12, B13, B14, B15, B16, B17
and B18 (B4 goes with B1). A2, B7 and B10 are described in §5 of the paper, as the method stands.

**Future improvement, in the design document only:** B5, B8, B11, B21, B22 and B23.

**No action:** B19 and B20.

### Implemented (2026-11-05)

Steps 1 to 5 are done; nothing is committed. `devtools::test()` gives 465 PASS, 0 FAIL; the Pima
`tests.R` passes 13 of 13; `devtools::check()` gives 0 errors, 0 warnings, 0 notes.

**What moved in Paper 1 (the hull, B10):**
- the TN target: nothing;
- the TP target, run 1: pair (2,3) becomes (3,4); search distance 1.50 to 1.84; counterfactual 1.42 to
  1.82. It now lowers every variable (*Age* 51 to 41.2) instead of raising *Age* to 59.6;
- the TP target, runs 4 and 5: no counterfactual now (figure 7(a) says so in its panel);
- run 7 (new): a counterfactual for the TN target with *Pregnancies* held fixed;
- the loan harness: unchanged.

**Sentences to look at again:**
- **A1:** the measured result and a filled-in sentence are under A1's follow-up.
- **A5:** Table 1's empty runs are now TP runs 3, 4, 5 and 7 and TN runs 3 and 4. The sentence itself
  holds.
- **Pima issue 11 (interpretation):** the "older patient" point no longer holds. Your choice (round
  5, in the chat): "Without the plausibility rule, the nearest counterfactual for the true positive
  raised Age from 51 to 59.6, outside the training data; inside the convex hull every counterfactual
  lowers Age, which is not actionable, so the constrained runs report none."
- **A1:** your choice (round 5) is the upper-bound sentence, recorded under A1's follow-up.
- The A7 numbers (TN, 0.14 to 0.42), B14's 27% and the B2 sentence are unchanged.

**Round 5 (2026-11-05), in the chat:**
- **The methodology note** has no approval boxes. It is organised as "Changes made" (history, core
  methodology, code improvements) and "Future work". The interpretation points sit in the sections
  they concern: issues 5 and 11 in §3.6, and the approximation gap as a third property in §3.12.
- **Outlier removal (the note's §3.11):** the code default of `bl_prepare_data(hull_fraction = )` is
  now 1, no removal. Each application keeps its own setting: 0.9 for the loan harness and vignette,
  1 for Pima. The help page explains how to use it. Tests: 465 PASS, 0 FAIL.

---

## Part A: the known issues

### A1. Measuring how far the 2D boundary lies from the full boundary

**Status:** Resolved: fix now (a), the local search only, with the A1.6 hull checks (plan 1)

**Today.**
- The searches find the decision boundary as a contour of the model's predictions on a grid in a 2D
  plane, not in the full p-dimensional space.
- The validity of a counterfactual is not in doubt: every local counterfactual is re-scored by the
  model and must cross the cutoff.
- Its distance is less certain. A counterfactual is the nearest boundary point *within the planes
  searched*, which may be further away than the nearest point of the full boundary.
- Paper 1's distances are therefore upper bounds of unknown tightness. Nobody has measured the gap.

**Where.** The global contour is built in `bl_build_grid()` (used by `bl_find_boundary()`). The local
contour is built per pair in the pair loop of `bl_find_local_cf()`.

**Effect on Paper 1.** None on the numbers. The paper can only say the distances are upper bounds.
Paper 1 also has no comparison with other methods (paper plan §9), so this is a likely reviewer
question.

**Options.**
- **(a) Measure it now, in the scripts only.** No package change.
  - For the 79 test rows of the four-feature model, find the smallest √d²_M in the full space for
    which the GAM crosses 0.5, within the training ranges. A multi-start constrained optimiser, for
    example `optim()` with a penalty or `nloptr`, would do this.
  - Compare it with the package's search point, taken before the integer rounding so that rounding
    does not confound the ratio.
  - Report the ratio of the package's distance to the full-space distance: its median, its range,
    and how often it is close to 1.
  - Moves Paper 1's numbers: no; it adds one result (a sentence or a small table in §8).
  - Effort: medium, about one session.
- **(b) Future improvement.** The paper states that the distances are upper bounds within the
  searched planes. Moves Paper 1's numbers: no. Effort: none.

**My recommendation: (a).**
- With p = 4 the full-space search is cheap.
- A measured ratio turns a caveat into a result, and goes some way to the missing comparison.
- One caution: the full-space optimum ignores the plausibility rules on the plane (the hull), so it
  can be nearer for reasons that are not a fault of the plane.

**Draft paper sentence.**
- If (b): "The counterfactual is the nearest boundary point within the planes that are searched, so
  its distance is an upper bound on the distance to the nearest point of the decision boundary in
  the full $p$-dimensional space. Measuring this gap, for example against a constrained search in
  the full space, is left for future work."
- If (a), the sentence reports the result: "On the Pima test data, the distance to the biplot
  counterfactual was [median ratio] times the distance found by a constrained search in the full
  four-dimensional space ([range])."

**Your reply**

- [ x] Fix now: option (a)
- [ ] Future improvement: option (b)
- [ ] No action
- [ ] Other (describe below)

Your solution:

>

Draft paper sentence: [ x] OK as drafted  [ ] Edit below

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: fix now, as an experiment in the scripts; the package does not
change. Six design points follow. My default is the first box in each: tick it, or tick another
and add a line.

**A1.1 Which searches are measured?**
- [ ] The local and the global search (my default).
  - The local search: the four-variable GAM, with each of the 79 test rows as the target, under the
    run 1 settings.
  - The global search: Phase 2's seven-variable GAM on the one global plane, for the same 79 rows.
  - The same optimiser serves both. The global gap is likely to be larger, since one plane is
    searched instead of six, and Phase 2's importance values rest on those distances.
  - The paper sentence then gives two ratios.
- [ ] The local search only, which gives the paper's counterfactuals.

**A1.2 What the full-space point must satisfy.**
- It must cross 0.5 from the target's side, and stay within the training ranges, as the package's
  feasibility filter does. No integer rounding, and no hull.
- The package's point is taken before the integer rounding: for the TN target, the search distance
  of 0.14, not the counterfactual's 0.42. The ratio then measures the 2D approximation alone.
- The package's point is also one of the optimiser's starts, so the ratio cannot fall below 1.

- [ ] Agree
- [ ] Change (write below)

**A1.3 The optimiser.** `nloptr`, which is installed; only the script depends on it. It minimises
d²_M subject to f(x) = 0.5, with the training ranges as bounds. Its starts are the package's point
and the boundary crossings found along 50 random directions from the target, with a fixed seed.

- [ ] Agree
- [ ] Change (write below)

**A1.4 Where the result goes.**
- A new script, `scripts/paper1-pima/05_full-space-gap.R`, run after `run_all.R`.
- Its output in `output/full_space_gap.csv`, and the summary in `numbers.rds`.
- A section on `review.qmd`, one sentence or a small table in the paper's §8, and the approximation
  entry of the final design document.

- [ ] Agree
- [ ] Change (write below)

**A1.5 The paper sentence.** The draft you approved names "the biplot counterfactual", but the
comparison uses the boundary point before rounding. Revised:

> On the Pima test data, the distance to the boundary point found on the biplot planes was [median
> ratio] times the distance found by a constrained search in the full four-dimensional space
> ([range]).

The numbers go in after the run. If the result is surprising, for example a median ratio well above
1, I bring it back to you before writing it.

- [ ] OK
- [ ] Edit below

**A1.6 A plausibility check in the same script (optional; linked to B10).**
- Whether each of the package's local boundary points lies inside the convex hull of the training
  data on its own rotated plane. This puts a number on the methodology note's §3.6, which today
  says "not measured", before your supervisor decides on it.
- Whether each full-space point lies inside the convex hull of the training data in the full space
  (with the `geometry` package). The ratio is then also reported for the rows whose full-space point
  is plausible, which answers the caution under my recommendation above.

- [ ] Add both (my default)
- [ ] Add the first only
- [ ] Leave out

Your replies:

>

**Round 3 (2026-11-05, in the chat).** The local search only; otherwise my defaults.
- A1.1: the local search only. The global (Phase 2) search is not measured, so the paper sentence
  gives one ratio.
- A1.2, A1.3 and A1.4: agreed as written.
- A1.5: the revised sentence, with the numbers filled in after the run.
- A1.6: both hull checks.

**Result (2026-11-05, `scripts/paper1-pima/05_full-space-gap.R`, 79 test rows, 50 random directions
per row).** Every row has both a package point and a full-space point, and no ratio falls below 1.

| | Median ratio | Quartiles | Range | At or below 1.05 |
|---|---|---|---|---|
| With the hull (the method) | 1.19 | 1.09 to 1.31 | 1.00 to 2.75 | 15% |
| Without the hull | 1.18 | 1.08 to 1.26 | 1.00 to 2.45 | 16% |

- By outcome group, with the hull: FN 1.13, FP 1.28, TN 1.16, TP 1.24.
- The TP target: 1.84 on the planes against 1.21 in the full space (ratio 1.52; 1.24 without the
  hull). The TN target: 0.140 against 0.134 (ratio 1.05).
- 59 of the 79 full-space points lie inside the convex hull of the training data in four dimensions;
  their median ratio is 1.18, so the gap is not explained by implausible full-space points.

**The paper sentence waits for you.** The median search distance on the biplot planes is about 19%
longer than the nearest boundary point in the full space, and for 15% of the rows it is within 5%.
The revised A1.5 sentence would read: "On the Pima test data, the distance to the boundary point found
on the biplot planes was a median 1.19 times the distance found by a constrained search in the full
four-dimensional space (quartiles 1.09 and 1.31)." Please confirm or edit it.

- [ ] OK, use this sentence
- [ ] Edit below

>

**Round 5 (2026-11-05, in the chat).** You chose the upper-bound wording: "The biplot search finds a
boundary point a median 19% further from the target than the nearest boundary point in the full
four-dimensional space, so the distances it reports are upper bounds (quartiles 9% and 31%
further)." Status: Resolved.

---

### A2. Keeping every pair's plane and boundary point

**Status:** Resolved: future improvement (b), sentence (ii)

**Today.**
- `bl_find_local_cf()` searches every pair of components, but keeps the plane, grid and boundary
  point of the winning pair only. The other pairs keep only their distances.
- So another pair's plane cannot be drawn, or its counterfactual reported, without a second search
  pinned to it. Run 5 (pair (1,2) only) exists for this reason.
- The `"fixed"` restore and the `"integer"` rounding are tried on the winning pair's point only. A
  target can therefore lose a counterfactual that another pair's point would have given (see B12).

**Where.** The per-pair distances are recorded in the pair loop of `bl_find_local_cf()`
(`all_distances_mahalanobis`). `best_result` overwrites the plane and the point in the same loop.

**Effect on Paper 1.** None on the numbers. Figure 6 already shows the search as a dot chart of the
per-pair distances.

**Options.**
- **(a) Store each pair's plane and boundary point.**
  - Keep `Vr_rot`, `tVr_rot`, `B_x`, `B_z`, `B_pred` and the distances for each pair, in a list
    returned when `keep_pairs = TRUE` (default `FALSE`).
  - Rebuild a pair's grid on demand for plotting.
  - Moves Paper 1's numbers: no. Effort: low to medium.
  - It enables B12, and a figure of the losing planes.
- **(b) Future improvement.** Moves Paper 1's numbers: no.

**My recommendation: (b).** Nothing in Paper 1 needs the other planes. Choose (a) if you want B12
(retry the restore on the next-best pair) for Paper 1.

**Draft paper sentence.** "Only the plane of the selected pair is retained in order to be visualised if required"

**Your reply**

- [ ] Fix now: option (a)
- [X ] Future improvement: option (b)
- [ ] No action
- [ ] Other (describe below)

Your solution:

>

Draft paper sentence: [x ] OK as drafted  [ ] Edit below

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: future improvement. Your sentence states what the package does
now. For the paper, I suggest the JCGS spelling ("visualized"). So that it also reads as future
work, I suggest a second sentence that names the improvement:

- [ ] (i) Your sentence, with the spelling only: "Only the plane of the selected pair is retained,
  so that it can be visualized if required."
- [ ] (ii) Your sentence plus the improvement (my recommendation, since it also covers B12): "Only
  the plane of the selected pair is retained, so that it can be visualized if required. Retaining
  the plane and boundary point of every pair would allow the alternative counterfactuals to be
  compared on their own planes, and the restore of a fixed variable to be retried on the next-best
  pair."
- [ ] Other (write below)

>

**Round 3 (2026-11-05, in the chat).** Sentence (ii).

---

### A3. One colour scheme for every biplot

**Status:** Resolved: fix now (a)

**Today.** `plot.bl_projection()` colours the outcome groups TP green, TN steel blue, FP orange and
FN red. Every other biplot uses the thesis scheme: TP red, TN blue, FP purple, FN orange. So red
means FN in one plot and TP in the others.

**Where.** The odd one out was the colour map set in `bl_build_projection()`, which
`plot.bl_projection()` uses. The thesis scheme was in `plot.bl_boundary()`, `bl_project_points()`
and `plot.bl_local_result()`. Since A3's fix all four use `.confusion_cols` (`R/utils.R`).

**Effect on Paper 1.** None. Paper 1 does not use `plot.bl_projection()`.

**Options.**
- **(a) One shared colour map** (a private helper) used by every plot method, on the thesis scheme.
  Moves Paper 1's numbers: no. Effort: low.
- **(b) Future improvement.**
- **(c) No action.**

**My recommendation: (a).** It is a small change, and the thesis scheme is the documented
convention.

**Draft paper sentence.** None; this is not a paper item.

**Your reply**

- [ X] Fix now: option (a)
- [ ] Future improvement: option (b)
- [ ] No action: option (c)
- [ ] Other (describe below)

Your solution:

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: fix now. No question.
- The three places that use the thesis scheme already agree exactly. They use R's named colours:
  TP `"red"`, TN `"blue"`, FP `"purple"` and FN `"orange"`.
- Without true labels, points are coloured `"blue"` or `"red"` by predicted class. A point with no
  label is `"grey50"`.
- One private helper will hold these colours, and `plot.bl_projection()` will use it instead of its
  green, steel blue, orange and red scheme.
- The green target and counterfactual markers do not change.

---

### A4. A configurable tolerance for fixed variables

**Status:** Resolved: future improvement (c), sentence agreed

**Today.**
- The local search accepts a `"fixed"` variable within ±0.5 of its observed value, in model units,
  for every variable. There is no argument to change it.
- Afterwards, the restore sets the variable back to its observed value exactly.
- The meaning of 0.5 depends on each variable's scale. In the raw-unit Pima model it is half a
  pregnancy or half a year. On the standardised loan harness it is half a training SD (±3,001 for
  `loan_amnt`).

**Where.** `.fixed_tol <- 0.5` (`R/local_cf.R`).

**Effect on Paper 1.**
- **Too tight leaves no candidates.** TP run 3 (*Pregnancies* and *Age* fixed) has no boundary point
  within the tolerance on any of the 6 pairs.
- **Too loose makes the restore fail more often** (A5).

**Options.**
- **(a) A `fixed_tol` argument to `bl_find_local_cf()`.** A single number, or a named vector per
  variable, in model units, with default 0.5. Moves Paper 1's numbers: no, at the default. Effort:
  low.
- **(b) As (a), with the tolerance in raw units or as a fraction of the within-group SD**
  (`sqrt(diag(W))`). The default must still equal 0.5 model units, or results move. Effort: low to
  medium.
- **(c) Future improvement.**

**My recommendation: (a).** The restore keeps every counterfactual exact whatever the tolerance, so
exposing it only gives the user control over how many candidates the search keeps.

**Draft paper sentence.**
- If (a) or (b): "The search accepts a fixed variable within a tolerance of its observed value, 0.5
  in the units of the model by default and adjustable per variable, and the reported counterfactual
  restores the observed value exactly."
- If (c): "The tolerance for a fixed variable is 0.5 in the units of the model, which is half a
  pregnancy or half a year for the Pima data. A tolerance set per variable, for example as a fraction
  of its within-group standard deviation, would suit variables measured on different scales."

**Your reply**

- [ ] Fix now: option (a)
- [ ] Fix now: option (b)
- [ X] Future improvement: option (c)
- [ ] No action
- [ ] Other (describe below)

Your solution:

>

Draft paper sentence: [ x] OK as drafted  [ ] Edit below

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: future improvement, with the (c) sentence. No question. The
sentence's first half describes the method as it is, so it goes in the paper's description of the
constraints. Its second half goes in the future-work paragraph on constraints (see G1).

---

### A5. A follow-up search when a fixed variable cannot be restored

**Status:** Resolved: future improvement (b), sentence as drafted

**Today.**
- `"fixed"` is strict (since 2026-11-04). When the search's point, with its fixed variables set
  back to their observed values, no longer crosses the boundary, no counterfactual is reported.
- This happens even when a counterfactual that holds the fixed variables exists.

**Where.** `.restore_fixed()` (`R/local_cf.R`). Also `documentation/fixed-constraint-restore.md`
§4 and §8 (first bullet).

**Effect on Paper 1.**
- The TN target's runs 3 and 4 find boundary points, but the restored points score 0.475 and 0.487,
  so neither run has a counterfactual.
- With TP run 3 (no boundary point, A4), three of the four runs with a fixed variable are empty.
  Only TP run 4 has a counterfactual (Pima issue 12).

**Options.**
- **(a) A follow-up search.**
  - From the restored point, extend the move in the free variables along the search's own direction
    until the model crosses 0.5, within `train_ranges` and the other constraints.
  - It is a 1-D, deterministic search. Its point lies off the plane, as a restored point already
    does.
  - Moves Paper 1's numbers: **yes**. TN runs 3 and 4 may gain counterfactuals, which changes
    Table 1, figure 7, the review page and the methodology note.
  - It is a new method change, so it needs your supervisor's approval.
  - Effort: medium.
- **(b) Future improvement.** The paper reports the empty runs as a result of the strict rule, as
  agreed for issue 12.

**My recommendation: (b), with a trade-off to weigh.**
- (b) keeps Paper 1 on the method your supervisor is reviewing now.
- The cost is that three of the constrained runs in Table 1 stay empty, which weakens the
  demonstration of the constraints.
- If you choose (a), A7 (a) costs little extra, because it uses the same search.

**Draft paper sentence.** "A fixed variable is held strictly: the counterfactual must keep its
observed value, and when the restored point no longer crosses the decision boundary, no
counterfactual is reported. A follow-up search in the free variables, with the fixed variables held
at their observed values, may recover a counterfactual in these cases."

**Your reply**

- [ ] Fix now: option (a)
- [X ] Future improvement: option (b)
- [ ] No action
- [ ] Other (describe below)

Your solution:

>

Draft paper sentence: [ ] OK as drafted  [ ] Edit below

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: future improvement. Two points to confirm:
- The sentence box is not ticked. Is the draft sentence OK?
- Table 1 then reports three runs with no counterfactual, each with its reason, as agreed for Pima
  issue 12: TP run 3 ("no boundary point"), and TN runs 3 and 4 ("the restore fails").

- [ ] Sentence OK as drafted, and Table 1 as described
- [ ] Edit below

>

**Round 3 (2026-11-05, in the chat).** The sentence is OK as drafted, and Table 1 is as described.

---

### A6. Searching in the actionable variables only

**Status:** Resolved: future improvement (b), sentence as drafted

**Today.** The rotation and the pair search use every variable, fixed ones included. The fixed
variables are held only by the tolerance during the search and by the restore afterwards. With half
the variables fixed, few points on a 2D plane hold them.

**Where.** `.bl_rotate()`, and the pair loop of `bl_find_local_cf()` (both in `R/local_cf.R`). Also
`documentation/fixed-constraint-restore.md` §8 (second bullet).

**Effect on Paper 1.** The same as A5: with two of the four variables fixed, only TP run 4 has a
counterfactual.

**Options.**
- **(a) Hold the fixed variables before the rotation.**
  - The local biplot and its boundary are then built in the subspace of the actionable variables.
  - Every boundary point holds the fixed variables exactly, so no restore is needed.
  - Moves Paper 1's numbers: **yes**, every fixed run, through a different projection for each set
    of constraints.
  - It is a change of method, and needs your supervisor's approval.
  - Effort: high.
- **(b) Future improvement.**

**My recommendation: (b).** It is the cleanest solution, but too large a change to make during
Paper 1.

**Draft paper sentence.** "The Procrustes rotation and the search over component pairs use all the
variables, including those that are held fixed. Building the local biplot in the subspace of the
actionable variables would ensure that every boundary point keeps the fixed variables at their
observed values, at the cost of a different projection for each set of constraints."

**Your reply**

- [ ] Fix now: option (a)
- [x ] Future improvement: option (b)
- [ ] No action
- [ ] Other (describe below)

Your solution:

>

Draft paper sentence: [ ] OK as drafted  [ ] Edit below

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: future improvement. The sentence box is not ticked.

- [ ] Sentence OK as drafted
- [ ] Edit below

>

**Round 3 (2026-11-05, in the chat).** The sentence is OK as drafted.

---

### A7. Letting the free variables compensate after integer rounding

**Status:** Resolved: future improvement (b), sentence as drafted

**Today.**
- `.round_integers()` tries only the whole values just below and just above each `"integer"`
  variable, with the other variables left where the search put them.
- When the nearer whole value does not cross the boundary, the further one is taken. This happens
  even when a small change in the free variables would let the nearer one cross.

**Where.** `.round_integers()` (`R/local_cf.R`). Also `documentation/fixed-constraint-restore.md`
§5.

**Effect on Paper 1.**
- The TN target's run 1 search point has *Pregnancies* 10.21 (√d²_M 0.14, prediction 0.500).
- *Pregnancies* 10 does not cross (0.487), so 11 is taken. √d²_M rises to 0.42 and the prediction
  to 0.548.
- Run 5 behaves the same way (√d²_M 0.17 to 0.44).
- The TN counterfactual's distance is therefore mostly the cost of a whole pregnancy (Pima issue
  13).

**Options.**
- **(a) Compensate.** When the nearer whole value does not cross, hold it and run the A5 follow-up
  search along the free variables. Keep the result if it is nearer than the other whole value.
  - Moves Paper 1's numbers: **yes**. The TN run 1 and run 5 distances (0.42 and 0.44) would
    probably fall.
  - It needs A5's search.
  - Effort: medium, or low if A5 (a) is chosen.
- **(b) Future improvement.** The paper explains that a count variable moves in whole steps, as
  agreed for issue 13.

**My recommendation: (b),** unless you choose A5 (a). The two then belong together.

**Draft paper sentence.** "A count variable can only change in whole steps, so the nearest valid
counterfactual may lie well past the decision boundary. For the true negative, rounding *Pregnancies*
from 10.21 to 11 triples the Mahalanobis distance, from 0.14 to 0.42; allowing the continuous
variables to compensate after the rounding may reduce this cost."

**Your reply**

- [ ] Fix now: option (a)
- [x ] Future improvement: option (b)
- [ ] No action
- [ ] Other (describe below)

Your solution:

>

Draft paper sentence: [ ] OK as drafted  [ ] Edit below

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: future improvement. The sentence box is not ticked. I checked
its numbers against `table1.csv`: the distance rises from 0.1405 to 0.4191, a factor of 2.98, so
"triples" holds.

- [ ] Sentence OK as drafted
- [ ] Edit below

>

**Round 3 (2026-11-05, in the chat).** The sentence is OK as drafted.

---

### A8. A README constraint example that finds a counterfactual

**Status:** Resolved: fix now (a), the run 4 constraints with their output (plan 2)

**Today.**
- The README's "Actionability constraints" example takes the first test patient the Pima GAM
  predicts as diabetic.
- It sets *Glucose* to `"decrease"`, and *Age* and *Pregnancies* to `"fixed"`.
- `bl_find_local_cf()` returns no counterfactual: the constraints eliminate every segment on all 15
  pairs (re-run on 2026-11-04).

**Where.** The example is in the README's "Actionability constraints" section, and its set-up (the
GAM and `bl_results_pima`) in its "Bringing your own model" section.

**Effect on Paper 1.** None on the numbers. But the paper will point readers and reviewers to the
package, and the README is the first thing they will run.

**Options.**
- **(a) Change the example so that it finds a counterfactual.** For example one fixed variable and
  one direction, chosen by running it, with its output shown. Moves Paper 1's numbers: no. Effort:
  low.
- **(b) Keep the example**, and add a sentence saying that it finds none and why. It then shows that
  constraints can rule out every counterfactual. Effort: low.

**My recommendation: (a),** before the release tag that the paper pins.

**Draft paper sentence.** None; this is not a paper item.

**Your reply**

- [a ] Fix now: option (a)
- [ ] Fix now: option (b)
- [ ] No action
- [ ] Other (describe below)

Your solution:

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: fix now (a). The README's Pima example has its own six-variable
GAM, not Paper 1's model, so its results will differ from Table 1.

**A8.1 Which constraints?**
- [ ] (i) Paper 1's run 4, on the same patient (my recommendation). *Pregnancies* is `"fixed"`,
  *Age* is `"increase"` (it cannot decrease), and *Glucose* is `"decrease"`.
  - It matches the paper's constrained example.
  - If it finds no counterfactual, I drop *Glucose*, then try the next patient the model predicts
    as diabetic. I will tell you which version was used.
- [ ] (ii) The current example with one constraint relaxed. *Glucose* is `"decrease"`, *Age* is
  `"fixed"`, and *Pregnancies* is `"integer"` instead of fixed. It also shows the integer
  constraint.
- [ ] Other (write below)

**A8.2 What the README shows.**
- [ ] The printed result (the Target and Counterfactual table), and one sentence saying that strict
  constraints can leave no counterfactual, with `blocking_constraint` reporting why. This is my
  default; the sentence keeps the point of option (b).
- [ ] The printed result only
- [ ] Neither: the code only, as now

>

**Round 3 (2026-11-05, in the chat).** A8.1 (i): the run 4 constraints on the same patient, with
the fallbacks as written. A8.2: the printed result, plus the sentence on `blocking_constraint`.

---

### A9. The loan harness's biplot title

**Status:** Resolved: fix now (a)

**Today.** The first-pass result in the loan harness is titled "Loan default (prior defaulters) --
XGB, CVA biplot". But Step 3 keeps only the applicants **without** a prior default on file.

**Where.** The title is in the `bl_build_result()` call of Steps 4–6 of
`development/03_loan_status_Boundary_Logic.R`, and the filter in its Step 3. The walkthrough quotes
the title in §6 of `development/reference/review_section4_to_6.md`.

**Effect on Paper 1.** None.

**Options.**
- **(a) Fix the title** to "Loan default (no prior default) -- XGB, CVA biplot", in the script and
  in the walkthrough. Effort: low.
- **(b) No action.**

**My recommendation: (a).**

**Draft paper sentence.** None; this is not a paper item.

**Your reply**

- [x ] Fix now: option (a)
- [ ] No action: option (b)
- [ ] Other (describe below)

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: fix now. The title becomes "Loan default (no prior default) --
XGB, CVA biplot", in the script and in the walkthrough. No question.

---

## Part B: the other open items

These items come from the rotation note (§4, §6), the methodology note (§3, Table 5),
`fixed-constraint-restore.md` (§8), CLAUDE.md (§9) and `progress.md`. I have pre-marked each with
my suggestion.

- **Your reply:** write `Now`, `Future` or `Drop` in the last column, or leave it blank to accept my
  suggestion.
- **Paper sentences:** I will draft them in the follow-up round, for the items that stay as future
  improvements and belong in the paper.

| # | Improvement area | Source | What it would do | Moves Paper 1's numbers? | My suggestion | Your reply |
|---|---|---|---|---|---|---|
| B1 | A minimal Givens rotation with a reproducible search over the planes that contain the target | Rotation note §5, §6.1; methodology note §3.12, Table 5 | Makes the local plane unique and reproducible, and searches beyond fixed pairs. The acceptance tests 4, 5 and 7 (§6.1) are incomplete. | Yes, every Phase 3 number | Future (after Paper 1) |Future (after Paper 1) |
| B2 | A cross-platform reproducibility check | Rotation note §4.1, §6.2 (1) | Re-run the Pima pipeline (`run_all.R`, `tests.R`) on Linux and macOS, for example through GitHub Actions, and compare `numbers.rds`. A check, not a change. This matters because of the plane's under-determination: a 5.0e-17 change in V moves the TP target from pair (2,3) to (1,3), and its distance from 1.42 to 1.70 (K3). | No; it tests whether they hold | Now, before submission |Future (after Paper 1) |
| B3 | Pair terminology | Rotation note §4.2 | Describe the search as planes through the target that are labelled by component pairs, not as eigenplanes. Only 21% of a selected (2,3) plane's second direction lies in the (2,3) eigenplane. | No; wording only | Now (paper wording; draft below) | Now (paper wording; draft below)|
| B4 | A frame convention for the local plot | Rotation note §4.3, §5.2 | Fixes the orientation of the axes within the local plane, which is now LAPACK's choice. | No numbers; the figures' orientation may change | Future (with B1) | Future (with B1) |
| B5 | A unit test that pins the CVA identity | Rotation note §6.2 (4) | A `testthat` test that d²_M = (n − g) · dist_z² for every pair, including two-group CVA. `17_3` §11 asserts it today; `tests/testthat/` does not. | No | Now |  Future (with B1)|
| B6 | Feasibility as a constraint inside the search | Rotation note §4.5, §5.4 | Seeks the nearest feasible point, instead of filtering contour vertices after contouring. On the loan data, `train_ranges` removes 18–56% of the candidates per pair. | Yes | Future | Future|
| B7 | Off-plane counterfactuals | Rotation note §4.6 | Marks the restored, rounded and sparse counterfactuals as off-plane on the plot, or recomputes `dist_*` for the reported point. | Only the printed distances | Future, plus a paper sentence on which distance is quoted (Table 1 reports the search's and the counterfactual's √d²_M) | Future, as recommended|
| B8 | Pair control | Rotation note §4.7 | Adds a `pairs =` argument to force any pair, and ties the `max_pairs` default to choose(p, 2). | No (p ≤ 7 already searches every pair) | Future | Future|
| B9 | A guard for a target at the centroid | Rotation note §4.8 | Checks what happens when the target is at the centre (the Procrustes matrix is zero), then warns. | No | Now | |
| B10 | A hull rebuilt on each rotated plane | Methodology note §3.6, Table 5 | Restores the thesis's plausibility rule (p.117). | Yes: some counterfactuals may move or disappear | Depends on your supervisor's decision on §3.6 | |
| B11 | A per-row feasibility flag in the global search, and the training-range filter as an option | Methodology note §3.3, Table 5 | Flags the global counterfactuals outside the training ranges (10 of 313 training and 5 of 79 test rows on the four-feature model). The option offers the thesis's rule (p.93). | No, at the default | Now |Future (with B1) |
| B12 | Retrying the fixed restore on the next-best pair | `fixed-constraint-restore.md` §8 | Rescues targets whose winning point fails the restore. | Yes, fixed runs | Future (needs A2) | Future|
| B13 | Integer rounding at other step sizes | `fixed-constraint-restore.md` §8 | A step size per variable in `set_filters()`, for example 0.5, or 100 for a loan amount. | No | Future | Future|
| B14 | A full Mahalanobis decomposition in Phase 2 | CLAUDE.md §9; methodology note §3.1 | A Shapley decomposition of the full distance. It removes the gap between the per-variable measure and the full distance: 27.4% on seven variables, 35.2% on four. | Yes: the importance values | Future |Future |
| B15 | Aggregated Shapley values | Thesis §6.2; Table 5 | Global importance built from local attributions. | No (a new output) | Future |Future |
| B16 | A CVA biplot for a singular W (p > n) | Table 5 | The generalised SVD of Ganey and Gardner-Lubbe (2026; package `wideRhino`). | No | Future | Future|
| B17 | LDA for more than two classes | CLAUDE.md §9; Table 5 | Multi-class models, with one boundary per pair of classes. | No | Future | Future, mention extension to Multi-class models being a key area of future research. Currently resolved by only targeting one of the multi classes as 1, the rest as 0, so counterfactual is out of the target class, not towards a specfic class|
| B18 | A different transformation per variable | CLAUDE.md §9; Table 5 | Display only. A non-affine transformation needs spline-calibrated axes. | No | Future | Future areas to incorporate spline axes - see Groenen, P., Le Roux, N.& Gardner-Lubbe, S. (2014). Spline-based nonlinear biplots. Advances in Data Analysis and Classification.https://doi.org/10.1007/s11634-014-0179-1|
| B19 | A non-interactive boundary-arrow layer | CLAUDE.md §9 | Arrows to the global counterfactuals for chosen rows, in one plot. | No | Future | Drop|
| B20 | Label direction per variable | CLAUDE.md §9 | A `label_dir_var` argument in `plot_biplotEZ()`. | No | Future | drop |
| B21 | Tests of computed values | CLAUDE.md §9 | Tests of the values of the local search result, the Shapley values on a real model, `bl_robustness()` and `bl_surrogate()`, and `bl_build_result()`. | No | Now, before the release | Now, before the release  |
| B22 | A vignette on the mathematical method | CLAUDE.md §9 | Derives the projection and the counterfactual geometry from first principles. | No | Future | Future|
| B23 | CRAN preparation | CLAUDE.md §9 | | No | Future |Future |
| B24 | Housekeeping before the release | `progress.md` | Rebuild the pkgdown site. Move the rotation note's citations from the 30 September 2025 draft to the final thesis pages. | No | Now, before the release | Now, before the release |

**Draft paper sentence for B3.** "Each pair of components labels one plane through the target
point; since the Procrustes rotation turns the plane, the selected plane need not coincide with the
plane spanned by the two components that name it."

**Your comments on Part B** (write the item number first, for example "B2: ..."):

>

#### Follow-up (Part B)

**Round 2 (2026-11-05).**

**Recorded without a question:**
- Future improvement: B1, B4, B6, B8, B12, B13, B14, B15, B16, B22 and B23.
- Future improvement, with your comment: B17 and B18 (sentences below).
- Future improvement, with a paper sentence on which distance is quoted: B7 (sentence below).
- No action: B19 and B20. They were removed from CLAUDE.md §9 on 2026-11-05.
- Now, before the release: B21 and B24 (one question each below).
- Now, as paper wording only: B3 (one question below).

**Questions:**

**B2. The cross-platform check (you: Future).** In K3, a change of 5.0e-17 in V moves the TP target
from pair (2,3) to pair (1,3), and its distance from 1.42 to 1.70. So a reader who re-runs the code
on macOS or Linux, or with another linear algebra library, may not reproduce Table 1. JCGS expects
reproducible code. Two practical facts:
- WSL has no Linux distribution installed on this machine, so a local Linux run is not available.
- `scripts/` is in `.gitignore`. GitHub Actions cannot run the Pima scripts until they are in a
  repository, for example the paper repository, as the supplement (task T4).

- [x ] (a) Future, with a mitigation (my recommendation if you keep Future).
  - The paper and the supplement state the platform (R 4.6.0 on Windows, with its linear algebra
    library), and the supplement includes `sessionInfo()`.
  - The paper says that the local plane is under-determined, and that the results reproduce on the
    same platform.
  - Draft sentence: "The Procrustes rotation determines the local plane only up to its second
    direction, which the linear algebra routine settles. The results are therefore reproducible on
    a given platform, but another platform may select a different plane, and hence a different
    counterfactual; the minimal rotation proposed as future work would make the plane unique."
- [ ] (b) Run the check before submission. Once the scripts are in the paper repository, one GitHub
  Actions run on Linux and macOS compares `numbers.rds`. About half a session.
- [ ] (c) Future, with no mitigation.

>

**B3. Pair terminology (Now).**
- [ x] The draft sentence above is OK as drafted
- [ ] Edit below

`?bl_find_local_cf` and the print notes speak of eigenvector pairs. Should the help pages use the
paper's wording? That is a sentence or two, with no code change.
- [ x] Yes, reword the help pages to match (my recommendation)
- [ ] No, the paper only

>

**B5. The CVA identity test (you: Future, with B1).** B5 does not depend on B1.
- The identity d²_M = (n − g) · dist_z² follows from biplotEZ's CVA scaling (V′WV = I/(n − g)).
- The per-pair rotation cancels, so the identity holds for the current rotation and for the
  minimal rotation alike.
- The test is about ten lines, and fits with the value tests of B21, which you marked Now.

- [x ] Fold into B21, now (my recommendation)
- [ ] Keep as Future, with B1

**B9. The centroid guard (blank, so Now by the form's rule).**
- It adds a warning to `.bl_rotate()` for a target at the centroid. There the Procrustes matrix is
  zero, so the plane is arbitrary.
- No result changes.
- You moved B5 and B11 to go with B1, so please confirm where B9 goes:

- [ x] Now (a warning only)
- [ ] Future, with B1

**B10. A hull on each rotated plane (blank: depends on your supervisor's decision on §3.6).**
- If the supervisor approves the departure, B10 becomes future work, with the sentence below.
- If the supervisor rejects it, B10 becomes Now. That moves Phase 3's results, so Table 1,
  figures 5 to 9 and the methodology note would be re-run.
- Until the supervisor decides, which way do you lean? A1.6 would give the supervisor a measured
  count of the package's points that fall outside the hull.

Draft sentence, if B10 stays future work: "The local search does not restrict the counterfactual
to the convex hull of the training data on the rotated plane, so a counterfactual may lie in a
region with little training data. Rebuilding the hull on each rotated plane would restore this
plausibility rule, at the cost of removing or moving some counterfactuals."

- [ x] I lean towards keeping the departure (B10 as future work)
- [ ] I lean towards restoring the hull (B10 now)
- [ ] No view until the supervisor decides
- [ ] Sentence OK as drafted, if B10 stays future work

>

**B11. A feasibility flag in the global search (you: Future, with B1).** B11 concerns Phase 2's
global search, not the rotation, so it does not depend on B1.
- On the seven-variable model that the paper uses for Phase 2, no global counterfactual falls
  outside the training ranges (0 of 313 training and 0 of 79 test rows).
- The 10 and the 5 are on the four-variable model, which the paper uses for Phase 3 only.
- So the paper needs no sentence on B11.

- [x ] Future on its own, not in the paper (my recommendation)
- [ ] Future, with B1
- [ ] Other (write below)

**B21. Tests of computed values (Now, before the release).** The plane is under-determined (B2), so
exact values of the local search may differ on another platform. Tests that pin them may then fail
on CRAN or in continuous integration.
- [x ] (a) My recommendation:
  - property tests everywhere: the counterfactual crosses the cutoff, meets every constraint and
    the training ranges, its distances recompute from its values, and the chosen pair minimises
    the distance;
  - exact values for everything that does not depend on the plane: the projection, the global
    search, the Shapley values for a given counterfactual, `bl_robustness()`, `bl_surrogate()` and
    `bl_build_result()`;
  - exact values of the local search on Windows only (`skip_on_os()`), with a comment pointing to
    B1.
- [ ] (b) Exact values everywhere.

The data: iris as two classes, as the existing tests use, plus one test on the bundled Pima file.

>

**B24. Housekeeping (Now, before the release).** It runs last, after the fixes and the supervisor's
decisions:
1. rebuild the pkgdown site;
2. move the rotation note's citations to the final thesis pages;
3. merge and tag, on your instruction.

Which version should the release that the paper pins carry? `DESCRIPTION` says 0.1.0.
- [x ] 0.2.0 (my recommendation: there are new arguments and behaviour since 0.1.0)
- [ ] 1.0.0
- [ ] Keep 0.1.0

**Draft paper sentences for the items that stay as future work.** They follow the JCGS rules. Tick
or edit each.

- **B1** (the minimal rotation): "The Procrustes rotation fixes only the direction of the target in
  the local plane, and the second direction is arbitrary. A minimal rotation through the target,
  combined with a reproducible search over the planes that contain the target, would make the plane
  unique and would search beyond fixed pairs of components."
  [x ] OK  [ ] Edit below
- **B6** (feasibility inside the search): "The constraints and the training ranges are applied as
  filters on the vertices of the contour, so the nearest feasible point may lie between two
  vertices. Treating feasibility as a constraint inside the search would find the nearest feasible
  counterfactual directly."
  [ x] OK  [ ] Edit below
- **B7** (the distances Table 1 quotes): "Table 1 reports two Mahalanobis distances: that of the
  boundary point found on the plane, and that of the reported counterfactual after the fixed
  variables are restored and the count variable is rounded. The reported counterfactual lies off
  the plane, so the local biplot shows only its projection onto the plane."
  [ ] OK  [ ] Edit below
- **B14** (the full Mahalanobis decomposition): "The per-variable importance standardizes each
  variable by its within-group standard deviation and ignores the correlations between the
  variables, so the squared per-variable measures do not sum to the squared Mahalanobis distance;
  on the seven-variable Pima model they exceed it by 27% in total. A Shapley decomposition of the
  full Mahalanobis distance would remove this gap, at a cost that grows exponentially with the
  number of variables."
  [x ] OK  [ ] Edit below
- **B15** (aggregated Shapley values): "Aggregating the local Shapley values over many observations
  would give a global measure of importance built from the local attributions."
  [x ] OK  [ ] Edit below
- **B16** (a singular W): "The canonical variate analysis biplot requires a nonsingular
  within-group covariance matrix, which fails when there are more variables than observations; the
  generalized singular value decomposition of Ganey and Gardner-Lubbe (2026) would extend the
  method to such data."
  [x ] OK  [ ] Edit below
- **B17** (multi-class models, from your comment): "The method is restricted to binary classifiers.
  A multi-class model can be explained one class at a time, by coding the class of interest as 1
  and the other classes as 0, so that the counterfactual moves the observation out of the class of
  interest rather than towards a specific other class. The canonical variate analysis biplot
  already separates more than two groups (here the four outcome groups), so extending the decision
  boundaries and the counterfactual search to multi-class models, with one boundary for each pair
  of classes, is a key area of future research."
  [ x] OK  [ ] Edit below
- **B18** (nonlinear transformations, from your comment): "The biplot axes are linear, so a variable
  can be displayed in its original units only when its transformation is affine. Spline-based
  nonlinear biplots (Groenen, le Roux, and Gardner-Lubbe 2015) would allow nonlinear
  transformations, such as a logarithm for a skewed variable, to be shown on calibrated nonlinear
  axes."
  [ x] OK  [ ] Edit below

  On the citation: the DOI you gave (10.1007/s11634-014-0179-1) was published online in 2014, but
  the journal issue is 2015. JCGS cites the issue year, so the reference reads: Groenen, P. J. F.,
  le Roux, N. J., and Gardner-Lubbe, S. (2015), "Spline-Based Nonlinear Biplots," *Advances in Data
  Analysis and Classification*, 9, 219–238.

Your edits:

>

**Round 3 (2026-11-05).** Recorded from your ticks, and from the chat for B7 and B10:
- B2: future, with the mitigation in (a); the draft sentence goes in the paper's first future-work
  paragraph (G1).
- B3: the sentence as drafted. The help pages are reworded to match (plan 2).
- B5: folded into B21, now (plan 3).
- B9: now, a warning only (plan 2).
- B10: open until your supervisor decides §3.6. You lean towards keeping the departure (future
  work). The sentence is OK as drafted, if it stays future work.
- B11: future on its own, not in the paper.
- B21: option (a) (plan 3).
- B24: version 0.2.0 (plan 4).
- The sentences for B1, B6, B7, B14, B15, B16, B17 and B18 are OK as drafted.

**Round 4 (2026-11-05).** Your instruction: "I don't want the tests. It must only work on my
computer. There can be a note in the paper that this may be a problem." Then: "add the convex hull
on the rotated plane as something that must be done now".
- **B21:** future, in the design document only. The tests are needed only if the package is used
  beyond your computer, so B21 goes with B23 (CRAN preparation).
- **B5:** future, with B21.
- **B2:** future, with this revised sentence for the paper: "The package was developed and run on
  one platform (R 4.6.0 on Windows). The Procrustes rotation determines the local plane only up to
  its second direction, which the linear algebra routine settles, so the results are reproducible
  on that platform, but another platform may select a different plane, and hence a different
  counterfactual; the minimal rotation proposed as future work would make the plane unique." The
  supplement includes `sessionInfo()`.
- **B10: now.** Your design choices (in the chat):
  - the hull on each rotated plane is built from the training data plus the target, as Appendix A
    of the thesis;
  - it is switchable: `hull = TRUE` is the default, and `hull = FALSE` gives the earlier results;
  - the rotation scripts (`17_1` to `17_3`) pin `hull = FALSE`, so the rotation note's
    measurements stay valid;
  - there is no enclosure rule for closed contours.

  The paper describes the hull in §5 as part of the method, not as future work; the draft future
  sentence above is no longer used. The methodology note's §3.6 records it as restored, decided by
  you on 2026-11-05.

  **The hull's effect, measured before the re-run (2026-11-05):**
  - over the 79 test rows (run 1), every row keeps a counterfactual: 68 are unchanged and 11 move,
    with a median distance 4.5% longer;
  - the TN target is unchanged in every run;
  - the TP target moves from pair (2,3) to pair (3,4), and its search distance from 1.50 to 1.84.
    Inside the hull, every TP counterfactual lowers *Age* from 51 to about 41, so any run that
    holds *Age* (runs 3 and 4) has no counterfactual. Run 5 (pair (1,2) only) has none either;
  - `hull = FALSE` reproduces the earlier results exactly.

  **Your decisions on this (in the chat):**
  - keep the hull on, and keep figure 7(a) on the TP target: the panel shows that no
    counterfactual exists inside the hull;
  - add a **run 7** to Table 1: *Pregnancies* fixed and *Glucose* unable to increase. For the TN
    target, a fixed variable then holds and a counterfactual exists. In the full run (2026-11-05)
    the search distance is 0.30 and the counterfactual's 0.29: *Pregnancies* is restored to 10,
    *Insulin* rises from 122 to 149.25, and *Glucose* (128.97) and *Age* (39.04) barely move. The
    TP target's run 7 has none (the restore fails).

---

## General follow-up

**G1. How the future work is grouped in the paper.** About 15 items stay in the paper as future
work. One sentence each would make a long list, so I suggest four short paragraphs in §9 (the
conclusion and further work), each built from the sentences agreed in this form:

1. **The local plane:** B1 (with B4) and B2.
2. **Constraints and the search:** A4, A5, A6, A7, B6, B12 and B13.
3. **Attribution:** B14 and B15.
4. **Scope:** B16, B17 and B18.

Some items describe the method as it is, so they sit elsewhere:
- A2, B3 and B7 go in §5 (local interpretation);
- B10 goes in §8 (properties and limitations);
- A1 is a result in §8.

- [ x] Four themed paragraphs, as above (my recommendation)
- [ ] One sentence per item
- [ ] Other (write below)

**G2. Items in the design document only, not in the paper.** These are software or documentation
items:
- A3, A8 and A9;
- B4 (folded into B1's sentence), B5, B8, B9 and B11;
- B19 to B24.

- [ x] Agree
- [ ] Move an item into the paper (write below)

>

**Round 4 (2026-11-05).** G1: B10 moves out of §8 (limitations) into §5, as part of the method.
G2: B21 joins the design-document-only items.

---

## Part C: for information, no reply needed

**Already resolved.** These items were fixed before or during the Task D list, and become design
decisions in the final document (the former `documentation/known-issues.md` §2):
- the deterministic nearest boundary vertex, which also removed the dependence on earlier calls;
- the cutoff passed to the grid;
- `rotate_deg` removed from the local plots;
- `xgboost::xgboost()` models accepted;
- empty CVA groups dropped;
- the thesis distance rule by projection (`distance = "auto"`, 2026-11-05);
- `used_in_sparse` for fixed variables;
- the stale-documentation sweep.

**Properties the paper states as limitations of the method, not as future work.** These come from
the rotation note §4.9:
- a 2D biplot is a slice of the decision surface, not a view of it;
- rotating to fit one target worsens the representation of the other points;
- the search finds the nearest boundary within the planes it draws, not in the full space (A1
  measures this);
- the CVA class means leave the plane;
- "optimal" refers only to the smallest distance.

If you think any of these should be future work instead, say so here:

>

**Round 2 (2026-11-05).** No comment, so Part C stands as written.
