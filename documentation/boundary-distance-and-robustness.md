# Distance to the Boundary and Robustness

**Package:** `boundarylogic`
**Location:** `documentation/boundary-distance-and-robustness.md`

This note covers how the package measures the distance from an observation to the decision boundary in Phase 2,
and what it builds on that distance:

- **the metric** $\mathbf{W}$ (the pooled within-class covariance) and its inverse;
- **the distance** from each observation to its nearest boundary point (`bl_find_boundary()`);
- **the per-variable importance** and the scalar robustness score (`plot.bl_boundary()`, `bl_robustness()`).

**Phase 3 is in a companion note.** The rotation of the local biplot and the choice of which plane to rotate,
including how the local search uses $\mathbf{W}$, are in
[`biplot-rotation/biplot_rotation_and_pair_selection.md`](biplot-rotation/biplot_rotation_and_pair_selection.md).

**Companion documents:**

- `../1 Foundation intro documents.txt`: methodological principles;
- `../2 implementation_summary.txt` §4.2.1, §4.5, §4.6 and §4.6.1: code mechanics;
- `../.claude/reference/mahalanobis_implementation_plan.md`: the approved plan (2026-05-20) that introduced
  $\mathbf{W}$ into the package.

---

## 1. The metric $\mathbf{W}$

### 1.1 What it is

`bl_build_projection()` computes a metric matrix once and stores three fields (`R/projection.R`,
`.compute_metric_inverse()`):

- `metric`;
- its inverse, `metric_inv`;
- a tag, `metric_type`.

All three propagate through `bl_assemble()` into `bl_result`. For training data $\mathbf{X}$ ($n\times p$) and a
class factor with $g$ levels:

$$
\mathbf{W} = \frac{1}{n-g}\sum_{k=1}^{g} (\mathbf{X}_k - \bar{\mathbf{x}}_k)'(\mathbf{X}_k - \bar{\mathbf{x}}_k)
$$

That is, centre each class on its own mean, sum the cross-products, and divide by the degrees of freedom. The
result is the unbiased within-class covariance under a common covariance across classes, the assumption behind LDA
and CVA. The code divides by `max(n - g, 1)`; the 1 only guards against a zero or negative divisor when there are
no more rows than classes.

**Which class factor**, in order of priority:

| `metric_type` | class factor | when |
|---|---|---|
| `"W_cva"` | the factor passed to `biplotEZ::CVA()`: by default the model's four outcomes (FN, FP, TN, TP; empty groups dropped), or a user-supplied `cva_classes` | `method = "CVA"` |
| `"W_binary"` | the binary 0/1 `class` column | `method = "PCA"` with a class column |
| `"Sigma"` | none — the total covariance `cov(X)` | no class information (rare) |

**$\mathbf{W}$ is computed on the unstandardised features**, i.e. the model units, whatever `standardise` is.

### 1.2 Its relation to the CVA biplot

biplotEZ (2.2) computes the same within-class matrix without the divisor, the sums of squares and products
$\mathbf{W}_{SSP}$, and returns it as `bp$Wmat`. It builds the CVA loadings as
$\mathbf{M} = \mathbf{W}_{SSP}^{-1/2}\mathbf{V}_{svd}$, so $\mathbf{M}'\mathbf{W}_{SSP}\mathbf{M} = \mathbf{I}$,
and in the package's scaling

$$
\mathbf{V}'\mathbf{W}\mathbf{V} = \mathbf{I}/(n-g).
$$

- **This is general.** It holds for any CVA biplot the package builds, including 2-group CVA, by construction.
- **The consequence:** distances measured with $\mathbf{W}$ are consistent with the CVA biplot's axes.
- **The PhD research code** used biplotEZ's `Wmat` directly. The package recomputes it with the divisor $n-g$.
  This scales every distance by the same constant and changes no comparison.
- **Why the package computes $\mathbf{W}$ rather than reading `bp$Wmat`:**
  - **PCA has no `Wmat`.** A biplotEZ PCA object has none, so `.compute_metric_inverse()` must compute the
    `"W_binary"` and `"Sigma"` metrics itself. Computing `"W_cva"` there too gives all three on one scale, in one
    code path.
  - **The divisor makes $\mathbf{W}$ a covariance.** $\sqrt{d_M^2}$ then counts within-class standard deviations,
    and $\sqrt{W_{jj}}$ is the within-class SD that §3 divides by. The methodology-changes note
    ([`methodology-changes/methodology-changes.qmd`](methodology-changes/methodology-changes.qmd), "The metric W as
    a covariance matrix") records this reason and its effect.
  - **Under CVA the two agree:** `bp$Wmat / (n - g)` equals `metric` to rounding (checked on iris with four
    groups, biplotEZ 2.2, 2026-11-06).

The rotation note (§3.4) proves the consequence that matters for the local search. Under CVA, for any displacement
lying in a biplot plane, rotated or not, the squared Mahalanobis distance is $(n-g)$ times the squared in-plane
Euclidean distance.

### 1.3 The inverse: Cholesky with a ridge fallback

For a positive-definite $\mathbf{M}$, the squared Mahalanobis distance of a displacement $v$ is
$d_M^2(v) = v'\mathbf{M}^{-1}v$.

- **When $\mathbf{M}$ is the identity**, it reduces to Euclidean distance.
- **When $\mathbf{M}$ is diagonal**, it reduces to per-feature standardised Euclidean distance.

**The Cholesky factor.** $\mathbf{W}$ is symmetric and positive-definite when it has full rank, which holds when
$n-g \ge p$. It then factorises uniquely as $\mathbf{W} = \mathbf{L}\mathbf{L}'$, with $\mathbf{L}$ lower
triangular with a positive diagonal. $\mathbf{L}^{-1}$ is the whitening transform:

$$
u = \mathbf{L}^{-1}v, \qquad \lVert u\rVert^2 = v'\mathbf{L}'^{-1}\mathbf{L}^{-1}v = v'\mathbf{W}^{-1}v.
$$

**The inversion.** The package inverts with `chol2inv(chol(W))` rather than `solve(W)`:

1. **Numerical stability.**
   - **Cholesky exploits the SPD structure** that LU decomposition does not, and it propagates less round-off.
   - **It avoids a non-positive $v'\mathbf{W}^{-1}v$** when $\mathbf{W}$ is near-singular.
2. **Singularity detection.**
   - **`chol()` fails loudly:** it errors when $\mathbf{W}$ is not positive-definite.
   - **The fallback it triggers** is a ridge, `chol2inv(chol(W + lambda * I))` with
     $\lambda = 10^{-6}\cdot\text{mean}(\text{diag}(\mathbf{W}))$.
   - **`solve()`, by contrast,** may return a meaningless inverse without warning.
3. **Speed.** It is roughly twice as fast for SPD matrices.

The factor $\mathbf{L}$ is not kept: only `metric` ($\mathbf{W}$) and `metric_inv` ($\mathbf{W}^{-1}$) are stored.
The same rule applies to any SPD inversion in the package (CLAUDE.md §5).

---

## 2. Distance to the nearest boundary — `bl_find_boundary()`

For every observation, `bl_find_boundary()` finds the nearest point on the opposite-class boundary contour of the
**global** biplot.

- **No rotation.** Phase 2 never calls `.bl_rotate()`.
- **The nearest vertex** is chosen by 2D Euclidean distance in the biplot plane. It is the exact minimum, so the
  result is deterministic (`nearest-vertex-tie-break.md`).
- **`train_ranges`** are checked *after* selection and reported, not used to discard points. Discarding would pick
  a farther, geometrically wrong point (`2 implementation_summary.txt` §6.11).
- **Closed contours.** A contour still closed after clipping to the training hull is kept only if it encloses
  *training* rows (`.encloses_data()`; the thesis's empty-region rule, §3.10 step 4, p.93). Until 2026-11-06 the
  rows being searched decided (as in the PhD code), so a row's counterfactual could depend on which other rows
  were searched with it (`documentation/algorithms.md`, finding F2). On the loan harness's 4-feature model the
  change keeps 18 contours instead of 12 and moves 346 of the 4,429 test rows' counterfactuals; on the Pima
  example and the 9-feature loan model nothing changes.

It returns three distances (`R/boundary.R`, the "Distances" block of `bl_find_boundary()`):

| field | definition |
|---|---|
| `dist_z` | Euclidean distance in the plane from the observation's projection `Z_obs` to the boundary point `B_z` |
| `dist_x` | Euclidean distance in X-space from the **observation** `x_obs` to the back-projected boundary point `B_x`, each feature divided by its total SD `X_sd` |
| `dist_mahalanobis` | the **squared** Mahalanobis distance of the same X-space move, $v'\mathbf{W}^{-1}v$ with $v = B_x - x_{obs}$ in model units and $\mathbf{W}^{-1}$ = `metric_inv` |

`dist_mahalanobis` is squared, as the field of the same name in `bl_find_local_cf()` is, so the two can be compared
directly; `sqrt()` gives the distance, and `print()` shows the mean of that. It is `NA` where no boundary point was
found and, with a warning, when `metric_inv` is missing (an older `bl_result`, §5).

**How they differ.**

- **`dist_x` and `dist_mahalanobis` include the off-plane residual.** In the global plane an observation is not on
  the plane, so both measure from the observation itself. `dist_x` scales each feature by its total SD and ignores
  correlations; `dist_mahalanobis` uses the within-class metric, correlations included.
- **`dist_z` measures only the in-plane displacement.**
- **Under CVA, `dist_z` is a Mahalanobis distance:** $(n-g)\cdot\texttt{dist\_z}^2$ is the squared Mahalanobis
  distance of that in-plane displacement (§1.2). It equals `dist_mahalanobis` only when the observation lies in the
  plane, which in the global biplot it rarely does.

**Why they differ, and how they are used** (`documentation/algorithms.md`, finding F3, decided 2026-11-06). Write
the observation's full coordinates as $(\mathbf{x}_i-\bar{\mathbf{x}})\mathbf{S}^{-1}\mathbf{V} = [\mathbf{z}_i,
\mathbf{w}_i]$: $\mathbf{z}_i$ holds the two plotted ones and $\mathbf{w}_i$ the other $p-2$. Then
$\mathbf{x}_i = \bar{\mathbf{x}} + \mathbf{z}_i\mathbf{V}_r^{-}\mathbf{S} + \mathbf{w}_i\mathbf{V}_{-r}^{-}\mathbf{S}$,
while the counterfactual is $\mathbf{x}^{cf}_i = \bar{\mathbf{x}} + \mathbf{b}_i\mathbf{V}_r^{-}\mathbf{S}$, a point of
the plane with $\mathbf{w} = \mathbf{0}$. So

$$
\mathbf{x}^{cf}_i - \mathbf{x}_i = (\mathbf{b}_i - \mathbf{z}_i)\mathbf{V}_r^{-}\mathbf{S} \;-\; \mathbf{w}_i\mathbf{V}_{-r}^{-}\mathbf{S}:
$$

the move along the plane to the boundary, plus the reset of the observation's off-plane part to the mean, which has
nothing to do with the boundary. The off-plane part is most of the global gap: 80% of its summed squared length on the
loan harness's 9-feature model, 59% on its 4-feature model and 78% on Pima's 7-feature model (within-class SD units;
the methodology note's F3 rows). Hence:

- **The variable importance and the robustness (§3) use only the in-plane move**, as the thesis does. They are not
  computed from `dist_x` or `dist_mahalanobis`, and their totals do not decompose them.
- **`dist_x` and `dist_mahalanobis` are reported only.** No global or local interpretation in the package uses them
  (the local search measures its own distances from the target, which lies on its rotated plane, §2.1). The Pima
  example of Paper 1 uses the global $\sqrt{d^2_M}$ to choose its two Phase 3 targets (the TP and TN test rows
  nearest their groups' median) and quotes its median.

### 2.1 Why the distance uses $\mathbf{W}$, not a rotated $\mathbf{W}$

The local search rotates the biplot for each eigenvector pair: pair $i$ has its own orthogonal rotation
$\mathbf{A}_i$ ($\mathbf{A}_i\mathbf{A}_i' = \mathbf{I}$). It is natural to ask whether a distance should allow for
that rotation, using a "rotated" metric. It need not, for three reasons.

1. **The distance is measured in X-space.** Both `bl_find_boundary()` and `bl_find_local_cf()` measure the move
   $v = B_x - x_{obs}$ between two points in the feature space. $\mathbf{W}$ is the within-class covariance of the
   features. It is a property of the data, not of any biplot plane, so the same $\mathbf{W}^{-1}$ measures every
   X-space move, whichever plane produced the boundary point.
2. **The rotation acts only on the Z coordinates.** $\mathbf{A}_i$ decides *which* plane is searched, and so which
   boundary point is found. It does not change how a move between two given X-space points is measured.
   `bl_find_boundary()` works in the global plane and has no rotation at all.
3. **A rotation would enter only if the distance were measured in Z.** For a move inside pair $i$'s plane, the
   metric in Z coordinates is the one $\mathbf{W}$ induces through that pair's back-projection.
   - **Under CVA it is $(n-g)\,\mathbf{I}$ for every pair:** $\mathbf{A}_i$ enters only as
     $\mathbf{A}_i\mathbf{A}_i' = \mathbf{I}$ and cancels. The rotation note §3.4 proves this, so a Z-space and an
     X-space measurement agree.
   - **Under PCA it is not a multiple of the identity,** so a Z-space distance would depend on the rotation.
     Measuring in X-space avoids that. (Under PCA `bl_find_local_cf()` compares the pairs by the Euclidean
     distance in X-space by default, `dist_x`, since 2026-11-05; this point applies when `"mahalanobis"` is
     passed.)

So "allowing for $\mathbf{A}$" is either automatic (CVA: it cancels) or the wrong thing to do (PCA: it would make
the distance depend on the plane chosen to display the move).

---

## 3. Per-variable distance and robustness — `plot.bl_boundary()`, `bl_robustness()`

### 3.1 The original method (thesis §4.2.2–4.2.3, pp.101–108)

**Each observation's change vector.**

- **On the plane:** the vector from observation $i$ to its nearest counterfactual is
  $\boldsymbol{\delta}_i = \mathbf{z}_{cf} - \mathbf{z}_i$. It is orthogonal to the local decision boundary.
- **In the original space:** it decomposes into per-variable components as
  $\mathbf{d}_{ij} = (\mathbf{z}_{cf} - \mathbf{z}_i)\mathbf{M}_2^{-}$.
- **Reading it:** the longer a variable's component, the more that variable contributes to the change needed to
  cross the boundary.

**Standardisation.** To compare variables, each component is divided by "the overall standard deviation of
that variable" (p.103). The thesis notes that outliers can distort that standard deviation.

**The per-variable total:**

$$
T_j = \sum_i \left|\frac{d_{ij}}{std_j}\right|, \qquad j = 1, \ldots, p.
$$

- **Interpreting it:** the variable with the smallest $T_j$ is the least important (compare Laugel et al., 2018).
- **The plot:** a jitter plot of the standardised components across all training points, coloured by outcome
  (Figure 4-5).

**Robustness** (§4.2.3) is the total distance of the training points to their counterfactuals.

- **The rule:** a larger total means a more robust model, one that needs larger input changes to alter its
  predictions.
- **The thesis's example:** a GBM with 1,000 trees totals 720, and one with 5,000 trees totals 450. The
  less-overfitted model is the more robust.

### 3.2 The current implementation

`plot.bl_boundary()` (`R/boundary_plot.R`, its "Compute signed standardised distance per variable" block; the same
code in `bl_robustness()`) builds the signed in-plane displacement per observation and feature:

```r
vec_to_boundary_Z  <- x$Z_obs - x$B_z
vec_to_boundary    <- vec_to_boundary_Z %*% tVr          # back to X-space
if (isTRUE(standardise)) vec_to_boundary <- sweep(vec_to_boundary, 2L, X_sd, "*")
denom <- if (distance == "mahalanobis") sqrt(diag(metric)) else X_sd
vec_to_boundary_sd <- sweep(vec_to_boundary, 2L, denom, "/")
sum_of_distance    <- colSums(abs(vec_to_boundary_sd), na.rm = TRUE)
```

- **What it plots:** the jitter or box plot shows `vec_to_boundary_sd`. The axis labels give each variable's total
  `sum_of_distance` ($T_j$).
- **The robustness score:** `bl_robustness()` returns the same totals and their sum, `robustness`. The score is
  computed on the rows the boundary was found on: by default `test_data`, where the thesis used the training data.
- **The denominator** is set by `distance`. The default, `"auto"` (since 2026-11-05), follows the thesis's rule by
  projection: the total SD for a PCA biplot, the within-class SD for a CVA biplot. Before, `"mahalanobis"` was the
  default for both.

| `distance` | denominator | |
|---|---|---|
| `"mahalanobis"` (CVA default) | $\sqrt{\mathbf{W}_{jj}}$, the within-class SD | the package's choice (2026-05) |
| `"euclidean"` (PCA default) | `X_sd`, the total SD | the thesis's rule |

### 3.3 Why the CVA default is the within-class SD

**The decomposition.** For each feature, the total variance splits into within-class and between-class parts:

$$
\sigma_j^2 = \mathbf{W}_{jj} + \mathbf{B}_{jj}.
$$

**Why the total SD under-weights the best features.**

- **What makes a good class separator:** a small $\mathbf{W}_{jj}$ (tight classes) and a large $\mathbf{B}_{jj}$
  (classes far apart).
- **What the total SD does with it:** the separator's total $\sigma_j$ is inflated by exactly the between-class
  signal the measure is meant to detect.
- **So:** dividing by the total SD down-weights the most informative features.

**The within-class SD measures what matters.** The variability that matters for classification is the variability
*within* a class, so the default divides by $\sqrt{\mathbf{W}_{jj}}$. The structure of the measure is unchanged:

- a signed $n\times p$ matrix;
- per-feature totals on the plot's axis;
- a scalar robustness.

Only the denominator differs from the thesis. Under PCA the default keeps the thesis's total SD (since
2026-11-05), so the departure applies to CVA biplots, where $\mathbf{W}$ is the metric the projection is built on.
Under PCA $\mathbf{W}$ comes from the 0/1 class (`W_binary`), which the PCA axes do not use.

---

## 4. Why diagonal, not full, Mahalanobis

### 4.1 The decomposition problem

The full squared Mahalanobis distance is a quadratic form:

$$
v'\mathbf{W}^{-1}v = \sum_j v_j^2\,(\mathbf{W}^{-1})_{jj} \;+\; 2\sum_{j<k} v_j v_k\,(\mathbf{W}^{-1})_{jk}.
$$

The diagonal terms can be attributed to features; the cross terms belong to *pairs* of features. There are three
options:

| option | per feature? | cross-correlation? | cost |
|---|---|---|---|
| **diagonal** (what the package uses) | yes, exact | ignored | $O(n)$ |
| **Cholesky-whitened components** | yes, additive | preserved | $O(n)$, but loses original-variable identity |
| **Shapley attribution** | yes, principled | distributed fairly | $O(n\cdot 2^p)$ |

**Whitened components lose the variable names.** They give additive contributions $u_k^2$ for
$u = \mathbf{L}^{-1}v$, but in the *whitened basis*. "Whitened direction 3 contributes X" cannot be read as
"`loan_amnt` contributes X", and naming variables is the purpose of `plot.bl_boundary()`.

**Shapley attribution is the theoretically correct answer.**

- **The game:** treat the features as players and $v'\mathbf{W}^{-1}v$ as the value function.
- **Its properties:** the attributions $\phi_j$ satisfy efficiency ($\sum_j\phi_j = v'\mathbf{W}^{-1}v$),
  symmetry, linearity and the null-player property.
- **The PhD identified it as the correct method.**
- **Why it is not used in Phase 2:** the cost is $O(n\cdot 2^p)$. Every observation contributes, and each costs
  $2^p$ subset evaluations.

Phase 3's `bl_shapley()` does use exact Shapley, for one observation at a time: $O(2^p)$, which is tractable.

**The diagonal is a deliberate trade-off:**

- the cost is linear;
- the interpretation by named variable is exact;
- the main correction, within-class rather than total noise, is kept;
- the cross-correlation effect is silently absorbed.

### 4.2 The Mahalanobis check

To make that loss visible, `plot.bl_boundary()` prints a check when `distance = "mahalanobis"`
(`R/boundary_plot.R`, the "Mahalanobis check" block of `plot.bl_boundary()`). It compares two squared distances,
each summed across the observations, both of the in-plane displacement (so its "full" value is not the sum of
`bl_find_boundary()`'s `dist_mahalanobis`, §2):

- **plotted**: $\sum_j v_j^2/\mathbf{W}_{jj}$, the squared length of the per-feature distances the plot shows
  ($v_j/\sqrt{\mathbf{W}_{jj}}$, the `vec_to_boundary_sd` matrix);
- **full**: $v'\mathbf{W}^{-1}v$, the full squared Mahalanobis distance;
- **difference**: full − plotted, the part the per-feature plot leaves out by ignoring the correlations in
  $\mathbf{W}$, also given as a percentage of full.

The harness's Step 7 (`development/reference/review_section7_to_8.md`, nine features):

```
Mahalanobis check (squared distances, summed across obs):
  plotted (per feature) = 9675.25 | full d_M^2 = 11754.70 | difference = 2079.45 (+17.7% of full)
```

For the four-feature model of Step 7b the difference is +15.9%.

**Reading it:**

- **The sign.** A positive difference means the plot *understates* the full distance: the correlations make the
  displacement longer than its per-feature parts suggest. A negative difference means correlated features are
  double-counted, so the plot *overstates* it.
- **The units.** The check is in squared units, because $d_M^2$ is a quadratic form. The plot's per-variable totals
  $T_j$ are $L_1$ sums of the same standardised components, so the percentage is a guide to how much the
  correlations matter, not an exact share of $T_j$.
- **When it warns.** If the absolute percentage is at least 25%, it adds a note pointing to
  `2 implementation_summary.txt` §4.6.1.
- **The numbers are returned** as `mahalanobis_check` (`plotted`, `full`, `difference`, `pct`) in
  `plot.bl_boundary()`'s invisible result; `NULL` for `"euclidean"`.

**Changed 2026-10-31.** The diagnostic's "diagonal" used to be $\sum_j v_j^2(\mathbf{W}^{-1})_{jj}$, the diagonal
of the *inverse*. That is not the plotted measure: for an SPD matrix $(\mathbf{W}^{-1})_{jj} \ge 1/\mathbf{W}_{jj}$,
with equality only when feature $j$ is uncorrelated with the others. On the same Step 7 run it reported −42.6%,
which mostly measured that mismatch.

---

## 5. Edge cases and safeguards

1. **Singular $\mathbf{W}$.** When $n-g < p$, `chol()` fails and the ridge fallback of §1.3 is used. `metric_type`
   is unchanged.
2. **PCA without class information.** The metric falls back to $\boldsymbol{\Sigma} = $ `cov(X)`: standard rather
   than within-class Mahalanobis.
3. **Older `bl_result` objects.**
   - **The cause:** objects built before 2026-05-20 have no `metric` / `metric_inv`.
   - **What `plot.bl_boundary()` / `bl_robustness()` do:** if `metric` is missing, they fall back to
     `"euclidean"`, with a warning to rebuild.
   - **What `bl_find_local_cf()` does:** the same, if `metric_inv` is missing.
   - **What `bl_find_boundary()` does:** returns `dist_mahalanobis` as `NA`, with a warning to rebuild.
4. **External targets.** No class information is needed when distances are computed: the metric was fixed when the
   projection was built.
5. **Reproducing the thesis.** `distance = "euclidean"` reproduces the thesis's Phase 2 measure; it is the default
   for a PCA biplot (`"auto"`, since 2026-11-05).

---

## 6. Code references

All paths are relative to this file.

- [`../R/projection.R`](../R/projection.R): `.compute_metric_inverse()` and its call in `bl_build_projection()`.
- [`../R/result.R`](../R/result.R): propagation through `bl_assemble()`.
- [`../R/boundary.R`](../R/boundary.R): `bl_find_boundary()`, with `dist_z`, `dist_x` and `dist_mahalanobis`.
- [`../R/boundary_plot.R`](../R/boundary_plot.R): `plot.bl_boundary()`, `bl_robustness()`, the denominator and the
  Mahalanobis check.
- `../2 implementation_summary.txt`:
  - §4.2.1: $\mathbf{W}$ and its inversion;
  - §4.5: the global search;
  - §4.6: robustness;
  - §4.6.1: Shapley versus the diagonal.
- `../development/reference/review_section7_to_8.md`: the Phase 2 walkthrough on the loan harness.
- **Thesis:** the final version of 27 February 2026 (`FINAL_ADRIAAN ROWAN_27022026.pdf`, in this folder;
  printed page = PDF page − 16).
  - §3.2: Euclidean and Mahalanobis distance;
  - §4.2.2–4.2.3: variable importance and robustness;
  - §6.4: calibration from the distance.

---

## 7. Future work

- **Shapley attribution of the full Mahalanobis distance in Phase 2.**
  - **What it gives:** the theoretically correct per-feature decomposition.
  - **Why it waits:** the cost is $O(n\cdot 2^p)$ (CLAUDE.md §9).
- **A Cholesky-whitened diagnostic plot.** It would show the whitened decomposition ($u_k^2$) across the training
  set, without claiming per-variable attribution. That is useful when the check's percentage is large; it has
  not been requested.
- **Calibrated probabilities from the boundary distance.** The thesis's §6.4 (pp.145–146) proposes Platt scaling on
  the change vector's Euclidean or Mahalanobis length. That would give a model-agnostic posterior class probability
  from the distance to the boundary.
