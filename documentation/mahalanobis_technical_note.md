# Technical Note: Mahalanobis Distance in `boundarylogic`

**Author:** Adriaan Rowan
**Date:** 2026-05-20
**Package:** `boundarylogic`
**Branch:** `method_developments`
**Location:** `documentation/mahalanobis_technical_note.md`
**Companion documents (paths relative to this file):**
- `../1 Foundation intro documents.txt` -- methodological principles
- `../2 implementation_summary.txt` Sections 4.2.1, 4.6.1, 4.8 -- code mechanics
- `../.claude/reference/mahalanobis_implementation_plan.md` -- the approved plan that was executed

---

## 1. What changed

Two distance measures in the package were upgraded from **per-feature Euclidean
(divided by total SD)** to **Mahalanobis-aware (using within-class covariance
`W`)**:

| Function | Before | After |
|---|---|---|
| `bl_find_local_cf()` -- *cross-pair* selector | Euclidean distance in rotated 2D Z-space | Squared Mahalanobis distance in original X-space using `W^{-1}` |
| `bl_robustness()` / `plot.bl_boundary()` -- per-feature standardisation denominator | Total per-feature SD (`X_sd`) | Within-class SD (`sqrt(diag(W))`) |

Both changes default to the new Mahalanobis-aware behaviour. The legacy
behaviour is available via `distance = "euclidean"`. The `W` matrix and its
Cholesky inverse are computed once in `bl_build_projection()`, stored on
`bl_projection$metric` / `metric_inv`, and propagated through `bl_assemble()`
into `bl_result`.

---

## 2. Why the change matters: the methodological problem

The package's two affected measures both standardise distances so different
features can be compared on a common scale. The denominator they used --
`X_sd`, the marginal per-feature standard deviation -- is **the wrong noise
model for a classification problem**.

### 2.1 The noise we care about is within-class noise

Consider a classification problem with two classes. For each feature `j`,
the total variance decomposes as:

```
total variance = within-class variance + between-class variance
sigma_j^2      = W[j,j]                 + B[j,j]
```

where `W` is the pooled within-class covariance and `B` is the between-class
covariance. `X_sd[j]^2 = sigma_j^2` -- the total.

A feature that is a **good class separator** has small `W[j,j]` (tight
clusters) and large `B[j,j]` (clusters far apart). Its total `sigma_j^2` is
*inflated* by between-class signal -- which is exactly the signal we want
to detect. Using `X_sd[j]` as a denominator therefore **down-weights the most
informative features** for distinguishing classes.

The classification "noise" is the variability *within* each class. Two
points in the same class should be considered "similar" regardless of how
far the class centroids are from each other. That is `W[j,j]`, not
`sigma_j^2`.

### 2.2 The per-eigenvector-pair problem in Phase 3

`bl_find_local_cf()` searches across multiple eigenvector pairs `(i, j)` to
find the rotated 2D plane in which the target observation is closest to the
decision boundary. Each pair gives a different SVD rotation of the loading
matrix, producing a different 2D Z-space slice.

Comparing distances *across pairs* using 2D Euclidean distance is
mathematically suspect:

- Each pair captures a different amount of total variance in its 2D slice.
- The CVA loading matrix `V` is **`W`-orthonormal** (`V^T W V = I`), not
  Euclidean-orthonormal. A 1-unit Euclidean step in Z-space therefore
  corresponds to a *different* feature-space change depending on the pair.
- The rotation `A` returned by `.bl_rotate()` is orthogonal in p-D space but
  the projected 2D slice `Vrho[, c(1,2)]` does not preserve Euclidean
  distances uniformly across pair choices.

The right comparison is a **pair-independent** metric -- something that
lives entirely in X-space and does not depend on which 2D rotation produced
the candidate counterfactual.

### 2.3 The PhD's original solution

In the original PhD research code (`scripts/1.3 Optimal Rotation.R`,
preserved on the `original_PhD_code` branch), the per-pair comparison uses
the squared Mahalanobis distance with within-class covariance:

```r
pc.distz[i] <- (v) %*% Wmat_inv %*% t(v)        # v = boundary - observed
```

where `Wmat` was computed in `scripts/1 Functions for biplot.R` as the
pooled within-class scatter matrix used by CVA:

```r
Wmat <- sum over classes k of crossprod(X_k - colMeans(X_k))
```

The current package delegates the CVA computation to `biplotEZ::CVA()`,
which does not expose `Wmat` on its return object. So the package -- before
this change -- had **no access to the metric `V` was orthonormal under**, and
therefore could not perform the correct pair comparison. It silently
substituted 2D Euclidean, which was both inconsistent across pairs and
inconsistent with the PhD's specification.

`1 Foundation intro documents.txt` Section 10.3 names this exact issue:

> *"CVA within-class scatter must use the full maths (`Wmat_inv`), not a
> simple matrix inverse."*

The implementation re-derives `W` from the same training data and the same
class factor used by CVA, then stores it.

---

## 3. Mathematical foundations

### 3.1 The squared Mahalanobis distance

For a positive-definite metric matrix `M` (p x p) and a displacement vector
`v` in R^p:

```
d_M^2(v) = v^T M^{-1} v
```

Properties:
- **Symmetric and non-negative**: `d_M^2(v) >= 0`, with equality iff `v = 0`.
- **Invariant to linear reparametrisation in the metric**: if `X' = A X` and
  `M' = A M A^T`, then `d_{M'}^2(A v) = d_M^2(v)`. So Mahalanobis distance
  is geometrically meaningful regardless of coordinate choice within the
  whitened space.
- **Reduces to Euclidean** when `M = I` (no covariance structure).
- **Reduces to per-feature standardised Euclidean** when `M = diag(sigma_j^2)`
  (no cross-correlations).

The classical Mahalanobis distance uses `M = Sigma` (total covariance). For
classification, the principled choice is `M = W` (within-class covariance) --
the metric that CVA itself optimises against.

### 3.2 Why Cholesky decomposition

`W` is symmetric positive-definite (whenever it has full rank, which holds
when `n - g >= p` -- training size exceeds features after class adjustment).
For SPD matrices, the **Cholesky decomposition** factorises:

```
W = L L^T
```

where `L` is lower triangular with positive diagonal entries. `L` is the
unique "matrix square root" of `W`. Geometrically, `L^{-1}` is the
**whitening transform** that maps the original (correlated) feature space
into an orthonormal frame where Euclidean distance equals Mahalanobis
distance:

```
u = L^{-1} v
||u||^2 = u^T u = v^T L^{-T} L^{-1} v = v^T (L L^T)^{-1} v = v^T W^{-1} v
```

The implementation inverts `W` using `chol2inv(chol(W))` rather than
`solve(W)` for three reasons:

1. **Numerical stability.** `solve()` uses LU decomposition with partial
   pivoting, which works for any invertible matrix but does not exploit
   SPD structure. `chol2inv(chol(.))` is the textbook stable inverse for
   SPD matrices: it propagates less round-off error and avoids producing a
   non-positive `v^T W^{-1} v` (which would silently corrupt distance
   measures) when `W` is near-singular.

2. **Singularity detection.** `chol()` errors immediately if `W` is not
   positive-definite, giving a clean trigger for a ridge-regularised
   fallback `W + lambda * I` (with `lambda = 10^-6 * mean(diag(W))`).
   `solve()` may return a meaningless inverse without warning.

3. **Speed.** `chol2inv()` is roughly 2x faster than `solve()` for SPD
   matrices.

The Cholesky factor `L` itself is not retained -- only `W` (as `metric`)
and `W^{-1}` (as `metric_inv`) are stored. Downstream code multiplies
vectors by `metric_inv` directly.

### 3.3 What `W` is, concretely

For training data `X` (n x p) with class factor `g` having g levels:

```
W = (1 / (n - g)) * sum over classes k of:  (X_k - colMeans(X_k))' (X_k - colMeans(X_k))
```

For each class, centre the within-class observations on the class mean,
take the outer product of those centred residuals, sum across classes,
divide by degrees of freedom. The result is an unbiased estimate of the
within-class covariance assuming all classes share a common covariance --
the same assumption underlying LDA and CVA.

In this package, the class factor is the same one passed to
`biplotEZ::CVA()` when building `V`. By default that is the confusion
labels (TP/TN/FP/FN, 4 levels) when a model is supplied; alternatively the
user can pass the binary class factor (2 levels). For PCA-based
projections, the package falls back to the binary 0/1 class, and finally
to total covariance `Sigma = cov(X)` when no class information is
available at all (rare in this package).

---

## 4. The two distance measures in detail

### 4.1 Phase 3 cross-pair selector (`bl_find_local_cf()`)

For each candidate eigenvector pair, the within-pair search identifies the
nearest valid boundary point `B_x` (back-projected to X-space). The
cross-pair comparison then uses:

```
v        = B_x - x_obs                              # 1 x p in X-space
dist_mah = (v) %*% metric_inv %*% (v)               # scalar squared d_M^2
```

The pair with the smallest `dist_mah` wins. Both `dist_z` (legacy 2D
Euclidean) and `dist_mah` are reported per pair so users can inspect both.

**What changes in practice.** For most observations the two selectors agree
because the local 2D geometry happens to be a reasonable proxy for the
global Mahalanobis geometry. They disagree when:

- One pair has a *short* Z-distance but a *long* X-space displacement
  (the 2D slice happens to be aligned with a low-variance direction in V's
  W-orthonormal frame). Euclidean would choose this pair; Mahalanobis
  correctly down-weights it.

- One pair has a *moderate* Z-distance but the displacement vector lies in
  a low-within-class-noise direction (a "tight" direction in feature
  space). Mahalanobis correctly elevates this pair as the more
  methodologically defensible counterfactual.

### 4.2 Phase 2 per-variable importance (`bl_robustness()`, `plot.bl_boundary()`)

Each observation has a boundary counterfactual `B_x`. The signed
displacement per feature is `v[j] = X_obs[j] - B_x[j]`. The previous code
standardised each `v[j]` by `X_sd[j]`:

```
vec_sd[i, j] = v[i, j] / X_sd[j]                # legacy
sum_of_distance[j] = sum_i |vec_sd[i, j]|       # per-feature total
robustness = sum_j sum_of_distance[j]           # scalar
```

The new default replaces `X_sd[j]` with `sqrt(W[j, j])` -- the **within-class
standard deviation** of feature `j`. The structure of the measure is
unchanged: it still produces a signed n x p matrix, per-feature totals on
the y-axis of `plot.bl_boundary()`, and a scalar `robustness`. Only the
denominator changes.

### 4.3 Why not full Mahalanobis with per-feature decomposition?

Full Mahalanobis `v^T W^{-1} v` is a quadratic form. Its decomposition is:

```
v^T W^{-1} v = sum_j v[j]^2 * W^{-1}[j, j]                   # diagonal terms
             + 2 * sum_{j < k} v[j] * v[k] * W^{-1}[j, k]   # cross-correlations
```

The diagonal terms can be attributed to features. The cross terms cannot --
they belong to *pairs* of features, not single features. Three options
exist:

| Option | Per-feature? | Cross-correlation? | Cost |
|---|---|---|---|
| **Diagonal Mahalanobis** (what we use) | yes, exact | ignored | `O(n)` |
| **Cholesky-whitened components** | yes, additive | preserved | `O(n)` but loses original-variable identity |
| **Shapley attribution** | yes, principled | distributed fairly | `O(n * 2^p)` |

Cholesky whitening produces additive per-component contributions
(`u_k^2` for `u = L^{-1} v`), but the components are in the *whitened
basis* -- linear combinations of the original features. We can no longer
say "feature `loan_amnt` contributes X units"; we can only say "whitened
direction 3 contributes X units". This defeats the purpose of
`plot.bl_boundary()`, which exists to surface *named-variable*
contributions.

**Shapley attribution** is the theoretically correct answer. Treat the
features as players in a cooperative game, with `v^T W^{-1} v` as the
value function. Compute each feature's Shapley value as its average
marginal contribution across all `2^p` feature subsets. The resulting
attributions `phi_j` satisfy:

- **Efficiency**: `sum_j phi_j = v^T W^{-1} v` (exact decomposition)
- **Symmetry**: features with identical impact get identical attributions
- **Linearity**: combining two value functions combines attributions linearly
- **Null player**: features with no impact get zero attribution

The PhD identified Shapley as the correct method. It is **not implemented
in Phase 2** because the cost is `O(n * 2^p)`: every observation contributes
to the global plot, and each observation's attribution costs `2^p` subset
evaluations. For the loan dataset (n ~ 30,000, p = 6) that is roughly 2
million model scorings just to draw one importance plot. For larger
feature sets (`p > 14`), Shapley becomes infeasible even with the
permutation approximation.

The diagonal approximation is a deliberate trade-off:

- **Cost**: `O(n)` model evaluations, which is the same as the legacy
  `X_sd` measure
- **Interpretation**: preserved exactly (per-feature, named variables)
- **Methodological correction**: captured (within-class denominator)
- **Information loss**: the cross-correlation effects of `W^{-1}` are
  absorbed silently

To make the silent loss visible, `plot.bl_boundary()` prints a
**cross-correlation diagnostic** whenever `distance = "mahalanobis"`:

```
Mahalanobis decomposition (sum across obs):
  diagonal = 12.34 | cross-correlation = 0.85 (+6.4% of total d_M^2)
```

A small cross-correlation percentage means the diagonal approximation
captures essentially all of the full Mahalanobis distance for this data.
A large percentage (|pct| >= 25%) triggers an explicit warning pointing
to this technical note and to `2 implementation_summary.txt` Section
4.6.1 for the theoretical alternatives.

Phase 3 (`bl_shapley()` for single observations) *does* use exact Shapley
because the cost is `O(2^p)` per observation, not `O(n * 2^p)`. For one
observation at `p = 14`, that is 16,384 subset evaluations -- a couple of
seconds, tractable. The same calculation across the full training set
would take 16384 * n seconds, hours to days.

---

## 5. Edge cases and safeguards

1. **Singular `W`.** When `n - g < p` (insufficient effective sample size),
   `W` is rank-deficient and `chol(W)` errors. The package applies a
   ridge-regularised fallback `chol2inv(chol(W + lambda * I))` with
   `lambda = 10^-6 * mean(diag(W))`. The `metric_type` attribute is
   unchanged; the regularisation is transparent to the user.

2. **PCA without class info.** Falls back to `Sigma = cov(X)`. Standard
   Mahalanobis (rather than within-class Mahalanobis), but still a
   principled improvement over `X_sd`-only standardisation.

3. **Older `bl_result` objects.** Before this change, `bl_result` had no
   `metric` field. The new code checks `is.null(bl_result$metric)` and
   falls back to Euclidean with a warning instructing the user to rebuild
   their `bl_result` with the current package version.

4. **External targets** (data frames not in train/test). `bl_target$x_obs`
   is in X-space; no class info is needed at distance-computation time
   because the metric was fixed when `bl_projection` was built.

5. **Backward compatibility.** Passing `distance = "euclidean"`
   reproduces the pre-change behaviour exactly. No existing tests broke.

---

## 6. References within the codebase

All paths below are relative to this file (`documentation/`).

- **Code:**
  - [../R/projection.R](../R/projection.R) -- `.compute_metric_inverse()`
    helper and integration into `bl_build_projection()`
  - [../R/result.R](../R/result.R) -- propagation through `bl_assemble()`
  - [../R/local_cf.R](../R/local_cf.R) -- `distance` parameter and per-pair
    Mahalanobis computation
  - [../R/boundary_plot.R](../R/boundary_plot.R) -- diagonal Mahalanobis
    denominator and cross-correlation diagnostic

- **Documentation:**
  - `../2 implementation_summary.txt` Section 4.2.1 (W metric and Cholesky
    inversion), Section 4.6.1 (Shapley vs diagonal trade-off), Section
    4.8 (cross-pair selector table)
  - `../CLAUDE.md` Section 6 (architectural facts), Section 9 (Shapley
    attribution as remaining deferred work)
  - `../.claude/reference/review_section7_to_8.md` Step 8 (Phase 2
    walkthrough)
  - `../.claude/reference/review_section9_to_15.md` Stage G (Phase 3
    walkthrough)
  - `../.claude/reference/mahalanobis_implementation_plan.md` -- the
    approved plan that drove this implementation

- **Original PhD research:**
  - `../1 Foundation intro documents.txt` Section 10.3 -- the original
    statement that within-class scatter must use the full `Wmat_inv`
    maths
  - `../scripts/1.3 Optimal Rotation.R` (on `original_PhD_code` branch) --
    the original `pc.distz` formula
  - `../scripts/1 Functions for biplot.R` (on `original_PhD_code` branch)
    -- the original `Wmat` construction

---

## 7. Future work explicitly deferred

- **Shapley attribution of full Mahalanobis distance in Phase 2.** The
  theoretically correct per-feature decomposition. Not implemented due to
  `O(n * 2^p)` cost. See `CLAUDE.md` Section 9 and Section 4.6.1 above.

- **Cholesky-whitened diagnostic plot.** A complementary plot showing the
  whitened-basis decomposition (`u_k^2` for `u = L^{-1} v`) could expose
  the structure of cross-correlations across the whole training set,
  without claiming per-original-variable attribution. Useful as a sanity
  check when the cross-correlation diagnostic shows a large percentage,
  but it has not been requested.

- **Mahalanobis as the within-pair selector** in `bl_find_local_cf()`.
  Currently the *within-pair* selection (nearest valid contour vertex)
  remains 2D Euclidean per the user's explicit decision. Could be
  revisited if the within-pair Euclidean choice ever turns out to
  conflict with the cross-pair Mahalanobis choice.
