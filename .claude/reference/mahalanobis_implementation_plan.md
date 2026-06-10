# Plan (IMPLEMENTED 2026-05-20): Mahalanobis distance for boundary logic

> **Status:** Approved and implemented. Plans A and B bundled in one commit.
> **Date implemented:** 2026-05-20
> **Branch:** `method_developments`
> **Companion document:** `documentation/mahalanobis_technical_note.md` --
> theoretical motivation and mathematical derivation in plain prose.
> **Test outcome:** 72 PASS, 0 FAIL (no regressions).
>
> This is the verbatim approved plan, preserved for future reference. The
> technical note linked above explains *why* the change was made; this file
> documents *what* was done and *how*.

---

# Plan: Mahalanobis distance for best-pair selection in `bl_find_local_cf()`

## Context

`bl_find_local_cf()` currently selects the best eigenvector pair by **Euclidean
distance in the rotated 2D Z-space** ([R/local_cf.R:497-503](R/local_cf.R)):

```r
dist_z <- sqrt(sum((Z_target - B_z_local)^2))
...
if (dist_z < best_dist) { best_dist <- dist_z; ... }
```

This is **not mathematically comparable across pairs**:

- Each pair `(i, j)` produces a different SVD rotation `Vrho = V %*% t(A)`. The
  rotation is orthogonal in p-D space but the 2D *slice* `Vrho[, c(1, 2)]` captures
  different amounts of total variance depending on which eigenvectors fed in.
- For CVA, `V` is **`W`-orthonormal**, not Euclidean-orthonormal (`V' W V = I`).
  Euclidean distances in Z-space therefore carry hidden scale distortions that
  differ from pair to pair.
- A 1-unit Euclidean step in the (1, 2) pair is not the same change in feature
  space as a 1-unit step in the (3, 4) pair.

The user's original research code (`1.3 Optimal Rotation.R`, preserved on the
`original_PhD_code` git branch -- see `1 Foundation intro documents.txt`
Section 9) uses the **squared Mahalanobis distance in the original X-space**
with within-class covariance `W`:

```r
pc.distz[i] <- (vec_to_boundary) %*% Wmat_inv %*% t(vec_to_boundary)
```

This is invariant to the eigenvector pair chosen because it lives entirely in
X-space using a fixed, projection-independent metric.

### Where W was originally calculated

In the PhD research scripts, `Wmat` is the pooled within-class scatter matrix
constructed at the same point as the CVA loading matrix:

- **Original source:** `scripts/1 Functions for biplot.R` (research-branch
  helper that builds the biplot inputs) computes
  `Wmat = sum over classes k of crossprod(X_k - colMeans(X_k))`, then
  `Wmat_inv = solve(Wmat / (n - g))` (g = number of CVA class levels).
- **Where it is consumed:** `1.3 Optimal Rotation.R` reuses `Wmat_inv` for the
  per-pair distance comparison (the `pc.distz` formula above), so the *same*
  matrix that defines the CVA eigenproblem also defines the comparison metric
  across rotated pairs. This is what makes the comparison internally
  consistent (V is `W`-orthonormal under exactly the same `W`).
- **Why it is not in the current package:** the package delegates the CVA
  computation entirely to `biplotEZ::CVA()` (`R/projection.R:177-180`).
  biplotEZ does not expose `Wmat` on its return object (confirmed by
  inspection), so the package never had a chance to inherit it. `bl_result`
  currently stores only `V`, `tV`, `X_center`, `X_sd`, `cva_classes` -- the
  raw within-class metric was dropped on the way through.
- **Authoritative caution:** `1 Foundation intro documents.txt` Section 10.3
  explicitly states "CVA within-class scatter must use the full maths
  (`Wmat_inv`), not a simple matrix inverse". The current Phase 3 Euclidean
  distance silently violates this principle. Plan A restores it.

The goal is to investigate mathematically sound options for replicating this
selection criterion inside `bl_find_local_cf()`.

---

## The math

For each candidate pair, after the within-pair search has identified the
nearest boundary point `B_x` (back-projected from the 2D rotated Z-space to
X-space, already computed at [R/local_cf.R:507-511](R/local_cf.R)), define
the displacement:

```
v = B_x - x_obs       # 1 x p row vector in X-space
```

The **squared Mahalanobis distance** with metric matrix `M` (positive-definite,
p x p) is:

```
d_M^2 = v %*% solve(M) %*% t(v)        # scalar
```

A pair is "best" when its `d_M^2` is minimum across all pairs that produced a
feasible counterfactual.

The choice of `M` is the central methodological decision.

---

## Options for the metric matrix `M`

### Option A -- Pooled within-class covariance `W` (CVA-native)

```
W = (1 / (n - g)) * sum over classes k of:  (X_k - mu_k)' (X_k - mu_k)
```

- `g` = number of class levels used
- This is the metric CVA optimises against (CVA solves the generalised eigenproblem
  `(B - lambda W) v = 0`)
- For the package's CVA path, the natural choice is the **same class factor used
  inside `bl_build_projection()`**, i.e. `cva_classes` (confusion labels TP/TN/FP/FN
  when `bl_model` was supplied, or binary 0/1 otherwise). This keeps the metric
  internally consistent with `V`.
- For PCA, `W` is not natural -- need a fallback.

**Pros:**
- Matches the user's existing research code exactly
- Mathematically aligned with CVA's optimisation criterion
- Pair-invariant and scale-invariant; comparable across pairs by construction
- Penalises moves along high within-class-noise directions (good: those moves are
  "easy" only because of noise spread, not because the model genuinely flips)

**Cons:**
- Requires class labels; not directly applicable to pure PCA without a model
- Singular when `n - g < p` (regularisation needed when p is large)
- For the loan dataset with `bl_model` supplied, confusion labels give `g = 4` so
  rank up to `n - 4` -- typically well-conditioned

### Option B -- Total data covariance Sigma (PCA-native)

```
Sigma = cov(X_train)
```

**Pros:**
- Always defined, no class info needed
- Works for any projection method
- Standard Mahalanobis distance, well known

**Cons:**
- Conflates within-class and between-class variance
- For classification problems this is the *wrong* noise model: a move along a
  high-between-class direction is exactly what we want to find, but Sigma
  *down-weights* such moves
- Less aligned with the spirit of the original research code

### Option C -- Identity (Euclidean in standardised X-space)

```
M = diag(X_sd^2)        # equivalent to scaling, no covariance term
```

- Reduces to plain Euclidean distance on standardised features
- Trivially invariant to pair choice (it does not depend on projection at all)
- Loses correlation structure -- a baseline rather than a methodological choice

**Use:** sanity-check baseline only.

### Option D -- Z-space Mahalanobis (per-pair rotated)

Use `cov(Z_train_rot)` per pair to whiten the 2D rotated Z-space, then take
Mahalanobis there. **Reject:** this is *less* sound than Option A. It just
re-derives Euclidean-on-whitened-Z, which is still pair-dependent because the
2D slice differs across pairs.

---

## Recommendation (confirmed with user)

**Implement Option A (W-based Mahalanobis) as the new selection criterion.**

User-confirmed choices:
- **Class factor for W:** CVA-consistent. Use the same `cva_classes` factor that
  was passed to `biplotEZ::CVA()` to build `V` (typically TP/TN/FP/FN, 4 classes
  when `bl_model` was supplied). This keeps the metric internally consistent
  with `V` -- which is `W`-orthonormal under exactly this `W`.
- **Default `distance` argument:** `"mahalanobis"`. This is a behaviour change;
  document in roxygen, `CLAUDE.md` Section 6, and progress notes. Users wanting
  the old behaviour pass `distance = "euclidean"`.

Fallback path for edge cases:

| Projection method | Class info available | Metric used                                                              |
|-------------------|----------------------|--------------------------------------------------------------------------|
| CVA               | `cva_classes` always present       | `W` from `cva_classes` (TP/TN/FP/FN when `bl_model` supplied; binary 0/1 otherwise) |
| PCA               | `train_data$class` present           | `W` from `class` (binary 0/1)                                             |
| PCA               | no class info                       | Sigma = `cov(X_train)` (Option B fallback)                                |

Rationale: in this package, classification class labels are almost always
available (`bl_data` requires them). PCA with no class labels is the exceptional
case; Sigma is the correct fallback there.

**Add a `distance` parameter to `bl_find_local_cf()`** with values
`c("mahalanobis", "euclidean")`, default `"mahalanobis"`. Setting
`distance = "euclidean"` preserves the current behaviour for users who want to
reproduce older results.

**Keep within-pair selection unchanged**: still use 2D Euclidean distance in
the rotated Z-space to pick the nearest valid contour vertex within a pair.
Only the *cross-pair* selection criterion changes. This matches the user's
explicit instruction.

---

## Files to modify

### 1. `R/projection.R` -- compute and store W (or Sigma)

In `bl_build_projection()`, after the existing W/Sigma-relevant data is
available (`X`, `cva_classes`, `train_data$class`), compute the metric matrix
once and store its inverse:

```r
metric_inv <- .compute_metric_inverse(
  X            = X,
  method       = method,
  cva_classes  = cva_classes,
  binary_class = if ("class" %in% names(train_data)) train_data[["class"]] else NULL
)
```

Return fields to add to the `bl_projection` object:
- `metric_inv` -- p x p matrix; the inverse of the chosen metric
- `metric_type` -- one of `"W_cva"`, `"W_binary"`, `"Sigma"`; for transparency

Add a private helper `.compute_metric_inverse()` (probably in `local_cf.R` or
a new `metric_utils.R`):

```r
.compute_metric_inverse <- function(X, method, cva_classes, binary_class) {
  # Pick class factor: CVA labels -> binary class -> NULL
  classes <- NULL
  type    <- NULL
  if (!is.null(cva_classes)) {
    classes <- as.factor(cva_classes); type <- "W_cva"
  } else if (!is.null(binary_class)) {
    classes <- as.factor(binary_class); type <- "W_binary"
  }

  if (is.null(classes)) {
    M <- stats::cov(X)
    type <- "Sigma"
  } else {
    n <- nrow(X); g <- nlevels(classes); p <- ncol(X)
    W <- matrix(0, p, p)
    for (lvl in levels(classes)) {
      Xk <- X[classes == lvl, , drop = FALSE]
      if (nrow(Xk) < 2L) next
      Xc <- sweep(Xk, 2L, colMeans(Xk))
      W  <- W + crossprod(Xc)
    }
    denom <- max(n - g, 1L)
    M <- W / denom
  }

  # Numerical inverse via Cholesky with ridge fallback if singular
  M_inv <- tryCatch(
    chol2inv(chol(M)),
    error = function(e) {
      lambda <- 1e-6 * mean(diag(M))
      chol2inv(chol(M + lambda * diag(nrow(M))))
    }
  )
  attr(M_inv, "metric_type") <- type
  M_inv
}
```

#### Why Cholesky decomposition for the inverse

**What it is.** The Cholesky decomposition factorises a symmetric
positive-definite matrix `M` into `M = L L^T`, where `L` is a unique lower
triangular matrix with positive diagonal entries. Geometrically, `L` is a
"square root" of `M`: it tells you how to whiten the data, since for any
vector `v`, `u = L^{-1} v` satisfies `u^T u = v^T M^{-1} v`. So `L^{-1}` maps
the original (correlated) feature space into an orthonormal space where
distances are exactly Mahalanobis distances under `M`.

**How it relates to `W`.** The within-class covariance matrix `W` is symmetric
by construction and positive-definite whenever it has full rank (which holds
when `n - g >= p`, i.e. there are enough training rows beyond the number of
classes to estimate p-dimensional covariance). Cholesky of `W` therefore
exists and is unique. R provides:

- `chol(M)` -- returns the upper-triangular factor `R` such that `M = R^T R`.
- `chol2inv(R)` -- computes `M^{-1}` directly from `R` by back-substitution.

**Why it is needed (vs `solve(M)`).** Three reasons, in priority order:

1. **Numerical stability.** `solve(M)` uses LU decomposition with partial
   pivoting; it works for any invertible square matrix but does not exploit
   the symmetric positive-definite structure of `W`. `chol2inv(chol(W))`
   is the textbook stable inverse for SPD matrices: it has better
   conditioning and avoids the propagated round-off error that LU can suffer
   when `W` is near-singular. For ill-conditioned `W` (small training sets,
   highly correlated features), `solve()` may return an inverse that
   produces a non-positive `v^T M^{-1} v`, which would silently corrupt the
   pair-selection criterion.

2. **Singularity detection.** `chol()` throws an error the moment `M` is not
   positive-definite -- giving the `tryCatch` block a clean trigger to apply
   the ridge fallback `M + lambda * I`. `solve()` does not detect this
   reliably; it may return a meaningless inverse without warning when the
   matrix is barely-singular.

3. **Speed.** `chol2inv()` is roughly 2x faster than `solve()` for SPD
   matrices because the triangular factor halves the work. Not the primary
   reason here, but a free benefit.

**Why the `chol(...)` call still appears even though we want only the
inverse.** `chol2inv()` accepts the Cholesky factor, not the matrix itself,
so the idiom `chol2inv(chol(M))` first factorises then inverts. We do not
need `L` separately for any other purpose in Plan A (`Plan B` Option D would
have used `L` directly, but is rejected -- see Plan B). For Plan A we only
ever multiply by `M^{-1}`, so retaining only the inverse is sufficient.

### 2. `R/result.R` -- propagate metric_inv into `bl_result`

In `bl_assemble()`, copy `metric_inv` and `metric_type` from `bl_projection`
into the `bl_result` object so downstream functions can read it without going
back to `bl_projection`.

### 3. `R/local_cf.R` -- consume `metric_inv` in `bl_find_local_cf()`

Add a `distance` parameter to `bl_find_local_cf()`:

```r
bl_find_local_cf <- function(bl_result, bl_target,
                             set_filters = NULL,
                             max_pairs   = 10L,
                             m           = 200L,
                             distance    = c("mahalanobis", "euclidean"),
                             verbose     = TRUE) {
  distance <- match.arg(distance)
  ...
}
```

Inside the per-pair loop ([R/local_cf.R:381-544](R/local_cf.R)), after the
nearest valid `B_z_local` and back-projected `B_x` are computed:

```r
# Existing: 2D Euclidean (kept for reporting and within-pair selection)
dist_z <- sqrt(sum((Z_target - B_z_local)^2))

# New: Mahalanobis in X-space
v        <- as.numeric(B_x[1L, var_names]) - x_obs   # 1 x p
M_inv    <- bl_result$metric_inv
dist_mah <- as.numeric(t(v) %*% M_inv %*% v)         # scalar (squared)

selector <- if (distance == "mahalanobis") dist_mah else dist_z

if (selector < best_selector) {
  best_selector <- selector
  ...  # store best_result with both distance fields
}
```

Note: `B_x` is currently only computed *inside* the `if (dist_z < best_dist)`
block. For Mahalanobis we need `B_x` for **every** pair (to compute the
selector). Move the back-projection of `B_z_local` to `B_x` outside the
`if`, into the body of each iteration. Cost: one back-projection (cheap) plus
no additional model scoring (model scoring happens lazily for `B_pred` only on
the winning pair).

Add two new fields to the returned object:
- `dist_mahalanobis` -- scalar, squared Mahalanobis distance at best pair
- `all_distances_mahalanobis` -- named numeric vector, one per pair tried
- (Existing `dist_z` and `all_distances` remain, populated for both modes)

### 4. Documentation updates

- `R/local_cf.R` roxygen for `bl_find_local_cf()`: document the new `distance`
  parameter and the methodological reasoning.
- `R/projection.R` roxygen for `bl_build_projection()`: document the new
  `metric_inv` / `metric_type` return fields.
- `2 implementation_summary.txt` Section 4.2 and 4.8: note the metric matrix is
  computed during projection and consumed during local CF selection. See the
  drafted subsection below for the Cholesky-decomposition explanation that
  must be added.
- `.claude/reference/review_section9_to_15.md` Stage G ("Update best result"):
  describe the new Mahalanobis selector.
- `CLAUDE.md` Section 6 ("Key Architectural Facts"): add bullet that
  `bl_find_local_cf()` selects the best pair via Mahalanobis distance in X-space
  by default; `dist_z` is reported but not the selector.

#### New subsection in `2 implementation_summary.txt`: the W metric matrix and Cholesky inversion

Add a new subsection to `2 implementation_summary.txt` (location: end of
Section 4.2 "Module B -- Projection and Biplot Construction", before Section
4.3) with the following content, drafted to match the existing prose style of
the document:

> **4.2.1 The within-class metric matrix W and Cholesky inversion**
>
> `bl_build_projection()` computes and stores a metric matrix `M` and its
> inverse `M_inv` for use by Phase 2 and Phase 3 distance measures. The
> matrix is:
>
> - `W` (pooled within-class covariance) when `cva_classes` or
>   `train_data$class` is available -- the standard CVA scatter matrix
>   `W = sum_k (X_k - mu_k)^T (X_k - mu_k) / (n - g)` where `g` is the
>   number of class levels and `mu_k` is the per-class feature mean.
> - `Sigma` (total covariance, `cov(X)`) as a fallback when no class
>   information is supplied (rare in this package; `bl_data` almost always
>   has classes).
>
> The metric is the same one CVA itself optimises against (V is
> `W`-orthonormal under exactly this `W`). Storing it makes the
> projection-to-distance pipeline internally consistent: distances measured
> in Phase 2 or Phase 3 use the same `W` that defined the biplot axes in
> Phase 1.
>
> **Cholesky decomposition for the inverse.** `W` is symmetric and
> positive-definite by construction, so it admits a unique factorisation
> `W = L L^T` where `L` is lower triangular with positive diagonal entries.
> This is the **Cholesky decomposition** -- the matrix equivalent of taking
> a square root. Geometrically, `L^{-1}` is the whitening transform that
> sends correlated data into an uncorrelated unit-variance frame, so that
> Euclidean distance in the whitened frame equals Mahalanobis distance
> under `W` in the original frame.
>
> The package inverts `W` using the idiom `chol2inv(chol(W))` rather than
> `solve(W)`, for three reasons:
>
> 1. **Numerical stability.** Cholesky exploits the symmetric
>    positive-definite structure that `solve()` cannot use; it propagates
>    less round-off error and avoids producing a non-positive
>    `v^T W^{-1} v` (which would silently corrupt distance measures) when
>    `W` is near-singular.
> 2. **Singularity detection.** `chol()` errors immediately if `W` is not
>    positive-definite -- giving a clean trigger to apply the ridge
>    fallback `W + lambda * I` (with `lambda = 10^-6 * mean(diag(W))`).
>    `solve()` would return a meaningless inverse without warning.
> 3. **Speed.** `chol2inv()` is roughly 2x faster than `solve()` for
>    SPD matrices; not the primary reason but a free benefit.
>
> The `chol(W)` factor itself is not retained -- only `W` (as `metric`) and
> `W^{-1}` (as `metric_inv`) are stored on `bl_projection`. The Cholesky
> step is purely an internal numerical detail of the inversion. Downstream
> code multiplies vectors by `metric_inv` directly.
>
> **Singular-W safeguard.** When `n - g < p` (more features than effective
> sample size after class adjustment), `W` is rank-deficient and Cholesky
> fails. The package applies a ridge-regularised fallback
> `chol2inv(chol(W + lambda * I))` so the distance measures remain
> well-defined. The `metric_type` attribute is unchanged; the
> regularisation is transparent to the user but logged.

### 5. CLAUDE.md Section 9 -- remove the Mahalanobis deferred item

The deferred reference to `memory/future_mahalanobis_distance.md` becomes
obsolete once this is implemented.

---

## Edge cases and safeguards

1. **Singular W** (`n - g < p`): caught by `tryCatch` in
   `.compute_metric_inverse()`; falls back to ridge-regularised inverse
   `(W + lambda*I)^-1` with `lambda = 1e-6 * mean(diag(W))`.

2. **External target with no class info**: `bl_target$x_obs` is in X-space; no
   class info is needed at distance-computation time -- the metric was already
   fixed when `bl_projection` was built.

3. **Numerical scale**: Mahalanobis distances can be tiny (W is in
   variance-units of the features). Report `sqrt(dist_mahalanobis)` to the user
   in `print()` and progress messages, with units noted as "Mahalanobis (W^{-1/2}
   feature-units)".

4. **Backward compatibility**: with `distance = "euclidean"`, behaviour is
   identical to current. No existing tests should break.

5. **`best_pair` may change**: under `distance = "mahalanobis"` the winning pair
   may differ from the Euclidean choice for many inputs. This is the desired
   behaviour. Document in the roxygen `@section`.

---

## Verification

```r
devtools::load_all()
devtools::test()
# Expect: 72 PASS, 0 FAIL (existing tests use default args; new default is
# "mahalanobis" so this is a behaviour-changing default. If any test relies on
# a specific best_pair, may need to pass distance = "euclidean" to that test.)

# --- Unit test for .compute_metric_inverse() ---
# - Verify it returns a positive-definite p x p matrix
# - Verify W = (1/(n-g)) * sum_k crossprod(X_k - mu_k) computed correctly
#   against a hand-rolled scalar example on iris

# --- Iris (CVA + GLM) side-by-side ---
bl_dat  <- bl_prepare_data(datasets::iris, class_col = "Species",
                            target_class = "versicolor")
bl_mod  <- bl_fit_model(bl_dat$train_data, bl_dat$var_names)
bl_proj <- bl_build_projection(bl_dat$train_data, bl_dat$var_names,
                                method = "CVA", bl_model = bl_mod)
bl_grid <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod)
bl_results <- bl_assemble(bl_dat, bl_model = bl_mod,
                           bl_projection = bl_proj, bl_grid = bl_grid)

# Check metric stored
stopifnot(!is.null(bl_results$metric_inv))
cat("Metric type:", attr(bl_results$metric_inv, "metric_type"), "\n")

# Compare Euclidean vs Mahalanobis selection
tgt   <- bl_select_target(bl_results, target = 1L)
bl_e  <- bl_find_local_cf(bl_results, tgt, max_pairs = 10L,
                           distance = "euclidean", verbose = FALSE)
bl_m  <- bl_find_local_cf(bl_results, tgt, max_pairs = 10L,
                           distance = "mahalanobis", verbose = FALSE)

cat("Euclidean   best pair:", bl_e$best_pair,
    "  dist_z =", bl_e$dist_z,
    "  dist_mah =", bl_e$dist_mahalanobis, "\n")
cat("Mahalanobis best pair:", bl_m$best_pair,
    "  dist_z =", bl_m$dist_z,
    "  dist_mah =", bl_m$dist_mahalanobis, "\n")

# Plot both for visual comparison
plot(bl_e); plot(bl_m)

# --- Loan dataset (XGB) ---
# Run scripts/03_loan_status_Boundary_Logic.R Steps 11-13 with both modes;
# confirm that downstream Phase 3 (bl_shapley, bl_find_sparse_cf) still
# functions with the Mahalanobis-selected pair.
```

---

## Out of scope (deferred)

- Mahalanobis-based **within-pair** vertex selection (currently 2D Euclidean
  on contour vertices). Per user instruction, keep within-pair selection
  Euclidean for now.
- Mahalanobis as an option for `bl_find_boundary()` (global Phase 2 search).
  That function operates in the original (unrotated) Z-space with no pair loop,
  so the comparability issue does not arise there.

---
---

# Plan B (separate): Mahalanobis-aware per-variable distance for `plot.bl_boundary()` and `bl_robustness()`

## Context

The Phase 2 distance-to-boundary "global importance measure" is computed in
[R/boundary_plot.R:77-87](R/boundary_plot.R) and
[R/boundary_plot.R:206-214](R/boundary_plot.R) as:

```r
vec_to_boundary_Z  <- Z_obs - B_z                         # n x 2 (Z-space)
vec_to_boundary    <- vec_to_boundary_Z %*% tVr           # n x p (X-space)
if (standardise) vec_to_boundary <- sweep(., 2, X_sd, "*")
vec_to_boundary_sd <- sweep(vec_to_boundary, 2, X_sd, "/")  # per-feature standardise
sum_of_distance    <- colSums(abs(vec_to_boundary_sd))    # per-variable importance
robustness         <- sum(sum_of_distance)                # scalar
```

The denominator `X_sd` (total per-feature standard deviation) is used to make
variables comparable across different scales. The per-variable totals appear
on the y-axis of `plot.bl_boundary()` and drive the variable ordering.

**The methodological issue** -- the same critique as Plan A:

- `X_sd` is the **total** SD per feature: `var_j = E[(X_j - E X_j)^2]`
- For a classification problem, the relevant noise is **within-class** SD:
  `sqrt(W[j, j])` where `W` is the pooled within-class covariance
- Total SD inflates the denominator for features that are good class separators
  (because their between-class variance is large). The per-variable importance
  measure therefore *under-weights* the features most useful for crossing the
  boundary -- exactly the opposite of what an "importance" measure should do.

**The constraint** (user's explicit requirement): keep the **per-variable
contribution structure**. A single scalar Mahalanobis distance loses the
decomposition that makes `plot.bl_boundary()` informative.

This makes full Mahalanobis distance (`v' W^{-1} v`) inadmissible as a direct
replacement, because the quadratic form has cross-terms that don't decompose
cleanly into "feature j's contribution".

The investigation: what variants of the Mahalanobis idea preserve a per-feature
breakdown?

---

## Options

### Option A -- Diagonal within-class SD (RECOMMENDED)

Replace `X_sd` with `sqrt(diag(W))`:

```r
w_sd       <- sqrt(diag(W))                       # length p
vec_per_sd <- sweep(vec_to_boundary, 2L, w_sd, "/")   # n x p
sum_of_distance <- colSums(abs(vec_per_sd))           # per-variable
robustness      <- sum(sum_of_distance)               # scalar
```

This is the *diagonal-only* part of the Mahalanobis metric -- equivalent to
assuming the off-diagonal correlations of `W` are zero for the purposes of
standardisation.

**Pros:**
- **Drop-in replacement**: identical structure to the current code, only the
  denominator changes. No restructuring of `plot.bl_boundary()` or
  `bl_robustness()`.
- **Per-variable interpretation preserved exactly**: "feature j moves
  `|v_ij| / sqrt(W[j,j])` within-class SDs to reach the boundary".
- **Fixes the methodological issue**: features that are good class separators
  (large between-class variance, small within-class variance) now get
  *amplified* in the importance measure -- correctly.
- **Reuses infrastructure from Plan A**: the same `W` matrix is computed once
  and stored in `bl_result`. `Plan A` uses `W^{-1}` (full inverse); `Plan B`
  uses `sqrt(diag(W))` -- both derived from the same `W`.
- **Aggregates cleanly**: `sum_of_distance` and `robustness` keep their
  meanings and dimensional units.

**Cons:**
- Ignores cross-feature correlations (same limitation as the current `X_sd`).
  The cross-correlation effects are absorbed silently; not made explicit.
- Not the *full* Mahalanobis distance; it is a diagonal approximation.

### Option B -- Quadratic decomposition with explicit cross-term

The full squared Mahalanobis distance decomposes as:

```
d_M^2(v) = v' W^{-1} v
        = sum_j v_j^2 * W^{-1}[j, j]                 # "purely feature j"
        + 2 * sum_{j<k} v_j * v_k * W^{-1}[j, k]     # cross-correlations
        = sum_j c_j + r                              # per-feature + single residual
```

Per-feature: `c_j = v_j^2 * [W^{-1}]_{jj}`.
Cross-term residual: `r = 2 * sum_{j<k} v_j v_k * [W^{-1}]_{jk}`.

**Pros:**
- **Mathematically exact**: `sum_j c_j + r == d_M^2`. No information lost.
- Reports both per-feature contributions and a single "interaction" scalar.
- Reveals when cross-correlations are large enough to matter.

**Cons:**
- The cross-term `r` can be positive or negative; not interpretable as
  "feature X's importance".
- `c_j` is in units of "squared distance", not "absolute distance" -- the
  signed visual cue of the current plot is lost (or has to be reintroduced as
  `sign(v_j) * sqrt(c_j)`, which loses additivity).
- Requires `W^{-1}`, not just `diag(W)` -- extra cost (already paid by Plan A).
- More to explain to users; less drop-in.

### Option C -- Shapley attribution of full Mahalanobis distance

Use Shapley values to attribute `d_M^2(v)` to each feature in a principled,
correlation-aware way. The `bl_shapley` framework (`R/shapley.R`) is already in
the package and can be adapted: the "game" becomes `d_M^2` as the value
function, with feature subset S having `v` values for `j in S` and zero
otherwise (or `v_j = 0` represents "feature not engaged in the move").

Per-feature Shapley `phi_j` satisfies the efficiency axiom: `sum_j phi_j ==
d_M^2(v)`.

**Pros:**
- **Most principled**: respects full covariance structure including cross-terms,
  while producing exactly per-feature contributions.
- Cross-correlations are *distributed* across features rather than lumped into
  a separate scalar.
- Reuses `bl_shapley` machinery for code coherence.

**Cons:**
- **Expensive**: exact computation requires 2^p subset evaluations per
  observation. For p <= 14 (the current `exact_max_vars` default), tractable
  but slower than Option A. Permutation approximation works for larger p.
- More complex to explain and document.
- Per-feature contributions are not always non-negative under Shapley
  (a feature may have negative attribution if its move "undoes" another's).
- Less directly comparable to the current `sum_of_distance` semantics.

### Option D -- Whitened component decomposition (Cholesky)

The Cholesky factor `L` of `W` (where `W = L L^T`, see Plan A's Cholesky
explanation block) gives a linear map `L^{-1}` that whitens the data: the
transformed displacement `u = L^{-1} v` lives in a coordinate system where
each axis is uncorrelated and has unit within-class variance. In that whitened
basis, the full squared Mahalanobis distance decomposes exactly as a sum of
squares:

```
d_M^2(v) = v^T W^{-1} v
       = v^T L^{-T} L^{-1} v
       = (L^{-1} v)^T (L^{-1} v)
       = u^T u
       = sum_k u_k^2
```

So each whitened component `u_k^2` is a clean additive contribution to the
total Mahalanobis distance. This is the mathematically cleanest decomposition.

**Reject for this task**: the components `u_k` correspond to the *whitened
basis* (rotated and rescaled mixtures of the original features), not the
original feature names. We can no longer tell the user "feature `loan_amnt`
contributes X units" -- only "whitened component 3 contributes X units". The
y-axis of `plot.bl_boundary()` would no longer carry interpretable variable
names, defeating the purpose of the per-variable importance measure.

---

## Recommendation

**Implement Option A** as the primary new measure -- the diagonal within-class
SD replacement of `X_sd`. It satisfies the user's constraint exactly
(per-variable contribution preserved), fixes the methodological issue
(within-class noise vs total variance), and reuses Plan A's `W` infrastructure
at zero marginal cost.

**Expose Option B as a secondary report**: include the cross-term residual `r`
in the printed console summary so users can see whether off-diagonal
correlations matter for their data. If `|r|` is small relative to
`sum_j c_j`, the diagonal approximation is essentially lossless.

---

## Implementation

### 1. Reuse Plan A's `bl_projection$metric_inv` and add `metric` matrix itself

Plan A computes `metric_inv` (i.e. `W^{-1}`). Plan B needs `diag(W)` (i.e. the
non-inverted diagonal). To support both cheaply, also store the metric matrix
itself, not just its inverse, in `bl_projection` / `bl_result`:

- `metric` -- p x p matrix `W` (or Sigma fallback)
- `metric_inv` -- p x p matrix `solve(W)` (Plan A)
- `metric_type` -- `"W_cva"`, `"W_binary"`, or `"Sigma"`

`diag(metric)` is then trivially available for Plan B.

### 2. Modify `R/boundary_plot.R`

Add a `distance` parameter to both functions, default `"mahalanobis"` (matches
Plan A's default convention):

```r
plot.bl_boundary <- function(x, type = c("jitter", "boxplot"),
                              distance = c("mahalanobis", "euclidean"),
                              ...) {
  distance <- match.arg(distance)
  ...
  denom <- if (distance == "mahalanobis") {
    sqrt(diag(bl_result$metric))           # within-class (or Sigma) SD per feature
  } else {
    bl_result$X_sd                          # current behaviour
  }
  vec_per_var <- sweep(vec_to_boundary, 2L, denom, "/")
  sum_of_distance <- colSums(abs(vec_per_var))
  ...
}
```

Same change to `bl_robustness()`.

### 3. Console output: report cross-term diagnostic (Option B as a check)

In `plot.bl_boundary()` console summary, when `distance = "mahalanobis"`:

```r
W_inv     <- bl_result$metric_inv
diag_sum  <- sum(vec_to_boundary^2 %*% diag(diag(W_inv)))   # diagonal d_M^2
total_sum <- sum(diag(vec_to_boundary %*% W_inv %*% t(vec_to_boundary)))
cross     <- total_sum - diag_sum
cat(sprintf(
  "Mahalanobis decomposition: diagonal = %.2f | cross-correlation = %.2f (%+.1f%%)\n",
  diag_sum, cross, 100 * cross / total_sum
))
```

A small cross-correlation percentage tells the user the diagonal approximation
is accurate. A large percentage signals to consider Shapley attribution.

### 4. Documentation

- `R/boundary_plot.R` roxygen: document the new `distance` parameter, default,
  and methodological reasoning. Note the cross-correlation diagnostic.
- `2 implementation_summary.txt` Section 4.6: describe both distance options
  and how to read the cross-correlation diagnostic.
- `.claude/reference/review_section7_to_8.md`: update Step 8 walkthrough.
- `CLAUDE.md` Section 6: add bullet that `bl_robustness()` and
  `plot.bl_boundary()` use within-class SD by default (Mahalanobis-diagonal).

#### New subsection in `2 implementation_summary.txt`: Shapley attribution -- theoretical correctness vs. computational cost

Add a new subsection to `2 implementation_summary.txt` (location: end of
Section 4.6 "Module E -- Distance and Robustness Analysis", before Section
4.7) with the following content, drafted to match the existing prose style of
the document:

> **4.6.1 Shapley attribution: theoretical correctness vs. computational cost**
>
> The per-variable importance measure produced by `bl_robustness()` and
> `plot.bl_boundary()` uses the *diagonal* of the within-class covariance
> matrix `W` to standardise per-feature distances to the boundary. This
> diagonal approximation preserves a clean per-feature interpretation but
> ignores cross-feature correlations encoded in the off-diagonal of `W^{-1}`.
>
> The **theoretically correct** decomposition of the full squared Mahalanobis
> distance `d_M^2(v) = v^T W^{-1} v` into per-feature contributions is the
> **Shapley value attribution** from cooperative game theory: treat the p
> features as players and `d_M^2` as the value function, then compute each
> feature's Shapley value as its average marginal contribution across all
> 2^p subsets. The resulting per-feature attributions `phi_j` satisfy the
> efficiency axiom `sum_j phi_j = d_M^2(v)`, distribute the cross-correlation
> structure of `W^{-1}` across features in a uniquely fair way, and respect
> linearity, symmetry, and null-player axioms.
>
> **This was identified in the PhD research as the methodologically correct
> per-feature attribution.** It is **not implemented in the package** because
> the computational cost is prohibitive for the global Phase 2 use case:
> exact Shapley requires 2^p subset evaluations per observation, and Phase 2
> distance-to-boundary measures aggregate across all training/test
> observations -- so the total cost is `O(n * 2^p)` model evaluations. For
> the loan dataset (`n approx 30,000`, `p = 6` after pruning), this is
> `n * 64 approx 2,000,000` model scorings just to build the importance plot.
> For larger feature sets (`p > 14`) the cost becomes infeasible even with
> the permutation approximation used elsewhere in the package
> (`R/shapley.R`).
>
> The diagonal approximation was chosen instead as a deliberate trade-off:
> linear cost (`O(n)` model evaluations), preserved per-variable
> interpretation, and the main methodological correction (within-class noise
> instead of total variance) achieved without the combinatorial expense. A
> cross-correlation diagnostic (printed by `plot.bl_boundary()`) reports
> what fraction of `d_M^2` is captured by the diagonal -- giving users a
> direct measure of how much information is foregone for their specific
> data. Where this fraction is small, the diagonal approximation is
> essentially lossless and the cost saving is "free".
>
> Phase 3 (`bl_shapley`) does use the exact Shapley formulation, but only
> for a single observation at a time -- the cost is `2^p` rather than
> `n * 2^p`, which is tractable for the local interpretation use case.

---

## Verification

```r
devtools::load_all()
devtools::test()

# Loan dataset -- compare per-variable importance under both measures
source("scripts/03_loan_status_Boundary_Logic.R")   # through bl_find_boundary()
old <- bl_robustness(bl_bnd)                              # legacy (X_sd)
new <- bl_robustness(bl_bnd, distance = "mahalanobis")    # new default

# Compare per-variable totals -- the ordering of variables in the plot may shift
old$sum_of_distance |> sort(decreasing = TRUE) |> head()
new$sum_of_distance |> sort(decreasing = TRUE) |> head()

# Visual side-by-side
plot(bl_bnd, distance = "euclidean")     # current behaviour
plot(bl_bnd, distance = "mahalanobis")   # new default

# Check cross-correlation diagnostic in console output
# A small percentage validates the diagonal approximation
```

Expected behavioural change: features that are good class separators (large
between-class variance) move *up* the importance ranking under Mahalanobis;
features that are noise-dominated move *down*.

---

## Synergy with Plan A and bundling decision

Plans A and B share infrastructure:

1. Both require computing `W` from `cva_classes` (or binary class, or Sigma
   fallback) in `bl_build_projection()`.
2. Both require propagating the `W` matrix and its inverse through
   `bl_assemble()` into `bl_result`.
3. Both add a `distance = c("mahalanobis", "euclidean")` parameter to a Phase 2
   or Phase 3 function with `"mahalanobis"` as the new default.
4. Both have the same methodological motivation: total SD is the wrong noise
   model for classification; within-class SD (Mahalanobis) is the right one.

**Recommended bundling**: implement Plans A and B as one commit. The shared
`bl_projection$metric`, `bl_projection$metric_inv`, and helper
`.compute_metric_inverse()` are written once and consumed by both. Splitting
them would require touching `bl_projection` twice and two rounds of
documentation updates.
