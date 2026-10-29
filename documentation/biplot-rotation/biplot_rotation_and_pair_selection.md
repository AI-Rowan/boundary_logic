# Local Biplot Rotation and Optimal Pair Selection

This note explains, in full technical depth, the two mechanisms `boundarylogic` uses to make a
2D biplot locally consistent with a model's prediction for one specific observation: (Part A) the
Procrustes-based rotation that forces a target point onto the biplot plane, and (Part B) the
search across candidate 2D projections that decides *which* plane to rotate. It reproduces the
governing mathematics from the original PhD thesis, maps every step onto the exact production
code, and closes with a dedicated section of known limitations that any user of the method should
understand before treating its output as more precise than it is.

Source material: PhD thesis, "Optimal rotation of the biplot" (§3.5–3.5.2); package source
`R/local_cf.R`; the internal technical note `documentation/mahalanobis_technical_note.md`; package
design documents
`1 Foundation intro documents.txt` and `2 implementation_summary.txt`. Full reference list at the
end. (The PhD thesis excerpt itself was a `.docx` in `documentation/`, read in full while writing
this note and since removed from the repo once its content was consolidated here.)

**Files in this folder.** Every measurement quoted in §5 and §6 is reproducible from the R scripts
sitting beside this note. `README.md` has the details; in brief:

| script | setting | backs |
|---|---|---|
| `17_1a_toy-2d-to-1d-rotation.R` | 5 points, 2 variables → 1D biplot | §2.5, §5.1, §5.8 |
| `17_1_loan-gam-pca-3d-3var-rotation.R` | 221 loans, 3 features, GAM — lossless | §5.1–§5.7 |
| `17_2_loan-gam-pca-4var-slices.R` | 221 loans, 4 features, GAM — lossy | §5.1, §5.4, §5.9, §6.1 |
| `17_3_minimal-rotation-verification.R` | both, console output only | §6 throughout |

Run `17_3` first: it asserts the geometry §6 depends on and reports the rest. Paths *inside* the
scripts point at each other and at this note; paths of the form `R/local_cf.R`,
`development/…` and `documentation/…` refer to the `boundarylogic` package repository.

**Section 5 (added 2026-09) collects measurements** taken while building the loan-data
worked examples. It quantifies the §4.1 caveat, records where §4.1's wording does not
generalise below $p=4$, and adds the interpretive point the published figures most need:
a 2D biplot is a **slice** of the decision surface, not a view of it. **§5.9 (added
2026-09-10) is the one to read first if you are about to quote a counterfactual:** the
reported distance depends on the session's call history, so "repeat the call and check
it is stable" is not sufficient advice.

**Section 6 (added 2026-09) is the recommendation** the §4 limitations point towards: a closed-form
deterministic replacement for the Procrustes rotation, and a reproducible search over the family of
planes that contain the target, with the measurements §4.1 asks for carried out.

---

## 1. Purpose and motivation

A 2D (or low-rank) biplot is built once, from the whole training set, using either PCA or CVA. It
is a *global* approximation: every point is projected onto the same fixed plane, and the decision
boundary drawn on that plane is the model's boundary as seen through that one projection.

For any individual observation, projection and prediction are, in general, **not the same
operation**. Projecting a point onto the biplot plane and reading off its position relative to the
drawn boundary can disagree with what the actual classifier `f(X)` says about the *original,
full-dimensional* point — because the plane is an approximation that only captures the two
directions of greatest global variance (PCA) or greatest class separation (CVA); it does not
generally contain any particular observation's exact location relative to the true, possibly
curved, decision surface in the full feature space.

This is not a hypothetical failure mode. The thesis illustrates it with the iris dataset: a
three-dimensional PCA biplot (three flower features, no dimensionality loss) correctly places a
target point in the "class 1" region, matching the true model prediction. Reducing to the
standard two-dimensional biplot (dropping the third, lowest-variance component) projects the same
point orthogonally onto the plane — and it lands in the "class 0" region. The model is right; its
2D visual approximation is wrong for this specific point, because the decision surface is
oriented obliquely relative to the biplot plane near that point.

For an explanation method whose whole purpose is to visualise *why* a specific observation was
classified the way it was, this is a serious problem: the picture can visibly contradict the
number. The rotation described in Part A exists to fix exactly this — for the one point you are
trying to explain.

---

## 2. Part A — Procrustes rotation of a given plane

### 2.1 The idea and the three-point construction

The rotation is deliberately framed as a **Procrustes problem**: the optimal one-to-one
transformation that maps points in one configuration onto points in another (Cox & Cox, 2008).
Ordinarily Procrustes analysis aligns two whole point clouds. Here it is specialised to a single
point: the target observation $x_{1\times p}$ (already centred, and — for a PCA biplot —
potentially scaled).

To pose the alignment as a Procrustes problem, the single point is embedded as **three points**,
symmetric about the origin:

$$
\ddot{\mathbf{X}}_{3\times p} =
\begin{bmatrix}
-x_1 & \cdots & -x_p \\
0    & \ddots & 0 \\
x_1  & \cdots & x_p
\end{bmatrix}
$$

The middle row is the origin (a fixed point under any rotation) and the two outer rows are $-x$
and $x$. Aligning this three-point configuration in the full $p$-dimensional space with its
counterpart in the reduced space is mathematically equivalent to aligning the single vector $x$
itself — the origin contributes nothing and $-x$ is redundant with $x$ — but casting it as three
points lets the rotation be solved with the standard Procrustes machinery.

### 2.2 Two projections of the same three points

$\ddot{\mathbf{X}}$ is projected two ways, using the same loading matrix $\mathbf{V}$ (the p×p
matrix from the SVD of the centred, and possibly scaled, training data $\mathbf{X}_{n\times p}$,
already used to build the ordinary biplot):

$$
\mathbf{Y} = \ddot{\mathbf{X}}\mathbf{V}
\qquad\text{(stays in the full $p$-dimensional PCA space)}
$$

$$
\mathbf{Z} = \ddot{\mathbf{X}}\mathbf{V}\mathbf{J}^{*} = \ddot{\mathbf{X}}\mathbf{V}_r
\qquad\text{(reduced to the $r$-dimensional biplot space, $\mathbf{V}_r$ = the columns of $\mathbf{V}$ used for the biplot)}
$$

If $r < p$, $\mathbf{Z}$ is zero-padded back out to $p$ columns so the two matrices are the same
shape. The thesis illustrates this for $r=1,\,p=2$:

$$
\mathbf{Z}_{3\times p} = \begin{bmatrix} -z_1 & 0 \\ 0 & 0 \\ z_1 & 0 \end{bmatrix}, \qquad
\mathbf{Y}_{3\times p} = \begin{bmatrix} -y_1 & -y_2 \\ 0 & 0 \\ y_1 & y_2 \end{bmatrix}
$$

$\mathbf{Y}$ is "where the target really sits" (full information); $\mathbf{Z}$ is "where the
biplot currently puts it" (reduced, lossy). The rotation's job is to find an orthogonal
transformation of $\mathbf{V}$ that makes these agree exactly for the target.

### 2.3 The Procrustes solution

$$
\mathbf{C} = \mathbf{Y}'\mathbf{Z}, \qquad
\mathbf{C} = \mathbf{U}\boldsymbol{\Lambda}\mathbf{V}^{*\prime}\ \text{(SVD)}, \qquad
\mathbf{A} = \mathbf{V}^{*}\mathbf{U}'
$$

($\mathbf{V}^{*}$ here is the SVD's right-singular-vector matrix — an artifact of the SVD
notation, unrelated to the loading matrix $\mathbf{V}$.) $\mathbf{A}$ is the optimal $p\times p$
orthogonal rotation matrix. The optimal dilation is

$$
\rho = \frac{\operatorname{tr}(\mathbf{Z}\mathbf{A}\mathbf{Y}')}{\operatorname{tr}(\mathbf{Z}'\mathbf{Z})}
$$

**The dilation is computed but deliberately not applied.** Using it would uniformly rescale
*every* point on the biplot plane in an attempt to recover magnitude information that dimension
reduction genuinely lost — a correction that has nothing to do with the target point and would
distort the axes for no benefit. Only the rotation $\mathbf{A}$ is used.

### 2.4 The new projection

$$
\mathbf{V}^{new} = \mathbf{V}\mathbf{A}'
$$

Biplots built from $\mathbf{V}^{new}$ are called **rotated component biplots (RC biplots)**. The
defining property — the reason for doing any of this — is that the target point now lies exactly
on the plane:

$$
\hat{\mathbf{x}}_i' = \mathbf{x}_i'\mathbf{V}_r^{new}(\mathbf{V}_r^{new})' = \mathbf{x}_i'
$$

and consequently the model's classification and the biplot-implied classification agree:

$$
f(\hat{\mathbf{x}}_i) = f(\mathbf{x}_i)
$$

The corresponding prediction (back-projection) matrix becomes $(\mathbf{V}^{new})' = \mathbf{A}\mathbf{V}'$,
and a point $z$ in rotated Z-space is mapped back to X-space via
$\hat{\mathbf{x}}' = \mathbf{z}'\mathbf{A}\mathbf{V}'\mathbf{J}'$.

### 2.5 Worked examples from the thesis

**5-point, 2D→1D toy example.** Five points in two dimensions are reduced to a one-dimensional
PCA biplot. Point 5 is the target. Before rotation its sample predictivity — how well the biplot
represents it — is 0.72; a 28° rotation with dilation $\rho = 1/1.13$ brings it to exactly 1
(perfect representation). The cost: the sum of squared residuals **for the dataset as a whole**
rises from 0.28 to 4 (out of a total sum of squares of 14) — every other point's approximation
gets slightly worse so that the target's becomes exact. Variable 2's axis, which happens to be
nearly parallel to the rotated plane, gains axis predictivity at the direct expense of variable 1,
which becomes nearly perpendicular to it.

**Iris 3D→2D example.** A three-variable PCA biplot, first shown at full rank (3D, no
approximation error at all) and then reduced to the usual 2D biplot. As described in §1, the
reduction misclassifies a specific target relative to the model's actual (correct) prediction,
because the point sits obliquely above the 2D plane. After rotating so the plane passes through
that point, the projected position moves into the correct (matching) coloured region, and the
projected/interpolated value becomes numerically identical to the true value — this is the
$\hat{\mathbf{x}}_i' = \mathbf{x}_i'$ identity from §2.4, demonstrated rather than just asserted.

### 2.6 The CVA variant (§3.5.2)

The same rotation machinery applies to a CVA transformation matrix. For CVA, the "loading matrix"
role is played by $\mathbf{M}$ (or the two-group optimised transform $\mathbf{M}^{opt}$), and the
rotated version is

$$
\mathbf{M}^{new} = \mathbf{M}^{opt}\mathbf{A}', \qquad
(\mathbf{M}^{new})^{-} = \mathbf{A}(\mathbf{M}^{opt})^{-1}
$$

This guarantees the target lands on the CVA biplot plane exactly as before, but it has a cost that
is specific to CVA and does not arise with PCA: **the class mean points are pulled off the
canonical plane by the rotation**, so their $(k-1)$-dimensional predictions — which CVA is
specifically constructed to represent exactly — stop being exact after rotation. PCA rotation has
no analogous "special point" that gets sacrificed; CVA rotation trades away class-mean fidelity to
buy target-point fidelity.

### 2.7 Code mapping — `.bl_rotate()`

`R/local_cf.R:45-66` is a direct, line-for-line port of the derivation above (originally adapted
from the research script `scripts/1.3 Optimal Rotation.R`, which itself calls this function
`Biplot_rotation()`):

```r
.bl_rotate <- function(x_target_st, V, proj_pair) {
  p    <- length(x_target_st)
  Vr   <- V[, proj_pair, drop = FALSE]           # p x 2

  Y    <- rbind(-x_target_st, rep(0, p), x_target_st)   # 3 x p
  YV   <- Y %*% V                                        # 3 x p
  YVr  <- Y %*% Vr                                       # 3 x 2

  # Pad YVr with zeros so it is 3 x p
  YVr_padded <- cbind(YVr, matrix(0, nrow = 3L, ncol = p - 2L))

  svd_res <- svd(t(YV) %*% YVr_padded)           # p x p SVD
  A       <- svd_res$v %*% t(svd_res$u)           # p x p rotation

  Vrho  <- V %*% t(A)                             # rotated full loading matrix
  tVrho <- solve(Vrho)

  list(
    Vr_rot  = Vrho[, c(1L, 2L), drop = FALSE],    # p x 2 -- SVD always puts target info in cols 1-2
    tVr_rot = tVrho[c(1L, 2L), , drop = FALSE]    # 2 x p
  )
}
```

Term-by-term correspondence with §2.1–2.4:

| Thesis | Code |
|---|---|
| $\ddot{\mathbf{X}}$ | `Y <- rbind(-x_target_st, rep(0, p), x_target_st)` |
| $\mathbf{Y} = \ddot{\mathbf{X}}\mathbf{V}$ | `YV <- Y %*% V` |
| $\mathbf{Z} = \ddot{\mathbf{X}}\mathbf{V}_r$, zero-padded | `YVr <- Y %*% Vr` then `YVr_padded` |
| $\mathbf{C} = \mathbf{Y}'\mathbf{Z}$ | `t(YV) %*% YVr_padded` |
| SVD $\mathbf{C} = \mathbf{U}\boldsymbol{\Lambda}\mathbf{V}^{*\prime}$ | `svd_res <- svd(...)`, giving `svd_res$u`, `svd_res$v` |
| $\mathbf{A} = \mathbf{V}^{*}\mathbf{U}'$ | `A <- svd_res$v %*% t(svd_res$u)` |
| $\mathbf{V}^{new} = \mathbf{V}\mathbf{A}'$ | `Vrho <- V %*% t(A)` |
| $(\mathbf{V}^{new})^{-1}$ / prediction matrix | `tVrho <- solve(Vrho)` |

Two points where the code makes an implementation choice not spelled out verbatim in the thesis
prose, both intentional and both consistent with it:

1. **No dilation step.** §2.3 above; the thesis explicitly says $\rho$ "is not required," and the
   code never computes or applies it.
2. **`solve(Vrho)`, not a transpose.** The thesis states the prediction matrix as
   $(\mathbf{V}^{new})' = \mathbf{A}\mathbf{V}'$ — a transpose, not a general inverse. That
   identity holds only when $\mathbf{V}^{new}$ is orthogonal, which is true in the PCA case
   because $\mathbf{V}$ itself is Euclidean-orthonormal and $\mathbf{A}$ is orthogonal, so their
   product is too. For a **CVA** projection, $\mathbf{V}$ is not Euclidean-orthonormal — it is
   $\mathbf{W}$-orthonormal ($\mathbf{V}'\mathbf{W}\mathbf{V} = \mathbf{I}$ for the within-class
   covariance $\mathbf{W}$; see `mahalanobis_technical_note.md` §2.2) — so $\mathbf{V}^{new}$ is
   not orthogonal in the ordinary sense and a transpose would not equal its inverse. The code uses
   the general matrix inverse `solve()` so both projection types are handled correctly by the same
   function.

One further departure from the thesis's general $r$: the package's biplot is **always
two-dimensional** ($r=2$ fixed, unlike the thesis's illustrative $r=1$ example). Because
`YVr_padded` therefore has non-zero entries only in its first two columns, the SVD construction
guarantees the target's rotated information concentrates in columns 1–2 of `Vrho` *regardless of
which eigenvector pair was used to build* `Vr` — which is why the function can simply always
return columns 1–2 (`Vrho[, c(1L, 2L)]`) rather than tracking where the pair's original indices
ended up.

---

## 3. Part B — Selecting which plane to rotate: the eigenvector-pair search

Part A rotates *whatever* 2D plane you hand it (via `proj_pair`, the two eigenvector/loading
columns chosen as $\mathbf{V}_r$). With $p$ features there are $\binom{p}{2}$ possible planes to
choose from before rotating. The thesis chapter given for this note describes the rotation of a
single, already-decided plane (its worked examples fix a specific PCA/CVA reduction and rotate
that one); it does not narrate a search over multiple candidate planes as part of this chapter's
prose.

That search is nonetheless a genuine part of the original PhD method, just documented in code
rather than in this chapter: `mahalanobis_technical_note.md` records that the original research
script `scripts/1.3 Optimal Rotation.R` already computed a per-pair distance
`pc.distz[i] <- (v) %*% Wmat_inv %*% t(v)` across candidate pairs, using the within-class scatter
matrix `Wmat` built in a companion script. The package's `bl_find_local_cf()` is a direct
generalisation and packaging of that mechanism into a documented, parameterised search.

### 3.1 Pair enumeration

`R/local_cf.R:18-32`:

```r
.generate_pairs <- function(max_dim, max_pairs) {
  pairs <- vector("list", max_pairs)
  count <- 0L
  for (j in 2L:max_dim) {
    for (i in 1L:(j - 1L)) {
      count <- count + 1L
      pairs[[count]] <- c(i, j)
      if (count >= max_pairs) {
        return(do.call(rbind, pairs[seq_len(count)]))
      }
    }
  }
  if (count == 0L) return(matrix(integer(0), nrow = 0L, ncol = 2L))
  do.call(rbind, pairs[seq_len(count)])
}
```

This produces pairs in canonical order $(1,2), (1,3), (2,3), (1,4), (2,4), (3,4), \ldots$, capped
at `max_pairs` (default `10L` in `bl_find_local_cf()`). The package design notes
(`1 Foundation intro documents.txt` §10.7) fix this as a hard constraint: never test more pairs
than $\binom{p}{2}$ actually exist.

### 3.2 The per-pair pipeline

Inside `bl_find_local_cf()` (`R/local_cf.R:432-622`), each candidate pair goes through the same
five-stage pipeline:

1. **Rotate** — call `.bl_rotate(x_target_st, V, pair)` (Part A above) to get `Vr_rot`/`tVr_rot`
   for this pair.
2. **Project** — the target and all training observations are projected into this pair's rotated
   Z-space (`Z_target`, `Z_train_rot`).
3. **Grid and score** — an $m\times m$ grid ($m=200$ by default) is built over the training data's
   rotated range (+10% padding), back-projected to X-space, and scored through the model in
   chunks of 50,000 rows.
4. **Extract and filter the boundary contour** — `grDevices::contourLines()` at levels
   $cutoff \pm b_{margin}$, then a four-stage filter cascade on the resulting vertices
   (`R/local_cf.R:490-541`):
   - (a) keep only contour segments on the class *opposite* the target's predicted class;
   - (b) **actionability** — `set_filters` constraints (`"decrease"`, `"increase"`, `"fixed"`,
     `c(min,max)`), applied first because they are typically the most selective, letting the loop
     short-circuit early (`if (!any(keep)) next`);
   - (c) **feasibility** — `train_ranges` (per-feature min/max observed in training data), applied
     only to the survivors of (b);
   - (d) **model-consistency re-score** — every surviving vertex is back-projected and re-scored
     through the actual model, and dropped if the model disagrees with which side of `cutoff` the
     contour level implied. This is the most expensive stage, which is exactly why (b) and (c) run
     first, to minimise how many vertices reach it.

   Note what is deliberately *not* applied here: the Z-space training-hull polygon used elsewhere
   in the package. Rotation invalidates that polygon's coordinates (it was built for the
   unrotated plane), so only `train_ranges` and `set_filters` constrain the local search
   (`1 Foundation intro documents.txt` §10.4).
5. **Within-pair nearest vertex** — among the vertices that survive all four filters, the nearest
   one to the target is found by plain 2D Euclidean distance in the rotated Z-space
   (`.nearest_idx_block()`). This distance metric is deliberately **basis-local**: within one
   fixed rotated plane, ordinary Euclidean distance is the natural and correct way to find "the
   nearest point on this plane."

### 3.3 Why a different metric is needed *across* pairs

Comparing the within-pair "nearest vertex" distances directly across different pairs by their raw
2D Z-space Euclidean value is not mathematically sound, for a reason specific to this package's
loading matrices (`mahalanobis_technical_note.md` §2.2):

- Each pair captures a different amount of the total variance/separation in its particular 2D
  slice.
- The loading matrix $\mathbf{V}$ is **$\mathbf{W}$-orthonormal**
  ($\mathbf{V}'\mathbf{W}\mathbf{V} = \mathbf{I}$ for the within-class covariance $\mathbf{W}$),
  not Euclidean-orthonormal. A one-unit step in rotated Z-space therefore corresponds to a
  *different* displacement in the original feature space depending on which pair produced that
  Z-space.
- The rotation $\mathbf{A}$ from Part A is orthogonal in the full $p$-dimensional space, but the
  2D slice it produces (`Vrho[, c(1,2)]`) does not preserve Euclidean distances uniformly across
  different pair choices.

So a metric is needed that lives entirely in X-space and does not care which 2D rotation produced
a given candidate. The chosen metric is the **squared Mahalanobis distance** using the pooled
within-class covariance $\mathbf{W}$ — the same matrix $\mathbf{V}$ is orthonormal under, and the
same one CVA itself optimises against:

$$
d_M^2(v) = v^{\mathsf T}\mathbf{W}^{-1}v, \qquad v = B_x - x_{obs}
$$

`R/local_cf.R:563-570`:

```r
dist_mah <- if (!is.null(metric_inv)) {
  v <- as.numeric(B_x[1L, var_names]) - x_obs
  as.numeric(t(v) %*% metric_inv %*% v)
} else {
  NA_real_
}
```

`metric_inv` ($\mathbf{W}^{-1}$) is computed once in `bl_build_projection()` via
`chol2inv(chol(W))` (with a ridge-regularised fallback `chol2inv(chol(W + \lambda I))`,
$\lambda = 10^{-6}\cdot\text{mean}(\text{diag}(\mathbf{W}))$, when $\mathbf{W}$ is singular — see
`mahalanobis_technical_note.md` §3.2 for why Cholesky rather than `solve()`), then stored on
`bl_result$metric_inv` and reused unchanged for every pair and every target.

This is not a novel addition on top of the original method — it is the **same criterion the
original PhD script already used** (`pc.distz <- v %*% Wmat_inv %*% t(v)`). The package's
`distance = "euclidean"` option (plain 2D Z-distance, used for cross-pair comparison) exists only
to reproduce a since-superseded intermediate state of the port — during the migration described in
`mahalanobis_technical_note.md`, the package had briefly regressed to Euclidean because
`biplotEZ::CVA()` does not expose the `Wmat` its own computation depends on internally, so the
package could not access the correct metric until it was reconstructed independently. Mahalanobis
is the default and the methodologically intended behaviour, matching the PhD original.

### 3.4 Winner selection

`R/local_cf.R:584-621` — a running-minimum loop over all tested pairs:

```r
selector <- if (distance == "mahalanobis") dist_mah else dist_z
if (is.finite(selector) && selector < best_selector) {
  best_selector <- selector
  ...
  best_result <- list(
    B_z = B_z_local, B_x = B_x, B_pred = B_pred[[1L]],
    dist_z = dist_z, dist_mahalanobis = dist_mah,
    Z_target = Z_target, best_pair = pair,
    Vr_rot = Vr_rot, tVr_rot = tVr_rot,
    Z_train_rot = Z_train_rot, Zgrid = Zgrid, grid_prob = grid_prob,
    col_value = col_value, ct_local = ct_local,
    xseq = xseq, yseq = yseq, min_val = min_val, max_val = max_val
  )
}
```

Every pair's `dist_z` and `dist_mahalanobis` are retained (`all_distances`,
`all_distances_mahalanobis`) even for pairs that don't win, so a user can inspect how close the
alternatives were. The winning pair's rotated loading matrices, grid, contour, and boundary point
are what get plotted and reported as the local counterfactual (`bl_local_result$best_pair`,
`$B_x`, `$B_z`, etc.).

---

## 4. Known issues and limitations

### 4.1 The rotation plane is under-determined (the core caveat)

This is the most significant limitation, and it is subtle enough to be easy to miss entirely.

**Why it happens.** $\ddot{\mathbf{X}} = [-x;\,0;\,x]$ has **rank 1** — all three rows are
multiples of the single vector $x$. Writing $a = (-1, 0, 1)'$ and $u = \mathbf{V}'x$, the
Procrustes cross-product becomes

$$
\mathbf{C} = \mathbf{Y}'\mathbf{Z}_{padded} = 2\,u\,w', \qquad w = (u_{[\text{pair}_1]}, u_{[\text{pair}_2]}, 0, \ldots, 0)
$$

— also rank 1. Minimising the Procrustes objective over orthogonal matrices then reduces to a
**single vector constraint**, $\mathbf{A}'u \propto w$, on a $p\times p$ orthogonal matrix. For
$p=4$ this pins down only 1 of 3 rotational degrees of freedom; the solution set is a
3-dimensional manifold of rotations that all achieve the identical Procrustes optimum. The SVD of
the rank-1 $\mathbf{C}$ has $p-1$ **exactly degenerate** zero singular values, and the
corresponding singular vectors — hence $\mathbf{A}$, hence the plane's second direction — are
whatever a particular LAPACK bidiagonalisation implementation happens to emit for a degenerate
subspace.

**Measured consequences** (loan-default GAM model, $p=4$, target row 49, `proj_pair = (2,3)` — the
pair the search actually selects for this target):

- Singular values of $\mathbf{C}$: `1.6274e+00, 6.7123e-17, 0.0000e+00, 0.0000e+00` — confirming
  rank 1 of 4, a 3-dimensional null space.
- Perturbing the input basis $\mathbf{V}$ by amounts at or below ordinary floating-point noise
  produces wildly different planes (first principal angle is always 0° — the target's own
  direction is always shared and stable; the *second* angle is the arbitrary part):

  | Perturbation to V | max abs change | 2nd principal angle |
  |---|---|---|
  | re-orthonormalise only | 4.44e-16 | 56.79° |
  | random, scale 1e-16 | 4.44e-16 | 12.16° |
  | random, scale 1e-15 | 1.44e-15 | 0.00° |
  | random, scale 1e-12 | 1.11e-12 | 35.26° |
  | random, scale 1e-08 | 2.24e-08 | 11.49° |
  | `svd(X_std)$v` vs biplotEZ's own `Lmat` (numerically the same matrix, an incidental real-world difference) | **4.65e-15** | **65.91°** |

- Only **46%** of the resulting plane's second direction actually lies inside
  $\text{span}(\mathbf{V}_{[,2]}, \mathbf{V}_{[,3]})$ — the eigenplane the pair label `(2,3)`
  claims to represent.
- The target's own coordinates are essentially identical in both planes regardless
  (`(0.764192, 1.346277, ...)` to 6 decimal places) — which is exactly why the discrepancy is easy
  to miss: the one thing you'd naturally check (does the target land correctly?) looks fine either
  way.

**Practical consequences:**

1. **Reproducibility risk.** The returned counterfactual depends on a LAPACK tie-break over an
   exactly degenerate subspace. A different BLAS/LAPACK implementation, platform, or R version may
   legitimately return a different plane — and therefore a different counterfactual — from
   identical data and identical code. For a method meant to support auditable, published
   explanations, this is a material risk to disclose.
2. **`proj_pair` is a weaker label than it appears.** The pair search (Part B) enumerates six
   named eigenvector pairs, but each candidate plane is only loosely tied to the pair that names
   it (46% overlap, in the measured case). It is more accurate to describe the search as *sampling
   six arbitrary planes through the target* than as *enumerating six eigenplanes*.
3. **The counterfactual is a local, not global, optimum within its own plane.** It is the nearest
   boundary point within one arbitrary plane through the target — not necessarily the nearest
   point within the *named* pair's actual eigenplane, and not the nearest boundary point in the
   full feature space (see §4.4).
4. Downstream artefacts built on `bl_local_result` — Shapley attribution and sparse
   counterfactuals — inherit whichever plane happened to be drawn.

**What is *not* affected**, stated explicitly so the finding stays in proportion: the target's
position is exact and stable regardless of the tie-break; the returned counterfactual is still a
genuine boundary point (every candidate is re-scored through the actual model in filter stage (d),
§3.2), just not a canonical or perfectly reproducible one; the cross-pair Mahalanobis selector
(§3.3) is unaffected because it is computed in X-space, a frame the rotation cannot distort; Phases
1 and 2 are untouched — `bl_build_projection()`, `bl_build_grid()`, `bl_find_boundary()`,
`bl_robustness()` and `bl_surrogate()` never call `.bl_rotate()`; and `bl_set_scaling()` is
display-only and provably cannot reach `V`, the grid, or the plane.

An explicit, deterministic alternative was considered and **deliberately not implemented**:
replace the rank-deficient Procrustes construction with a minimal (Givens-style) rotation acting
only in the plane spanned by the target $x$ and its shadow $x_P = PP'x$ ($P$ = the named pair's
eigenvectors), which would make `proj_pair` mean exactly what it says and be fully reproducible. It
was not actioned because it would change every Phase 3 counterfactual the package has ever
produced, and the original PhD method's validated behaviour is the current, under-determined
construction — a change of this kind needs its own prototyping and measurement (reproducibility,
variance retained, and effect on Mahalanobis distances) before being adopted, not a unilateral
substitution.

**§6.2 works this alternative out** — the rotation in closed form, its four properties with proofs,
and the measurements this paragraph asks for. Measured headline: under the minimal rotation the
plane's first principal angle to the eigenplane its label names is **0 exactly** for every pair,
where the current construction is up to 55.62° off.

### 4.2 Within-pair vs. cross-pair metric mismatch

Within one rotated plane, the nearest boundary vertex is chosen by plain 2D Euclidean distance
(§3.2, stage 5). Across planes, the winner is chosen by X-space Mahalanobis distance (§3.3). These
are not necessarily wrong *together* — they operate at different scopes, and within a single fixed
plane Euclidean is the natural local metric — but it is an inconsistency worth naming for anyone
auditing the method end to end: the thing that decides "which point is nearest" is not the same
thing that decides "which plane is best." This is recorded as an open question, not yet resolved,
in the package's own documentation. **§6.3 recommends resolving it** by using squared Mahalanobis
for the within-plane step as well, which costs no additional model evaluations.

### 4.3 Global fidelity is traded for local accuracy

Directly from the thesis's own worked example (§2.5): rotating the biplot to fit one target
increases the total sum-of-squares residual for **every other point** on the biplot (0.28 → 4
against a total of 14, in the toy example), and improves one variable's axis predictivity strictly
at the expense of others. This is not a defect — it is the explicit, acknowledged mechanism of the
method — but it means a rotated (RC) biplot must never be read as a global summary of the model or
the data. It is only valid as a *local* picture, centred on the one target it was rotated for.

### 4.4 A rotated plane does not guarantee the globally nearest boundary

The thesis is explicit about the scope of what rotation fixes: it resolves the *classification
consistency* problem (does the model and the picture agree about this one point?), not the
*nearest-boundary-in-the-full-feature-space* problem. There is no guarantee that the nearest
decision boundary to the target, in the full $p$-dimensional space, lies anywhere near the rotated
2D plane. If multiple decision-boundary surfaces exist in the full space, this method can only
report on the one visible within whichever plane it happened to draw.

**This limitation is retained.** §6 improves the choice of plane *within* the biplot framework and
deliberately does not address the full-space problem — see §6.6.

### 4.5 CVA-specific consequence

Specific to §2.6: rotating a CVA transformation pulls the group mean points off the canonical
plane, so their exact $(k-1)$-dimensional predictions — a property CVA is otherwise built to
guarantee — are lost after rotation. This cost has no analogue in PCA-based rotation.

### 4.6 Adjacent mechanisms not to confuse with this one

- The diagonal-only Mahalanobis approximation used elsewhere in the package (Phase 2's
  `bl_robustness()`/`plot.bl_boundary()`) is a related but separate simplification — full
  Shapley-based decomposition of the Mahalanobis distance would be the theoretically exact
  per-feature attribution but costs $O(n\cdot 2^p)$ and is deliberately not implemented at that
  scale (it *is* implemented, tractably, for single-observation Phase 3 explanations via
  `bl_shapley()`, at $O(2^p)$).
- `rotate_deg` (an argument on the local-biplot plotting function) is a purely **cosmetic**
  post-hoc 2D rotation applied on top of the already-SVD-rotated plane, for display convenience
  only. It has nothing to do with the Procrustes mechanism described here and should not be
  confused with it.

---

## 5. Measured findings from the loan-data worked examples (2026-09)

Everything in this section is a **measurement**, not theory, taken from the scripts in this folder,
built against `boundarylogic` 0.1.0 under R 4.6.0:

| script | setting | role |
|---|---|---|
| `17_1a_toy-2d-to-1d-rotation.R` | 5 points, 2 variables → **1D** biplot (the thesis's toy) | the rotation at minimum scale |
| `17_1_loan-gam-pca-3d-3var-rotation.R` | 221 loans, **3** features → 3 PCs, GAM | the rotation with **no** dimension loss |
| `17_2_loan-gam-pca-4var-slices.R` | 221 loans, **4** features → 3 PCs, GAM | rotation **plus** truncation; two planes on one surface |
| `17_3_minimal-rotation-verification.R` | both settings above | the acceptance harness for §6's recommendation |

The 3-variable script is the workhorse: because a 3×3 PCA basis is a lossless rotation, anything
that changes between its views is attributable to the rotation alone. Unless stated otherwise the
numbers below are from it, with target row 6 and the planes taken from the package's own
`Vr_rot`.

### 5.1 The dimensional ladder — how much freedom the rotation actually has

§4.1 establishes that the Procrustes matrix is rank 1 and the plane therefore under-determined.
Counting the degrees of freedom makes the severity explicit and shows where it starts:

| setting | subspace DOF | minus "contains the target" | free |
|---|---|---|---|
| $p=2,\ r=1$ (the thesis's toy) | 1 | −1 | **0 — unique** |
| $p=3,\ r=2$ | 2 | −1 | 1 |
| $p=4,\ r=2$ | $2(4-2)=4$ | −1 | 3 (matches §4.1's "3-dimensional manifold") |

**Correction — see §6.1.** The $p=4$ row counts free *rotations*, not free *planes*. Containing a
fixed line is two conditions on $\operatorname{Gr}(2,4)$, not one, so the family of admissible
planes has dimension **2** — measured numerically in `17_2` §8a. The two counts differ by the
in-plane spin, the same free rotation §5.7 records as an entry-wise difference of 1.61 at a plane
difference of 2.2e-16. The $p=2$ and $p=3$ rows are correct as written.

**The toy example in §2.5 is the only case in the family where the rotation is determined.** That
is exactly why it teaches well — and exactly why the intuition it builds is wrong from $p=3$
onwards. Any exposition that introduces the method through the 2D→1D picture should say so.

Confirmed in the toy: rotating onto **PC1** and rotating onto **PC2** produce the *identical*
rotation matrix $\mathbf{A}$, the identical retained axis, and identical adequacy / RSS / sample
predictivity. Only $\rho$ differs (1.134058 vs 2.120244), and $\rho$ is deliberately discarded
(§2.3). The reason is geometric: with $r=1$ the retained axis is forced to be the target's own
direction — verified $|\cos| = 1.000000000000$ for all five points — and in 2D only one line
through the origin contains the target. There is nothing to choose.

**Caution on a §4.1 sentence.** "The first principal angle is always 0° — the target's own
direction is always shared and stable" is a $p=4$ statement and does **not** transfer to $p=3$.
In 3-space *any* two 2-planes intersect in a line, so a first principal angle of 0° is forced by
dimension and carries no information; three random plane pairs give `0.00, 83.34` / `0.00, 47.79`
/ `0.00, 87.81`. Nor is the shared line the target's direction there: for row 6 only 0.4435 of
that direction lies in the PC1–PC2 plane (which is exactly its unrotated sample predictivity).

### 5.2 The pair label is weak — now quantified

For row 6, rotating onto each named pair:

| named pair | $\theta$ | 2nd principal angle vs `span(V[,pair])` | 2nd direction inside the named eigenplane |
|---|---|---|---|
| (1,2) | 48.25° | 86.26° | **56.08%** |
| (1,3) | 3.03° | 3.32° | 99.94% |
| (2,3) | 41.59° | 86.97° | **44.35%** |

Worse, per column: **for all three pairs alike, rotated column 1 comes out closest to PC3 and
column 2 closest to PC1.** Three different labels, the same resulting axes — a blunter
demonstration of §2.7 than the algebra alone.

**And $\theta$ does not predict the overlap.** Across all 663 (row, pair) combinations on this
data: Pearson $r = +0.33$, Spearman $+0.14$, and among the 74 cases with $\theta < 5°$ the overlap
spans the **full range 0% to 100%**. Row 6's tidy pairing of $\theta = 3.03°$ with 99.94% is a
coincidence; row 137 has $\theta = 1.49°$ with **0.1%**. So there is no cheap screen for §4.1's
problem — the overlap has to be computed, and it cannot be inferred from how far the plane had to
move.

§6.2 removes the problem rather than screening for it: under the minimal rotation the plane shares a
direction with its named eigenplane exactly, and its only departure is $\theta$ itself.

### 5.3 Reproducibility is worse than §4.1 records

§4.1 anticipates that a different BLAS/LAPACK, platform or R version may return a different plane.
Measured here: **`best_pair` can vary between identical runs in a single session.** Row 102 flipped
between (1,2) and (1,3) across five consecutive calls with the same objects, its (1,3) distance
moving between 0.232776 and 0.228868. Rows 161 and 140 kept a stable winner with wobbling
distances; rows 6 and 104 were bit-identical every run.

*Practical consequence for the posts:* if a figure or a quoted counterfactual must be reproducible,
choose a target whose `best_pair` is bit-stable and say that you checked, by repeating the call.

§6.2 removes the cause: with one determinate plane per label there is nothing left to tie-break, and
the returned basis was measured bit-identical across five identical calls.

**And this is not the whole of it — see §5.9.** The advice above (repeat the call and check) turns
out to be necessary but not sufficient, because the reported distance also depends on *what was
called before it* in the session, and repeating a call is one of the things that changes it.

### 5.4 The 2D biplot is a **slice** of the decision surface, not a view of it

This is the most load-bearing interpretive point found, and the one most likely to mislead a reader
of the published figures.

- A **3D** figure of a 3-variable model shows the whole $p\approx0.5$ surface. Because the basis is
  a lossless rotation, every rotated view holds the same surface, the same cloud and the same
  points, re-expressed — rotating only moves the camera.
- A **2D biplot** shows where the model crosses 0.5 **within the drawn plane**, with the off-plane
  coordinate pinned at 0. Two planes through the same target cut one surface along two genuinely
  different curves.

So `plot(bl_local)` for two different pairs can look completely different while the 3D picture
looks unchanged. Neither plot is wrong; they are different cross-sections.

Measured for row 6 (pinned (1,2) versus the selected (1,3)):

| | points | model $p$ along the curve | $\max|p-0.5|$ |
|---|---|---|---|
| pinned (1,2) | 700 | [0.49816, 0.50202] | 2.02e-03 |
| selected (1,3) | 498 | [0.49812, 0.50198] | 1.98e-03 |

Both are genuine boundary curves, yet nearest-neighbour distance between them is median **2.30**,
max **3.85** against a data span of **6.71** — they separate by **57% of the data's extent**. The
two planes are **90.00°** apart here.

**They nonetheless share exactly one point, and this is forced.** Both planes contain the origin (a
rotation has no translation) and the target; two distinct planes meet in exactly one line;
therefore their intersection is $\operatorname{span}(x_t)$, the target's own ray. Verified
$|\cos(n_B \times n_C,\ x_t)| = 1.000000000000$. The surface crosses that ray once for row 6, at
$t=0.41223$, and **both curves pass within 0.0036 and 0.0010 of it** against a contour grid step of
0.0363.

*Why it is not cosmetic:* the counterfactuals differ — $B_x = (+0.759, +0.062, +0.165)$ versus
$(+0.861, +0.381, +0.084)$, **0.345 apart**, with $\sqrt{\text{Mahalanobis}}$ 1.2527 versus 1.1455.
Different slice → different nearest boundary point → different advice. This is the concrete case
for why the pair search exists, and why §4.1's under-determination has consequences rather than
being a curiosity.

### 5.5 A plane's tilt needs two angles, and tilt is not a selection criterion

Relative to the unrotated PC1–PC2 plane, a rotated plane's orientation needs **two** numbers, not
one: writing it as $s_3 = s_1 s_{1\text{ax}} + s_2 s_{2\text{ax}}$, the component tilts are
$\arctan s_1$ along PC1 and $\arctan s_2$ along PC2, and the dihedral angle is their resultant,
$\tan^2\theta = s_1^2 + s_2^2$ (verified to 6 dp). The single dihedral hides the direction of lean:
for row 6 the pinned (1,2) plane is $(-48.26°, +2.04°)$ — almost a pure PC1 tilt — while the
selected (1,3) plane is $(+65.85°, +88.64°)$, leaning substantially both ways.

"Contains the target" is again one equation, $s_1 a_1 + s_2 a_2 = a_3$, so the admissible tilts
form a **line** in the $(s_1,s_2)$ plane — a 1-parameter family. For row 6 that family admits
dihedrals from **48.25° to 89.05°**; the pinned (1,2) plane came out at 48.28° (the shallow end)
and the selected (1,3) plane at 88.64° (the steep end).

**Terminology discipline, worth carrying into the posts:** "optimal" and "best" should refer only
to what `bl_find_local_cf()` measures — the smallest (squared) Mahalanobis distance from target to
counterfactual. How steeply a plane is tilted is *not* a criterion the package optimises, and here
the two run in opposite directions. A shallower plane is easier to read; it is not a better answer.

### 5.6 What the rotation costs, per variable

`adequacy` and `axis predictivity` index the **rows** of $\mathbf{V}$ — they are per *variable*, not
per axis. §4.3 records that global fidelity is traded for local accuracy; the per-variable view
shows the trade is a **redistribution**, and that *which* variable is sacrificed depends on where
the rotated plane lands. Row 6, $\mathbf{J} = \operatorname{diag}(1,1,0)$:

| | adequacy (income / score / rate) | RSS (of 660) | target sample predictivity | mean of the other 220 |
|---|---|---|---|---|
| unrotated | 0.5023 / 0.4978 / 0.9999 | **168.656** | 0.4435 | 0.7161 |
| rotated → (1,2) | 0.9988 / **0.0442** / 0.9570 | 229.465 | **1.000000** | 0.6710 |
| rotated → (1,3) | 0.9671 / 0.9609 / **0.0720** | 213.624 | **1.000000** | 0.6758 |

The two rotations sacrifice **different** variables — `credit_score` in one, `loan_int_rate` in the
other. The target reaches sample predictivity exactly 1 in both (the §2.4 identity
$\hat{x}_i = x_i$, measured at 2.2e-15), and the unrotated RSS equals the dropped component's
variance $d_3^2$ exactly, which is a useful arithmetic check on any reimplementation.

### 5.7 Practical notes for reproducing the package's objects

- **Use the package's own objects rather than re-deriving them.** `bl_local$Vr_rot` is the rotated
  basis the local biplot actually uses; `bl_local$ct_local` and `bl_result$biplot_grid$ct` are its
  boundary contours, already in the plane's 2D coordinates, and plot verbatim.
- **Rotating from `svd(X)$v` instead of `bl_result$V` lands on a different plane.** The two matrices
  agree to ~1.4e-15, and that is enough — the completion is degenerate (§4.1). Rotating from
  `bl_result$V` reproduces the package's plane to 2.2e-16.
- Comparing planes: compare **projection matrices** $\mathbf{P} = \mathbf{V}_r\mathbf{V}_r'$, never
  entry-wise. Even with the plane matching exactly, the in-plane basis can differ by a large
  rotation (measured 1.61 entry-wise at 2.2e-16 plane difference).
- **`max_pairs = 1L` pins the search to pair (1,2)**, because `.generate_pairs()` emits
  colexicographically. This is the only supported way to force a specific plane, and it makes a
  pinned-versus-selected comparison possible.
- **Any replication of the pair search must include the `train_ranges` filter** (§3.2 stage (c)).
  It is not a formality: measured removal of boundary candidates was 31% / 18% / 56% by pair for
  row 6 and 56% / 51% / 49% for row 199. Omitting it changes which candidate is nearest and hence
  which pair wins — a replication without it agreed with the package on only about 3 rows in 13.
- `all_distances_mahalanobis` stores the **squared** distance (§3.3 states the metric; the field is
  easy to misread). The package's own summary prints its square root.

### 5.8 Caveat on the §2.5 worked-example numbers

`17_1a_toy-2d-to-1d-rotation.R` reproduces the thesis's rotation angle (**28.14°** — measured
28.1411°) and dilation ($\rho = 1.134 \approx 1/1.13$) exactly, but **not** its predictivity or
residual figures: the thesis reports sample predictivity 0.72 → 1 and RSS 0.28 → 4 of 14, where its
five points give **0.7776 → 1** and **4.307 → 4.927 of 11.4**. The data behind §2.5 is therefore
probably not the data in that script, and the two sets of numbers should not be quoted together.

(On the residual figure: RSS here is the **trace** of $\mathbf{X}'\mathbf{X}$, not the sum of its
entries. The distinction matters — summing every entry is not a sum of squares and makes the
rotation appear to *improve* global fit, which is the opposite of the trade-off the method actually
makes. With the trace, the unrotated RSS equals $d_2^2$ exactly and the rotated RSS is correctly
larger, as §4.3 requires. `17_1a` and `17_1` both compute it that way.)

### 5.9 The full search depends on what was called before it (2026-09-10)

Found while building `17_3` and not previously recorded. It is a stronger reproducibility problem
than §5.3's, because it identifies a **trigger** rather than an intermittent flip.

For one fixed target ($p=4$, row 49), `bl_find_local_cf()`'s reported distance depends on the
session's call history. Measured, in a single R session, on the same `bl_result` and the same
target object (`17_3` §8b):

| call | | `dist_mahalanobis` |
|---|---|---|
| 1 | full search (first call in the session) | **1.147879** |
| 2 | `max_pairs = 1L` | 1.208670 |
| 3 | full search | **1.153789** |
| 4 | full search | 1.153789 |
| 5 | `max_pairs = 1L` | 1.202898 |
| 6 | full search | 1.153789 |
| 7 | full search | 1.147879 |

The selected pair is **(2,3) throughout** — this is not §5.3's pair flip. What varies is the
distance, across two values for the full search (1.147879 / 1.153789) and two for the pinned search
(1.208670 / 1.202898). The exact sequence is not itself stable run to run; the *set* of values is.

**What it is not.** 200 unrelated `svd()` calls interleaved between two full searches leave the
result unchanged, and so does building a second `bl_result` in between. So this is not ambient
floating-point state and not §4.1's input-perturbation sensitivity.

**Not only earlier `bl_find_local_cf()` calls (added 2026-10-17, measured 2026-10-06).** The first
version of this section said the effect was specific to earlier `bl_find_local_cf()` calls on the
same objects. That is too narrow. In the *Understanding Boundaries* project's slide script for a
4-variable GAM (row 49, `k = 8`, `max_pairs = 1L`, `distance = "euclidean"`), three Phase 1–2
calls on the same `bl_result` before the first `bl_find_local_cf()` also moved the answer. The
calls were `plot(bl_result)`, `bl_find_boundary()` and `bl_surrogate()`:

| order | distance | `loan_amnt` change at the counterfactual |
|---|---|---|
| `bl_find_local_cf()` first on a fresh `bl_result` | 0.6567 | +2,101 |
| after `plot()`, `bl_find_boundary()`, `bl_surrogate()` | 0.6567 | +2,146 |

Two things differ from the table above. The **distance is unchanged** to four decimals and the
**point moves**, so this looks like a choice between near-tied boundary points rather than a
different plane. And **the plane is pinned** (`max_pairs = 1L`), so pair selection cannot be
involved. The three calls were not separated, so which of them is the trigger is not known. §4.1
notes that Phases 1 and 2 never call `.bl_rotate()`, and that still holds, so any shared state
lies outside the rotation itself. The question for the package is therefore wider: does *any*
exported function modify the `bl_result` it is given, or cache something on it, that
`bl_find_local_cf()` later reads?

**Mechanism not established.** The installed package ships as a lazy-load database, so the source
was not available to inspect from this project. This is a measurement and a bug report, not a
diagnosis; it is recorded here so the package project can look at whether `bl_find_local_cf()`
mutates the `bl_result` or `bl_target` it is passed, or caches anything across calls. Since
2026-10-06, the same question applies to `plot.bl_result()`, `bl_find_boundary()` and
`bl_surrogate()`.

**Practical consequence.** §5.3's advice — repeat the call and check stability — is necessary but
not sufficient, because *repeating* a call is one of the things that changes the answer. Take a
quoted counterfactual from a **first** call on a freshly built `bl_result`, with no other function
called on that object beforehand, not even `plot()`, and record that this is what you did. For a published figure, rebuild the result object rather than reusing one that has
already been searched.

---

## 6. Recommendation: a deterministic rotation and a reproducible plane search

Sections 4 and 5 diagnose two faults and stop there. §4.1 sketches the fix in one paragraph —
replace the rank-deficient Procrustes construction with "a minimal (Givens-style) rotation acting
only in the plane spanned by the target $x$ and its shadow $x_P = PP'x$" — and declines it pending
"prototyping and measurement (reproducibility, variance retained, and effect on Mahalanobis
distances)."

This section works that proposal out: the construction in closed form, a convention for the
in-plane frame, a replacement for the six-label sample, and the measurements §4.1 asks for. Every
number below is reproduced by **`17_3_minimal-rotation-verification.R`** in this folder, on the loan
data at $p=3$ target row 6 and $p=4$ target row 49 — the same objects scripts `17_1` and `17_2` use.
That script asserts the geometry and reports the timings, so a drift fails loudly rather than
quietly changing a quoted figure. Run it before relying on anything here.

**Scope.** Everything here stays inside the 2D biplot framework. The rotation's job is to put the
target *on* the drawn plane so the local reading is exact; the search's job is to find the nearest
boundary *available within that framework*. Searching the full $p$-dimensional space for a nearest
boundary point is **out of scope** — §4.4's limitation is retained, and §6.6 restates it.

### 6.1 What exactness actually requires, and how much freedom that leaves

The biplot's reconstruction is the orthogonal projector onto the plane's span, so

$$
\hat{x} = \mathbf{V}_r(\mathbf{V}_r'\mathbf{V}_r)^{-1}\mathbf{V}_r'x = x
\quad\Longleftrightarrow\quad
x \in \operatorname{span}(\mathbf{V}_r).
$$

That is a condition on the **plane**, not on the frame, and it is independent of how
$\mathbf{V}_r$ is scaled or oriented — so it reads identically for PCA and CVA. Nothing beyond
"the plane contains the target" is needed for the exact local interpretation. Consequently the
admissible family is

$$
\mathcal{A} \;=\; \{\,2\text{-planes } P \subset \mathbb{R}^p \;:\; x_t \in P \,\},
$$

and since every such plane is $\operatorname{span}(x_t, d)$ for a direction $d$ in the
$(p-1)$-dimensional space orthogonal to $x_t$, taken modulo sign,
$\mathcal{A} \cong \mathbb{P}^{p-2}$ — dimension $p-2$:

| $p$ | $\dim\mathcal{A}$ | shape |
|---|---|---|
| 2 | 0 | a single plane — determined, which is why the thesis's toy (§2.5) teaches so cleanly |
| 3 | 1 | a pencil of planes about the axis $\operatorname{span}(x_t)$ — one angle |
| 4 | 2 | a two-parameter family |

**This corrects §5.1's table.** That table reports "free = 3" at $p=4$, citing §4.1's
"3-dimensional manifold of rotations". The 3 is right for **rotations** — those satisfying
$\mathbf{A}'u \propto w$ form a coset of $\operatorname{Stab}(u) \cong SO(p-1)$, and
$\dim SO(4) - \dim SO(3) = 6 - 3 = 3$ — and wrong for **planes**: containing a fixed line is *two*
conditions on $\operatorname{Gr}(2,4)$, not one, leaving $4 - 2 = 2$. Script `17_2` §8a measures
the plane family's tangent dimension numerically and gets **2**. The two counts differ by the
in-plane spin, the same free rotation §5.7 measures as an entry-wise difference of 1.61 at a plane
difference of 2.2e-16. At $p=3$ the table's "$-1$" happens to be right (planes through a line in
3-space form a 1-parameter pencil), which is why the error never surfaced in `17_1`.

So the plane choice is entirely free with respect to exactness. The two recommendations below are
about how to spend that freedom: **deterministically** (§6.2) and **well** (§6.3).

### 6.2 Part 1 — a minimal Givens rotation, in closed form

Given a named eigenplane $P_0 = \operatorname{span}(\mathbf{V}_{[,a]}, \mathbf{V}_{[,b]})$, write
$\mathbf{V}_{r0} = \mathbf{V}_{[,(a,b)]}$ and the target's shadow
$x_P = \mathbf{V}_{r0}\mathbf{V}_{r0}'x_t$. For $0 < \lVert x_P\rVert < \lVert x_t\rVert$:

$$
u = \frac{x_P}{\lVert x_P\rVert}, \qquad
r = x_t - x_P, \qquad
w = \frac{r}{\lVert r\rVert}, \qquad
\cos\theta = \frac{\lVert x_P\rVert}{\lVert x_t\rVert}, \qquad
\sin\theta = \frac{\lVert r\rVert}{\lVert x_t\rVert}
$$

$$
\boxed{\;\mathbf{R} = \mathbf{I} + (\cos\theta - 1)(uu' + ww') + \sin\theta\,(wu' - uw'),
\qquad \mathbf{V}_r = \mathbf{R}\,\mathbf{V}_{r0}\;}
$$

$\{u,w\}$ is orthonormal ($r \perp P_0$ and $u \in P_0$), and $\mathbf{R}$ is a rotation by
$\theta$ inside $\operatorname{span}(u,w)$ and the identity on its complement. Note
$\lVert x_P\rVert^2 + \lVert r\rVert^2 = \lVert x_t\rVert^2$, so $\cos^2\theta + \sin^2\theta = 1$
as written.

**Properties.**

| property | proof |
|---|---|
| **Exactness**: $x_t \in \operatorname{span}(\mathbf{V}_r)$ | $\mathbf{R}u = \cos\theta\,u + \sin\theta\,w = x_t/\lVert x_t\rVert$, so $\mathbf{R}'x_t = \lVert x_t\rVert u \in P_0$, i.e. $x_t \in \mathbf{R}P_0$ |
| **Deterministic**: no SVD, no degenerate subspace, no tie-break | every quantity is a closed-form expression in $x_t$ and $\mathbf{V}_{[,(a,b)]}$ |
| **Minimal**: the smallest rotation achieving exactness | $\theta$ *is* the angle between $x_t$ and $P_0$; any $\mathbf{S}$ with $x_t \in \mathbf{S}P_0$ needs $\mathbf{S}'x_t \in P_0$ with $\lVert \mathbf{S}'x_t\rVert = \lVert x_t\rVert$, so it must move $x_t$'s direction by at least $\theta$ |
| **`proj_pair` means what it says**: principal angles $(\mathbf{R}P_0, P_0) = (0,\theta)$ | the direction of $P_0$ orthogonal to $u$ lies in $\operatorname{span}(u,w)^\perp$, which $\mathbf{R}$ fixes pointwise — so the two planes share it |
| **$O(p^2)$** — cheaper than the current $O(p^3)$ SVD | matrix assembly only |

**Measured, $p=4$, target row 49.** The last property is the sharp one, because it is exactly what
the current construction gets wrong. `ang1` is the first principal angle against the eigenplane the
pair *names*:

| pair | $\theta$ | `ang1` minimal | `ang2` minimal | exactness residual | `ang1` current |
|---|---|---|---|---|---|
| (1,2) | 38.4522° | **0.00e+00** | 38.4522° | 1.40e-15 | 1.1111° |
| (1,3) | 34.9160° | **0.00e+00** | 34.9160° | 1.99e-15 | 5.1471° |
| (2,3) | 70.1512° | **0.00e+00** | 70.1512° | 1.08e-15 | 55.6173° |
| (1,4) | 19.8488° | **0.00e+00** | 19.8488° | 4.19e-15 | 17.2008° |
| (2,4) | 55.0840° | **0.00e+00** | 55.0840° | 2.69e-15 | 5.4437° |
| (3,4) | 51.5478° | **0.00e+00** | 51.5478° | 2.34e-15 | 51.3428° |

Under the minimal rotation the plane shares a direction with its named eigenplane exactly, and its
only departure is the $\theta$ it was forced to move. Under the current construction it shares
nothing — up to 55.62° off for pair (2,3) — which is the mechanism behind §5.2's 21–70% overlap
figures and behind three different labels yielding the same nearest axes. **Determinism**: the
returned basis was bit-identical across five identical calls.

**Edge cases**, each needing a documented convention:

- $\lVert x_P\rVert = 0$ (target orthogonal to $P_0$; $\theta = 90°$): $u$ is undefined. Take
  $u = \mathbf{V}_{[,a]}$, so the first named axis is the one rotated onto the target.
- $\lVert x_P\rVert = \lVert x_t\rVert$ (target already in $P_0$): $\theta = 0$,
  $\mathbf{R} = \mathbf{I}$, no rotation. Guard the $0/0$ in $w$.
- $x_t = 0$ (target at the centroid): every plane contains it, no rotation is defined or needed —
  return the unrotated eigenplane and warn.

**The in-plane frame, fixed by convention.** The plane does not determine an ordered orthonormal
frame, and the picture's appearance depends on which one is chosen (§5.7). Nothing *measurable*
does, so it is safe to fix. Recommended:
$\mathbf{V}_r = (\mathbf{R}\mathbf{V}_{[,a]},\ \mathbf{R}\mathbf{V}_{[,b]})$ — the named axes
carried through the rotation. RC1 and RC2 then have a defensible reading for the first time
("PC $a$ and PC $b$, tilted by $\theta$"), and §2.7's "columns 1–2 of $\mathbf{V}_\rho$ regardless
of which pair was used" convention, with the axis-identity confusion it creates, is no longer
needed. A documented alternative is the target-aligned frame
$(x_t/\lVert x_t\rVert,\ \cdot)$, which puts every target on the horizontal axis and makes figures
comparable across targets, at the cost of loading arrows that move between targets.

**PCA versus CVA.** §3.3 describes the loading matrix as $\mathbf{W}$-orthonormal
($\mathbf{V}'\mathbf{W}\mathbf{V} = \mathbf{I}$). That is a CVA property: for `method = "PCA"`,
`17_1` and `17_2` both assert $\mathbf{V}'\mathbf{V} = \mathbf{I}$ to 1e-10 and pass. The
construction above is written for a Euclidean-orthonormal $\mathbf{V}$. Under CVA, carry it out in
$\mathbf{W}$-whitened coordinates ($\mathbf{W}^{-1} = \mathbf{L}\mathbf{L}'$, already available as
`bl_result$metric_inv`) and map back; the exactness condition of §6.1 is basis-independent, so it
survives the change of metric. §4.5's CVA-specific cost is unaffected either way.

**The thesis's worked example is preserved.** At $p=2,\ r=1$ the admissible family is a single
line, so the minimal rotation and Procrustes must agree. Verified on `17_1a`'s five-point
data (target row 5): Procrustes **28.1411°**, minimal rotation **28.1411°**, agreeing to 6.75e-14
and reproducing the thesis's 28.14° (§2.5). The change therefore does not invalidate the
theoretical exposition — it only removes freedom that exists from $p=3$ upward.

### 6.3 Part 2 — searching the admissible family, not six labels

§6.2 alone makes reproducibility a non-issue: each of the $\binom{p}{2}$ labels now maps to one
determinate plane, so §5.3's flip between (1,2) and (1,3) across identical calls cannot occur —
there is nothing left to tie-break. But the labels still sample only $\binom{p}{2}$ points of a
$(p-2)$-dimensional continuum — three at $p=3$, six at $p=4$ — so "best of the labels" is not best
of anything well defined. Since $\mathcal{A}$ has an explicit parameterisation, search it.

**Parameterisation.** Let $q = x_t/\lVert x_t\rVert$ and let
$\mathbf{U} \in \mathbb{R}^{p\times(p-1)}$ be a deterministic orthonormal basis of $q^\perp$ — the
Householder reflector carrying $q \to e_1$, minus its first column. Then every admissible plane is

$$
P(d) = \operatorname{span}(q,\ \mathbf{U}d), \qquad d \in \mathbb{R}^{p-1},\ \lVert d\rVert = 1,\ d \equiv -d.
$$

| $p$ | search space | recommended search |
|---|---|---|
| 3 | one angle $\varphi \in [0,\pi)$ | **exhaustive** 1° sweep, then golden-section refinement |
| 4 | hemisphere of $S^2$ | deterministic Fibonacci lattice, then local refinement from the best few |
| $\ge 5$ | $\mathbb{P}^{p-2}$ | the $\binom{p}{2}$ minimal-rotation planes as seeds plus local refinement; **document as a local optimum, not exhaustive** |

No RNG anywhere: every candidate set is a deterministic function of $p$ and the resolution.

**Objective, in two stages.** Stage 1 ranks planes cheaply; stage 2 runs the existing §3.2 pipeline
($m\times m$ grid → `contourLines` → constraints → nearest vertex) on the best few, so the reported
boundary and every plotting object keeps the package's established definition.

Stage 1 is a **ray scan**: in the candidate plane, cast $K$ rays from the target at equal angles
and bisect each to its **first** cutoff crossing, taking the smallest squared Mahalanobis distance
over rays. This is the right objective for a nearest-point search — the nearest in-plane boundary
point is reached by some ray from $x_t$, and along that ray it must be the *first* crossing, since
an earlier crossing would be nearer.

**Measured, $p=3$, target row 6, exhaustive 1° sweep of all 180 planes.**

| | squared $d_M$ | $d_M$ |
|---|---|---|
| sweep minimum (at $\varphi = 134°$) | **1.237877** | **1.1126** |
| sweep median | 1.490022 | 1.2207 |
| sweep maximum | 1.579964 | 1.2570 |
| minimal rotation, label (1,2) | 1.580762 | 1.2573 |
| minimal rotation, label (1,3) | 1.293857 | 1.1375 |
| minimal rotation, label (2,3) | 1.568419 | 1.2524 |
| **the package, as it stands** | 1.312277 | 1.1456 |

Three things this shows.

1. **The plane choice matters, and the search wins.** The sweep improves on the package by **5.7%
   in squared distance (2.9% in distance)**, and monotonicity holds as it must — every plane the
   current method draws lies in $\mathcal{A}$, so the sweep's optimum cannot be worse.
2. **Pinning to (1,2) lands at the worst plane in the family for this target.** Its 1.580762 is
   within 0.05% of the sweep's *maximum* of 1.579964. A pinned pair is not a neutral default.
3. **The curve has two local minima.** Near-ties are real, so they should be reported rather than
   silently resolved.

At $p=4$, target row 49, over a 200-point Fibonacci lattice: best label (1,2) **1.147994**, package
**1.147879**, family best **1.097519** — a **4.4%** improvement in squared distance.

**One caveat on those comparisons.** The sweep's figures come from the stage-1 ray scan and the
package's from its own $200\times200$ grid, so the two are measured with different discretisations
and the percentages are indicative rather than exact. The ray scan's 10° angular resolution can
only *miss* the nearest point, so it should if anything understate the improvement — but a
like-for-like number requires running stage 2 (the package's own pipeline) on the sweep's winning
plane, which is the gap recorded against test 4 in §6.5.

**Cost — measured, and it is a cost, not a saving.** A
plane-at-a-time sweep is **6.3×** the package's time (5.75 s vs 0.91 s at $p=3$), even though the
prediction counts are comparable (129,600 vs 120,000). The cost is per-call overhead, not
arithmetic: the package's grid does 40,000 rows in one vectorised call per pair, while a per-plane
ray scan does 36 rows per call. Batching the scan **across planes** — all $180 \times 36 = 6{,}480$
rows in one call per march/bisect step, 20 calls for the entire sweep — brings it to **1.65×**
(0.81 s vs 0.49 s). So the honest statement is: an exhaustive sweep at $p=3$ costs roughly
**1.6× the current search when batched, and about 6× if implemented naïvely.** That is affordable
for a per-target explanation, and the sweep is embarrassingly parallel besides, but it is a real
cost and should not be advertised as a saving.

**Two further fixes, both nearly free.**

1. **Feasibility as a constraint, not a post-hoc filter.** §5.7 measured `train_ranges` removing
   18–56% of boundary candidates *after* contouring. Filtering afterwards means the reported point
   is "nearest among survivors" and no better feasible point was ever sought — a bias, not just a
   loss of candidates. Restrict the grid and contour extraction to the feasible region instead.
   The same argument applies to `set_filters` actionability, which becomes a box or half-space
   constraint on the search rather than filter stage (b).
2. **One metric throughout**, resolving §4.2. Use squared Mahalanobis for the within-plane
   nearest-vertex step as well as across planes, so "nearest" means one thing end to end. The
   vertices are already back-projected to X-space for filter stage (d), so this costs no additional
   model evaluations.

**Reporting.** Replace `all_distances_mahalanobis` — six arbitrary planes — with the objective over
the searched family, plus the runner-up local minima, so near-ties are visible instead of silently
resolved. At $p=3$ that is a distance-versus-$\varphi$ **curve**, a genuinely informative
per-target diagnostic: for row 6 it shows the reported distance moving by a factor of 1.13 across
the family. Keep $\theta$ and the nearest named pair in the output for continuity with `best_pair`.

### 6.4 Adoption path

§4.1's own bar applies: this changes every Phase 3 counterfactual the package has produced, so it
should not be substituted unilaterally. Land it as a new option — `rotation = c("procrustes",
"minimal")` and `plane_search = c("pairs", "family")` — with the defaults unchanged, and switch the
defaults only once §6.5's tests have been run and reported across the full training set rather than
the two targets measured here.

### 6.5 Acceptance tests

These are the measurements §4.1 asks for. All of them are run by
`17_3_minimal-rotation-verification.R`; the "runs in" column names the section of that script, so a
claim can be traced to the code that produced it. Tests 1, 2, 3 and 6 are asserted there — the
script fails if they stop holding. Tests 4, 5 and 7 remain partial, and the table says so rather
than rounding them up.

| # | test | runs in | status |
|---|---|---|---|
| 1 | **Exactness** — $x_t \in \operatorname{span}(\mathbf{V}_r)$ to ~1e-12 for every candidate | `17_3` §2 | ✅ asserted; residuals 1.1e-15 to 4.2e-15, all six pairs |
| 2 | **Determinism** — identical across ≥5 calls and across BLAS/LAPACK builds | `17_3` §3 | ✅ asserted bit-identical across 5 calls; cross-BLAS **not yet run** |
| 3 | **Minimality / label integrity** — principal angles $(0,\theta)$ exactly | `17_3` §2 | ✅ asserted; `ang1` = 0.00e+00 and `ang2` $=\theta$, all six pairs |
| 4 | **Monotonicity** — family optimum $\le$ current, per target | `17_3` §7, §8 | ◐ asserted and holding at both targets, but the family figure comes from the stage-1 ray scan and the package figure from its own 200×200 grid. Needs stage 2 re-run on the winning plane for a like-for-like number, and a sweep over all rows to quantify what the label sample costs |
| 5 | **Variance retained** — adequacy, axis predictivity, RSS at $\mathbf{J}=\operatorname{diag}(1,1,0,\ldots)$, minimal versus current | — | ⬜ **not run.** Reuse `17_1` §13 / `17_2` §11's diagnostic block |
| 6 | **Thesis reproduction** — §2.5's 28.14° at $p=2, r=1$ | `17_3` §10 | ✅ asserted; 28.1411°, matching Procrustes to 6.75e-14 |
| 7 | **Cost** | `17_3` §9 | ◐ measured at $p=3$ and reported, not asserted — wall clock is machine-dependent. $p\ge4$ **not yet timed** |

One further check the script performs that is not on §4.1's list, because it is a property of the
*current* implementation rather than the proposal: `17_3` §8b measures the order dependence
recorded in §5.9 below.

### 6.6 What this does *not* fix

- **§4.4 stands.** This finds the best plane *within the biplot framework* — the smallest in-plane
  distance to the in-plane slice of the boundary. It is **not** the nearest boundary in the full
  $p$-dimensional feature space, and the improved search must not be read as implying otherwise.
- **§5.4 stands, unchanged.** The drawn curve is still a cross-section of the decision surface, not
  a view of it. Choosing the plane better does not make the slice representative.
- **§4.3 and §5.6 stand.** Global fidelity is still traded for local accuracy. The minimal rotation
  should *reduce* that trade for a given eigenplane, since it is the smallest rotation achieving
  exactness — but that is test 5, and it has not yet been measured.
- **§4.1 point 4 stands.** Shapley attribution and sparse counterfactuals still inherit whichever
  plane is drawn. The search makes that choice principled and reproducible; it does not remove the
  dependence. Sparse counterfactuals mix coordinates of the target and the counterfactual, so they
  generally leave the plane and remain approximate in the picture — no worse than at present, but
  not fixed here.

---

## 7. References

- PhD thesis, "Optimal rotation of the biplot," §3.5–3.5.2 (theoretical source for Part A).
- `R/local_cf.R` — `.generate_pairs()` (pair enumeration), `.bl_rotate()` (Procrustes rotation),
  `bl_find_local_cf()` (the full per-pair search and selection pipeline).
- `documentation/mahalanobis_technical_note.md` — full derivation of the $\mathbf{W}$ metric, the
  Cholesky-based inversion, and the cross-pair selection criterion behind §3.3.
- `2 implementation_summary.txt` §4.2, §4.2.1, §4.8 — code-mechanics description of the rotation
  and the local counterfactual search.
- `1 Foundation intro documents.txt` §10.3, §10.4, §10.7 — original design constraints (CVA
  scatter-matrix requirement, hull-polygon invalidation under rotation, the `max_pairs` cap).
- Original PhD research code (branch `original_PhD_code`): `scripts/1.3 Optimal Rotation.R`
  (`Biplot_rotation()`, the `pc.distz` cross-pair criterion) and `scripts/1 Functions for biplot.R`
  (`Wmat` construction).
