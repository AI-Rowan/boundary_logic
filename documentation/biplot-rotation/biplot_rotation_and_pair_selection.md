# Local Biplot Rotation and Selection of the Plane

`boundarylogic` explains one observation by drawing a 2D biplot that passes through it, then reading off the
nearest decision boundary on that biplot as the counterfactual. Two mechanisms make this work:

- a **rotation** that forces the target onto the biplot plane;
- a **search** that chooses *which* plane to rotate, so that the counterfactual is as near as the framework allows.

This note covers both, in five parts:

- §2: the original PhD method, from the thesis and from the research code;
- §3: the method as the package implements it, mapped line by line onto the code;
- §4: the current problems;
- §5: the proposed solution;
- §6: future developments.

Appendix A holds the measurements behind §3–§5, and Appendix B the practical notes for reproducing the package's
objects.

**Companion note.** [`../boundary-distance-and-robustness.md`](../boundary-distance-and-robustness.md) covers the
rest of the distance machinery:

- the metric $\mathbf{W}$ (the pooled within-class covariance), its definition and its Cholesky inverse;
- Phase 2's distance to the boundary;
- the per-variable robustness measures.

This note uses $\mathbf{W}$ only to choose the plane under CVA, or under PCA with `distance = "mahalanobis"`
(§3.4).

**Sources:**

- **The thesis:** the final version of 27 February 2026 (`1 FINAL_ADRIAAN ROWAN_27022026.pdf`, kept outside the
  repository in `Explain XAI code/reference-materials/`), cited by section and printed page (printed page = PDF
  page − 16). Until 2026-11-05 this note cited the 30 September 2025 draft
  (`../Adriaan Rowan - 13551361 - Thesis - 30 Sept 2025.pdf`), whose section numbers and pages differ.
- **The research code:** branch `origin/original_PhD_code`, cited as `scripts/<file>:<line>`.
- **Package source:** `R/local_cf.R` and `R/projection.R`.
- **Package design documents:** `1 Foundation intro documents.txt` and `2 implementation_summary.txt`.

**Files in this folder.** Every measurement in this note is reproduced by the R scripts beside it; `README.md` has
the details. In brief:

| script | setting | backs |
|---|---|---|
| `17_1a_toy-2d-to-1d-rotation.R` | 5 points, 2 variables → 1D biplot | §2.3, A.4, A.8 |
| `17_1_loan-gam-pca-3d-3var-rotation.R` | 221 loans, 3 features, GAM — lossless | §3.4 (PCA), §4.2, A.2–A.7, Appendix B |
| `17_2_loan-gam-pca-4var-slices.R` | 221 loans, 4 features, GAM — lossy | §4.1, §5.1, A.1, A.4, A.5 |
| `17_3_minimal-rotation-verification.R` | both of the above, plus iris; console output only | §3.4 (CVA), §5, §6.1, A.2, A.9, A.10 |

**Run `17_3` first.**

- **What it checks:** it asserts the geometry §5 depends on, the package's determinism (A.2), and the CVA identity
  of §3.4. It reports the rest.
- **Paths:** paths *inside* the scripts point at each other and at this note. Paths of the form `R/local_cf.R`,
  `development/…` and `documentation/…` refer to the package repository.
- **When measured:** all numbers were measured on 2026-10-30, and the iris and loan figures in §3.4 on 2026-10-31,
  against the package on branch `paper_nr1_changes` (R 4.6.0).
- **Which selector:** the measurements use the Mahalanobis pair selector, the package's default for both projections
  until 2026-11-05. Since then the default follows the thesis by projection (§3.4, §4.4), and the scripts pin
  `distance = "mahalanobis"` so that every number here is reproduced unchanged (re-run on 2026-11-05).
- **Which search:** since 2026-11-05 `bl_find_local_cf()` keeps only boundary points inside the convex hull of the
  training data and the target on each rotated plane (`hull = TRUE`, `documentation/local-convex-hull.md`). The
  measurements here concern the rotation and the pair search before that filter, so the scripts also pin
  `hull = FALSE`, the search they measured.

---

## 1. Purpose and motivation

A 2D biplot is built once, from the whole training set, using either PCA or CVA. It is a *global* approximation:

- every point is projected onto the same fixed plane;
- the decision boundary drawn on that plane is the model's boundary as seen through that one projection.

For any individual observation, projection and prediction are, in general, **not the same operation**. Projecting a
point onto the biplot plane and reading off its position relative to the drawn boundary can disagree with what the
classifier `f(X)` says about the *original, full-dimensional* point. The reason is that the plane is an
approximation:

- it captures only the two directions of greatest global variance (PCA) or greatest class separation (CVA);
- it does not generally contain any particular observation's exact location relative to the true, possibly curved,
  decision surface in the full feature space.

The thesis illustrates the failure with the iris dataset.

- **In 3D:** a three-variable PCA biplot (no dimensionality loss) correctly places a target point in the "class 1"
  region, matching the true model prediction.
- **In 2D:** reducing to the standard two-dimensional biplot (dropping the third, lowest-variance component)
  projects the same point orthogonally onto the plane, and it lands in the "class 0" region.

The model is right; its 2D visual approximation is wrong for this point, because near it the decision surface is
oblique to the biplot plane.

For an explanation method whose purpose is to visualise *why* a specific observation was classified the way it was,
the picture can then visibly contradict the number. The rotation (§2.1) fixes exactly this, for the one point being
explained. The search (§2.4) then chooses among the planes the rotation can produce.

---

## 2. The original PhD method

The thesis defines the method in three places:

- **the rotation:** §3.5–3.5.2 (pp.52–61);
- **the choice of projection and distance:** §3.8.3–3.9.3 (pp.82–93);
- **the search over component pairs:** §4.3.1 (pp.113–118).

The research code that produced the thesis's results is on the branch `origin/original_PhD_code` (§2.6).

### 2.1 The rotated component biplot (thesis §3.5)

**The three-point construction.** The rotation is framed as a **Procrustes problem**: the optimal one-to-one
transformation that maps points in one configuration onto points in another (Cox & Cox, 2008). Ordinarily
Procrustes analysis aligns two whole point clouds. Here it is specialised to a single point, the target observation
$x_{1\times p}$, which is centred and, for a PCA biplot, possibly scaled. The single point is embedded as **three
points**, symmetric about the origin:

$$
\ddot{\mathbf{X}}_{3\times p} =
\begin{bmatrix}
-x_1 & \cdots & -x_p \\
0    & \ddots & 0 \\
x_1  & \cdots & x_p
\end{bmatrix}
$$

- The middle row is the origin, a fixed point under any rotation.
- The two outer rows are $-x$ and $x$.

Aligning this configuration is equivalent to aligning the vector $x$ itself, because the origin contributes nothing
and $-x$ is redundant with $x$. Casting it as three points lets the standard Procrustes machinery solve it.

**Two projections of the same three points.** $\ddot{\mathbf{X}}$ is projected two ways, using the loading matrix
$\mathbf{V}$ (the $p\times p$ matrix from the SVD of the centred, possibly scaled, training data):

$$
\mathbf{Y} = \ddot{\mathbf{X}}\mathbf{V}
\qquad\text{(the full $p$-dimensional space)}
$$

$$
\mathbf{Z} = \ddot{\mathbf{X}}\mathbf{V}\mathbf{J}^{*} = \ddot{\mathbf{X}}\mathbf{V}_r
\qquad\text{(the $r$-dimensional biplot space, $\mathbf{V}_r$ = the columns of $\mathbf{V}$ used for the biplot)}
$$

If $r < p$, $\mathbf{Z}$ is zero-padded to $p$ columns. For $r=1,\,p=2$:

$$
\mathbf{Z}_{3\times p} = \begin{bmatrix} -z_1 & 0 \\ 0 & 0 \\ z_1 & 0 \end{bmatrix}, \qquad
\mathbf{Y}_{3\times p} = \begin{bmatrix} -y_1 & -y_2 \\ 0 & 0 \\ y_1 & y_2 \end{bmatrix}
$$

$\mathbf{Y}$ is where the target really sits; $\mathbf{Z}$ is where the biplot currently puts it. The rotation finds
an orthogonal transformation of $\mathbf{V}$ that makes the two agree for the target.

**The Procrustes solution.**

$$
\mathbf{C} = \mathbf{Y}'\mathbf{Z}, \qquad
\mathbf{C} = \mathbf{U}\boldsymbol{\Lambda}\mathbf{V}^{*\prime}\ \text{(SVD)}, \qquad
\mathbf{A} = \mathbf{V}^{*}\mathbf{U}'
$$

Here $\mathbf{V}^{*}$ is the SVD's right-singular-vector matrix, unrelated to the loading matrix $\mathbf{V}$.
$\mathbf{A}$ is the optimal $p\times p$ orthogonal rotation. The optimal dilation is

$$
\rho = \frac{\operatorname{tr}(\mathbf{Z}\mathbf{A}\mathbf{Y}')}{\operatorname{tr}(\mathbf{Z}'\mathbf{Z})}
$$

and the thesis states that it "is not required": it would uniformly rescale every point on the plane, to recover
magnitude that dimension reduction genuinely lost. Only $\mathbf{A}$ is used.

**The new projection.**

$$
\mathbf{V}^{new} = \mathbf{V}\mathbf{A}'
$$

The thesis calls biplots built from $\mathbf{V}^{new}$ **rotated component biplots (RC biplots)**. Their defining
property is that the target lies exactly on the plane,

$$
\hat{\mathbf{x}}_i' = \mathbf{x}_i'\mathbf{V}_r^{new}(\mathbf{V}_r^{new})' = \mathbf{x}_i'
$$

so the model's classification and the biplot's agree: $f(\hat{\mathbf{x}}_i) = f(\mathbf{x}_i)$.

- **The prediction (back-projection) matrix** is $(\mathbf{V}^{new})' = \mathbf{A}\mathbf{V}'$.
- **A point $z$ in rotated Z-space** maps back to X-space as
  $\hat{\mathbf{x}}' = \mathbf{z}'\mathbf{A}\mathbf{V}'\mathbf{J}'$.

### 2.2 The CVA variant (thesis §3.5.2, p.61)

The same machinery applies to the CVA transformation matrix $\mathbf{M}$, or to the two-group optimised
$\mathbf{M}^{opt}$:

$$
\mathbf{M}^{new} = \mathbf{M}^{opt}\mathbf{A}', \qquad
(\mathbf{M}^{new})^{-} = \mathbf{A}(\mathbf{M}^{opt})^{-1}
$$

**The general rule.** "The only change required between using either the CVA and PCA biplots is replacing the
matrices $\mathbf{V}$ and $\mathbf{V}'$ with $\mathbf{M}$ and $\mathbf{M}^{-1}$" (§3.9.3, p.92).

**The cost specific to CVA.**

- Rotation pulls the **class mean points** off the canonical plane.
- Their $(k-1)$-dimensional predictions, which CVA is built to represent exactly, stop being exact.

PCA has no analogous special point. CVA rotation trades class-mean fidelity for target-point fidelity.

### 2.3 Worked examples (thesis §3.5–3.5.1)

**5-point, 2D→1D toy.**

- **The set-up:** five points in two dimensions are reduced to a one-dimensional PCA biplot, with point 5 as the
  target.
- **The rotation:** before rotation the target's sample predictivity is 0.72. A 28° rotation with dilation
  $\rho = 1/1.13$ brings it to exactly 1.
- **The cost:** the sum of squared residuals for the dataset rises from 0.28 to 4, out of a total of 14. Variable
  2's axis gains axis predictivity at the expense of variable 1.

`17_1a` reproduces the angle and the dilation exactly, but not the predictivity or residual figures (A.8). The
thesis's numbers and the script's should not be quoted together.

**Iris 3D→2D.**

- A three-variable PCA biplot is shown first at full rank, then reduced to 2D.
- The reduction misclassifies a target relative to the model's prediction.
- After the rotation, the projected point moves into the matching region, and its interpolated value equals the
  true value: the identity of §2.1, demonstrated.

### 2.4 Searching all component pairs (thesis §3.9.3 p.93, §4.3.1 pp.113–118)

The rotation turns *a given* plane. The thesis's local explanation tries every plane that two of the components
can span:

- **Every pair is tried.** "the projection matrix is recalculated for all possible two-dimensional biplots formed
  from pairs of components". That is $p!/\bigl(2!(p-2)!\bigr)$ combinations, each selected by a $p\times 2$ matrix
  $\mathbf{J}_i$ of zeros with ones in the rows of the chosen components (p.117).
- **Each is rotated** so its plane passes through the target. "The new rotated component biplot is created using
  the first two components of $\mathbf{M}^{new}$" (p.115).
- **Constraints** (p.92, p.117):
  - **actionability:** for the Pima data, Pregnancies and Age cannot decrease;
  - **feasibility:** each feature stays within its training range;
  - **plausibility:** the counterfactual "must be selected from a local decision boundary that is within the
    convex hull of the training data", and boundary areas that contain no data points are discarded.
- **Within a plane**, the nearest boundary point to the target is found on the biplot plane (p.86).
- **Across planes**, "the nearest counterfactual is then selected as the counterfactual that corresponds to the
  shortest change vector measured as the Mahalanobis distance of the change in the original $p$-dimensional space"
  (pp.117–118):

$$
\mathbf{d}_i = (\mathbf{z}_{i,cf} - \mathbf{z}_i)(\mathbf{M}_i^{new})^{-}, \qquad
\hat{\mathbf{x}}_{cf}:\ \min_i\ \mathbf{d}_i'\mathbf{W}^{-1}\mathbf{d}_i,
$$

  with $\mathbf{W}$ "the total within groups sums-of-squares-and-products (SSP)".
- **If no plane yields a constrained counterfactual,** "it is assumed that no constraint counterfactual exists for
  this data point" (p.118).

**The Pima example** (p.118, Figure 4.11):

- seven variables give 21 biplots;
- 11 of them yield constrained counterfactuals;
- the nearest comes from components 1 and 7.

The thesis's conclusion lists the search as a contribution: it "systematically constructs biplots using all pairs
of components, thereby approximating different slices of the decision surface" (§4.4, p.125).

### 2.5 The distance criterion in the thesis

The thesis fixes the cross-plane criterion **by projection**:

- **p.86 (§3.9.1):** "the nearest counterfactual … is selected as the counterfactual that corresponds to the
  shortest change vector $\mathbf{d}_i$ measured in the full dimensional space using either **Euclidean distance
  for PCA biplots and Mahalanobis distance for CVA biplots**", with
  $\mathbf{d}_i = (\mathbf{z}_{i,cf}-\mathbf{z}_{i,0})\mathbf{V}_2'$ or
  $(\mathbf{z}_{i,cf}-\mathbf{z}_{i,0})\mathbf{M}_2^{-}$.
- **p.83 (§3.8.3) and p.91 (§3.9.3):** "In the CVA biplot … the Euclidean distances on the biplot plane
  correspond to Mahalanobis distances in the original feature space." By contrast, a standardised PCA biplot
  "ignores the correlation" between features.
- **p.92 (§3.9.3):** PCA is for the initial data investigation. "The CVA biplot is then used to generate local and
  global interpretations, since it ensures that the estimated counterfactual reflects the shortest Mahalanobis
  distance in the original feature space", with the groups defined by the model's four outcomes
  (TP, TN, FP, FN).
- **p.124 (§4.4):** the conclusion describes a "dual approach": Mahalanobis distance in the feature space for
  selection, and the $L_2$ distance on the plane for consistency with the picture.

**The two halves of the "dual approach" are one quantity under CVA.** The thesis's own p.91 statement says so, and
§3.4 proves it for every rotated plane. The in-plane $L_2$ step therefore already finds the Mahalanobis-nearest
point, and the cross-plane comparison ranks planes by the same measure. **Under PCA** the thesis's rule is
Euclidean throughout; it is the package's default since 2026-11-05 (§3.4).

### 2.6 The research code (`origin/original_PhD_code`)

**`Biplot_rotation()`** (`scripts/1.3 Optimal Rotation.R:32-81`):

- **Target vector:** centred; scaled only for PCA with `standardise = TRUE`. CVA forces no scaling (42–47).
- **Y:** `Y <- rbind(-X.st[tdp,], rep(0, num_vars), X.st[tdp,])` (50–52).
- **Rotation:** `YVr` is the pair's two columns, zero-padded. Then `A <- svd(t(YV) %*% YVr)` with
  `A = v %*% t(u)` (59–66).
- **ρ:** computed as a "Rho diagnostic" (68–69), never applied, and not returned.
- **Output:** the full `Vrho <- V %*% t(A)` and `tVrho <- solve(Vrho)`. Callers then use columns 1–2
  (`4.3:91`, `4.3:215`).

**The pair search** (`scripts/4.3 Different Eigenvalues.R`). No other script loops over pairs. In order:

1. **Enumerate:** `n_contrib <- t(combn(num_vars, 2))`, i.e. all $\binom{p}{2}$ pairs with no cap (49–50).
2. **Data:** the target followed by the **test** data, `rbind(biplot_dataT[tdp, ], biplot_dataT)` (27).
3. **Rotate:** each pair is rotated (77–95).
4. **Grid and hull:** a new grid is built in the rotated plane. The **convex hull is rebuilt in that plane** from
   the rotated target and test data, and contours are clipped to it (`1.1:131-142`, `1.4:121-150`).
5. **Contour rule:** closed contours are kept only if they enclose the target (`1.4:130-141`).
6. **Filters:** contour points are back-projected, then two filters run:
   - `set_filters` (range expressions only; there are no fixed or direction constraints);
   - `train_ranges` (`1.4:159-169`).
7. **Re-score:** each contour keeps only the points whose prediction is on its side (`1.4:174-184`).
8. **Nearest point per contour:** by 2D squared Euclidean distance with `max.col(-d2)`, whose default
   `ties.method = "random"` breaks ties at random (`1.4:40-45`). Then the nearest opposite-class contour is taken
   (`1.4:237-266`).
9. **Distance per pair** (162–166):

   ```r
   vec_to_boundary_Z <- as.matrix(biplot_search_ev$B_boundary_z_value - biplot_grid_rotate$Z.st[1, , drop = FALSE])
   vec_to_boundary   <- vec_to_boundary_Z %*% tVr
   pc.distz[i] <- (as.matrix(vec_to_boundary, nrow = 1)) %*% Wmat_inv %*% t(as.matrix(vec_to_boundary, nrow = 1))
   ```

   The in-plane difference is back-projected by the rotated prediction matrix, and $\mathbf{W}^{-1}$ is applied in
   X-space. A Euclidean alternative is commented out (167).
10. **Select:** the argmin of `pc.distz` after dropping zeros (178–186). The winning pair is then rotated and
    searched again.

**`Wmat`.** It is biplotEZ's `bp$Wmat` (`scripts/3.1 Biplot search - Train data.r:183-189`):

- the within-class SSP with **no divisor**;
- on the four-outcome factor (TP/TN/FP/FN, `pred_col`);
- from the same `CVA()` call that gives `Vmat = Lmat`.

The search runs under **CVA only**. A PCA biplotEZ object has no `Wmat`, so `solve(Wmat)` errors. Since biplotEZ's
CVA loadings satisfy $\mathbf{M}'\mathbf{W}\mathbf{M}=\mathbf{I}$, `pc.distz` is exactly $\lVert\Delta z\rVert^2$
(§3.4). The research code's Mahalanobis selection was always the in-plane distance selection.

**Defects in the research code**, recorded so that the package's behaviour is not mistaken for a change of
method:

- **The B = 0 fallback.** When no contour survives a pair, `Biplot_boundary_search()` returns
  `B_boundary_z_value = 0` (`1.4:100-113`). The pair's "distance" is then the target's distance from the origin,
  and the pair is not excluded.
- **Inconsistent plot data.** The in-loop plot coordinates come from uncentred data (`4.3:99`); the final step's
  come from `bp$X` (`4.3:225`).
- **The default flow does not rotate.** `3.2` sets `Rotate <- F` (`3.2:90`). The Shapley script `4.4` uses that
  unrotated result by default; the pair-search result is only a commented-out alternative (`4.4:25-29`).

---

## 3. The current implementation

### 3.1 The rotation — `.bl_rotate()`

`R/local_cf.R:45-66` is a line-for-line port of §2.1 (from the research code's `Biplot_rotation()`):

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

Three implementation choices, all consistent with the thesis:

1. **No dilation.** ρ is never computed (§2.1).
2. **`solve(Vrho)`, not a transpose.**
   - **The thesis's form:** the prediction matrix $(\mathbf{V}^{new})' = \mathbf{A}\mathbf{V}'$ is a transpose.
   - **When that is valid:** only when $\mathbf{V}^{new}$ is orthogonal. That is true under PCA, where $\mathbf{V}$
     is Euclidean-orthonormal.
   - **Under CVA:** $\mathbf{V}$ is $\mathbf{W}$-orthonormal ($\mathbf{V}'\mathbf{W}\mathbf{V} =
     \mathbf{I}/(n-g)$; §3.4), so a transpose would not be its inverse.
   - The general inverse handles both projections.
3. **Columns 1–2 are returned directly.**
   - **Why columns 1–2:** the package's biplot is always two-dimensional. `YVr_padded` is non-zero only in its
     first two columns, so the target's rotated coordinates land in columns 1–2 of `Vrho` whichever pair built
     `Vr`.
   - **Equivalent to the original:** the research code returns the full matrix and its callers take columns 1–2
     (§2.6).

### 3.2 Pair enumeration

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

**Order and cap:**

- Pairs come in colexicographic order, $(1,2), (1,3), (2,3), (1,4), (2,4), (3,4), \ldots$.
- They are capped at `max_pairs`, default `21L` in `bl_find_local_cf()` (10 until 2026-10-31).
- The design notes (`1 Foundation intro documents.txt` §10.7) also cap the count at the $\binom{p}{2}$ pairs that
  exist.

**The cap is not in the thesis or the research code**, which try every pair. The default covers every pair up to
$p=7$ ($\binom{7}{2}=21$), so the Pima example's search matches the thesis's. It first leaves pairs out at $p=8$
($\binom{8}{2}=28$) (§4.7).

### 3.3 The per-pair pipeline

Inside `bl_find_local_cf()` (`R/local_cf.R:590`; the pair loop is lines 650–840; line numbers as of 2026-11-04), each pair goes through five
stages:

1. **Rotate.** `.bl_rotate(x_target_st, V, pair)` (§3.1) gives `Vr_rot` / `tVr_rot`.
2. **Project.** The target and all training observations go into the pair's rotated Z-space (`Z_target`,
   `Z_train_rot`).
3. **Grid and score.**
   - An $m\times m$ grid ($m=200$ by default) spans the training data's rotated range, plus 10% padding.
   - It is back-projected to X-space and scored through the model in chunks of 50,000 rows.
4. **Extract and filter the boundary contour.** `grDevices::contourLines()` at $cutoff \pm b_{margin}$, then four
   filters on the vertices (`R/local_cf.R:708-759`):
   - (a) **opposing class:** only contour segments on the class *opposite* the target's prediction;
   - (b) **actionability:** the `set_filters` constraints (`"decrease"`, `"increase"`, `"fixed"`, `c(min,max)`;
     `"integer"` is not applied here, but to the counterfactual after the search, §3.5).
     They run first because they are typically the most selective, so the loop can stop early
     (`if (!any(keep)) next`);
   - (c) **feasibility:** `train_ranges` (per-feature training min/max), applied to the survivors of (b);
   - (d) **model re-score:** every surviving vertex is back-projected and re-scored, and dropped if the model
     disagrees with the contour's side of `cutoff`. This is the most expensive filter, which is why (b) and (c)
     run first.

   **The hull, rebuilt per plane (since 2026-11-05).** The global biplot's polygon was built for the unrotated
   plane, so rotation invalidates its coordinates and the search never uses it. With `hull = TRUE` (the default)
   a convex hull of the training data and the target is built on each rotated plane, and only the contour
   vertices inside it are kept, before filters (b)-(d) (`../local-convex-hull.md`). Until 2026-11-05 no hull
   was applied, and that is the search every measurement in this note used: the scripts pin `hull = FALSE`.
5. **Within-pair nearest vertex.** Among the survivors, the nearest to the target is found by 2D Euclidean distance
   in the rotated plane (`.nearest_idx_block()`). It is the exact minimum, taking the first vertex on an exact tie,
   so the step is deterministic (`../nearest-vertex-tie-break.md`).

### 3.4 The selector and its metric

The winning pair is the one with the smallest distance in X-space, by `distance`. The default, `"auto"` (since
2026-11-05), follows the thesis rule (§2.5): the **squared Mahalanobis distance** under CVA and the **Euclidean
distance on the total-SD scale** under PCA (the pair loop of `bl_find_local_cf()`, `R/local_cf.R`):

$$
d_M^2(v) = v^{\mathsf T}\mathbf{W}^{-1}v, \qquad d_x(v) = \Big(\sum_j (v_j / s_j)^2\Big)^{1/2}, \qquad v = B_x - x_{obs}
$$

where $s_j$ is the total standard deviation of feature $j$ (`X_sd`), as in `bl_find_boundary()`'s `dist_x`.

```r
dist_mah <- if (!is.null(metric_inv)) {
  v <- as.numeric(B_x[1L, var_names]) - x_obs
  as.numeric(t(v) %*% metric_inv %*% v)
} else {
  NA_real_
}
dist_x <- sqrt(.cf_sq_distance(
  matrix(as.numeric(B_x[1L, var_names]) - x_obs, nrow = 1L),
  "euclidean", metric_inv, X_sd))
```

**The metric.** `metric_inv` is $\mathbf{W}^{-1}$, computed once in `bl_build_projection()` and reused for every
pair and every target.

- **What $\mathbf{W}$ is:** the pooled within-class covariance, divisor $n-g$, on the unstandardised features.
- **Which class factor:** the CVA groups (`metric_type = "W_cva"`), else the binary class (`"W_binary"`, the PCA
  case).
- **How it is inverted:** by Cholesky, with a ridge fallback.

The companion note has the definition and the inversion. `distance = "euclidean"` selects by `dist_x`; until
2026-11-05 it selected by the in-plane distance `dist_z`. All three distances are stored for every pair
(`all_distances`, `all_distances_x`, `all_distances_mahalanobis`).

#### Under CVA the in-plane distance *is* the Mahalanobis distance, for every pair

**What differs between pairs.** Each pair has its own rotation $\mathbf{A}_i$, so its own prediction matrix
$(\mathbf{M}\mathbf{A}_i')^{-1} = \mathbf{A}_i\mathbf{M}^{-1}$. The same in-plane step $\Delta z$ therefore maps to
a different X-space step $\Delta x$ in every pair, and the Euclidean length of $\Delta x$ differs too. What does
*not* differ is its Mahalanobis length.

**The algebra.**

1. **biplotEZ's construction.** biplotEZ (2.2) builds the CVA loadings as
   $\mathbf{M} = \mathbf{W}_{SSP}^{-1/2}\mathbf{V}_{svd}$, with $\mathbf{W}_{SSP}$ the within-class SSP (it is
   returned as `bp$Wmat`). So $\mathbf{M}'\mathbf{W}_{SSP}\mathbf{M} = \mathbf{I}$, and
   $\mathbf{M}^{-1} = \mathbf{M}'\mathbf{W}_{SSP}$.
2. **The step in X-space.** Let $\mathbf{J}$ select columns 1–2. In pair $i$, an in-plane step maps to
   $\Delta x = \Delta z\,\mathbf{J}'\mathbf{A}_i\mathbf{M}'\mathbf{W}_{SSP}$.
3. **Its Mahalanobis length:**

$$
\Delta x\,\mathbf{W}_{SSP}^{-1}\,\Delta x'
= \Delta z\,\mathbf{J}'\mathbf{A}_i(\mathbf{M}'\mathbf{W}_{SSP}\mathbf{M})\mathbf{A}_i'\mathbf{J}\,\Delta z'
= \Delta z\,\mathbf{J}'\mathbf{A}_i\mathbf{A}_i'\mathbf{J}\,\Delta z'
= \lVert\Delta z\rVert^2 .
$$

$\mathbf{A}_i$ enters only as $\mathbf{A}_i\mathbf{A}_i' = \mathbf{I}$ and cancels, whatever the pair. The package's
metric is $\mathbf{W} = \mathbf{W}_{SSP}/(n-g)$, so for the package

$$
d_M^2 = (n-g)\cdot \texttt{dist\_z}^2 \quad\text{for every pair.}
$$

**The geometric picture.**

- $\mathbf{M}$ whitens the within-class scatter: in canonical coordinates $u = x\mathbf{M}$, Mahalanobis distance
  in X is ordinary Euclidean distance.
- $\mathbf{A}_i$ is a rotation *inside* that space, and each pair's plane is a differently tilted flat slice of it.
- Distance measured within a flat slice equals distance in the whole space, for every slice.

Different pairs give different slices and so different boundary points and different $\Delta z$, but each $\Delta z$
is measured in the same metric. **The two selectors therefore rank the pairs identically.** The within-pair
Euclidean step (§3.3, stage 5) already finds the Mahalanobis-nearest vertex.

**Conditions.** The identity needs three things:

1. **CVA:** under PCA, $\mathbf{V}'\mathbf{V}=\mathbf{I}$ instead (below).
2. **The CVA's own $\mathbf{W}$:** the package builds `metric` from the same class factor it passes to `CVA()`.
3. **A displacement in the drawn plane:** true for the search's point, because the rotation puts the target in the
   plane and the boundary vertex is on it. False for a restored "fixed" or a rounded "integer" counterfactual
   (§3.5, §4.6).

**It also holds for 2-group CVA.** biplotEZ then adds the second axis by `CVAlowdim(low.dim = "sample.opt")`,
which builds `cbind(M[, 1], M[, -1] %*% <orthogonal matrix>)` and so keeps $\mathbf{M}'\mathbf{W}\mathbf{M} =
\mathbf{I}$. The other option, `"Bhattacharyya.dist"`, would not; the package does not use it.

**Measured.**

| case | $\max\lvert\mathbf{V}'\mathbf{W}\mathbf{V}(n-g)-\mathbf{I}\rvert$ | $d_M^2/\texttt{dist\_z}^2$, every pair | source |
|---|---|---|---|
| iris, CVA on TP/TN/FP/FN ($n-g = 146$) | 1.29e-14 | 146 (6 of 6 pairs) | `17_3` §11 |
| iris, 2-group CVA ($n-2 = 148$) | 6.13e-15 | 148 (6 of 6) | `17_3` §11 |
| loan harness, `bl_results_v2` ($n = 15{,}996$, $g = 4$) | 1.01e-13 | 15,992 in the constrained, free and external searches (6 of 6 each) | harness driver |
| loan harness, 2-group CVA | 8.81e-14 | 15,994 | harness driver |

The same in-plane step $\Delta z = (0.7, -0.4)$, sent through each pair's prediction matrix (iris, target 1;
`17_3` §11):

| pair | prediction matrix − pair (1,2)'s (max abs) | $\lVert\Delta x\rVert^2$ (Euclidean) | $\Delta x\,\mathbf{W}^{-1}\Delta x'$ | ratio to $\lVert\Delta z\rVert^2$ |
|---|---|---|---|---|
| (1,2) | 0 | 12.7235 | 94.9 | 146 |
| (1,3) | 5.434 | 15.3715 | 94.9 | 146 |
| (2,3) | 9.878 | 85.0205 | 94.9 | 146 |
| (1,4) | 16.561 | 200.3843 | 94.9 | 146 |
| (2,4) | 20.649 | 202.9754 | 94.9 | 146 |
| (3,4) | 13.099 | 239.3651 | 94.9 | 146 |

In the package's own search on that target, the Mahalanobis selector and the smallest in-plane distance both pick
pair (3,4); on the loan harness's free search, both pick (3,4). The thesis states the identity in prose (p.91, §2.5).

#### Under PCA they are different metrics

- **The in-plane distance:** $\mathbf{V}$ is orthonormal in the (standardised) feature space, and so is every
  rotated basis. For in-plane points, `dist_z` is the Euclidean distance between the standardised feature
  vectors. That is comparable across pairs; it is the thesis's PCA criterion (§2.5). The X-space distance `dist_x`
  divides each feature by the same total SD that standardised the PCA, so for the search's point
  `dist_x` = `dist_z` (`tests/testthat/test-distance_by_method.R`); `dist_x` also measures the off-plane
  `"integer"` candidates (§3.5).
- **The Mahalanobis selector:** it uses $\mathbf{W}$ = `W_binary`, the within-class covariance by the 0/1 class on
  unstandardised features. That is a different metric, and the ratio varies by pair: 1.30 to 1.58 on iris
  (`17_3` §11). It was the PCA default until 2026-11-05.
- **How often they disagree:**
  - over the 221 rows of `17_1`'s PCA set-up, the two selectors pick different pairs on 51 (23%);
  - in the standardised-iris fixture of `tests/testthat/test-boundary.R` (PCA, target 1), Z-distance would pick
    (1,4) and Mahalanobis (3,4).

**So, by projection** (since 2026-11-05):

| `distance` | CVA | PCA |
|---|---|---|
| `"auto"` (default) | `"mahalanobis"` | `"euclidean"` |
| `"mahalanobis"` | the thesis's criterion; ranks the pairs as the in-plane distance does | `W_binary` Mahalanobis; the PCA default until 2026-11-05, not the thesis's rule |
| `"euclidean"` | Euclidean in X-space on the total-SD scale; a different metric, not the thesis's CVA rule | Euclidean in X-space on the total-SD scale, equal to the in-plane distance at the search point; the thesis's rule (p.86) |

Until 2026-11-05 `"euclidean"` was the in-plane distance itself, identical to `"mahalanobis"` under CVA.

### 3.5 Winner selection, the "fixed" restore and the "integer" rounding

The pair loop of `bl_find_local_cf()` keeps a running minimum over the tested pairs:

```r
selector <- if (distance == "mahalanobis") dist_mah else dist_x
if (is.finite(selector) && selector < best_selector) {
  best_selector <- selector
  ...
  best_result <- list(
    B_z = B_z_local, B_x = B_x, B_pred = B_pred[[1L]],
    dist_z = dist_z, dist_x = dist_x, dist_mahalanobis = dist_mah,
    Z_target = Z_target, best_pair = pair,
    Vr_rot = Vr_rot, tVr_rot = tVr_rot,
    Z_train_rot = Z_train_rot, Zgrid = Zgrid, grid_prob = grid_prob,
    col_value = col_value, ct_local = ct_local,
    xseq = xseq, yseq = yseq, min_val = min_val, max_val = max_val
  )
}
```

**What is reported:**

- **Every pair's distances** are kept, winners or not: `dist_z`, `dist_x` and `dist_mahalanobis`
  (`all_distances`, `all_distances_x`, `all_distances_mahalanobis`).
- **The winning pair's objects** are what the result reports and plots: its rotated loading matrices, grid,
  contour and boundary point (`best_pair`, `B_x`, `B_z`, …).
- **No pair can be a winner by accident.** A pair with no surviving vertex is skipped (`next`); the research
  code's B = 0 defect (§2.6) has no counterpart.

**After the search: the "fixed" restore** (the post-search block of `bl_find_local_cf()` and `.restore_fixed()`, `R/local_cf.R`).

- **During the search:** a `"fixed"` constraint is a tolerance of 0.5 model units, because a vertex of the 2D
  contour rarely holds a feature at exactly its observed value.
- **After the winner is chosen:** each fixed feature the search moved is set back to its observed value. If that
  point still crosses the boundary and is within `train_ranges`, it becomes the counterfactual (`B_x`, `B_pred`).
  Otherwise no counterfactual is reported (strict, since 2026-11-04). `fixed_check` records the outcome.
- **The restored point is generally not on the plane:**
  - `B_z` is its projection, which is what the local plot draws;
  - `dist_z`, `dist_x`, `dist_mahalanobis` and the per-pair tables still describe the search's point.
- Details: [`../fixed-constraint-restore.md`](../fixed-constraint-restore.md) §1–§4.

**Then the "integer" rounding** (`.round_integers()`, since 2026-11-04).

- **What it does:** for each `"integer"` feature that is not fixed and has moved off a whole number, it tries the
  whole values (raw units) just below and above. It keeps the combinations that still cross the boundary and
  meet every constraint and `train_ranges`, and takes the nearest in X-space by the same distance as the pair
  search (`.cf_sq_distance()`; until 2026-11-05 the Euclidean option measured the projections on the plane).
- **Its result:** the chosen point replaces `B_x`, `B_z` (its projection, `.project_local()`) and `B_pred`. With
  no valid combination no counterfactual is reported. `integer_check` records the outcome.
- **Like the restored point, it lies off the plane,** and the distances still describe the search's point.
- Details: [`../fixed-constraint-restore.md`](../fixed-constraint-restore.md) §5.

### 3.6 How the implementation differs from the original

| aspect | thesis / research code | package |
|---|---|---|
| projection for the local search | CVA (TP/TN/FP/FN); the research code cannot run it under PCA | PCA or CVA |
| rotation | Procrustes, no dilation, columns 1–2 | the same (§3.1) |
| pairs tried | all $\binom{p}{2}$ | up to `max_pairs = 21L`: all pairs for $p \le 7$ (§3.2) |
| grid built over | the target and the test data (research code) | the training data's rotated range |
| plausibility | the convex hull, rebuilt in each rotated plane; boundary areas without data discarded; closed contours must enclose the target | the convex hull of the training data and the target, rebuilt in each rotated plane (`hull = TRUE`, since 2026-11-05; before, no hull, §3.3); no enclosure rule |
| actionability | range expressions (`set_filters`) | `"decrease"`, `"increase"`, `"fixed"` (with the post-search restore), `c(min,max)`, `"integer"` (rounded after the search) |
| feasibility | `train_ranges` | `train_ranges`, filter (c) |
| nearest vertex in a plane | 2D Euclidean, random tie-break | 2D Euclidean, exact minimum (§3.3) |
| cross-plane criterion | $\mathbf{d}'\mathbf{W}_{SSP}^{-1}\mathbf{d}$ under CVA; Euclidean under PCA (thesis only) | by `distance = "auto"`: $v'\mathbf{W}^{-1}v$ with $\mathbf{W} = \mathbf{W}_{SSP}/(n-g)$ under CVA; total-SD Euclidean in X-space under PCA (since 2026-11-05; before, Mahalanobis for both) (§3.4) |
| a pair with no boundary point | B = 0, kept with a spurious distance | skipped |
| the explanation used downstream | the research code's default flow is unrotated (§2.6) | always the rotated winner |

The divisor $n-g$ scales every distance by the same constant, so it never changes the winner. Under CVA the package
selects exactly as the thesis does, and under PCA it does so by default since 2026-11-05 (§4.4).

---

## 4. Current problems

Ranked by consequence. Each states the problem and why it happens, with a headline number; Appendix A has the
evidence.

### 4.1 The plane is under-determined, and cross-platform reproducibility is untested

**Why it happens.**

- $\ddot{\mathbf{X}} = [-x;\,0;\,x]$ has **rank 1**: all three rows are multiples of $x$.
- Write $\ddot{\mathbf{X}} = a\,x'$ with $a = (-1, 0, 1)'$, and $u = \mathbf{V}'x$. The Procrustes cross-product
  is then $\mathbf{C} = 2\,u\,w'$, with $w = (u_{[\text{pair}_1]}, u_{[\text{pair}_2]}, 0, \ldots, 0)$. That is
  also rank 1.
- So the Procrustes optimum is **a single vector constraint**, $\mathbf{A}'u \propto w$, on a $p\times p$
  orthogonal matrix.
- At $p=4$ the optimal rotations form a 3-dimensional family, and the planes they produce a 2-dimensional one
  (§5.1).
- The SVD of $\mathbf{C}$ has $p-1$ **exactly degenerate** zero singular values. The plane's second direction is
  whatever LAPACK emits for that degenerate subspace.

**Consequences.**

- **Inputs that differ only in their last bits give different planes.** Rotating from `svd(X_std)$v` instead of
  `bl_result$V`, a difference of 4.65e-15, gives planes **65.91°** apart. The counterfactual moves with the plane
  (A.1).
- **An instance in the package itself.**
  - **What happened:** reordering the CVA group levels (TP/TN/FP/FN), a change with no mathematical effect,
    perturbed `V` enough to change every Phase 3 result in the loan harness. Target 13's constrained
    counterfactual moved from 0.547 to 0.509.
  - **What the package does now:** it pins the level order (CLAUDE.md §6).
- **Within one platform the search is deterministic** (A.2):
  - `.bl_rotate()` is bit-identical for identical inputs;
  - the nearest vertex is the exact minimum.
- **Across platforms it is untested.** A different BLAS/LAPACK, OS or R version may legitimately return a different
  plane, and so a different counterfactual, from identical data and code. For a method meant to support auditable,
  published explanations, this is a material risk to disclose.

**What is *not* affected**, to keep this in proportion:

- the target's position is exact whichever plane is drawn;
- the counterfactual is a genuine boundary point, re-scored through the model (filter (d)), just not a canonical
  one;
- the cross-pair distance is measured in X-space, which the rotation cannot distort;
- Phases 1 and 2 never call `.bl_rotate()`: `bl_build_projection()`, `bl_build_grid()`, `bl_find_boundary()`,
  `bl_robustness()` and `bl_surrogate()`;
- `bl_set_scaling()` is display-only and cannot reach `V`, the grid or the plane.

The research code has the same construction (§2.6): this is inherited from the original method, not a porting
error.

### 4.2 The pair labels are weak

- **The naming problem.** The search enumerates named eigenvector pairs, but each rotated plane is only loosely tied
  to the pair that names it. At $p=4$, only **21%** of the selected (2,3) plane's second direction lies in the
  (2,3) eigenplane (`17_2` §9).
- **The same axes from different labels.** At $p=3$, all three labels give rotated axes closest to the same
  principal components (A.3).
- **No cheap screen exists.** The rotation angle $\theta$ does not predict the overlap: among 74 cases with
  $\theta<5°$, the overlap spans 0% to 100% (A.3).

So "Best pair: (3,4)" in the output, and the wording in any paper, imply an eigenplane that is not what is drawn. It
is more accurate to describe the search as *sampling $\binom{p}{2}$ planes through the target* than as
*enumerating eigenplanes*. Since 2026-11-05 the package's help pages and messages say "component pair" and state
that the selected plane need not coincide with the plane of the two components that name it; Paper 1 uses the
same wording.

### 4.3 The local plot's orientation is arbitrary

- **Where the frame comes from.** The plane does not fix an ordered orthonormal frame within it. The local plot
  uses columns 1–2 of `Vrho`, which is LAPACK's completion.
- **The scale of it.** Two bases can span the same plane to 2.2e-16 and still differ entry-wise by 1.61 (Appendix
  B). The picture's orientation and axis layout are therefore incidental.
- **No user correction is possible.** `rotate_deg` has been removed from the local plots, since the orientation
  belongs to the method. §5.2 proposes a frame convention that fixes it.

### 4.4 The metric: one selector under CVA; the PCA default (resolved 2026-11-05)

- **Under CVA there was one selector, presented as two.**
  - Until 2026-11-05 `distance = "euclidean"` selected by the in-plane distance, which already measures the
    Mahalanobis distance, so the option could not change the result (§3.4).
  - Since then `"euclidean"` is the X-space Euclidean distance on the total-SD scale: under CVA a different metric,
    and not the thesis's rule. The default, `"auto"`, is Mahalanobis under CVA.
  - Before 2026-10-31 the package's documentation claimed the opposite (that Z-distances are not comparable across
    pairs).
- **Under PCA the default departed from the thesis (resolved 2026-11-05).**
  - The thesis's rule is standardised Euclidean (p.86). Until 2026-11-05 the default was `W_binary` Mahalanobis,
    which selects a different pair on 23% of `17_1`'s rows (§3.4).
  - That default also mixed metrics: within a plane the nearest vertex is chosen by standardised Euclidean
    distance, across planes by `W_binary` Mahalanobis. The point judged "nearest" was not nearest in the measure
    that judged the "best" plane.

**Decision (2026-10-31): the PCA default is deferred.**

- **Paper 1 compares counterfactuals by Mahalanobis distance only**, under CVA, where it equals the in-plane
  distance (§3.4). The PCA departure does not touch its results.
- **The PCA default (`W_binary`) stays as it is** until the rotation problem (§4.1) is fully addressed (§5). The
  choice of metric under PCA (§5.4) is taken then, not before.

**Resolution (2026-11-05, on the user's instruction).** The default became `distance = "auto"`, the thesis's rule by
projection, before the rotation problem was addressed:

- **Under PCA:** the pairs are compared by the Euclidean distance in X-space, each feature divided by its total SD.
  It equals the in-plane distance at the search point, so one metric now acts within and across planes.
- **Under CVA:** Mahalanobis, as before.
- **Phase 2's variable importance** follows the same rule (total SD under PCA).

Paper 1 (CVA) is unchanged. The measurements in this note use the Mahalanobis selector, which `17_1`–`17_3` now
pin.

### 4.5 Constraints are applied after contouring

`train_ranges` and `set_filters` filter contour vertices *after* the contour is extracted (filters (b)–(c)).

- **The bias.** The reported point is "nearest among the survivors", and no nearer feasible point is ever sought.
- **The scale.** `train_ranges` alone removes 18–56% of boundary candidates per pair (Appendix B).
- **It changes the winner.** A replication of the pair search without that filter agrees with the package on only
  76 of 221 rows.

### 4.6 Off-plane counterfactuals

Three counterfactuals leave the drawn plane:

- **The restored "fixed" counterfactual** (§3.5), which resets coordinates to the target's.
- **The rounded "integer" counterfactual** (§3.5), which moves coordinates onto whole values.
- **The sparse counterfactual** (`bl_find_sparse_cf()`), which mixes coordinates of the target and the
  counterfactual.

**What follows:**

- **On the plot:** all three are drawn at their projections, which need not lie on the drawn boundary.
- **In the output:** for a restored or rounded counterfactual, `dist_z`, `dist_x` and `dist_mahalanobis` describe
  the search's point, not the reported one.
- **For the identity:** the §3.4 identity does not apply to them.

### 4.7 Only pair (1,2) can be forced, and `max_pairs` caps the search

- **Forcing a pair.** `max_pairs = 1L` pins the search to (1,2), because pairs are emitted colexicographically.
  There is no way to force any other pair.
- **Capping the search.**
  - The default `max_pairs = 21L` tries every pair up to $p=7$. At the Pima example's $p=7$ it tries all 21 pairs,
    including the thesis's winning pair (1,7), which is 16th in the package's order.
  - From $p=8$ ($\binom{8}{2}=28$) it tries fewer pairs than the thesis: every pair involving the highest components
    is left out.
  - Until 2026-10-31 the default was `10L`, which left pairs out from $p=6$ and never reached (1,7).

### 4.8 A target at the centroid

**Verified on 2026-11-05** (the README's Pima CVA model, with the target at `X_center`):

- if $x_t = 0$, the Procrustes matrix is zero, so the SVD returns the identity and $\mathbf{A} = \mathbf{I}$;
- every pair label then returns the unrotated plane of components (1,2): the plane for pair (3,4) equals
  `V[, 1:2]`;
- the search still runs and can find a counterfactual.

**The guard (since 2026-11-05):** `bl_find_local_cf()` warns once, before the pair loop, when
`sqrt(sum(x_target_st^2)) < sqrt(.Machine$double.eps)`: the target lies at the centre of the biplot, so the
Procrustes rotation is undefined and the local plane is arbitrary. `.bl_rotate()` is unchanged, and so are the
results. §5.2's edge cases give the convention the proposed minimal rotation would follow.

### 4.9 Inherent limitations — documented, not to be "fixed"

- **Global fidelity is traded for local accuracy.** Rotating to fit one target worsens every other point's
  representation. In the thesis's toy, the residual rises from 0.28 to 4 of 14 (§2.3).
  - The trade is a **redistribution** across variables, and which variable is sacrificed depends on where the plane
    lands (A.7).
  - A rotated biplot is a local picture only.
- **It is not the nearest boundary in the full space.**
  - Rotation resolves *classification consistency* for the target. It does not find the nearest boundary in
    $p$-dimensional space.
  - The search finds the nearest boundary *within the planes it draws*.
- **CVA class means leave the plane** (§2.2).
- **A 2D biplot is a slice of the decision surface, not a view of it.**
  - Two planes through the same target cut one surface along two different curves. For `17_1`'s row 6 those curves
    are 57% of the data's extent apart, with counterfactuals 0.345 apart (A.5).
  - Neither plot is wrong; they are different cross-sections.
- **Tilt is not a criterion.**
  - "Optimal" and "best" refer only to the smallest distance from target to counterfactual.
  - A shallower plane is easier to read, not a better answer (A.6).
- **Adjacent mechanisms are not this one.**
  - `rotate_deg` turns the *global* biplot for display only (`plot_biplotEZ()`, `plot.bl_surrogate()`). The local
    plots have no `rotate_deg` and stop if one is passed.
  - Phase 2's per-feature distances are a separate measure; see the companion note.

---

## 5. Proposed solution

**This is a proposal; nothing in the shipped package works this way.** It changes every Phase 3 counterfactual, so
it needs an explicit decision (CLAUDE.md §9) and the adoption path of §5.6.

**What it does:**

- replaces the rank-deficient Procrustes construction with a closed-form minimal rotation;
- replaces the $\binom{p}{2}$-label sample with a reproducible search over every plane that contains the target;
- uses one metric end to end;
- treats feasibility as a constraint.

**Scope.** Everything stays inside the 2D biplot framework.

- **The rotation's job:** put the target *on* the drawn plane, so the local reading is exact.
- **The search's job:** find the nearest boundary *available within that framework*.
- **Out of scope:** searching the full $p$-dimensional space (§4.9, §6.3).

### 5.1 Exactness requires only that the plane contains the target

The biplot's reconstruction is the orthogonal projector onto the plane, so

$$
\hat{x} = \mathbf{V}_r(\mathbf{V}_r'\mathbf{V}_r)^{-1}\mathbf{V}_r'x = x
\quad\Longleftrightarrow\quad
x \in \operatorname{span}(\mathbf{V}_r).
$$

**What this condition is.**

- It is a condition on the **plane**, not on the frame.
- It does not depend on how $\mathbf{V}_r$ is scaled or oriented, so it reads the same for PCA and CVA.
- Nothing beyond "the plane contains the target" is needed for the exact local interpretation.

The admissible family is

$$
\mathcal{A} \;=\; \{\,2\text{-planes } P \subset \mathbb{R}^p \;:\; x_t \in P \,\}.
$$

Every such plane is $\operatorname{span}(x_t, d)$ for a direction $d$ in the $(p-1)$-dimensional space orthogonal
to $x_t$, taken modulo sign. So $\mathcal{A} \cong \mathbb{P}^{p-2}$, of dimension $p-2$:

| $p$ | $\dim\mathcal{A}$ | shape |
|---|---|---|
| 2 | 0 | a single plane — determined, which is why the thesis's toy teaches so cleanly |
| 3 | 1 | a pencil of planes about the axis $\operatorname{span}(x_t)$ — one angle |
| 4 | 2 | a two-parameter family |

**Two counts, both right.**

- **Rotations.** The rotations with $\mathbf{A}'u \propto w$ form a coset of $\operatorname{Stab}(u) \cong
  SO(p-1)$, of dimension $(p-1)(p-2)/2$: 3 at $p=4$. That is §4.1's count.
- **Planes.** Containing a fixed line is *two* conditions on $\operatorname{Gr}(2,4)$, not one, leaving
  $4-2=2$. `17_2` §8a measures the family's tangent dimension numerically as **2**.
- **The difference** is the in-plane spin of §4.3.

**The toy is the exception, not the rule.** The thesis's 2D→1D toy is the only case in the family where the
rotation is determined (A.4). That is why it teaches cleanly, and why the intuition it builds fails from $p=3$
onwards.

The plane choice is entirely free with respect to exactness. §5.2 spends that freedom **deterministically**, §5.3
spends it **well**.

### 5.2 A minimal Givens rotation, in closed form

Given a named eigenplane $P_0 = \operatorname{span}(\mathbf{V}_{[,a]}, \mathbf{V}_{[,b]})$, write
$\mathbf{V}_{r0} = \mathbf{V}_{[,(a,b)]}$ and the target's shadow $x_P = \mathbf{V}_{r0}\mathbf{V}_{r0}'x_t$. For
$0 < \lVert x_P\rVert < \lVert x_t\rVert$:

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

- $\{u,w\}$ is orthonormal, because $r \perp P_0$ and $u \in P_0$.
- $\mathbf{R}$ rotates by $\theta$ inside $\operatorname{span}(u,w)$ and is the identity on its complement.
- $\lVert x_P\rVert^2 + \lVert r\rVert^2 = \lVert x_t\rVert^2$, so $\cos^2\theta + \sin^2\theta = 1$ as written.

| property | proof |
|---|---|
| **Exactness**: $x_t \in \operatorname{span}(\mathbf{V}_r)$ | $\mathbf{R}u = \cos\theta\,u + \sin\theta\,w = x_t/\lVert x_t\rVert$, so $\mathbf{R}'x_t = \lVert x_t\rVert u \in P_0$, i.e. $x_t \in \mathbf{R}P_0$ |
| **Deterministic**: no SVD, no degenerate subspace, no tie-break | every quantity is a closed-form expression in $x_t$ and $\mathbf{V}_{[,(a,b)]}$ |
| **Minimal**: the smallest rotation achieving exactness | $\theta$ *is* the angle between $x_t$ and $P_0$; any $\mathbf{S}$ with $x_t \in \mathbf{S}P_0$ needs $\mathbf{S}'x_t \in P_0$ with $\lVert \mathbf{S}'x_t\rVert = \lVert x_t\rVert$, so it must move $x_t$'s direction by at least $\theta$ |
| **`proj_pair` means what it says**: principal angles $(\mathbf{R}P_0, P_0) = (0,\theta)$ | the direction of $P_0$ orthogonal to $u$ lies in $\operatorname{span}(u,w)^\perp$, which $\mathbf{R}$ fixes pointwise — so the two planes share it |
| **$O(p^2)$** — cheaper than the current $O(p^3)$ SVD | matrix assembly only |

**Measured at $p=4$** (target row 49; A.9):

- the plane's first principal angle to the eigenplane its label names is **0 exactly** for every pair, where the
  current construction is up to 55.62° off;
- the basis is bit-identical across five calls.

**Edge cases**, each with a convention:

- $\lVert x_P\rVert = 0$ (target orthogonal to $P_0$, $\theta = 90°$): $u$ is undefined. Take
  $u = \mathbf{V}_{[,a]}$.
- $\lVert x_P\rVert = \lVert x_t\rVert$ (target already in $P_0$): $\theta = 0$ and $\mathbf{R} = \mathbf{I}$.
  Guard the $0/0$ in $w$.
- $x_t = 0$ (target at the centroid): every plane contains it. Return the unrotated eigenplane and warn (§4.8).

**The in-plane frame, fixed by convention** (resolving §4.3). Nothing *measurable* depends on the frame, so it is
safe to fix.

- **Recommended:** $\mathbf{V}_r = (\mathbf{R}\mathbf{V}_{[,a]},\ \mathbf{R}\mathbf{V}_{[,b]})$, the named axes
  carried through the rotation.
  - RC1 and RC2 then read as "PC $a$ and PC $b$, tilted by $\theta$".
  - The "columns 1–2 regardless of pair" convention of §3.1 is no longer needed.
- **Alternative:** the target-aligned frame $(x_t/\lVert x_t\rVert,\ \cdot)$.
  - It puts every target on the horizontal axis, making figures comparable across targets.
  - The cost is loading arrows that move between targets.

**PCA versus CVA.** The construction is written for a Euclidean-orthonormal $\mathbf{V}$, as under PCA (`17_1` and
`17_2` assert $\mathbf{V}'\mathbf{V} = \mathbf{I}$ to 1e-10).

- Under CVA, carry it out in $\mathbf{W}$-whitened coordinates (the canonical space of §3.4) and map back.
- The exactness condition of §5.1 is basis-independent, so it survives the change of metric.
- The CVA cost of §2.2 is unaffected.

**The thesis's worked example is preserved.**

- At $p=2,\ r=1$ the family is a single line, so the minimal rotation and Procrustes must agree.
- On `17_1a`'s five points (target row 5) both give **28.1411°**, agreeing to 6.75e-14 and reproducing the
  thesis's 28.14°.

### 5.3 Searching the admissible family, not $\binom{p}{2}$ labels

**Why labels are not enough.** §5.2 makes each label map to one determinate plane. The labels still sample only
$\binom{p}{2}$ points of a $(p-2)$-dimensional continuum, so "best of the labels" is not best of anything well
defined.

**This is the thesis's own future work.** Its §6.3 (p.144), "Grand Tour of data instead of component selection",
proposes extending the search beyond fixed pairs to a smooth sequence of planes. Each frame would be rotated to place
the target on the biplot, and the candidates scored by Mahalanobis distance. $\mathcal{A}$ has an explicit
parameterisation, so it can be searched directly.

**Parameterisation.**

- Let $q = x_t/\lVert x_t\rVert$.
- Let $\mathbf{U} \in \mathbb{R}^{p\times(p-1)}$ be a deterministic orthonormal basis of $q^\perp$: the
  Householder reflector carrying $q \to e_1$, minus its first column.

Then every admissible plane is

$$
P(d) = \operatorname{span}(q,\ \mathbf{U}d), \qquad d \in \mathbb{R}^{p-1},\ \lVert d\rVert = 1,\ d \equiv -d.
$$

| $p$ | search space | recommended search |
|---|---|---|
| 3 | one angle $\varphi \in [0,\pi)$ | **exhaustive** 1° sweep, then golden-section refinement |
| 4 | hemisphere of $S^2$ | deterministic Fibonacci lattice, then local refinement from the best few |
| $\ge 5$ | $\mathbb{P}^{p-2}$ | the $\binom{p}{2}$ minimal-rotation planes as seeds plus local refinement; **document as a local optimum, not exhaustive** |

There is no RNG anywhere: every candidate set is a deterministic function of $p$ and the resolution.

**The objective, in two stages.**

- **Stage 1** ranks planes cheaply by a **ray scan**:
  - cast $K$ rays from the target at equal angles within the plane;
  - bisect each to its **first** cutoff crossing;
  - take the smallest squared Mahalanobis distance over the rays.

  The nearest in-plane boundary point is reached by some ray, and along that ray it is the first crossing.
- **Stage 2** runs the existing §3.3 pipeline on the best few planes, so the reported boundary and every plotting
  object keep the package's definition.

**Measured** (A.10):

| target | improvement in squared distance | where it comes from |
|---|---|---|
| $p=3$, row 6 | **5.7%** over the package | an exhaustive 1° sweep of 180 planes |
| $p=4$, row 49 | **4.4%** | a 200-point lattice |

At $p=3$, pinning to (1,2) lands within 0.05% of the *worst* plane in the family.

- **Cost:** about **1.5×** the current search when batched across planes, and an order of magnitude more when
  implemented naïvely. That is affordable per target, but a cost, not a saving.
- **A caveat:** stage 1 and the package use different discretisations, so the percentages are indicative until
  stage 2 is run on the winning plane (§6.1, test 4).

### 5.4 One metric end to end, and feasibility as a constraint

1. **One metric, chosen by projection** (resolving §4.4; **the default was adopted on 2026-11-05**: the thesis's
   standardised Euclidean under PCA, Mahalanobis under CVA).
   - **Under CVA** there is already one metric, because the in-plane distance is Mahalanobis (§3.4). The
     `distance` option can be documented as having no effect, or removed for CVA.
   - **Under PCA**, choose one metric and use it both within and across planes:
     - either the thesis's standardised Euclidean, which needs no change within the plane;
     - or `W_binary` Mahalanobis within the plane as well. The vertices are already back-projected for filter (d),
       so this costs no extra model evaluations.

     Which one was a methodological decision for the user, who chose the thesis's rule (2026-11-05).
2. **Feasibility as a constraint, not a filter** (resolving §4.5).
   - Restrict the grid and contour extraction to the feasible region, so the search seeks the nearest *feasible*
     point.
   - `set_filters` actionability likewise becomes a box or half-space constraint rather than filter (b).

### 5.5 Reporting

**What changes in the output:**

- **The distance table.** Replace `all_distances_mahalanobis` (one value per sampled label) with the objective over
  the searched family, plus the runner-up local minima, so near-ties are visible instead of silently resolved.
- **At $p=3$ this is a curve.** Distance against $\varphi$ is a genuinely informative per-target diagnostic: for
  row 6 the reported distance moves by a factor of 1.13 across the family.
- **Continuity.** Keep $\theta$ and the nearest named pair in the output, for continuity with `best_pair`.

### 5.6 Adoption path

**Land it as new options, defaults unchanged:**

- `rotation = c("procrustes", "minimal")`;
- `plane_search = c("pairs", "family")`.

**Switch the defaults** only once the §6.1 acceptance tests have been run and reported across a full training set,
not just the two targets measured here.

---

## 6. Future developments

### 6.1 Acceptance tests for §5

All are run by `17_3`; the "runs in" column names its section. Tests 1, 2, 3 and 6 are asserted, so the script fails
if they stop holding.

| # | test | runs in | status |
|---|---|---|---|
| 1 | **Exactness** — $x_t \in \operatorname{span}(\mathbf{V}_r)$ to ~1e-12 for every candidate | `17_3` §2 | ✅ asserted; residuals 1.1e-15 to 4.2e-15, all six pairs |
| 2 | **Determinism** — identical across ≥5 calls and across BLAS/LAPACK builds | `17_3` §3 | ✅ asserted bit-identical across 5 calls; cross-BLAS **not yet run** |
| 3 | **Minimality / label integrity** — principal angles $(0,\theta)$ exactly | `17_3` §2 | ✅ asserted; `ang1` = 0.00e+00 and `ang2` $=\theta$, all six pairs |
| 4 | **Monotonicity** — family optimum $\le$ current, per target | `17_3` §7, §8 | ◐ asserted and holding at both targets, but the family figure comes from the stage-1 ray scan and the package figure from its own 200×200 grid. Needs stage 2 re-run on the winning plane for a like-for-like number, and a sweep over all rows to quantify what the label sample costs |
| 5 | **Variance retained** — adequacy, axis predictivity, RSS at $\mathbf{J}=\operatorname{diag}(1,1,0,\ldots)$, minimal versus current | — | ⬜ **not run.** Reuse `17_1` §13 / `17_2` §11's diagnostic block |
| 6 | **Thesis reproduction** — the 28.14° at $p=2, r=1$ | `17_3` §10 | ✅ asserted; 28.1411°, matching Procrustes to 6.75e-14 |
| 7 | **Cost** | `17_3` §9 | ◐ measured at $p=3$ and reported, not asserted — wall clock is machine-dependent. $p\ge4$ **not yet timed** |

**Properties of the current package that `17_3` also asserts** (not on this list):

- §5 and §8b: within-platform determinism (A.2);
- §11: the CVA identity (§3.4).

### 6.2 Further work on the current method

Each item needs an explicit instruction (CLAUDE.md §9). Ranked as in §4:

1. **Cross-platform reproducibility** (§4.1). Run `17_3` and the harness on another OS or BLAS (macOS, Linux,
   OpenBLAS) and compare. This matters for Paper 1, where a reviewer on another platform may not reproduce the
   quoted counterfactuals.
2. **The pair terminology** (§4.2). Decide how to describe a pair in output and in papers; report $\theta$ or the
   overlap.
3. **A frame convention** for the local plot (§4.3, §5.2).
4. **The metric under PCA** (§4.4, §5.4). **Resolved 2026-11-05:** the thesis's Euclidean rule is the PCA default
   (`distance = "auto"`). Still to do: add a test that pins the CVA identity: the ratio is $n-g$ for every pair, including 2-group
   CVA. `17_3` §11 asserts it today; `tests/testthat/` does not.
5. **Feasibility as a constraint** (§4.5, §5.4).
6. **Off-plane counterfactuals** (§4.6). Mark them on the plot, or recompute `dist_z` / `dist_mahalanobis` for the
   reported point.
7. **Pair control** (§4.7):
   - a `pairs =` argument to force any pair;
   - a `max_pairs` default tied to $\binom{p}{2}$, so the default search covers every pair for any $p$, as the
     thesis's does. The default of 21 already does so up to $p=7$.
8. **A centroid guard** (§4.8). Done on 2026-11-05: verified, and `bl_find_local_cf()` now warns.

### 6.3 What the proposal does not fix

- **§4.9 stands.** The search finds the best plane *within the biplot framework*, not the nearest boundary in the
  full $p$-dimensional space. The drawn curve is still a slice of the decision surface. Global fidelity is still
  traded for local accuracy:
  - the minimal rotation should *reduce* that trade for a given eigenplane, because it is the smallest rotation
    achieving exactness;
  - that is test 5, and it is not yet measured.
- **Downstream objects still inherit the plane.** Shapley attribution and sparse counterfactuals inherit whichever
  plane is drawn. The search makes that choice principled and reproducible; it does not remove the dependence.
- **Off-plane points are unchanged.** Sparse, restored and rounded counterfactuals still leave the plane (§4.6).

---

## Appendix A — Measurements

Unless stated otherwise:

- **`17_1`** is the workhorse: target row 6, with the planes taken from the package's own `Vr_rot`. A 3×3 PCA basis
  is a lossless rotation, so anything that changes between its views is the rotation alone.
- **`17_2`** is lossy: four features into three PCs, target row 49.

### A.1 The plane under perturbation (`17_2`; supports §4.1)

**Setting:** loan-default GAM model, PCA, $p=4$, target row 49; the search selects pair (2,3).

- **The rank.** The singular values of $\mathbf{C}$ for pair (2,3) are `1.63e+00, 6.71e-17, 0, 0`: rank 1 of 4,
  with a 3-dimensional null space. All six pairs are rank 1 (`17_2` §8).
- **Identical inputs give an identical plane** (A.2).
- **Inputs that differ in their last bits give different planes.**
  - The first principal angle between two such planes is always 0°, because both contain the target.
  - The *second* angle is the arbitrary part:

  | Input change (pair rotated onto) | max abs change | 2nd principal angle | source |
  |---|---|---|---|
  | re-orthonormalise V only, (1,2) | 5.55e-16 | 24.09° | `17_2` §8b |
  | random noise on V, scale 1e-16, (1,2) | 5.55e-16 | 24.09° | `17_2` §8b |
  | random noise on V, scale 1e-15, (1,2) | 1.44e-15 | 71.57° | `17_2` §8b |
  | random noise on V, scale 1e-12, (1,2) | 1.61e-12 | 24.09° | `17_2` §8b |
  | random noise on V, scale 1e-08, (1,2) | 1.01e-08 | 36.87° | `17_2` §8b |
  | `svd(X_std)$v` instead of `bl_result$V` (numerically the same matrix), (2,3) | **4.65e-15** | **65.91°** | `17_2` §7b |
  | target standardised as a row of `scale()` instead of the package's `(x - X_center)/X_sd`, (1,2) | 3.33e-16 | 29.81° | `17_2` §7b |

- **The counterfactual moves with the plane.** At pair (1,2), noise of 1e-15 on V moves the nearest boundary point
  by 0.138 in feature space (√$d_M$ 1.1262 → 1.1088; `17_2` §8c).
- **The label overlap.** Only **21%** of the selected (2,3) plane's second direction (taken orthogonal to the
  target) lies inside $\text{span}(\mathbf{V}_{[,2]}, \mathbf{V}_{[,3]})$ (`17_2` §9).
- **Why it is easy to miss.** The target lies in every one of these planes (distance ~1e-15; `17_2` §14). The one
  thing a reader would naturally check, whether the target lands correctly, looks fine either way.

### A.2 Within-platform determinism (supports §4.1)

For identical inputs:

- `.bl_rotate()` returns a bit-identical basis;
- the nearest boundary vertex is the exact minimum.

So a repeated call, or the same call after any other package function, returns the identical counterfactual.

| target | calls | result |
|---|---|---|
| `17_1` set-up, row 102 | 5 identical full searches | `best_pair` (1,2) every call; squared $d_M$ (1,2) 0.229143, (1,3) 0.232776, (2,3) 0.232243 |
| `17_1` set-up, rows 6, 104, 140, 161 | 5 identical full searches each | one result per row |
| `17_2` set-up, row 49 | 7 mixed full and pinned searches (`17_3` §8b) | full 1.147879 and pinned 1.208670 on every call, in any order |

Asserted by:

- `17_3` §5 and §8b;
- `tests/testthat/test-boundary.R`, which checks for identical results under any RNG state and after
  `bl_find_boundary()` / `bl_surrogate()`.

Before 2026-10-30, repeated calls could differ through a random tie-break in the nearest-vertex step
(`../nearest-vertex-tie-break.md`).

### A.3 The pair labels at $p=3$ (`17_1` §12; supports §4.2)

**Setting:** row 6, rotating onto each named pair with `17_1`'s own copy of `.bl_rotate()` from `svd(X_std)$v`.

| named pair | $\theta$ | 2nd principal angle vs `span(V[,pair])` | 2nd direction inside the named eigenplane |
|---|---|---|---|
| (1,2) | 48.25° | 86.26° | **56.08%** |
| (1,3) | 3.03° | 3.32° | 99.94% |
| (2,3) | 41.59° | 86.97° | **44.35%** |

**These planes versus the package's.** They come from `svd(X_std)$v`, which is 1.4e-15 from the package's
`bl_result$V`.

- For (1,3) the package draws the same plane.
- For (1,2) its plane is a different valid completion, **86.42°** away (`17_1` §15c). A.6 measures that plane's
  tilt as 48.28°.

The two figures for "the (1,2) plane" are §4.1 in action, not a contradiction.

**Per column:** for all three pairs alike, rotated column 1 comes out closest to PC3 and column 2 closest to PC1.
Three different labels give the same resulting axes.

**$\theta$ does not predict the overlap.** Across all 663 (row, pair) combinations:

- Pearson $r = +0.33$, Spearman $+0.14$;
- among the 74 cases with $\theta < 5°$, the overlap spans the full range, 0% to 100%;
- row 6's pairing of $\theta = 3.03°$ with 99.94% is a coincidence: row 137 has $\theta = 1.49°$ with **0.1%**.

### A.4 The dimensional ladder and the toy (`17_1a`, `17_1`, `17_2`; supports §5.1)

| setting | free rotations | admissible planes (§5.1) |
|---|---|---|
| $p=2,\ r=1$ (the thesis's toy) | 0 — unique | 0 — unique |
| $p=3,\ r=2$ | 1 | 1 |
| $p=4,\ r=2$ | 3 | 2 (`17_2` §8a) |

**The toy.** Rotating onto **PC1** and onto **PC2** gives:

- the *identical* rotation matrix $\mathbf{A}$, retained axis, adequacy, RSS and sample predictivity;
- only $\rho$ differs (1.134058 versus 2.120244), and $\rho$ is discarded.

With $r=1$ the retained axis is forced to be the target's own direction ($|\cos| = 1.000000000000$ for all five
points), and in 2D only one line through the origin contains the target.

**At $p=3$ a first principal angle of 0° carries no information.**

- In 3-space any two 2-planes meet in a line. Three random plane pairs give `0.00, 83.34` / `0.00, 47.79` /
  `0.00, 87.81`.
- Nor is the shared line the target's direction: for row 6, only 0.4435 of that direction lies in the PC1–PC2
  plane, which is its unrotated sample predictivity.

### A.5 A slice, not a view (`17_1`; supports §4.9)

- **In 3D:** a figure of a 3-variable model shows the whole $p\approx0.5$ surface. Every rotated view holds the same
  surface, cloud and points; rotating only moves the camera.
- **In 2D:** a biplot shows where the model crosses 0.5 **within the drawn plane**, with the off-plane coordinate
  pinned at 0.

Row 6, pinned (1,2) versus the selected (1,3):

| | points | model $p$ along the curve | $\max|p-0.5|$ |
|---|---|---|---|
| pinned (1,2) | 700 | [0.49816, 0.50202] | 2.02e-03 |
| selected (1,3) | 498 | [0.49812, 0.50198] | 1.98e-03 |

**Two genuine boundary curves, far apart.**

- Both are genuine boundary curves, yet the nearest-neighbour distance between them is median **2.30**, max
  **3.85**, against a data span of **6.71**. They separate by 57% of the data's extent.
- The planes are **90.00°** apart.

**They share exactly one point, and this is forced.**

- Both planes contain the origin and the target. Two distinct planes meet in one line, so their intersection is the
  target's ray, $\operatorname{span}(x_t)$; verified $|\cos(n_B \times n_C,\ x_t)| = 1.000000000000$.
- The surface crosses that ray at $t=0.41223$. Both curves pass within 0.0036 and 0.0010 of that crossing, against
  a contour grid step of 0.0363.

**The counterfactuals differ.**

- $B_x = (+0.759, +0.062, +0.165)$ versus $(+0.861, +0.381, +0.084)$: **0.345 apart**.
- $\sqrt{d_M}$ is 1.2527 versus 1.1455.

A different slice gives a different nearest boundary point and so different advice.

### A.6 Tilt needs two angles, and is not a criterion (`17_1`; supports §4.9)

**Two numbers, not one.** Relative to the unrotated PC1–PC2 plane, a rotated plane's orientation needs **two**
numbers.

- Write it as $s_3 = s_1 s_{1\text{ax}} + s_2 s_{2\text{ax}}$.
- The component tilts are $\arctan s_1$ along PC1 and $\arctan s_2$ along PC2.
- The dihedral angle is their resultant, $\tan^2\theta = s_1^2 + s_2^2$ (verified to 6 dp).

**For row 6:**

- the pinned (1,2) plane is $(-48.26°, +2.04°)$, almost a pure PC1 tilt;
- the selected (1,3) plane is $(+65.85°, +88.64°)$.

**The admissible tilts form a line.** "Contains the target" is one equation, $s_1 a_1 + s_2 a_2 = a_3$, so the
admissible tilts form a **line** in the $(s_1, s_2)$ plane.

- For row 6 that line admits dihedrals from **48.25° to 89.05°**.
- The pinned plane came out at 48.28° (the shallow end), the selected plane at 88.64° (the steep end).
- The criterion and the tilt run in opposite directions here.

### A.7 What the rotation costs, per variable (`17_1`; supports §4.9)

**What is measured.** `adequacy` and axis predictivity index the **rows** of $\mathbf{V}$: they are per *variable*.
Row 6, $\mathbf{J} = \operatorname{diag}(1,1,0)$:

| | adequacy (income / score / rate) | RSS (of 660) | target sample predictivity | mean of the other 220 |
|---|---|---|---|---|
| unrotated | 0.5023 / 0.4978 / 0.9999 | **168.656** | 0.4435 | 0.7161 |
| rotated → (1,2) | 0.9988 / **0.0442** / 0.9570 | 229.465 | **1.000000** | 0.6710 |
| rotated → (1,3) | 0.9671 / 0.9609 / **0.0720** | 213.624 | **1.000000** | 0.6758 |

**What it shows:**

- **Different sacrifices.** The two rotations sacrifice different variables: `credit_score` in one,
  `loan_int_rate` in the other.
- **Exactness.** The target reaches sample predictivity exactly 1 in both, measured at 2.2e-15.
- **An arithmetic check.** The unrotated RSS equals the dropped component's variance $d_3^2$ exactly, which is a
  useful check on any reimplementation.

### A.8 The thesis toy's numbers (`17_1a`; supports §2.3)

**What reproduces:**

- **The angle and dilation.** `17_1a` reproduces the thesis's rotation angle (**28.14°**; measured 28.1411°) and
  dilation ($\rho = 1.134 \approx 1/1.13$) exactly.

**What does not:**

- **Predictivity and residuals.** The thesis reports sample predictivity 0.72 → 1 and RSS 0.28 → 4 of 14. The
  script's five points give **0.7776 → 1** and **4.307 → 4.927 of 11.4**.

The data behind the thesis's example is probably not the script's.

**On the residual figure:**

- RSS here is the **trace** of $\mathbf{X}'\mathbf{X}$, not the sum of its entries.
- Summing every entry is not a sum of squares, and it makes the rotation appear to *improve* global fit.
- With the trace, the unrotated RSS equals $d_2^2$ exactly and the rotated RSS is larger, as §4.9 requires.

### A.9 The minimal rotation's label integrity (`17_3` §2; supports §5.2)

**Setting:** $p=4$, target row 49. `ang1` is the first principal angle against the eigenplane the pair *names*.

| pair | $\theta$ | `ang1` minimal | `ang2` minimal | exactness residual | `ang1` current |
|---|---|---|---|---|---|
| (1,2) | 38.4522° | **0.00e+00** | 38.4522° | 1.40e-15 | 1.1111° |
| (1,3) | 34.9160° | **0.00e+00** | 34.9160° | 1.99e-15 | 5.1471° |
| (2,3) | 70.1512° | **0.00e+00** | 70.1512° | 1.08e-15 | 55.6173° |
| (1,4) | 19.8488° | **0.00e+00** | 19.8488° | 4.19e-15 | 17.2008° |
| (2,4) | 55.0840° | **0.00e+00** | 55.0840° | 2.69e-15 | 5.4437° |
| (3,4) | 51.5478° | **0.00e+00** | 51.5478° | 2.34e-15 | 51.3428° |

**The two constructions compared:**

- **Under the minimal rotation**, the plane shares a direction with its named eigenplane exactly, and its only
  departure is $\theta$.
- **Under the current construction**, it shares nothing. That is the mechanism behind the 0.3–70% overlaps `17_2`
  §9 measures at $p=4$, and behind three labels giving the same axes (A.3).

The returned basis was bit-identical across five calls (`17_3` §3).

### A.10 The family search (`17_3` §6–§9; supports §5.3)

**$p=3$, target row 6, exhaustive 1° sweep of all 180 planes:**

| | squared $d_M$ | $d_M$ |
|---|---|---|
| sweep minimum (at $\varphi = 134°$) | **1.237877** | **1.1126** |
| sweep median | 1.490022 | 1.2207 |
| sweep maximum | 1.579964 | 1.2570 |
| minimal rotation, label (1,2) | 1.580762 | 1.2573 |
| minimal rotation, label (1,3) | 1.293857 | 1.1375 |
| minimal rotation, label (2,3) | 1.568419 | 1.2524 |
| **the package, as it stands** | 1.312277 | 1.1456 |

**Three findings:**

1. **The search wins.** It improves on the package by **5.7% in squared distance** (2.9% in distance).
   Monotonicity holds as it must: every plane the current method draws lies in $\mathcal{A}$.
2. **Pinning to (1,2) lands at the worst plane** for this target: 1.580762, within 0.05% of the sweep's maximum of
   1.579964. A pinned pair is not a neutral default.
3. **The curve has two local minima.** Near-ties are real, and should be reported rather than silently resolved.

**At $p=4$**, target row 49, over a 200-point Fibonacci lattice:

- best label (1,2): **1.147994**;
- the package: **1.147879**;
- the family's best: **1.097519**, a **4.4%** improvement in squared distance.

**Discretisation.**

- **What differs:** the sweep's figures come from the stage-1 ray scan, the package's from its 200×200 grid.
- **Which way it biases:** the ray scan's 10° angular resolution can only *miss* the nearest point, so if anything
  it understates the improvement.
- **What is still needed:** a like-for-like figure needs stage 2 on the winning plane (§6.1, test 4).

**Cost.**

- **Plane at a time:** at $p=3$ the sweep takes **6–12×** the package's search time, though the prediction counts
  are comparable (129,600 vs 120,000).
- **Why:** the cost is per-call overhead. The package's grid scores 40,000 rows in one call per pair; a per-plane
  ray scan scores 36 rows per call.
- **Batched across planes:** all 180 × 36 = 6,480 rows go in one call per march/bisect step, 20 calls in all. That
  brings it to **1.4–1.7×**.
- **The recorded ratios** (`17_3` §9): 6.3× / 1.65× on 2026-09-04, and 12.1× / 1.44–1.63× on 2026-10-30. They move
  with machine load; the ratio, not the seconds, is the point.

---

## Appendix B — Reproducing the package's objects

- **Use the package's own objects rather than re-deriving them.**
  - `bl_local$Vr_rot` is the rotated basis the local biplot uses.
  - `bl_local$ct_local` (the local plane) and `bl_result$biplot_grid$ct` (the global biplot) are the boundary
    contours, already in their plane's 2D coordinates. They plot verbatim.
- **To reproduce the package's plane, use its exact inputs:** `bl_result$V`, and the target standardised as the
  package does, `(x - X_center) / X_sd`.
  - Rotating from `svd(X)$v` (1.4e-15 away), or from a target computed another way (3e-16 away), can land on a
    different plane (§4.1; `17_1` §15c, `17_2` §7b).
  - With identical inputs, the mirror reproduces `.bl_rotate()` exactly (`17_1` §15a, `17_2` §7a).
- **Compare planes by their projection matrices** $\mathbf{P} = \mathbf{V}_r\mathbf{V}_r'$, never entry-wise. Even
  with the plane matching exactly, the in-plane basis can differ by a large rotation: 1.61 entry-wise at a plane
  difference of 2.2e-16 (§4.3).
- **`max_pairs = 1L` pins the search to pair (1,2)**, because `.generate_pairs()` emits colexicographically. It is
  the only supported way to force a plane, and it makes a pinned-versus-selected comparison possible (§4.7).
- **Any replication of the pair search must include the `train_ranges` filter** (§3.3, filter (c)).
  - **Its effect:** it removes 31% / 18% / 56% of boundary candidates by pair for row 6, and 56% / 51% / 49% for
    row 199.
  - **Without it:** a replication (`17_1`'s `select_pair()`) agrees with the package on 76 of 221 rows, and on 84
    when it also rotates from `bl_result$V` (§4.5).
- **`all_distances_mahalanobis` stores the squared distance.** The package's summary prints its square root.

---

## References

- **The thesis:** Rowan, A. (2026), PhD thesis, final version of 27 February 2026 (`1 FINAL_ADRIAAN ROWAN_27022026.pdf`):
  - §3.5–3.5.2 (pp.52–61): the rotation and its CVA variant;
  - §3.8.3, §3.9.1 and §3.9.3 (pp.82–93): PCA versus CVA, the distance criterion, and the procedure;
  - §4.3.1 (pp.113–118): the search over all component pairs (the Pima example);
  - §4.4 (pp.124–125): the conclusions;
  - §6.3 (p.144): the Grand Tour as future work.
- **Package source:**
  - `R/local_cf.R`: `.generate_pairs()`, `.bl_rotate()`, `bl_find_local_cf()`, `.restore_fixed()`,
    `.round_integers()`, `.project_local()`;
  - `R/projection.R`: `.compute_metric_inverse()`, `bl_build_projection()`.
- **Companion and related notes:**
  - [`../boundary-distance-and-robustness.md`](../boundary-distance-and-robustness.md): $\mathbf{W}$, its Cholesky
    inverse, and Phase 2's distance and robustness measures;
  - [`../nearest-vertex-tie-break.md`](../nearest-vertex-tie-break.md): why the nearest vertex is the exact minimum;
  - [`../fixed-constraint-restore.md`](../fixed-constraint-restore.md): the post-search "fixed" restore and "integer" rounding (§3.5).
- **Design documents:**
  - `2 implementation_summary.txt` §4.2, §4.2.1 and §4.8: code mechanics;
  - `1 Foundation intro documents.txt` §10.3, §10.4 and §10.7: the CVA scatter-matrix requirement, hull
    invalidation under rotation, and the `max_pairs` cap.
- **The research code** (branch `origin/original_PhD_code`), `scripts/`:
  - `1.3 Optimal Rotation.R` (`Biplot_rotation()`);
  - `4.3 Different Eigenvalues.R` (the pair search and `pc.distz`);
  - `1.4 Biplot Boundary Search_optimized.R` (the in-plane search);
  - `3.1 Biplot search - Train data.r` (`Wmat`);
  - `1 Functions for biplot.R` (the CVA call).
- **biplotEZ** 2.2: `CVA.biplot()` and `CVAlowdim()`, read from the installed package.
- Cox & Cox (2008), on Procrustes analysis, as cited in the thesis §3.5 (full reference in the thesis's
  bibliography).
