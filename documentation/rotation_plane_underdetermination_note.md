# Technical note: the local biplot rotation plane is under-determined

**Status:** Finding recorded. **No code changed.** Decision deferred to the author.
**Date:** 2026-07-15
**Applies to:** `.bl_rotate()` (`R/local_cf.R:45`), consumed by `bl_find_local_cf()`. Phase 3 only.
**Origin:** inherited from `origin/original_PhD_code:scripts/1.3 Optimal Rotation.R`
(`Biplot_rotation()`). **The package is a faithful port; this is a property of the original
method, not a porting defect.**

---

## 1. Summary

`.bl_rotate(x_target_st, V, proj_pair)` rotates the loading matrix so the target observation
lies in the local biplot plane. It succeeds at that. But the plane it produces is **not
determined by its inputs**: only the target's direction is pinned, while the plane's *second*
direction is fixed by an arbitrary tie-break inside LAPACK's SVD.

Consequences, all measured (§4):

- Two bases differing by **4.6e-15** produce planes **65.9 degrees** apart.
- Only **46%** of the resulting plane's second direction lies inside the `proj_pair` eigenplane
  it is named after.
- The counterfactual returned by `bl_find_local_cf()` therefore depends on a floating-point
  tie-break, and could differ across BLAS implementations, machines, or R versions.

What is **not** affected is important and is set out in §6: the target's own coordinates are
identical regardless, and the counterfactual remains a genuine boundary point.

---

## 2. What the rotation is for

Phase 3 explains a single observation by rotating the biplot so the target lies in the drawn
plane, then finding the nearest decision-boundary point *within that plane*. The original
research script states the intent in its header:

> Computes a rotation matrix A that aligns 2D biplot plane spanned by `proj`
> (typically PCs 1 & 2) with a selected target point (tdp)

So the intended object is **the `proj` eigenplane, rotated to contain the target**. Two
properties are wanted: (a) the target lies in the plane, and (b) the plane is still recognisably
the one `proj` names. The implementation delivers (a) but not (b).

## 3. The construction and why it is under-determined

The original (`1.3 Optimal Rotation.R`), reproduced line-for-line by `.bl_rotate()`:

```r
Y   <- rbind(-X.st[tdp, ], rep(0, num_vars), X.st[tdp, ])   # 3 x p
YV  <- Y %*% V
YVr <- Y %*% Vr
if (num_vars > r) YVr <- cbind(YVr, matrix(0, ncol = num_vars - r, nrow = 3))
svdyx <- svd(t(YV) %*% YVr)
A     <- svdyx$v %*% t(svdyx$u)
Vrho  <- V %*% t(A)
```

**`Y` has rank 1.** Its three rows are `-x`, `0`, `x` — all multiples of a single vector. The
zero row contributes nothing; `-x` is redundant with `x`. So despite presenting as a three-point
Procrustes problem, this is effectively a **one-vector alignment**.

Write `a = (-1, 0, 1)'`, `u = V'x`, and `w = (u[p1], u[p2], 0, ..., 0)`. Then
`Y = a x'`, `YV = a u'`, `YVr_padded = a w'`, and

```
M = t(YV) %*% YVr_padded = 2 * u w'      -- rank 1
```

Minimising the Procrustes objective over orthogonal `R` reduces to

```
|| a u' R - a w' ||_F^2  =  2 || R'u - w ||^2
```

and since `||R'u|| = ||u||` is fixed, the minimiser is characterised **solely** by

```
R'u = ||u|| * v1 ,      v1 = w / ||w||
```

That is a **single vector constraint on a p x p orthogonal matrix**. For `p = 4` the solution
set is a manifold of dimension 3 (an `SO(3)`'s worth of rotations), *all achieving the identical
optimum* of the script's own `rho` diagnostic. The SVD of a rank-1 `M` has `p - 1` exactly
degenerate zero singular values, and the corresponding singular vectors — hence `A` — are
whatever LAPACK's bidiagonalisation happens to emit.

### 3.1 What the arbitrary part actually controls

Decompose the coordinate basis. With `v1perp` the unit vector in `span(e1, e2)` orthogonal to
`v1`, and `U2`, `V2` the (arbitrary) degenerate completions:

```
R e1 = a1 * u1 + b1 * q ,   R e2 = a2 * u1 + b2 * q ,   q = U2 V2' v1perp
```

where `a_i = v1'e_i` and `b_i = v1perp'e_i`. Therefore

```
span(R e1, R e2) = span(u1, q)
```

The plane always contains the target direction `u1` — that part is genuine and reliable. Its
second direction `q` is a function of the arbitrary completion and **is not determined by the
maths**.

## 4. Measurements

Loan data, GAM model, `p = 4`, target row 49, `proj_pair = c(2, 3)` (the pair
`bl_find_local_cf()` actually selects for this target).

**The SVD is exactly degenerate.** Singular values of `M`:

```
1.6274e+00,  6.7123e-17,  0.0000e+00,  0.0000e+00      -> rank 1 of 4, null space dim 3
```

**The plane is discontinuous in its inputs.** Perturb `V`, re-orthonormalise, re-rotate, and
measure principal angles against the unperturbed plane (first angle is always 0.00 — the target
direction is always shared):

| perturbation to `V` | max abs dV | principal angles (deg) |
|---|---|---|
| re-orthonormalise only | 4.44e-16 | 0.00, **56.79** |
| 1e-16 | 4.44e-16 | 0.00, **12.16** |
| 1e-15 | 1.44e-15 | 0.00, 0.00 |
| 1e-12 | 1.11e-12 | 0.00, **35.26** |
| 1e-08 | 2.24e-08 | 0.00, **11.49** |
| `svd(Xstd)$v` vs biplotEZ `Lmat` (a real, incidental difference) | **4.65e-15** | 0.00, **65.91** |

The last row is not hypothetical: it is the difference between a hand-rolled `svd(Xstd)$v` and
the package's `bp$Lmat` — numerically the same matrix — yet the planes end up 66 degrees apart.

**The plane is not the one `proj_pair` names.** Fraction of each plane's second direction lying
inside `span(V[,2], V[,3])`, and the principal angles against the canonical minimal rotation of
that eigenplane onto the target:

| plane | fraction inside `span(V2,V3)` | angles vs canonical | SS retained |
|---|---|---|---|
| canonical (eigenplane rotated minimally onto the target) | 1.000 | -- | 0.5489 |
| **package (`.bl_rotate`)** | **0.461** | **0.00, 81.17** | 0.5743 |
| hand-rolled `svd(Xstd)$v` | 0.985 | 0.00, 29.58 | 0.5585 |

**The target's coordinates are identical regardless** — which is exactly why this is easy to
miss:

```
package    target coords = (0.764192, 1.346277, -1.75e-16, -1.27e-16)
hand-rolled target coords = (0.764192, 1.346277,  1.39e-17, -2.15e-16)
```

Same point, to 6 decimal places, in two planes 66 degrees apart.

## 5. Consequences

1. **Reproducibility.** The returned counterfactual depends on a LAPACK tie-break over an
   exactly degenerate subspace. A different BLAS (OpenBLAS vs MKL vs reference), platform, or
   LAPACK version may legitimately return a different plane and therefore a different
   counterfactual, from identical data and identical code. For a method intended to support
   published, auditable explanations this is the material risk.
2. **`proj_pair` is a weaker label than it appears.** `bl_find_local_cf()` enumerates six pairs
   and reports a winner (`best_pair`), but each candidate plane is only loosely tied to the pair
   that names it (46% here). The pair search is better described as sampling six planes through
   the target than as enumerating the six eigenplanes.
3. **The counterfactual is a local, not a global, minimum.** It is the nearest boundary point
   within one arbitrary plane through the target, not the nearest within the pair's eigenplane
   nor in the full space.
4. **Downstream artefacts inherit the plane**: `bl_shapley()` and `bl_find_sparse_cf()` build on
   `bl_local`, so they inherit whichever plane was drawn.

## 6. What is NOT affected

Deliberately recorded, to keep the finding in proportion:

- **The target's position is exact and stable** (§4). Any statement about where the target sits
  is unaffected.
- **The counterfactual is still a genuine boundary point.** `bl_find_local_cf()` re-scores every
  back-projected candidate through the model (filter 4, `R/local_cf.R:531-538`) and keeps only
  vertices that actually cross the cutoff. An arbitrary plane yields a *valid* counterfactual,
  just not a canonical or reproducible one.
- **The cross-pair selector is basis-independent.** `dist_mah = v' metric_inv v` with
  `v = B_x - x_obs` is computed in X-space (`R/local_cf.R:565-570`), so the comparison between
  pairs is measured in a frame the rotation cannot distort.
- **Phases 1 and 2 are untouched.** `bl_build_projection()`, `bl_build_grid()`,
  `bl_find_boundary()`, `bl_robustness()` and `bl_surrogate()` never call `.bl_rotate()`.
- **`bl_set_scaling()` is display-only** and provably cannot reach `V`, the grid, or the plane.

## 7. Option considered but NOT actioned

**Explicit minimal (Givens) rotation.** Replace the rank-deficient Procrustes with the rotation
acting in `span(x, x_P)` that carries the target's shadow `x_P` onto the target `x`, where
`x_P = P P' x` and `P = V[, proj_pair]`. Because the direction `q*` in `P` orthogonal to `x_P`
is orthogonal to `x` as well, it is fixed by that rotation, giving

```
R_min(P) = span(x_hat, q*)
```

— the eigenplane rigidly rotated onto the target. This is deterministic, reproducible, distorts
the plane least, and makes `proj_pair` mean what it says. It has a genuine degenerate case the
current construction lacks: when the target is orthogonal to `P` (`||x_P|| = 0`) the minimal
rotation is undefined and should raise, rather than silently return an arbitrary plane.

**Not done, by decision (2026-07-15):** it would move every Phase 3 counterfactual, so it is not
covered by the "outputs proven equivalent" rule in CLAUDE.md §4, and the current behaviour is
inherited from the validated PhD method. Any future attempt should first prototype outside the
package and measure: reproducibility under perturbation, SS retained, and whether counterfactual
Mahalanobis distances improve or worsen across a range of targets.

## 8. Reproducing these numbers

Not committed as a script (the analysis was one-off). To reconstruct: build a `bl_result` on the
loan data exactly as `scripts/18_loan-glm-4var-pcabiplot.R` does, then for target row 49 and
`proj_pair = c(2, 3)` compare `bl_local$Vr_rot` against `.bl_rotate()` applied to
`svd(scale(train))$v`, using principal angles

```r
prin_angles <- function(A, B) {
  qa <- qr.Q(qr(A)); qb <- qr.Q(qr(B))
  acos(pmin(1, pmax(-1, svd(t(qa) %*% qb)$d))) * 180 / pi
}
```

and inspect `svd(t(YV) %*% YVr_padded)$d` to see the rank-1 degeneracy directly. Per CLAUDE.md
§5, write such probes to a `.R` file and run `Rscript file.R` — never `Rscript -e`.

## 9. Related

- `.claude/plans/compare_3d_vs_2d_biplot_target49.md` — the investigation that surfaced this,
  including the curve-derivation findings for `scripts/17_loan-gam-pca-3d-4var-surface.R`.
- `2 implementation_summary.txt` §4.2 — describes `.bl_rotate()`; its phrase "aligns the target
  with the chosen projection plane" is true of the target but not of the plane.
- `documentation/mahalanobis_technical_note.md` — the `W` metric used by the (basis-independent)
  cross-pair selector.
