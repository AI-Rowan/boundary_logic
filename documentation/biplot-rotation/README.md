# Biplot rotation and pair selection — note and worked examples

This folder is a self-contained bundle: the technical note on the `boundarylogic`
local-biplot rotation and the selection of the plane, plus the R scripts that
reproduce its measurements. The note runs from the original PhD method (§2)
through the current implementation (§3), its problems (§4) and a proposed solution
(§5) to future developments (§6). Its measurements are in Appendix A. Phase 2's
distance and robustness measures are in a separate note,
`../boundary-distance-and-robustness.md`.

This repo is the bundle's home: edit the note and the scripts here. All four
scripts run end to end against the package, every assertion passes, and the
note takes its numbers from their console output (the few it takes from elsewhere
name their source). The scripts pin two settings of `bl_find_local_cf()`, so
re-running `17_1`-`17_3` reproduces every number the note takes from them:

- `distance = "mahalanobis"`, the selector the note measures. The package's
  default, `distance = "auto"`, is Euclidean in X-space under PCA and Mahalanobis
  under CVA (note §3.4).
- `set_filters(tgt, hull = FALSE)`, the search without the convex hull that
  `bl_find_local_cf()` builds on each rotated plane by default
  (`../local-convex-hull.md`).

## Contents

| file | what it demonstrates | backs |
|---|---|---|
| `biplot_rotation_and_pair_selection.md` | the note itself — the original method, the code mapping, the problems, the proposal, and future developments | — |
| `17_1a_toy-2d-to-1d-rotation.R` | the rotation at minimum scale: 5 points, 2 variables → a **1D** biplot. The only setting where the rotation is fully determined | §2.3, A.4, A.8 |
| `17_1_loan-gam-pca-3d-3var-rotation.R` | 221 loans, **3** features → 3 PCs, GAM. A 3×3 PCA basis is a lossless rotation, so anything that differs between views is the rotation alone | §3.4 (PCA), §4.2, A.2–A.7, Appendix B |
| `17_2_loan-gam-pca-4var-slices.R` | 221 loans, **4** features → 3 PCs, GAM. Rotation *plus* truncation; two planes on one surface, where a 3D picture of 4D data misleads, and the perturbation measurements behind §4.1 | §4.1, §5.1, A.1, A.4, A.5 |
| `17_3_minimal-rotation-verification.R` | the acceptance harness for the §5 proposal — the closed-form minimal rotation, the label-integrity table, determinism, the plane dimension, the family search, and the cost — plus asserted checks that the package's own local search is deterministic and that under CVA its in-plane distance is the Mahalanobis distance for every pair (iris, §11) | §3.4 (CVA), **§5 throughout**, A.2, A.9, A.10; the acceptance tests of `../future-improvements-review.md` §3.4 |

## Running them

```
Rscript 17_3_minimal-rotation-verification.R
```

All four are self-contained: they read `loan_data.csv` from the `boundarylogic`
package they load, so there is no working directory to set and no data to
supply. Requirements are `boundarylogic` (0.2.0 or later), `pkgload` (installed
with `devtools`), `mgcv`, and — for `17_1`/`17_2` — `plotly` and `reshape2`.

- `17_3` is **console output only** and is the one to run first. It asserts the
  geometry and reports the timings, and finishes in about a minute. If an
  assertion fails, something has drifted: each asserted value carries its
  expected number in a comment so the failure is diagnosable.
- `17_1` and `17_2` build interactive **plotly** figures and are best run a
  figure at a time in RStudio. They also run end to end through `Rscript` (set
  `options(browser = function(...) invisible(NULL))` so the auto-printed figures do
  not open a browser); their console output carries the numbers the note quotes.
- `17_1a` draws a base-R six-panel figure on load.

**Which package they load.** `17_1`, `17_2` and `17_3` load the package
themselves, in their "Load boundarylogic" block:

- **Inside the `boundary_logic` repo** they load its working tree with
  `pkgload::load_all()`, so the numbers describe the code you are looking at. They
  find the repo by looking up from their own folder (when run with `Rscript` or
  `source()`) and from the working directory (an open `boundary_logic.Rproj` sets it
  to the repo root), and print `boundarylogic 0.2.0 (working tree)`.
- **Elsewhere** they use the installed package, and stop with a message unless it is
  0.2.0 or later, because they call `set_filters(hull = )`. To update the
  installed copy, run `devtools::install()` from the repo.
- **A development copy already loaded** (`devtools::load_all()`) is kept.

## Two things to know before relying on the numbers

**The note's §5 is a proposal, not implemented behaviour.** It proposes
replacing the rank-1 Procrustes construction with a closed-form minimal Givens
rotation, and replacing the component-pair sample with a search over the
family of planes containing the target. Nothing in the package works that way.
Adoption waits on measurement (§5.6), and `17_3` is that measurement.

**Not every acceptance test is complete.** The tests and their status are in
`../future-improvements-review.md` §3.4. Still open:

- test 4 at full scale — the family optimum against the current search, like for
  like (stage 2 on the winning plane) and over all rows, not two targets;
- test 5 — variance retained (adequacy, axis predictivity, RSS) under the minimal
  rotation versus the current one;
- test 7 at $p \ge 4$ — cost;
- test 2 across BLAS/LAPACK builds — determinism, which matters precisely because
  the current construction's plane depends on the exact input bits (a tie-break
  over a degenerate subspace, note §4.1). Within one platform the package's
  results are reproducible (note A.2).

`17_3` provides the harness for these; filling them in is separate work.

## Why `A3 Rotated Biplot.R` is not here

`17_1a` replaces it. It uses A3's five points and panel layout, names the axes RC,
and computes the residual as the **trace** of $\mathbf{X}'\mathbf{X}$. A3 sums every
entry of $\mathbf{X}'\mathbf{X}$ instead, which makes the rotation appear to
*improve* global fit rather than trade global fidelity for local accuracy.
`17_1a`'s header credits A3, so the name appears in its comments; the file itself
is left out so that a script with a wrong RSS does not sit beside its
replacement.
