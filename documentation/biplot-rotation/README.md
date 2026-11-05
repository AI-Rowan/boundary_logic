# Biplot rotation and pair selection — note and worked examples

This folder is a self-contained bundle: the technical note on the `boundarylogic`
local-biplot rotation and the selection of the plane, plus the R scripts that
produced every measurement in it. The note runs from the original PhD method (§2)
through the current implementation (§3), its problems (§4) and a proposed solution
(§5) to future developments (§6). Its measurements are in Appendix A. Phase 2's
distance and robustness measures are in a separate note,
`../boundary-distance-and-robustness.md`.

**Originated 2026-09-10** in the `Explain XAI code` project's sandbox
(`R/In-development/`). This repo is now this bundle's home — edit the note and
scripts here. **Last re-measured 2026-10-30**: all four scripts run end to end
against the package on branch `paper_nr1_changes`, every assertion passes, and the
note's numbers are taken from those runs. `17_3` §11 (the CVA identity) was added
on 2026-10-31. **2026-11-05:** the package's default pair selector became
`distance = "auto"` (Euclidean in X-space under PCA, Mahalanobis under CVA). The
scripts pin `distance = "mahalanobis"`, the selector the note measured, and
re-running `17_1`-`17_3` reproduces every number in the note.

## Contents

| file | what it demonstrates | backs |
|---|---|---|
| `biplot_rotation_and_pair_selection.md` | the note itself — the original method, the code mapping, the problems, the proposal, and future developments | — |
| `17_1a_toy-2d-to-1d-rotation.R` | the rotation at minimum scale: 5 points, 2 variables → a **1D** biplot. The only setting where the rotation is fully determined | §2.3, A.4, A.8 |
| `17_1_loan-gam-pca-3d-3var-rotation.R` | 221 loans, **3** features → 3 PCs, GAM. A 3×3 PCA basis is a lossless rotation, so anything that differs between views is the rotation alone | §3.4 (PCA), §4.2, A.2–A.7, Appendix B |
| `17_2_loan-gam-pca-4var-slices.R` | 221 loans, **4** features → 3 PCs, GAM. Rotation *plus* truncation; two planes on one surface, where a 3D picture of 4D data misleads, and the perturbation measurements behind §4.1 | §4.1, §5.1, A.1, A.4, A.5 |
| `17_3_minimal-rotation-verification.R` | the acceptance harness for the §5 proposal — the closed-form minimal rotation, the label-integrity table, determinism, the plane dimension, the family search, and the cost — plus asserted checks that the package's own local search is deterministic and that under CVA its in-plane distance is the Mahalanobis distance for every pair (iris, §11) | §3.4 (CVA), **§5 throughout**, §6.1, A.2, A.9, A.10 |

## Running them

```
Rscript 17_3_minimal-rotation-verification.R
```

All four are self-contained: they read `loan_data.csv` from the installed
`boundarylogic` package, so there is no working directory to set and no data to
supply. Requirements are `boundarylogic`, `mgcv`, and — for `17_1`/`17_2` —
`plotly` and `reshape2`.

- `17_3` is **console output only** and is the one to run first. It asserts the
  geometry and reports the timings, and finishes in about a minute. If an
  assertion fails, something has drifted: each asserted value carries its
  expected number in a comment so the failure is diagnosable.
- `17_1` and `17_2` build interactive **plotly** figures and are best run a
  figure at a time in RStudio. They also run end to end through `Rscript` (set
  `options(browser = function(...) invisible(NULL))` so the auto-printed figures do
  not open a browser); their console output carries the numbers the note quotes.
- `17_1a` draws a base-R six-panel figure on load.

**Running these in the `boundary_logic` package repo** (rather than against whatever
build is currently installed): run `devtools::load_all()` from the repo root before
sourcing a script, or `devtools::install()` the branch under test first, so the
numbers describe the working tree you are actually looking at.

## Two things to know before relying on the numbers

**The note's §5 is a proposal, not implemented behaviour.** It proposes
replacing the rank-1 Procrustes construction with a closed-form minimal Givens
rotation, and replacing the eigenvector-pair sample with a search over the
family of planes containing the target. Nothing in the shipped package works that
way yet. Adoption waits on measurement (§5.6), and `17_3` is that measurement.

**Not everything in §6.1's acceptance table is verified.** Still open:

- test 5 — variance retained (adequacy, axis predictivity, RSS) under the minimal
  rotation versus the current one;
- test 7 at $p \ge 4$ — cost;
- cross-BLAS/LAPACK determinism, which matters precisely because the current
  construction's plane depends on the exact input bits (a tie-break over a
  degenerate subspace, note §4.1). Within one platform the package's results are
  reproducible (note A.2).

`17_3` provides the harness for these; filling them in is separate work.

## Why `A3 Rotated Biplot.R` is not here

`17_1a` supersedes it. It took A3's five points and panel layout, renamed the axes
to RC, and fixed A3's residual calculation — which summed every entry of
$\mathbf{X}'\mathbf{X}$ where the **trace** is meant, making the rotation appear to
*improve* global fit rather than trade global fidelity for local accuracy.
`17_1a`'s header still credits A3 as its predecessor, so the name appears in the
comments; the file itself is deliberately absent rather than shipping a superseded
script with a known-wrong RSS beside its replacement.
