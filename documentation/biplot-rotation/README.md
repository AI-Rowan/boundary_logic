# Biplot rotation and pair selection — note and worked examples

This folder is a self-contained bundle: the technical note on the `boundarylogic`
local-biplot rotation, plus the R scripts that produced every measurement in it.
It is meant to be read in place or copied wholesale into the `boundary_logic`
package project, where §6 of the note is a proposal for that codebase to act on.

**Originated 2026-09-10** in the `Explain XAI code` project's sandbox
(`R/In-development/`). This repo is now this bundle's home — edit the note and
scripts here.

## Contents

| file | what it demonstrates | backs |
|---|---|---|
| `biplot_rotation_and_pair_selection.md` | the note itself — the maths, the code mapping, the limitations, and the recommendation | — |
| `17_1a_toy-2d-to-1d-rotation.R` | the rotation at minimum scale: 5 points, 2 variables → a **1D** biplot. The only setting where the rotation is fully determined | §2.5, §5.1, §5.8 |
| `17_1_loan-gam-pca-3d-3var-rotation.R` | 221 loans, **3** features → 3 PCs, GAM. A 3×3 PCA basis is a lossless rotation, so anything that differs between views is the rotation alone | §5.1–§5.7 |
| `17_2_loan-gam-pca-4var-slices.R` | 221 loans, **4** features → 3 PCs, GAM. Rotation *plus* truncation; two planes on one surface, and where a 3D picture of 4D data misleads | §5.1, §5.4, §6.1 |
| `17_3_minimal-rotation-verification.R` | the acceptance harness for the §6 recommendation — the closed-form minimal rotation, the label-integrity table, determinism, the plane-dimension correction, the family search, and the cost | **§6 throughout** |

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
  figure at a time in RStudio rather than through `Rscript`.
- `17_1a` draws a base-R six-panel figure on load.

**Running these in the `boundary_logic` package repo** (rather than against whatever
build is currently installed): run `devtools::load_all()` from the repo root before
sourcing a script, or `devtools::install()` the branch under test first, so the
numbers describe the working tree you are actually looking at.

## Two things to know before relying on the numbers

**The note's §6 is a recommendation, not implemented behaviour.** It proposes
replacing the rank-1 Procrustes construction with a closed-form minimal Givens
rotation, and replacing the six-eigenvector-pair sample with a search over the
family of planes containing the target. Nothing in the shipped package works that
way yet. §4.1 declined the change pending "prototyping and measurement" —
`17_3` is that measurement.

**Not everything in §6.5's acceptance table is verified.** Still open:

- test 5 — variance retained (adequacy, axis predictivity, RSS) under the minimal
  rotation versus the current one;
- test 7 at $p \ge 4$ — cost;
- cross-BLAS/LAPACK determinism, which matters precisely because the current
  construction's plane depends on a tie-break over a degenerate subspace.

`17_3` provides the harness for these; filling them in is separate work.

## Why `A3 Rotated Biplot.R` is not here

`17_1a` supersedes it. It took A3's five points and panel layout, renamed the axes
to RC, and fixed A3's residual calculation — which summed every entry of
$\mathbf{X}'\mathbf{X}$ where the **trace** is meant, making the rotation appear to
*improve* global fit rather than trade global fidelity for local accuracy.
`17_1a`'s header still credits A3 as its predecessor, so the name appears in the
comments; the file itself is deliberately absent rather than shipping a superseded
script with a known-wrong RSS beside its replacement.
