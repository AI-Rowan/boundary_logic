# Technical Note: A Deterministic Nearest Boundary Vertex

**Author:** Adriaan Rowan
**Date:** 2026-10-30
**Package:** `boundarylogic`
**Branch:** `paper_nr1_changes`
**Location:** `documentation/nearest-vertex-tie-break.md`
**Companion documents (paths relative to this file):**
- `../2 implementation_summary.txt` Section 4.5 (`bl_find_boundary()`): code mechanics
- `biplot-rotation/biplot_rotation_and_pair_selection.md` Appendix A.2: the
  reproducibility findings this note explains
- `../tests/testthat/test-boundary.R`: the tests that pin the behaviour
- `../.claude/reference/nearest-vertex-tie-break-plan.md`: the approved plan that was executed

---

## 1. What changed

All three boundary searches find their counterfactual as the contour vertex nearest
to an observation in the 2D biplot: `bl_find_boundary()` (global), `bl_surrogate()` and
`bl_find_local_cf()` (local, in each rotated plane). They share one helper,
`.nearest_idx_block()` in `R/boundary.R`, and one line in it changed:

| | Line | Vertex returned |
|---|---|---|
| Before | `min_loc <- max.col(-d2)` | A random draw from every vertex "tied" with the nearest, within a tolerance (Section 2) |
| After | `min_loc <- max.col(-d2, ties.method = "first")` | The exact nearest vertex, or the first one on an exact tie |

With the change, the three searches are deterministic. The same call on the same
objects returns the same counterfactual every time, whatever ran before it in the
session.

## 2. The mechanism

`d2` holds squared distances from each observation (rows) to each contour vertex
(columns), so `max.col(-d2)` was meant to return the column of the row minimum. The
PhD code that the helper was ported from comments the line "index of min per row". But
`max.col()`'s default is `ties.method = "random"`, and R documents it as follows
(`?max.col`):

> When `ties.method = "random"`, as per default, ties are broken at random. In this
> case, the determination of a tie assumes that the entries are probabilities: there is
> a relative tolerance of 1e-5, relative to the largest (in magnitude, omitting
> infinity) entry in the row.

Two features of that rule matter here.

1. **The tolerance is relative to the farthest vertex.**
   - The entries are negated squared distances, so the largest magnitude belongs to the
     contour vertex farthest from the observation.
   - Every vertex whose squared distance lies within 1e-5 × (farthest squared distance)
     of the nearest counts as tied.
   - The band therefore widens with the size of the contour set and with the
     observation's distance from it. That is why the loan harness's external applicant,
     about 116 units from the boundary, was the local result most affected (Section 5).
2. **The draw uses the session's random-number stream.** The pick therefore depends on
   the state of that stream, which every earlier call that uses random numbers moves.
   That includes:
   - the searches themselves, since all three call this helper;
   - `set.seed()` inside `bl_prepare_data()`;
   - `bl_shapley(seed = )`;
   - any user code.

   This explains both reproducibility findings the rotation note once recorded (it now
   records the deterministic behaviour, in its Appendix A.2):
   - **Call order:** the reported distance depended on what was called before.
   - **Repeat calls:** `best_pair` flipped across identical calls when a near-tie decided
     which pair was closest.

It also explains why two fresh R sessions running the same script agreed with each
other. `bl_prepare_data()` calls `set.seed()`, so the same script makes the same draws.
The `bl_wrap_data()` path seeds nothing, which is why the §5.9 sequence differed from
run to run.

The two tied vertices are almost equidistant from the observation **in the plane**,
but they can be quite different counterfactuals. Each back-projects to a different
feature vector, and the model scores that vector, not the plane. In the iris example
below, the two candidates in pair (1,3) are 1.82300 from the target in the plane (equal
to five decimals), yet their Mahalanobis distances in feature space are 2.0812 and 2.0852.

## 3. How it was diagnosed

The test case is the standardised-iris pipeline from `test-boundary.R`: GLM, PCA,
`m = 60`, `b_margin = 0.01`, target 1. On it, pair (1,3) returned one of two
counterfactuals.

| Check | Result |
|---|---|
| `V`, `X_center`, model coefficients, `metric_inv` across 6 rebuilds | bit-identical |
| Rotation (`.bl_rotate()`), its SVD, the grid, model scores, contours, across 12 calls | bit-identical |
| 10 calls, no seed | 2 distinct results |
| `set.seed(42)` before each call | 1 distinct result |
| Same seed, plus one `runif(1)` before the call | result changes |
| Same seed, plus 200 `svd()` calls before the call | result unchanged |
| Tie set at the nearest-vertex step, pair (1,3) | 2 distinct vertices inside the band (band 2.66e-4 in squared distance; farthest squared distance 26.57) |
| `ties.method = "first"`, 10 calls, no seed | 1 distinct result |

- **Inputs and computation:** every stage before the nearest-vertex step is
  bit-identical across calls.
- **Floating point ruled out:** BLAS is R's 200 KB reference `Rblas.dll`, which is
  deterministic. That excludes multithreaded floating-point noise.
- **No mutation:** nothing modifies or caches anything on `bl_result` or `bl_target`.

The only state that changed between calls was the random-number stream, and the only
consumer of it on the search path was `max.col()`.

## 4. Origin

The helper is a direct port of the PhD research code. The line is identical in
`scripts/1.4 Biplot Boundary Search_optimized.R:45` on branch `original_PhD_code`, and
the PhD local workflow (`4.3 Different Eigenvalues.R`) uses the same boundary search.
So the random pick is part of the method as it was run for the thesis.
- **Which thesis numbers carry it:** any counterfactual, distance, robustness score or
  surrogate accuracy produced by these searches.
- **How large the effect is:** Section 5 gives the scale on the loan data.
- **What this means for Task E:** this is a bug fix rather than a change of method,
  since the thesis describes the *nearest* boundary point. It belongs in the
  methodology-changes note for the supervisor.

## 5. Measured impact

### 5.1 Loan harness, before and after

Each column is a full run of `development/03_loan_status_Boundary_Logic.R`, and the two
console logs were diffed. Everything not listed is identical.

| Output | Before | After |
|---|---|---|
| Phase 1 (data, model, projection, grid) | — | unchanged |
| Robustness total, v1 / v2 model | 9203.63 / 6813.76 | 9203.29 / 6813.80 |
| Per-variable totals | — | change in the 2nd decimal; **ranking unchanged**, so the pruned feature set and `bl_results_v2` are identical |
| Surrogate accuracy vs model, v1 / v2 | 0.8334 / 0.8843 | 0.8337 / 0.8845 |
| Target 13, constrained and unconstrained local CF, Shapley, sparse CF | — | unchanged |
| External applicant: CF prediction / sparse CF prediction | 0.488 / 0.427 | 0.368 / 0.347 |
| External applicant: √d²_M / best pair | 14638.44 / (1,2) | 14638.39 / (1,2) |

The "before" column is one particular sequence of draws. Section 5.2 measures the spread
over draws.

### 5.2 How much the random pick could move a result

The old rule was run under 10 seeds and compared with the exact nearest vertex.

| Search | Result |
|---|---|
| Global, v1 model, 4,429 test rows | 51 rows (1.15%) get a different counterfactual under at least one seed (23–31 per seed) |
| Global, v2 model, 4,429 test rows | 133 rows (3.0%) (59–79 per seed) |
| Size of those differences (both models) | plane distance changes by < 1e-4, the point in the plane by < 0.004, but the counterfactual's feature values by up to **0.36 training SDs** and its model score by up to **0.18** |
| Local, target 13 (constrained and unconstrained) | identical under every seed |
| Local, external applicant | 4 distinct counterfactuals over 10 seeds; CF prediction 0.368 or 0.488 |

Over 62,024 nearest-vertex queries (both global searches and the three local
calls):
- 492 (0.79%) had more than one distinct vertex inside the old tie band;
- **none** had an exact tie between distinct vertices. `ties.method = "first"` therefore
  never falls back on vertex order in practice.

### 5.3 The rotation note's findings

| Finding | Before | After |
|---|---|---|
| call order, `17_3` §8b, row 49: full search / pinned search | 1.147879 or 1.153789 / 1.208670 or 1.202898 | 1.147879 / 1.208670, every call |
| repeat calls, `17_1` set-up, row 102, five calls | `best_pair` (1,2), (1,3), (1,3), (1,2), (1,3) | (1,2) every call |

For row 102, the exact nearest vertex in pair (1,3) has a squared Mahalanobis distance
of 0.232776, which loses to (1,2)'s 0.229143. A near-tied vertex at 0.228868 sometimes
won. The rotation note had attributed the flip to the rotation's plane
under-determination (its §4.1). That attribution was wrong: `.bl_rotate()` is
bit-identical across calls.

## 6. Why `ties.method = "first"`

- **It is the exact minimum.** R documents `"first"` as equal to
  `unname(apply(m, 1, which.max))`, with no tolerance. The searches promise the nearest
  boundary point, and this returns it.
- **It is deterministic and leaves the random-number stream alone.** A search no longer
  changes the user's RNG state, and the user's RNG state no longer changes a search.
- **Exact ties are resolved by vertex order,** which is reproducible. Section 5.2 found
  none in practice.
- **Across blocks it stays consistent.** The helper processes long contours in blocks,
  and keeps a block's minimum only if it is strictly smaller (`<`). An exact tie across
  blocks therefore keeps the earlier block, which is the same "first" rule.

Rejected alternatives:
- **Seeding inside the helper.** This makes the pick reproducible, but it is still not
  the nearest vertex, and it would reset the caller's random-number stream on every call.
- **A deterministic rule with a tolerance** (for example, the first vertex within the
  band). This is reproducible, but again not the nearest, and the band has no meaning
  for distances.

## 7. What this does not fix

- **The rotation's plane under-determination (rotation note §4.1)** is unchanged. The
  local plane is still LAPACK's arbitrary completion of a rank-1 problem.
  - Within one session it is bit-identical, so it is reproducible.
  - It can change with a different BLAS/LAPACK build, platform or R version, or with any
    numerical perturbation of the inputs. A known example: reordering the CVA group
    levels moved every Phase 3 number (CLAUDE.md §6, 2026-10-30).
  - The rotation note's §5 proposal addresses this, not the problem in this note.
- **A 2D search is still a 2D search.** The nearest vertex in the plane is not the
  nearest boundary point in the full feature space (rotation note §4.9).

## 8. Records corrected

- **Rotation note:**
  - §5.3 and §5.9: "Explained and fixed" paragraphs, with the before-and-after numbers
    in Section 5.3 above;
  - §6.3: no longer claims §6.2 removes §5.3's cause;
  - §5: a dated line saying its search numbers predate the fix.
- **`17_1`, `17_2`, `17_3`:** comments and printed text that repeated the old
  attribution or the "first call on a fresh object" advice.
- **CLAUDE.md §6:** the "distance can differ between repeated calls" warning is replaced
  by the deterministic rule.
- **`review_section7_to_8.md`:** robustness figures from a fresh harness run.
- **The Task C plan:** Task D's "best-pair instability" and "call-history dependence"
  items are now marked fixed.

## 9. Reproducing

`tests/testthat/test-boundary.R` pins the behaviour:
- the true nearest vertex inside the old band, under 20 seeds;
- equality with `which.min`;
- first-on-exact-tie;
- no use of random numbers;
- identical `bl_find_local_cf()`, `bl_find_boundary()` and `bl_surrogate()` results under
  any RNG state and after other calls.

With the old line put back, 5 of its 6 tests fail. The rule itself, on one row:

```r
d2 <- matrix(c(1 + 1e-4, 1, 2e4), nrow = 1)   # the nearest is column 2
max.col(-d2)                                  # 1 or 2, at random: both within 1e-5 * 2e4
max.col(-d2, ties.method = "first")           # 2
```
