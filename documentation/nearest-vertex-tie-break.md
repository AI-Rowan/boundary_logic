# Technical Note: A Deterministic Nearest Boundary Vertex

Date created: 2026-10-30
Date last edited: 2026-11-16

**Author:** Adriaan Rowan
**Package:** `boundarylogic`
**Location:** `documentation/nearest-vertex-tie-break.md`
**Companion documents (paths relative to this file):**
- `../2 implementation_summary.txt` Section 4.5 (`bl_find_boundary()`): code mechanics
- `biplot-rotation/biplot_rotation_and_pair_selection.md` Appendix A.2: the within-platform
  determinism of the local search, which this rule provides for the nearest-vertex step
- `../tests/testthat/test-boundary.R`: the tests that pin the behaviour
- `methodology-changes/methodology-changes.md` §2 and §4 ("Nearest boundary vertex"): the history
  of this step

---

## 1. The rule

All three boundary searches find their counterfactual as the contour vertex nearest to an
observation in the 2D biplot: `bl_find_boundary()` (global), `bl_surrogate()` (the surrogate
model's class for each row) and `bl_find_local_cf()` (local, in each rotated plane). They share one
helper, `.nearest_idx_block()` in `R/boundary.R`, which takes the nearest vertex with

```r
min_loc <- max.col(-d2, ties.method = "first")
```

`d2` holds the squared distances from each observation (rows) to each contour vertex (columns), so
the line returns the column of the row minimum: the exact nearest vertex, or the first one on an
exact tie.

The three searches are therefore deterministic. The same call on the same objects returns the same
counterfactual every time, whatever ran before it in the session. They use no random numbers, so a
search neither depends on nor changes the session's random-number state.

## 2. What R's default tie-break would do

`max.col()`'s default is `ties.method = "random"`, and R documents it as follows (`?max.col`):

> When `ties.method = "random"`, as per default, ties are broken at random. In this
> case, the determination of a tie assumes that the entries are probabilities: there is
> a relative tolerance of 1e-5, relative to the largest (in magnitude, omitting
> infinity) entry in the row.

Two features of that rule would matter here.

1. **The tolerance is relative to the farthest vertex.**
   - The entries are negated squared distances, so the largest magnitude belongs to the contour
     vertex farthest from the observation (among the vertices the helper compares at once).
   - Every vertex whose squared distance lies within about 1e-5 × (farthest squared distance) of
     the nearest would count as tied.
   - The band therefore widens with the extent of the contour set and with the observation's
     distance from it: an observation far from the boundary has the widest band.
2. **The draw uses the session's random-number stream** (`max.col()` draws only when a row has a
   tie). The pick would then depend on the state of that stream, which every earlier call that uses
   random numbers moves. That includes:
   - the searches themselves, since all three call this helper;
   - `set.seed()` inside `bl_prepare_data()`;
   - `bl_shapley(seed = )`;
   - any user code.

   Under that rule:
   - **Call order would matter:** a counterfactual would depend on what was called before it.
   - **Repeat calls could differ:** when a near-tie decides which component pair is closest,
     `best_pair` could change between identical calls of `bl_find_local_cf()`.
   - **Two fresh R sessions running the same script would still agree,** because
     `bl_prepare_data()` calls `set.seed()`, so the same script makes the same draws. A workflow
     through `bl_wrap_data()`, which seeds nothing, would differ from run to run.

The vertices inside the band are almost equidistant from the observation **in the plane**, but they
can be quite different counterfactuals. Each back-projects to a different feature vector, and the
model scores that vector, not the plane. Section 5 measures how different on the loan data.

## 3. The thesis's rule

The thesis selects the nearest counterfactual (§3.10, step 7, p.93) and states that the biplot
method is deterministic: it always generates the same counterfactuals from the same projection
matrices (§3.9.3, "Stability", p.92). Its code in Appendix A takes the nearest boundary point with
`which.min` (A3, pp.184–185; A6, p.206), which returns the first minimum, with no tolerance.

R documents `ties.method = "first"` as equal to `unname(apply(m, 1, which.max))`. On `-d2` the
helper therefore returns what `which.min` returns on `d2`: the package's rule is the thesis's rule,
not a change of method.

## 4. What the nearest vertex decides

The nearest vertex is each search's counterfactual, so it decides every number built on one:
- **Phase 2:** `bl_find_boundary()`'s `B_z`, `B_x`, `B_pred` and its distances `dist_z`, `dist_x`
  and `dist_mahalanobis`, and through them the per-variable importance and the robustness of
  `plot.bl_boundary()` and `bl_robustness()`;
- **the surrogate:** `bl_surrogate()`'s class for each row, and its accuracy;
- **Phase 3:** each component pair's boundary point and distance, the best pair and the
  counterfactual, and with them the Shapley values and the sparse counterfactual.

A rule that changed the vertex could move any of these. Section 5 measures how far.

## 5. What a random tie-break would change on the loan harness

The measurement patches `.nearest_idx_block()` in a driver (CLAUDE.md §5, "To measure one rule
alone"): a copy of its body with only `ties.method = "first"` changed to `"random"`. The driver runs
the loan harness, `development/03_loan_status_Boundary_Logic.R`, with the exact rule, then repeats
the harness's searches with the random rule, twice without setting a seed and once under each of
`set.seed(1)` to `set.seed(10)`, and compares every result with the exact rule's. A copy with
`"first"` left in place reproduces the exact results, and the helper is restored afterwards.

### 5.1 The local searches

| Search | Exact rule | Random rule, 12 runs |
|---|---|---|
| Target 13, constrained (Step 11) | pair (3,4), √d²_M 0.6934, prediction 0.693 | identical in all 12 |
| Target 13, unconstrained (Step 15) | pair (3,4), √d²_M 0.5554, prediction 0.542 | identical in all 12 |
| External applicant (Step 16) | pair (1,2), √d²_M 2.8241, prediction 0.541 | identical in all 12 |

Every pair's distance is identical as well, so here the random rule picks the exact nearest vertex
each time. These three targets lie close to their boundary, where the band is narrow (Section 2); a
target far from the boundary could still change.

### 5.2 The global search and the surrogate

| | v1 model (9 features) | v2 model (4 features) |
|---|---|---|
| Test rows with a different counterfactual under at least one of the 10 seeds | 51 of 4,429 (1.15%) | 128 of 4,429 (2.89%) |
| Rows that change under a single seed | 23–31 | 53–74 |
| Largest change in the plane distance | 9.6e-5 | 4.5e-5 |
| Largest move of the point in the plane | 0.0037 | 0.0027 |
| Largest change in a feature of the counterfactual | 0.36 training SDs | 0.24 training SDs |
| Largest change in the counterfactual's model score | 0.18 | 0.08 |
| Robustness: exact rule | 9203.29 | 6790.20 |
| Robustness: random rule, 10 seeds | 9203.04 to 9204.09 | 6790.06 to 6790.66 |
| Variable ranking (the pruning order) | unchanged under all 10 seeds | unchanged under all 10 seeds |
| Surrogate accuracy vs model: exact rule | 0.8337 | 0.8845 |
| Surrogate accuracy vs model: random rule, 5 seeds | 0.8327 to 0.8337 | 0.8841 to 0.8850 |

The vertices the random rule could pick are almost equidistant from the observation in the plane
(the distance changes by less than 1e-4), yet they can be far apart as counterfactuals: up to 0.36
training standard deviations in a feature and 0.18 in the model's score. The totals built from
them move only slightly, and the ranking of the variables does not change.

### 5.3 How often a query has a choice

The searches above make 88,598 nearest-vertex queries (the two global searches and the three local
ones; a query is one observation against one block of contour vertices). Of these:
- 503 (0.57%) have more than one distinct vertex inside the random rule's band;
- **none** has an exact tie between distinct vertices, so `ties.method = "first"` never falls back
  on vertex order here.

## 6. Why `ties.method = "first"`

- **It is the exact minimum.** It is `which.min` on the squared distances (Section 3), with no
  tolerance. The searches promise the nearest boundary point, and this returns it.
- **It is deterministic and leaves the random-number stream alone.** A search does not change the
  user's random-number state, and the user's random-number state does not change a search.
- **Exact ties are resolved by vertex order,** which is reproducible. Section 5.3 found none.
- **Across blocks it stays consistent.** The helper processes long contours in blocks, and keeps a
  block's minimum only if it is strictly smaller (`<`). An exact tie across blocks therefore keeps
  the earlier block, which is the same "first" rule.

Two other rules would be reproducible but would not return the nearest vertex:
- **Seeding inside the helper** would make the random pick reproducible, but the pick would still
  not be the nearest vertex, and it would reset the caller's random-number stream on every call.
- **A deterministic rule with a tolerance** (for example, the first vertex within the band) would
  be reproducible, but again not the nearest, and the band has no meaning for distances.

**Never take the nearest point on distances with a rule that breaks ties at random,** such as
`max.col()` with its default `ties.method = "random"`: it would bring back the behaviour of
Sections 2 and 5.

## 7. What the exact rule does not settle

- **The rotation's plane under-determination (rotation note §4.1).** The local plane is LAPACK's
  completion of a rank-1 problem.
  - Within one session it is bit-identical, so it is reproducible.
  - It can change with a different BLAS/LAPACK build, platform or R version, or with any numerical
    perturbation of the inputs, such as a different order of the CVA groups (methodology note
    §3.12).
  - The rotation note's §5 proposal addresses this, not the problem in this note.
- **A 2D search is still a 2D search.** The nearest vertex in the plane is not the nearest
  boundary point in the full feature space (rotation note §4.9).

## 8. Reproducing

`tests/testthat/test-boundary.R` pins the behaviour in six tests:
- the true nearest vertex inside the random rule's band, under 20 seeds;
- equality with `which.min` on every row;
- the first nearest vertex on an exact tie;
- no use of random numbers;
- identical `bl_find_local_cf()` results under any random-number state and after other searches;
- identical `bl_find_boundary()` and `bl_surrogate()` results under any random-number state.

With `ties.method = "random"` in the helper, 5 of the 6 would fail. Only the `which.min` test would
pass, since its random points have no near-ties.

The rule itself, on one row:

```r
d2 <- matrix(c(1 + 1e-4, 1, 2e4), nrow = 1)   # the nearest is column 2
max.col(-d2)                                  # 1 or 2, at random: both within 1e-5 * 2e4
max.col(-d2, ties.method = "first")           # 2
```

To measure the random rule on other data, as in Section 5, swap the helper in the namespace and
restore it afterwards:

```r
ns  <- asNamespace("boundarylogic")
old <- get(".nearest_idx_block", envir = ns)
src <- sub('ties.method = "first"', 'ties.method = "random"', deparse(old), fixed = TRUE)
new <- eval(parse(text = src)); environment(new) <- ns
unlockBinding(".nearest_idx_block", ns); assign(".nearest_idx_block", new, envir = ns)
# ... run the searches under several seeds, then restore:
assign(".nearest_idx_block", old, envir = ns); lockBinding(".nearest_idx_block", ns)
```
