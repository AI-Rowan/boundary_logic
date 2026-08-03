# Code Walkthrough: Step 7 of `03_loan_status_Boundary_Logic.R`

## Context

This document continues from `review_section4_to_6.md`. At this point `bl_results` is the fully assembled Phase 1 anchor object and `test_pts` is the projected test-set `bl_points` object. Step 7 forms Phase 2 (first pass): finding the global boundary counterfactuals and visualising distance to boundary in one combined block. Step 8 (surrogate model) is independent and not covered here.

---

## Step 7 — `bl_find_boundary()` → `bl_bnd`

**File:** `R/boundary.R`

**Call:**
```r
bl_bnd <- bl_find_boundary(bl_results)
print(bl_bnd)
plot(bl_results, points = test_pts)   # or plot_biplotEZ(...)
```

### What `bl_find_boundary()` does

For every observation in `bl_results$test_data`, it finds the nearest point on the decision boundary (in Z-space), then back-projects it to original feature space to produce a **counterfactual** — the minimal perturbation needed to flip the model's prediction.

The function works in nine stages:

---

#### Stage 1 — Unpack and project observations to Z-space

Parameters extracted from `bl_results`:

| Parameter | Source field | Purpose |
|---|---|---|
| `V`, `tV` | `bl_results$V`, `bl_results$tV` | Full loading matrix + its inverse |
| `Vr`, `tVr` | `V[, proj_dims]`, `tV[proj_dims, ]` | Reduced 2-column/2-row versions for the biplot plane |
| `X_center`, `X_sd` | `bl_results$X_center/sd` | Training means and SDs for standardisation |
| `standardise` | `bl_results$standardise` | `FALSE` for CVA (as in this script) |
| `cutoff` | `0.5` | Decision threshold |
| `ct` | `bl_results$biplot_grid$ct` | Pre-computed contour line segments at `cutoff ± b_margin` |
| `polygon` | `bl_results$polygon` | Convex hull of training data in Z-space |
| `train_ranges` | `bl_results$train_ranges` | Per-feature `c(min, max)` from training data |

Input data defaults to `bl_results$test_data` (all rows, no `tdp` index).

**Projection formula:**
```
X_centered = X_test - X_center        # subtract training column means
Z_obs      = X_centered %*% Vr        # project to 2D Z-space  (n × 2)
```
(No division by `X_sd` because CVA uses `standardise = FALSE`.)

Each observation is also scored through the XGBoost model via `.pred_function()` to produce `pred_obs` (probabilities) and `pred_class_obs` (0/1 at cutoff 0.5).

---

#### Stage 2 — Clip contour segments to convex hull

The pre-computed contours in `ct` extend across the entire grid. Only points inside the training convex hull polygon are meaningful.

For each raw contour segment:
```
Mi_clipped = rows of Mi where sp::point.in.polygon(Mi, polygon) > 0
```
This is done via `.poly_clip()` (defined in `R/boundary.R`), which uses `sp::point.in.polygon()` to test each contour vertex against the polygon.

Result: `z_boundaries_list` — list of clipped contour matrices; `z_boundary_type` — the probability level of each contour (either `cutoff - b_margin` = 0.49 or `cutoff + b_margin` = 0.51, because `b_margin = 0.01`).

---

#### Stage 3 — Prune: keep only contours that enclose observations

For each **closed** contour segment (where first row == last row, tested via `.is_closed()`):
- Use `sp::point.in.polygon()` to check whether any observation falls inside the closed contour ring
- Keep the contour only if at least one observation lies inside it
- Open segments (arc-shaped, not closed) are always kept

If no contours survive pruning, all original contours are retained as a fallback.

**Why:** A closed contour that encloses no observations is irrelevant — no observation would ever project to it as its nearest boundary.

---

#### Stage 4 — Consistency pruning: drop boundary points where the model disagrees

Each surviving contour segment is back-projected to X-space point by point:
```
Bx = (Mi %*% tVr) + X_center          # Z → X
pred_bnd = .pred_function(model, "XGB", Bx)
```

A contour at probability level `p_level` is consistent only if back-projected X-space predictions agree:
- Contour with `p_level < 0.5` (class-0 side): keep rows where `pred_bnd < 0.5`
- Contour with `p_level >= 0.5` (class-1 side): keep rows where `pred_bnd >= 0.5`

Rows that fail this check are dropped from the contour. If an entire contour empties, it is removed.

**Why this matters:** The biplot is a 2D linearisation of a p-dimensional space. When the model is nonlinear (e.g. XGBoost), the contour extracted from the grid may include points that, when back-projected to the full feature space, are actually on the wrong side of the model's true decision surface. This pruning corrects for that.

---

#### Stage 5 — Find nearest boundary point per observation (geometry)

Helper: `.nearest_idx_block(Z_obs, Mi, block = 5000L)`

For each retained contour segment `Mi` (a matrix of k boundary vertices in Z-space):
- Chunk the n observations into blocks of 5 000 rows
- For each block, compute pairwise squared Euclidean distances to all k contour vertices: `d2 = (z_i - m_j)^2 summed over dimensions`
- Find the column minimum per row (nearest vertex) via `max.col(-d2)`

This gives, for each observation, the index of the nearest vertex on each contour — collected into `B_boundary_z` (list of n × 2 matrices, one per contour).

**Memory design:** Computing the full n × k distance matrix at once would be expensive for large n. Chunking keeps the working matrix at `5000 × k` regardless of n.

---

#### Stage 6 — Select the opposing-side boundary per observation

For a class-1 observation (predicted to default), the counterfactual is the nearest class-0 boundary (the nearest point where the model predicts non-default). For a class-0 observation, it's the nearest class-1 boundary.

**Algorithm:**
1. Classify contours: `boundary_class_1 = (z_boundary_type >= 0.5)`
2. For each observation j, compute squared distances to all contour candidates
3. Take `which.min(d2)` separately within class-1 contours and class-0 contours → `pos1`, `pos0`
4. Select the opposing side: class-1 predictions use `pos0` (nearest class-0 contour); class-0 predictions use `pos1`
5. Store selected boundary coordinate as `B_z[j, ]`

If no opposing-side contour exists for an observation, `B_z[j, ]` is set to `NA`.

---

#### Stage 7 — Back-project boundary to X-space

```
B_x = (B_z %*% tVr) + X_center       # n × p data frame
```
Column names set to `var_names`. This is the counterfactual in original feature units — "if this applicant had these feature values, the model would predict the opposite class."

---

#### Stage 8 — Post-selection feasibility check (informational only)

After selecting the nearest boundary, `train_ranges` are used to check whether the back-projected counterfactual falls within the observed range of each training feature. If any counterfactuals are out of range, a `message()` is printed (not a warning or error).

**Why post-selection (not pre-filter):** Applying `train_ranges` *before* finding the nearest boundary point would filter out boundary segments that may be the geometrically correct nearest point. The correct approach is to find the true nearest point first, then report feasibility separately.

---

#### Stage 9 — Score and measure

```
B_pred  = .pred_function(model, "XGB", B_x)          # predicted probability at counterfactual
dist_z  = sqrt(rowSums((Z_obs - B_z)^2))              # Euclidean distance in Z-space
diff_sd = (x_obs_mat - B_x_mat) / X_sd               # per-feature diff in SD units
dist_x  = sqrt(rowSums(diff_sd^2))                    # Euclidean distance in standardised X-space
```

`B_pred` values should cluster near 0.5 (the cutoff). Larger deviations indicate the back-projection landed in a region where the nonlinear model's prediction does not match the contour's probability level.

---

### Returned object — `bl_bnd` (S3 class `"bl_boundary"`)

| Field | Type | Dimensions | Content |
|---|---|---|---|
| `B_z` | numeric matrix | n × 2 | Z-space coordinates of the selected boundary point per observation. NA rows = no boundary found. |
| `B_x` | data.frame | n × p | Back-projected counterfactual in original feature space. Column names = `var_names`. |
| `B_pred` | numeric vector | n | Model probability at the counterfactual; should be ≈ 0.5. NA if prediction failed. |
| `dist_z` | numeric vector | n | Euclidean distance (obs → boundary) in Z-space. |
| `dist_x` | numeric vector | n | Euclidean distance (obs → boundary) in standardised X-space. |
| `Z_obs` | numeric matrix | n × 2 | Z-space coordinates of each input observation. |
| `x_obs` | data.frame | n × p | Original feature values of each input observation. |
| `pred_obs` | numeric vector | n | Model probability at the observation (floor-rounded to 3 d.p.). |
| `class_obs` | integer vector | n | True loan status (0/1), or NA if absent. |
| `z_boundaries_list` | list | nr_boundaries elements | Clipped + pruned contour matrices; each element is a k × 2 matrix. |
| `z_boundary_type` | numeric vector | nr_boundaries | Probability level of each surviving contour (0.49 or 0.51 here). |
| `nr_boundaries` | integer | 1 | Count of surviving contour segments. |
| `position_store` | integer matrix | n × 2 | Contour index used per observation; col 1 = class-1 contour, col 2 = class-0. |
| `bl_result` | S3 object | — | Reference to parent `bl_results` (used by downstream methods). |

---

### `print(bl_bnd)` — console summary

**File:** `R/boundary.R` (S3 method `print.bl_boundary`)

Prints to console:
```
<bl_boundary>
  Observations   : [n test rows]
  Boundaries used: [nr_boundaries]
  Dist Z (mean)  : [mean of dist_z, 4 d.p.]
  Dist X (mean)  : [mean of dist_x, 4 d.p.]
  B_pred range   : [min(B_pred)] – [max(B_pred)]
```

`B_pred range` is the diagnostic check: values should be close to 0.5. Wide range or values far from 0.5 indicate consistency issues.

---

### `plot(bl_results, points = test_pts)` (or `plot_biplotEZ(bl_results, points = test_pts)`)

`plot_biplotEZ()` no longer accepts a `boundary =` parameter. The biplot renders without
arrow overlays. To inspect boundary counterfactuals for individual observations, call
`bl_pick_point()` on the active plot after rendering:

```r
plot(bl_results, points = test_pts)   # or plot_biplotEZ(...)
bl_pick_point(bl_results, bl_boundary = bl_bnd)
```

Each click identifies the nearest training observation and draws its counterfactual on the
plot: an `×` cross at `B_z` (the counterfactual position in Z-space) and an arrow from
the observation to `B_z`. If no boundary was found for that observation (`B_z = NA`), a
message is printed. Press Escape to finish clicking.

**Why interactive instead of global:** The global overlay drew n arrows at once, causing
zero-length arrow warnings from R's graphics device (when an observation is exactly on the
boundary) and visual clutter at large n. The per-click approach fires at most one arrow per
interaction, in interactive mode where the warning is visible and actionable.

**Complete layer order:**

| Layer | Content |
|---|---|
| 1 | Grey axes + variable labels (biplotEZ) |
| 2 | 200 × 200 probability grid (blue→white→red) |
| 3 | Test data points coloured by confusion category (from `test_pts`) |
| 4 | Darker axes redrawn on top |
| 5 | Decision boundary contour lines |

**Per-click additions (via `bl_pick_point(bl_results, bl_boundary = bl_bnd)`):**

| Element | Content |
|---|---|
| Yellow circle + row label | Picked observation highlighted on the plot |
| `×` cross | Counterfactual position `B_z` in Z-space |
| Grey arrow | Direction from observation to counterfactual |

---

### `plot(bl_bnd)` — distance-to-boundary visualisation

**File:** `R/boundary_plot.R` (S3 method `plot.bl_boundary`)

Both calls render the same data in two different chart styles. The core computation is identical.

---

### Core computation (shared by both plot types)

**Step 1 — Back-project the Z-space distance vector to X-space**

The distance vector in Z-space (from observation to boundary) is mapped back to understand which features contribute most:

```
vec_to_boundary_Z  = Z_obs - B_z            # n × 2 distance vector in Z-space
vec_to_boundary    = vec_to_boundary_Z %*% tVr   # n × p — back-project to X-space
                                            # tVr is the 2 × p inverse-loading slice
```

Because CVA has `standardise = FALSE`, no `X_sd` multiplication is needed here.

Then standardise by dividing by a per-feature denominator so all features are on the same scale. The denominator depends on the `distance` argument:
```
# distance = "mahalanobis" (default, new):
denom              = sqrt(diag(bl_results$metric))    # within-class SD per feature
vec_to_boundary_sd = vec_to_boundary / denom          # n x p signed standardised distances

# distance = "euclidean" (legacy):
vec_to_boundary_sd = vec_to_boundary / X_sd           # n x p, total SD denominator
```

Each cell `[i, j]` is the **signed standardised distance** from observation `i` to its counterfactual in the direction of feature `j`. The sign indicates which side of the boundary the observation is on for that feature.

**Why within-class SD (the new default).** `X_sd` mixes within-class noise with between-class signal. For a classification problem the relevant noise model is within-class, so dividing by `sqrt(diag(W))` correctly amplifies features that are good class separators (small within-class variance, large between-class variance). Phase 2 distance plots under the new default emphasise the features that *actually matter* for the boundary, not those that simply happen to have a wide marginal spread.

A cross-correlation diagnostic is printed to the console when `distance = "mahalanobis"`, comparing the diagonal d_M^2 contribution to the full d_M^2: small cross-correlation percentage means the diagonal approximation is essentially lossless. See `2 implementation_summary.txt` Section 4.6.1 for why the full Shapley decomposition is theoretically more correct but is not used (`O(n * 2^p)` cost).

**Step 2 — Compute per-variable importance**

```
sum_of_distance = colSums(abs(vec_to_boundary_sd))   # p-vector
```

Each element is the sum of absolute standardised distances across all n observations for one feature. This is the **variable-level importance proxy**: a higher value means the population collectively has more "room to move" in that feature's direction before hitting the boundary.

**Step 3 — Sort ascending (least important at bottom of y-axis)**

```
ord       = order(sum_of_distance)       # ascending
y_labels  = paste0(var_sorted, " : ", round(sd_sorted, 0))
```

The Y-axis label format `"variable_name : total"` shows both the feature name and its aggregated distance. Variables with the smallest total appear at the bottom, the largest at the top.

**Step 4 — Confusion category colours**

Because `bl_bnd$class_obs` contains true labels (the test set has a `"class"` column):
- `TP` (predicted 1, actual 1) → red
- `TN` (predicted 0, actual 0) → blue
- `FP` (predicted 1, actual 0) → purple
- `FN` (predicted 0, actual 1) → orange

**Step 5 — Reshape to long format for ggplot**

The n × p matrix is stacked into a three-column data frame:
- `values`: signed standardised distance
- `Variable`: feature name (as factor, ordered by importance)
- `grp`: confusion category

---

### `plot(bl_bnd)` — Jitter plot (default `type = "jitter"`)

**What is drawn:**

- **Vertical line at x = 0**: represents the boundary itself (zero distance)
- **One point per observation per variable**: position on x = `vec_to_boundary_sd[i, j]`; slight vertical jitter (`height = 0.25`) prevents overplotting
- **Colour**: confusion category (TP/TN/FP/FN)
- **Y-axis**: variables sorted ascending by total distance, with total shown in label

**Reading the chart:**
- Points to the **right of 0** (positive): the observation's feature value is higher than its counterfactual's value for that feature
- Points to the **left of 0** (negative): the observation's feature value is lower than its counterfactual's
- **Spread** along the x-axis for a variable shows how varied the distance-to-boundary is across observations
- **Variables near the top** (largest total) are most important for distinguishing the boundary — the model uses those features most to draw the decision boundary

**Console output after plot:**
```
Robustness (total distance): [sum of all sum_of_distance values]

Per-variable totals (descending):
loan_int_rate       : [value]
loan_percent_income : [value]
...
```

---

### `plot(bl_bnd, type = "boxplot")` — Boxplot by confusion group

**What is drawn:**

Same axes and y-axis labels as the jitter plot. For each variable, one box per confusion group (TP, TN, FP, FN), side by side (`position_dodge(width = 0.7)`):
- **Box**: interquartile range (25th–75th percentile) of distances for that group
- **Middle line**: median
- **Whiskers**: data extremes or IQR ± 1.5 (default ggplot whisker rule)
- **Outlier dots**: observations beyond the whiskers
- **Fill colour**: confusion category (semi-transparent, `alpha = 0.4`)

**Reading the chart:**
- **FN boxes far from 0 in a feature**: false negatives (borrowers who will default but are predicted safe) differ most from the boundary in that feature — that feature would need to change most for them to be flagged
- **TP boxes close to 0**: true positives are near the boundary; a small change in features would flip their prediction
- **Separation between TP and TN boxes** shows whether the two "correct" groups sit on opposite sides of the boundary or both far from it
- Boxplots reveal distributional shape (skew, spread, outliers) that the jitter plot obscures when n is large

---

### Return value of `plot.bl_boundary()`

Returns a list invisibly:

| Field | Type | Content |
|---|---|---|
| `plot` | ggplot object | The rendered chart (can be re-displayed with `print(result$plot)`) |
| `sum_of_distance` | named numeric vector (length p) | Per-feature total absolute standardised distance |
| `robustness` | numeric scalar | Sum of `sum_of_distance` across all features |
| `vec_to_boundary_sd` | numeric matrix (n × p) | Full signed standardised distance matrix |

---

## Complete Object Flow: Step 7

```
bl_results  [bl_result]           ← Phase 1 anchor (from Steps 4-6)
test_pts    [bl_points]           ← test data projected to Z-space
     │
     │  bl_find_boundary(bl_results)
     │
     │  [internal chain]:
     │    project test_data → Z_obs (n × 2)
     │    score obs → pred_obs, pred_class_obs
     │    clip ct segments to polygon → z_boundaries_list
     │    prune closed contours by obs membership
     │    consistency pruning: back-project + re-score each segment
     │    .nearest_idx_block() → nearest vertex on each contour per obs
     │    select opposing-side contour per obs
     │    back-project B_z → B_x (n × p)
     │    score B_x → B_pred; measure dist_z, dist_x
     ▼
bl_bnd  [bl_boundary]
  ├── B_z          (n × 2, boundary in Z-space)
  ├── B_x          (n × p, counterfactual in X-space)
  ├── B_pred       (predicted probability at counterfactual, ≈ 0.5)
  ├── dist_z       (Z-space distance per obs)
  ├── dist_x       (standardised X-space distance per obs)
  ├── Z_obs        (observation positions in Z-space)
  ├── x_obs        (observation feature values)
  ├── pred_obs     (predicted probability at obs)
  ├── class_obs    (true labels)
  └── bl_result    (reference to bl_results)
     │
     ├── print(bl_bnd)
     │     └─→ console: n, boundaries, mean distances, B_pred range
     │
     ├── plot(bl_results, points=test_pts)   [or plot_biplotEZ(...)]
     │     └─→ clean biplot (no boundary overlay)
     │
     ├── bl_pick_point(bl_results, bl_boundary=bl_bnd)    [interactive]
     │     └─→ per-click: yellow circle + ×-cross at B_z + arrow from Z_obs → B_z
     │
     ├── plot(bl_bnd)                     [type = "jitter"]
     │     ├── compute vec_to_boundary_sd  (n × p signed distances)
     │     ├── sort variables by sum_of_distance
     │     ├── render jitter plot (one dot per obs per variable)
     │     └─→ console: robustness scalar + per-variable totals
     │
     └── plot(bl_bnd, type = "boxplot")   [type = "boxplot"]
           ├── same vec_to_boundary_sd computation
           └─→ box-and-whisker plot per variable × confusion group
```

---

## Key Methodological Points

1. **Z-space search, X-space result.** The nearest boundary search is entirely in 2D Z-space (fast, geometric). The result is then back-projected to p-dimensional X-space to produce an actionable counterfactual. This two-step design is what makes the method scalable.

2. **Consistency pruning is the accuracy guard.** Because the biplot is a linear 2D projection of a nonlinear model's decision surface, contour segments can be misleading. The back-project-and-re-score step ensures that only contour vertices that are genuinely near the true decision boundary (in X-space) are used.

3. **Opposing-side selection defines the counterfactual direction.** A class-1 observation's counterfactual must lie on the class-0 boundary — the nearest point where the prediction flips. This is what makes `dist_x` a meaningful measure of "how far is this applicant from becoming safe?"

4. **`sum_of_distance` is a variable importance proxy for the boundary.** Unlike model-level feature importance (e.g. XGBoost's `gain`), `sum_of_distance` tells you which features matter for reaching the boundary from where the population currently sits — a local, population-anchored importance measure.

5. **`bl_result` stored inside `bl_bnd`.** Every downstream method (`plot.bl_boundary`, `bl_robustness`) re-extracts `V`, `tV`, `X_sd`, etc. from `bl_bnd$bl_result`. This means the boundary object is self-contained and fully reproducible without requiring `bl_results` to be in scope.
