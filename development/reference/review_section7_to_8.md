# Code Walkthrough: Step 7 and variable pruning of `03_loan_status_Boundary_Logic.R`

## Context

This document continues from `review_section4_to_6.md`. At this point `bl_results` is the fully assembled Phase 1 anchor object and `test_pts` is the projected test-set `bl_points` object. Step 7 forms Phase 2 (first pass): finding the global boundary counterfactuals and visualising distance to boundary in one combined block.

The variable-pruning block that follows Step 7 uses its per-variable distances to choose the four features of the reduced model `bl_results_v2`, which `review_section9_to_15.md` uses throughout. Steps 8 and 8b (the surrogate models for the full and the reduced model) run together after Step 7b. They are not used directly in the workflow and are not covered here.

---

## Step 7 — `bl_find_boundary()` → `bl_bnd`

**File:** `R/boundary.R`

**Call:**
```r
bl_bnd <- bl_find_boundary(bl_results)
hist(bl_bnd$B_pred)                  # diagnostic: B_pred should cluster near 0.5

plot(bl_results, points = test_pts)  # or plot_biplotEZ(...)

plot(bl_bnd)
plot(bl_bnd, type = "boxplot")
```

`print(bl_bnd)` is not called in Step 7; the reduced-model pass (Step 7b) prints it. In the harness run `bl_find_boundary()` also reports: "27 identified counterfactual(s) are outside the training feature ranges, but are retained." For the reduced model of Step 7b the same message reports 97.

**Parameters (`bl_find_boundary()`):**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `data` | `NULL` | default | Rows to search; `NULL` uses `bl_results$test_data`, and `train_data` gives a training-data search. `bl_pick_point()` defaults to `train_data`, so it needs `data = bl_results$test_data` to show these counterfactuals |
| `tdp` | `NULL` | default | Row indices of `data` to process; `NULL` processes every row |

Full argument list: `?bl_find_boundary`.

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
| `metric_inv` | `bl_results$metric_inv` | The inverse within-class metric `W^-1`, for `dist_mahalanobis` (Stage 9) |

Input data defaults to `bl_results$test_data` (all rows, no `tdp` index).

**Projection formula:**
```
X_centered = X_test - X_center        # subtract training column means
Z_obs      = X_centered %*% Vr        # project to 2D Z-space  (n × 2)
```
(No division by `X_sd` because CVA uses `standardise = FALSE`.)

Each observation is also scored through the XGBoost model via `.pred_function()` (model type `"custom"`, which calls the stored `predict_fn`) to produce `pred_obs` (probabilities) and `pred_class_obs` (0/1 at cutoff 0.5).

---

#### Stage 2 — Clip contour segments to convex hull

The pre-computed contours in `ct` extend across the entire grid. Only points inside the training convex hull polygon are meaningful.

**The polygon** is `bl_results$polygon`: the convex hull of the projected training data in the final biplot plane.
- **Where it is built:** `.build_hull_polygon()` (`R/hull_utils.R`) builds it inside `bl_build_grid()` (Section 6.2 of `review_section4_to_6.md`).
- **Why clip to it:** a boundary point outside the hull lies where the model has seen no training data, so a counterfactual there would rest on extrapolation.
- **Graphics device:** `.build_hull_polygon()` computes the hull by drawing it with `aplpack::plothulls()` on a temporary off-screen `png()`. It closes only that device and re-activates the caller's device, so building the result leaves the user's active plot where it was.

For each raw contour segment:
```
Mi_clipped = rows of Mi where sp::point.in.polygon(Mi, polygon) > 0
```
This is done via `.poly_clip()` (defined in `R/hull_utils.R`), which uses `sp::point.in.polygon()` to test each contour vertex against the polygon.

Result: `z_boundaries_list` — list of clipped contour matrices; `z_boundary_type` — the probability level of each contour (either `cutoff - b_margin` = 0.499 or `cutoff + b_margin` = 0.501, because `b_margin = 0.001`).

---

#### Stage 3 — Prune: keep only contours that enclose training data

For each **closed** contour segment (where first row == last row after clipping, tested via `.is_closed()`), `.encloses_data()` (`R/hull_utils.R`):
- projects `bl_result$train_data` onto the plane, as the search projects `data`;
- uses `sp::point.in.polygon()` to check whether any **training** row falls inside or on the closed contour ring;
- keeps the contour only if at least one training row lies inside it.

Open segments (arc-shaped, not closed) are always kept. If no contours survive pruning, all original contours are retained as a fallback, with a warning.

**Why:** A closed contour that encloses no training data marks a region the data do not support: the thesis discards "any local decision boundary area that does not contain data points" (§3.10 step 4, p.93). The training data decide, not the rows being searched, so a row's counterfactual does not depend on which other rows are searched with it (`tdp`, `data`); `bl_find_local_cf()` applies the same rule on each rotated plane. In this harness the rule keeps 18 contours on the reduced model (Step 7b) and 2 on the full model (Step 7).

---

#### Stage 4 — Consistency pruning: drop boundary points where the model disagrees

Each surviving contour segment is back-projected to X-space point by point:
```
Bx = (Mi %*% tVr) + X_center          # Z → X
pred_bnd = .pred_function(model, model_type, Bx)   # model_type "custom" here
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
- Find the column minimum per row (nearest vertex) via `max.col(-d2, ties.method = "first")`: the exact nearest vertex, the first one on an exact tie. The result is deterministic. R's default, `"random"`, would treat every vertex within 1e-5 of the farthest vertex's squared distance as tied and draw one at random (see `documentation/nearest-vertex-tie-break.md`).

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

**Why post-selection (not pre-filter):** the counterfactual is a point on the biplot plane mapped back to the p features through the two plotted components (`B_z %*% tVr`), so every component that is not plotted is set to its mean (a rank-2 reconstruction). Its feature values can fall outside the observed ranges even when the point lies well inside the convex hull, most often under PCA (the case that prompted the change, `ea6a12f`) but also under CVA. Applying `train_ranges` *before* finding the nearest boundary point would then discard the geometrically correct nearest point and select a farther one, changing the distances the variable importance and robustness aggregate. So the search finds the true nearest point first and reports feasibility separately, as a count only (no per-row flag). The thesis text filters by the ranges (§3.10 step 4, p.93); the departure is recorded for the supervisor in `documentation/methodology-changes/` §3.3.

---

#### Stage 9 — Score and measure

```
B_pred  = .pred_function(model, model_type, B_x)     # predicted probability at counterfactual
dist_z  = sqrt(rowSums((Z_obs - B_z)^2))              # Euclidean distance in Z-space
diff_sd = (x_obs_mat - B_x_mat) / X_sd               # per-feature diff in SD units
dist_x  = sqrt(rowSums(diff_sd^2))                    # Euclidean distance in standardised X-space
v       = B_x_mat - x_obs_mat                         # the move to the boundary, model units
dist_mahalanobis = rowSums((v %*% metric_inv) * v)    # squared Mahalanobis distance, v' W^-1 v
```

`dist_mahalanobis` is squared, as the field of the same name in `bl_find_local_cf()` is; `sqrt()` gives the distance. It uses `W`, the within-class scatter of the features, because it measures an X-space move; `documentation/boundary-distance-and-robustness.md` explains why no rotation enters it.

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
| `dist_mahalanobis` | numeric vector | n | Squared Mahalanobis distance (obs → boundary) in X-space, `v' W^-1 v`. NA where no boundary was found, or (with a warning) when `bl_result$metric_inv` is missing. |
| `Z_obs` | numeric matrix | n × 2 | Z-space coordinates of each input observation. |
| `x_obs` | data.frame | n × p | Original feature values of each input observation. |
| `pred_obs` | numeric vector | n | Model probability at the observation (floor-rounded to 3 d.p.). |
| `class_obs` | integer vector | n | True loan status (0/1), or NA if absent. |
| `z_boundaries_list` | list | nr_boundaries elements | Clipped + pruned contour matrices; each element is a k × 2 matrix. |
| `z_boundary_type` | numeric vector | nr_boundaries | Probability level of each surviving contour (0.499 or 0.501 here). |
| `nr_boundaries` | integer | 1 | Count of surviving contour segments. |
| `position_store` | integer matrix | n × 2 | Contour index used per observation; col 1 = class-1 contour, col 2 = class-0. |
| `bl_result` | S3 object | — | Reference to parent `bl_results` (used by downstream methods). |

---

### `print(bl_bnd)` — console summary

**File:** `R/boundary.R` (S3 method `print.bl_boundary`)

Prints to console. Step 7b's `print(bl_bnd_v2)` in the harness run:
```
<bl_boundary>
  Observations   : 4429
  Boundaries used: 18
  Dist Z (mean)  : 0.0095
  Dist X (mean)  : 1.8120
  Dist M (mean)  : 1.9545 (sqrt of dist_mahalanobis)
  B_pred range   : [0.2570, 0.8490]
```

The lines are the number of rows, the surviving contour segments, the means of `dist_z`, `dist_x` and `sqrt(dist_mahalanobis)`, and the range of `B_pred`, each to 4 d.p. The `Dist M` line is left out when `dist_mahalanobis` is all NA.

`Dist X` and `Dist M` are reported only. They measure from the observation itself to its counterfactual, so they include the part of the observation off the biplot plane, which the counterfactual resets to the mean; the variable importance and robustness of `plot(bl_bnd)` use only the in-plane move (`documentation/boundary-distance-and-robustness.md` §2).

`B_pred range` is the diagnostic check: values should be close to 0.5. Wide range or values far from 0.5 indicate consistency issues.

---

### `plot(bl_results, points = test_pts)` (or `plot_biplotEZ(bl_results, points = test_pts)`)

`plot_biplotEZ()` draws no boundary arrows. To inspect boundary counterfactuals for
individual observations, call `bl_pick_point()` on the active plot after rendering:

```r
plot(bl_results, points = test_pts)   # or plot_biplotEZ(...)
bl_pick_point(bl_results, data = bl_results$test_data, bl_boundary = bl_bnd)
```

`bl_bnd` holds one row per `bl_results$test_data` row (the `bl_find_boundary()` default), so
`bl_pick_point()` must search the same rows: `data = bl_results$test_data`. Its own default is
`train_data`, which would pair clicks with the wrong counterfactuals (it warns when the row
counts differ). Each click identifies the nearest test observation and draws its
counterfactual on the plot: an `×` cross at `B_z` (the counterfactual position in Z-space) and
an arrow from the observation to `B_z`. If no boundary was found for that observation
(`B_z = NA`), a message is printed. Press Escape to finish clicking.

**Why interactive, one observation at a time:** drawing all n arrows at once clutters the plot
at large n and raises a zero-length arrow warning for every observation that lies exactly on
the boundary. Drawing one arrow per click keeps the plot readable and confines any such warning
to the observation being inspected.

**Complete layer order:**

| Layer | Content |
|---|---|
| 1 | Grey axes + variable labels (biplotEZ) |
| 2 | Probability grid (blue→white→red), the 200 × 200 grid trimmed to the hull |
| 3 | Test data points coloured by confusion category (from `test_pts`) |
| 4 | Darker axes redrawn on top |
| 5 | Decision boundary contour lines |

**Per-click additions (via `bl_pick_point(bl_results, data = bl_results$test_data, bl_boundary = bl_bnd)`):**

| Element | Content |
|---|---|
| Yellow circle + row label | Picked observation highlighted on the plot |
| `×` cross | Counterfactual position `B_z` in Z-space |
| Grey arrow | Direction from observation to counterfactual |

---

### `plot(bl_bnd)` — distance-to-boundary visualisation

**File:** `R/boundary_plot.R` (S3 method `plot.bl_boundary`)

Both calls render the same data in two different chart styles. The core computation is identical.

**Parameters:**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `type` | `"jitter"` | `"jitter"`, then `"boxplot"` | Chart style; the computation below is the same for both |
| `distance` | `"auto"` | default | Per-feature denominator, by the thesis rule: under CVA (as here) `"mahalanobis"`, `sqrt(diag(W))` (within-class SD), which also prints the Mahalanobis check; under PCA `"euclidean"`, `X_sd` (total SD) |

Full argument list: `?plot.bl_boundary`.

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
# distance = "mahalanobis" (the default under CVA, used here):
denom              = sqrt(diag(bl_results$metric))    # within-class SD per feature
vec_to_boundary_sd = vec_to_boundary / denom          # n x p signed standardised distances

# distance = "euclidean" (the default under PCA; the thesis's total SD):
vec_to_boundary_sd = vec_to_boundary / X_sd           # n x p, total SD denominator
```

Each cell `[i, j]` is the **signed standardised distance** from observation `i` to its counterfactual in the direction of feature `j`. The sign indicates which side of the boundary the observation is on for that feature.

**Why within-class SD (the CVA default).** `X_sd` mixes within-class noise with between-class signal. For a classification problem the relevant noise model is within-class, so dividing by `sqrt(diag(W))` correctly amplifies features that are good class separators (small within-class variance, large between-class variance). Phase 2 distance plots under the default emphasise the features that *actually matter* for the boundary, not those that simply happen to have a wide marginal spread.

A Mahalanobis check is printed to the console when `distance = "mahalanobis"`. It compares the squared length of the plotted distances, `sum(vec_to_boundary_sd^2)` (that is, the sum of v_j^2 / W_jj), with the full squared Mahalanobis distance v' W^-1 v, both summed across observations. Their difference is what the per-feature plot leaves out by ignoring the correlations in W: positive means the plot understates the distance, negative that it overstates it. See `2 implementation_summary.txt` Section 4.6.1 for why the full Shapley decomposition is theoretically more correct but is not used (`O(n * 2^p)` cost).

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

**Console output after plot** (each `plot(bl_bnd)` call prints this; from the harness run):
```
Distance metric: mahalanobis
Robustness (total distance): 9203.29

Per-variable totals (descending):
             loan_int_rate        loan_percent_income 
                   3929.71                    1570.78 
                 loan_amnt              person_income 
                   1301.68                     777.71 
              credit_score                 person_age 
                    556.06                     390.66 
            person_emp_exp cb_person_cred_hist_length 
                    323.83                     305.42 
             person_gender 
                     47.43 

Mahalanobis check (squared distances, summed across obs):
  plotted (per feature) = 9675.25 | full d_M^2 = 11754.70 | difference = 2079.45 (+17.7% of full)
```

The check appears only for `distance = "mahalanobis"`. A note is added only when the difference is at least 25% of the full distance in absolute value. Here it is +17.7% of the full distance: the correlations in W make the full squared distance 21.5% larger than the plotted per-feature one (11754.70 against 9675.25), so the per-feature importance is a fair summary. For the four-feature model of Step 7b the check gives `plotted (per feature) = 7968.57 | full d_M^2 = 9486.73 | difference = 1518.16 (+16.0% of full)`.

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
| `distance` | character | The denominator used, after `"auto"` is resolved: `"mahalanobis"` (within-class SD) or `"euclidean"` (total SD) |
| `mahalanobis_check` | list or `NULL` | `plotted`, `full`, `difference`, `pct`: the console check's numbers; `NULL` under `"euclidean"` |

---

## Variable pruning — `bl_robustness()` → `feature_cols_v2`

**File:** `R/boundary_plot.R` (`bl_robustness()`)

The block after Step 7 uses the global counterfactuals in `bl_bnd` to choose the features of a smaller model. Its output, `feature_cols_v2`, drives the refit (Steps 4b–6b), and the refit's `bl_results_v2` is the anchor for every step in `review_section9_to_15.md`.

**Call:**
```r
rob     <- bl_robustness(bl_bnd)
var_imp <- sort(rob$sum_of_distance)     # ascending: weakest variable first
print(round(var_imp, 2))

cat(sprintf("Overall robustness (not used here): %.2f\n", rob$robustness))

vars_to_drop <- c("person_gender", "person_emp_exp", "cb_person_cred_hist_length",
                  "person_income", "person_age")
feature_cols_v2 <- setdiff(bl_dat$var_names, vars_to_drop)
cat("Retained features:", paste(feature_cols_v2, collapse = ", "), "\n")
```

**Parameters:**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `distance` | `"auto"` | default | Resolves to `"mahalanobis"` under CVA: per-feature denominator `sqrt(diag(W))`, the same as `plot(bl_bnd)`, so the totals match Step 7's console output |

Full argument list: `?bl_robustness`.

**What it does:** `bl_robustness()` repeats the per-variable part of the `plot(bl_bnd)` computation without drawing a plot (see "Core computation", Steps 1 and 2, above):
1. Back-projects each observation's move to its boundary point, `Z_obs - B_z`, into X-space with `tVr`.
2. Divides each feature by its within-class SD, `sqrt(diag(W))`.
3. Sums the absolute values over the 4,429 test observations, giving one total per feature: `sum_of_distance`.
4. Adds the per-feature totals into one overall score: `robustness`.

`plot(bl_bnd)` prints the same totals and returns them invisibly; `bl_robustness()` returns them as a value without a plot.

**Returned object — `rob` (plain list):**

| Field | Type | Content | Used here |
|---|---|---|---|
| `sum_of_distance` | named numeric (length 9) | Per-variable totals: the components | Yes: ranks the features for pruning |
| `robustness` | numeric scalar | Overall score, `sum(sum_of_distance)` | No: printed for reference only |
| `distance` | character | `"mahalanobis"` | No |

**Console output** (from the harness run):
```
             person_gender cb_person_cred_hist_length 
                     47.43                     305.42 
            person_emp_exp                 person_age 
                    323.83                     390.66 
              credit_score              person_income 
                    556.06                     777.71 
                 loan_amnt        loan_percent_income 
                   1301.68                    1570.78 
             loan_int_rate 
                   3929.71 
Overall robustness (not used here): 9203.29
Retained features: loan_amnt, loan_int_rate, loan_percent_income, credit_score 
```

The overall score equals the "Robustness (total distance)" line that `plot(bl_bnd)` prints in Step 7. It is one number for the whole test population, so it cannot rank features, and this example does not use it.

**How the variables are chosen.** The pruning uses the components, the per-variable totals extracted from `bl_bnd`. A larger total marks a feature that matters more for reaching the boundary ("Core computation", Step 2). The analyst reads `var_imp` and lists the features to drop in `vars_to_drop`:

| Feature | Total | Decision |
|---|---|---|
| `person_gender` | 47.43 | Dropped: the weakest |
| `cb_person_cred_hist_length` | 305.42 | Dropped |
| `person_emp_exp` | 323.83 | Dropped |
| `person_age` | 390.66 | Dropped |
| `credit_score` | 556.06 | **Kept:** actionable |
| `person_income` | 777.71 | **Dropped:** redundant, although it outranks `credit_score` |
| `loan_amnt` | 1301.68 | Kept |
| `loan_percent_income` | 1570.78 | Kept |
| `loan_int_rate` | 3929.71 | Kept: the strongest |

The ranking guides the choice, and two judgements override it:
- **`person_income` is redundant.** `loan_percent_income` is `loan_amnt / person_income` rounded to two decimals (it matches in 44,992 of the 45,000 raw rows, and differs by at most 0.01; the harness driver prints both). With both loan variables kept, income adds little.
- **`credit_score` is actionable.** The applicant can act on it: Step 10 constrains it to `c("increase", "integer")`, a whole-number score that may only rise, and the local counterfactual in `review_section9_to_15.md` raises it.

**What it feeds.** `feature_cols_v2` keeps the `var_names` order: `loan_amnt`, `loan_int_rate`, `loan_percent_income`, `credit_score`. The refit, Steps 4b–6b, rebuilds the pipeline on these four features:
- `bl_prepare_data()`, then z-scoring fitted on the training split and recorded with `bl_set_scaling()` (`review_section4_to_6.md`, "Optional: raw-unit biplot axes via `bl_set_scaling()`").
- An XGBoost model wrapped through the direct `"XGB"` path.
- `bl_build_result(method = "CVA", b_margin = 0.01)`, which returns `bl_results_v2`.

In the harness run the reduced model's training accuracy is 0.8161 (Gini 0.7747), against 0.8280 (Gini 0.8093) for the full model. Step 7b repeats Step 7 on `bl_results_v2`, and Phase 3 (`review_section9_to_15.md`) works only with `bl_results_v2`.

---

## Complete Object Flow: Step 7 and variable pruning

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
     │    prune closed contours by training-row membership
     │    consistency pruning: back-project + re-score each segment
     │    .nearest_idx_block() → nearest vertex on each contour per obs
     │    select opposing-side contour per obs
     │    back-project B_z → B_x (n × p)
     │    score B_x → B_pred; measure dist_z, dist_x, dist_mahalanobis
     ▼
bl_bnd  [bl_boundary]
  ├── B_z          (n × 2, boundary in Z-space)
  ├── B_x          (n × p, counterfactual in X-space)
  ├── B_pred       (predicted probability at counterfactual, ≈ 0.5)
  ├── dist_z       (Z-space distance per obs)
  ├── dist_x       (standardised X-space distance per obs)
  ├── dist_mahalanobis (squared Mahalanobis distance per obs)
  ├── Z_obs        (observation positions in Z-space)
  ├── x_obs        (observation feature values)
  ├── pred_obs     (predicted probability at obs)
  ├── class_obs    (true labels)
  ├── z_boundaries_list, z_boundary_type, nr_boundaries  (surviving contours)
  ├── position_store  (contour index used per obs)
  └── bl_result    (reference to bl_results)
     │
     ├── hist(bl_bnd$B_pred)
     │     └─→ diagnostic: B_pred should cluster near 0.5
     │
     ├── plot(bl_results, points=test_pts)   [or plot_biplotEZ(...)]
     │     └─→ clean biplot (no boundary overlay)
     │
     ├── bl_pick_point(bl_results, data=bl_results$test_data, bl_boundary=bl_bnd)    [interactive]
     │     └─→ per-click: yellow circle + ×-cross at B_z + arrow from Z_obs → B_z
     │
     ├── plot(bl_bnd)                     [type = "jitter"]
     │     ├── compute vec_to_boundary_sd  (n × p signed distances)
     │     ├── sort variables by sum_of_distance
     │     ├── render jitter plot (one dot per obs per variable)
     │     └─→ console: metric, robustness scalar, per-variable totals,
     │                  Mahalanobis check (+ note when |difference| >= 25%)
     │
     ├── plot(bl_bnd, type = "boxplot")   [type = "boxplot"]
     │     ├── same vec_to_boundary_sd computation
     │     └─→ box-and-whisker plot per variable × confusion group
     │
     │  Variable pruning: bl_robustness(bl_bnd)   [same totals as plot(bl_bnd), no plot]
     ▼
rob  [list]
  ├── sum_of_distance  (per-variable totals: used to rank features)
  ├── robustness       (overall score 9203.29: printed, not used)
  └── distance         "mahalanobis"
     │
     │  var_imp <- sort(rob$sum_of_distance)
     │  vars_to_drop (analyst's choice) → setdiff(bl_dat$var_names, vars_to_drop)
     ▼
feature_cols_v2  (loan_amnt, loan_int_rate, loan_percent_income, credit_score)
     │
     │  Steps 4b–6b: bl_prepare_data() → z-score + bl_set_scaling() → XGB → bl_build_result()
     ▼
bl_results_v2  [bl_result]   → Step 7b, then Phase 3 (review_section9_to_15.md)
```

---

## Key Methodological Points

1. **Z-space search, X-space result.** The nearest boundary search is entirely in 2D Z-space (fast, geometric). The result is then back-projected to p-dimensional X-space to produce an actionable counterfactual. This two-step design is what makes the method scalable.

2. **Consistency pruning is the accuracy guard.** Because the biplot is a linear 2D projection of a nonlinear model's decision surface, contour segments can be misleading. The back-project-and-re-score step ensures that only contour vertices that are genuinely near the true decision boundary (in X-space) are used.

3. **Opposing-side selection defines the counterfactual direction.** A class-1 observation's counterfactual must lie on the class-0 boundary — the nearest point where the prediction flips. This is what makes `dist_x` a meaningful measure of "how far is this applicant from becoming safe?"

4. **`sum_of_distance` is a variable importance proxy for the boundary.** Unlike model-level feature importance (e.g. XGBoost's `gain`), `sum_of_distance` tells you which features matter for reaching the boundary from where the population currently sits — a local, population-anchored importance measure.

5. **`bl_result` stored inside `bl_bnd`.** Every downstream method (`plot.bl_boundary`, `bl_robustness`) re-extracts `V`, `tV`, `X_sd`, etc. from `bl_bnd$bl_result`. This means the boundary object is self-contained and fully reproducible without requiring `bl_results` to be in scope.

6. **Pruning reads the components, not the overall score.** `rob$robustness` adds the per-variable totals into one number for the whole test population, so it cannot rank features; the pruning ranks them by `rob$sum_of_distance`. The ranking guides the choice, and domain knowledge settles it: `person_income` is dropped as redundant and `credit_score` is kept as actionable, against their ranks.
