# Code Walkthrough: Steps 9–15 of `03_loan_status_Boundary_Logic.R`

## Context

This document continues from `review_section7_to_8.md`. At this point `bl_results_v2` is the Phase 1 anchor built on the reduced feature set (after variable pruning), and `test_pts_v2` is the projected test-set `bl_points` object. Steps 9–15 form Phase 3: single-observation local interpretation — selecting one target, finding its local counterfactual, attributing feature contributions, and building the sparsest possible flip.

Note: Phase 3 uses `bl_results_v2` (6-feature reduced model), not `bl_results` (the original 10-feature model).

---

## Step 9 — `bl_predict()` + `bl_select_target()` + `plot()` (or `plot_biplotEZ()`) with target

### 9a — `bl_predict()` → `pred_summary`

**File:** `R/project_points.R`

**Call:**
```r
pred_summary <- bl_predict(bl_results_v2)
```

**What it does:**

1. Defaults `data` to `bl_results_v2$test_data`.
2. Calls `bl_project_points(data, bl_results_v2)` internally — projects test data to Z-space and scores through XGBoost (see `review_section4_to_6.md`).
3. Assembles a flat data frame, one row per test observation:

**Returned object — `pred_summary` (plain `data.frame`):**

| Column | Content |
|---|---|
| `row` | 1-based integer row index into test data |
| `pred_prob` | Predicted default probability (rounded to 3 d.p.) |
| `pred_class` | Predicted class: 0 or 1 |
| `true_class` | True loan status (0/1) — included because test data has a `"class"` column |
| `confusion` | `"TP"`, `"TN"`, `"FP"`, or `"FN"` |
| `person_age` … `credit_score` | All feature columns (in `var_names` order) |

**Purpose:** The analyst reads this table to identify a target row — typically a False Negative (predicted 0, true 1: the model clears an applicant who will actually default).

**Where `bl_predict()` is called from:**

No internal callers — `bl_predict()` is purely a user-facing convenience wrapper. Nothing in the package consumes its output; it exists solely for analyst inspection at the "pick a target" step.

External callers — the canonical workflow pattern (per README and `1 Foundation intro documents.txt`) places `bl_predict()` immediately before `bl_select_target()`:

```
... → bl_assemble() → bl_predict() → bl_select_target() → set_filters() → bl_find_local_cf() → ...
```

Used this way in scripts 00, 01, 03, 05, 06 and both vignettes. The analyst reads the returned `pred_summary` data frame, picks an interesting row (often a False Negative), and passes that row index to `bl_select_target()`. One non-target use exists: `scripts/02_contour_inspection.R:60` calls `bl_predict(bl_results, data = X_orig)` to score back-projected contour points for an inspection table.

`bl_project_points()` vs `bl_predict()` — same projection + scoring under the hood; the only difference is the return shape. Use `bl_project_points()` when you want a `bl_points` object for plotting; use `bl_predict()` when you want a tidy data frame for table-style inspection. (`bl_predict()` literally calls `bl_project_points()` at `R/project_points.R:229`.)

---

### 9b — `bl_select_target()` → `tgt`

**File:** `R/local_cf.R`

**Call:**
```r
tdp <- 2
tgt <- bl_select_target(bl_results_v2, target = tdp)
```

`tdp = 2` means row 2 of `bl_results_v2$test_data`.

**What it does:**

1. Extracts `x_obs = test_data[2, var_names]` — a single-row data frame of feature values.
2. **Projects to Z-space** using the same formula as everywhere else:
   ```
   X_centered = x_obs - X_center          # subtract training means
   Z_obs      = X_centered %*% V[, proj_dims]   # project to 2D
   ```
   (No `X_sd` scaling — CVA, `standardise = FALSE`)
3. **Scores through XGBoost** via `.pred_function()`: `pred_prob` scalar, `pred_class` integer (0 or 1).

**Alternative call — external data frame:**
```r
tgt_ext <- bl_select_target(bl_results_v2, target = new_applicant)
```
When `target` is a data frame instead of an integer, `x_obs` is taken directly from it and `row_id` is set to `NA_integer_` (external point, no true label available).

**Returned object — `tgt` (S3 class `"bl_target"`):**

| Field | Type | Content |
|---|---|---|
| `x_obs` | data.frame (1 row × p cols) | Feature values of the selected observation |
| `z_obs` | matrix (1 × 2) | Z-space coordinates, columns `"x"` and `"y"` |
| `pred_prob` | numeric scalar | Predicted default probability |
| `pred_class` | integer (0 or 1) | Predicted class at cutoff 0.5 |
| `row_id` | integer | Row index in test data (2 here); `NA_integer_` for external |

**`print(tgt)` console output:**
```
-- bl_target --
  Row ID         : 2
  Predicted prob : 0.3120  (for example)
  Predicted class: 0
  Feature values :
  person_age  loan_amnt  loan_int_rate  ...
      28.0000  15000.0000      11.5000  ...
```

---

### 9c — `plot()` with `target_point` and biplot label customisation

```r
plot(
  bl_results_v2,
  points       = test_points,          # bl_project_points(test_data[tdp, ], ...)
  target_point = target_value,         # bl_results_v2$test_data[tdp, ]
  target_label = tdp,
  label_dir         = "Paral",
  label_offset_var  = c("person_age", "loan_amnt", "loan_int_rate",
                        "loan_percent_income", "credit_score"),
  label_offset_dist = c(0, 0, 0.5, 0, 0)
)
```

Note: `test_points` is the projection of just the target row (`bl_results_v2$test_data[tdp, ]`), used so the highlight plot shows only the target observation rather than a wider subset -- per the script's Step 9 comment, "Can be empty, or only the target data point to avoid unnecessary information on the biplot." The `label_*` arguments are the same cosmetic biplot-label controls documented in `review_section4_to_6.md` Step 7, passed here via `plot()`'s `...` argument to `plot_biplotEZ()`.

If the points layer should be suppressed entirely -- the "can be empty" case in the script's Step 9 comment -- pass `plot_points = FALSE` to `plot()`/`plot_biplotEZ()`. This skips the data-points layer regardless of what `points` resolves to, including the `points = NULL` default (which would otherwise auto-project and draw the full training set).

Same 5-layer biplot as before (grid, test points, axes, contour), plus:

**Layer 6 — Target circle:**
- `target_point = target_value` (i.e. `bl_results_v2$test_data[tdp, ]`): data frame row; `plot_biplotEZ()` calls `as.numeric(target_point[var_names])` internally, so no `unlist()` is needed
- Inside `plot_biplotEZ()`, this is centered (`- X_center`) and projected to Z-space: `target_z = target_st %*% V[, proj_dims]`
- Rendered as a large filled circle (`pch = 21`, `cex = 1.8`) with the confusion category colour as background (or red/blue if true class unknown)
- `target_label = 2` is displayed inside the circle

---

## Step 10 — `set_filters()` → `flt`

**File:** `R/local_cf.R`

**Call:**
```r
flt <- set_filters(
  tgt,
  person_age    = "fixed",
  loan_int_rate = "increase",
  credit_score  = "increase"
)
```

**What it does:**

1. Validates `tgt` is a `"bl_target"` object.
2. Validates all named arguments are feature names in `tgt$x_obs`.
3. Validates each constraint value is one of `"decrease"`, `"increase"`, `"fixed"`, or a `c(min, max)` numeric pair.
4. Stores everything.

**Constraint semantics (applied in X-space after back-projecting contour points):**

| Value | Applied as | Effect on search |
|---|---|---|
| `"decrease"` | `B_x[nm] <= x_obs[nm]` | Counterfactual must borrow less than observed |
| `"increase"` | `B_x[nm] >= x_obs[nm]` | Counterfactual must be greater than observed |
| `"fixed"` | `abs(B_x[nm] - x_obs[nm]) <= 0.5` | Counterfactual value within ±0.5 of observed; in sparse CF, always reverts to observed |
| `c(min, max)` | `B_x[nm] >= min & B_x[nm] <= max` | Counterfactual must lie within absolute range |

No constraint = unconstrained (all contour points accepted for that feature).

**Returned object — `flt` (S3 class `"bl_filters"`):**

| Field | Type | Content |
|---|---|---|
| `constraints` | named list | Each element is a constraint value (`"decrease"`, `"increase"`, `"fixed"`, or `c(min, max)`) |
| `bl_target` | `bl_target` object | The `tgt` passed in — carried through so `bl_find_sparse_cf()` can access it |

**`print(flt)` console output:**
```
-- bl_filters --
  person_age           : fixed
  loan_int_rate        : increase
  credit_score         : increase
```

---

## Step 11 — `bl_find_local_cf()` → `bl_local`

**File:** `R/local_cf.R`

**Call:**
```r
bl_local <- bl_find_local_cf(
  bl_result   = bl_results_v2,
  set_filters = flt,
  bl_target   = tgt
)
```

Parameters:

| Parameter | Default | Description |
|---|---|---|
| `bl_result` | — | Phase 1 anchor |
| `bl_target` | — | Single-observation target |
| `set_filters` | `NULL` | Actionability constraints; `NULL` = unconstrained |
| `max_pairs` | `10L` | How many eigenvector pairs to try |
| `m` | `200L` | Grid resolution per pair |
| `verbose` | `TRUE` | Print progress messages |

**Why this function exists:** The global `bl_find_boundary()` uses the original biplot projection (fixed `proj_dims = c(1, 2)`). For a single target, we can do better: rotate the entire biplot so the target observation lies along the first axis of the projected plane, then search across multiple eigenvector pairs to find whichever rotation brings the boundary closest. This local rotation maximises the chance of finding a valid counterfactual.

**`bl_local$bl_counterfactual` (companion to `bl_local$bl_target`):** the result carries a
`"bl_counterfactual"` object packaging the boundary point (`x_cf = B_x` model units, `z_cf`,
`pred_prob = B_pred`, `pred_class`, `scaling`), `NULL` when no solution was found. Like
`bl_target`, its print method shows the counterfactual in original (raw) units when a scaling
was recorded via `bl_set_scaling()` (reuses `.scale_to_raw()`). `print(bl_local)` shows the
Target and Counterfactual blocks together. The raw counterfactual values are also available in
the sparse table (`B_x` column, Step 14).

---

### The algorithm: 10-pair loop

**Pre-loop: standardise target and training data**
```
x_target_st  = (x_obs - X_center) / sv    # sv = 1s for CVA (standardise = FALSE)
X_train_st   = scale(train_data[, var_names], center = X_center, scale = FALSE)
```

**Generate pairs:** `.generate_pairs(p, max_pairs)` produces up to 10 eigenvector index pairs in canonical order: (1,2), (1,3), (2,3), (1,4), (2,4), (3,4), … Each pair defines a 2D projection plane from the full p-dimensional loading matrix.

**For each pair `(i, j)`, the loop does:**

#### Stage A — SVD rotation via `.bl_rotate()`

**File:** `R/local_cf.R` (private helper)

The rotation aligns the biplot so the target observation lies along the first axis of the chosen eigenvector pair. This is done via SVD:

```
Y = rbind(-x_target_st, rep(0, p), x_target_st)   # 3 × p: [-target, 0, target]
YV  = Y %*% V                                      # 3 × p in full loading space
YVr = Y %*% Vr                                     # 3 × 2 in pair subspace

# Pad YVr to 3 × p so shapes match
YVr_padded = cbind(YVr, zeros)

# SVD of cross-covariance: finds the rotation that best maps YVr_padded onto YV
svd_res = svd(t(YV) %*% YVr_padded)
A       = svd_res$v %*% t(svd_res$u)   # p × p orthogonal rotation matrix

Vrho  = V %*% t(A)        # Rotated full loading matrix
tVrho = solve(Vrho)        # Its inverse

Vr_rot  = Vrho[, c(1, 2)]    # p x 2 -- always cols 1-2; SVD concentrates target info here
tVr_rot = tVrho[c(1, 2), ]   # 2 x p inverse (for back-projection)
```

**Why:** By constructing Y from `[-target, 0, target]` and finding the rotation that maps the pair-space projection of Y onto its full-space projection, we ensure the target observation aligns with the first axis of the rotated biplot. This places the target at a known position and makes boundary search more reliable.

**Key invariant:** `YVr_padded` is non-zero only in its first two columns, so the SVD rotation `A` always concentrates the target's information in columns 1 and 2 of `Vrho = V %*% t(A)` -- regardless of which `proj_pair` was used to build `Vr`. The return therefore always selects `Vrho[, c(1, 2)]`, not `Vrho[, proj_pair]`. Selecting `proj_pair` columns (the former bug) would return near-zero columns for any pair other than `c(1, 2)`, placing the target at approximately the biplot origin in `plot(bl_local)`.

#### Stage B — Build m×m grid in rotated Z-space

```
Z_target    = x_target_st %*% Vr_rot       # 1 × 2
Z_train_rot = X_train_st %*% Vr_rot        # n × 2

z_range = range(Z_train_rot)
z_pad   = diff(z_range) * 0.10             # 10% padding
xseq    = seq(z_range[1]-pad, z_range[2]+pad, length.out = m)
Zgrid   = expand.grid(xseq, xseq)          # m² × 2
```

The 10% padding is needed because decision boundary contours frequently run near the *edge* of the training data cloud — that is where the class transition occurs. Without padding the grid would be flush with the outermost training observations, and contour lines that fall at or beyond the data range would be clipped, leaving no valid contour vertex for that pair. Adding 10% of the total range on each side gives boundary segments room to exist slightly outside the training envelope. The same rationale applies to the global grid built in `bl_build_grid()`.

#### Stage C — Back-project grid and score

```
Xgrid = Zgrid %*% tVr_rot + X_center      # m² × p in original feature space
grid_prob = .pred_function(model, "XGB", Xgrid)   # chunked, 50 000 rows at a time
```

#### Stage D — Extract contour lines

```
ct_local = contourLines(xseq, yseq, matrix(grid_prob, ncol=m),
                         levels = c(cutoff - b_margin, cutoff + b_margin))
```
`b_margin = 0.01` (set directly via the `b_margin` parameter), so contours at 0.49 and 0.51.

#### Stage E — Filter contour segments (three successive filters)

For each contour segment in `ct_local`:

1. **Opposing-class filter:** Keep only segments where `seg_class != target_pred_class`. A class-1 target needs class-0 contours (where it would flip to "safe").

2. **Back-project + feasibility:** Back-project segment to X-space, apply `train_ranges` filter, then apply `set_filters` actionability constraints via `.apply_actionability()`:
   - `"decrease"`: `B_x[nm] <= x_obs[nm]`
   - `"increase"`: `B_x[nm] >= x_obs[nm]`
   - `"fixed"`: `abs(B_x[nm] - x_obs[nm]) <= 0.5`
   - `c(min, max)`: range check

3. **Model consistency:** Re-score filtered back-projected points; keep only rows where prediction agrees with contour's probability level.

If all segments are eliminated, the pair is skipped.

#### Stage F — Find nearest surviving point

All surviving contour rows are combined into `all_Mi`. The nearest to `Z_target` is found via `.nearest_idx_block(Z_target, all_Mi, block=1)`. Distance = `sqrt(sum((Z_target - B_z_local)^2))`.

#### Stage G — Update best result

For every pair (whether or not it wins), the back-projection from `B_z_local` to `B_x` is computed and the **squared Mahalanobis distance in X-space** is calculated:

```
v        = B_x - x_obs                              # 1 x p in X-space
dist_mah = as.numeric(t(v) %*% metric_inv %*% v)    # scalar squared Mahalanobis
```

where `metric_inv = W^{-1}` is the inverse of the pooled within-class covariance matrix stored on `bl_result$metric_inv` (added in 2026-05-20). The Cholesky decomposition is used for the inversion -- see `2 implementation_summary.txt` Section 4.2.1 for the maths and singularity safeguards.

The pair selector then depends on the `distance` argument of `bl_find_local_cf()`:

- `distance = "mahalanobis"` (default): pair selected by minimum `dist_mah`. Pair-invariant; matches the PhD's `pc.distz <- v %*% Wmat_inv %*% t(v)` formula.
- `distance = "euclidean"` (legacy): pair selected by minimum `dist_z`. Not comparable across pairs because each pair captures a different amount of variance in its 2D slice.

If the chosen selector for this pair beats the running `best_selector`, save everything for this pair as the new best result. The loop continues through all configured pairs, keeping the global minimum. Both `dist_z` and `dist_mah` are reported for the winning pair, plus per-pair tables `all_distances` and `all_distances_mahalanobis`.

---

**Returned object — `bl_local` (S3 class `"bl_local_result"`):**

| Field | Type | Content |
|---|---|---|
| `B_z` | matrix (1 × 2) | Nearest boundary point in rotated Z-space |
| `B_x` | data.frame (1 row × p) | Back-projected counterfactual in X-space |
| `B_pred` | numeric scalar | Model probability at counterfactual (should be ≈ 0.5) |
| `dist_z` | numeric scalar | Z-space distance from target to boundary |
| `Z_target` | matrix (1 × 2) | Target position in rotated Z-space |
| `best_pair` | integer (length 2) | Winning eigenvector pair, e.g. `c(1L, 2L)` |
| `Vr_rot` | matrix (p × 2) | Rotated loading matrix for best pair |
| `tVr_rot` | matrix (2 × p) | Inverse loading matrix for best pair |
| `Z_train_rot` | matrix (n × 2) | Training data in rotated Z-space |
| `Zgrid` | matrix (m² × 2) | Grid coordinates in rotated Z-space |
| `grid_prob` | numeric (m²) | Model probabilities at grid points |
| `col_value` | character (m²) | Colour per grid point (blue→red palette) |
| `ct_local` | list | Raw contour list from `contourLines()` |
| `xseq`, `yseq` | numeric (m) | Grid axis sequences |
| `min_val`, `max_val` | numeric | Grid axis limits |
| `all_distances` | named numeric | Best Z-distance per pair tried (NA if no solution) |
| `solution_found` | logical | `TRUE` if at least one valid boundary was found |
| `blocking_constraint` | character or NULL | Reason message if no solution found |
| `bl_target` | `bl_target` | The `tgt` object |
| `bl_result` | `bl_result` | Reference to `bl_results_v2` |
| `set_filters` | `bl_filters` or NULL | The `flt` object (needed by `bl_find_sparse_cf()`) |

**`print(bl_local)` console output:**
```
-- bl_local_result --
  Solution found : TRUE
  Best pair      : (1, 2)
  Distance (Z)   : 0.1342
  Boundary pred  : 0.4900

  Target:
-- bl_target --
  Row ID         : 2
  ...

  All Z-distances by pair:
  pair   distance
  (1,2)   0.1342
  (1,3)       NA
  ...
```

---

## Step 12 — `plot(bl_local)` → `plot.bl_local_result()`

**File:** `R/local_cf.R`

Default: `plot_points = FALSE` (training points hidden), `no_grid = FALSE`, `no_contour = FALSE`.
`rotate_deg = 0` and `label_cex = 1` are now available. `label_offset_var` accepts variable names or integer indices. `label_dir = "Paral"` (border-adaptive) is the new default.

`new_title` (default `NA`) overrides the auto-generated `"Local biplot -- target N [pair (i,j) | p=..., class ...]"` title; `NA` keeps the auto-generated title.

**Raw-unit axis relabel:** if `bl_result$scaling` is set (via `bl_set_scaling()` — see
`review_section4_to_6.md`), the rotated local biplot's axes are relabelled into original feature
units. After the SVD rotation has patched `Z`/`Lmat`/`ax.one.unit` (and any extra `rotate_deg`),
`plot.bl_local_result()` calls `.bl_rescale_biplot_axes(biplot_plot, bl_result$scaling, var_names)`.
The relabel transforms only the biplotEZ `$means`/`$sd`, which are **rotation-invariant**, so it
composes cleanly with the local rotation and leaves the plotted geometry untouched. The same
applies to `plot.bl_surrogate()` and, by pass-through, `plot.bl_sparse_result()` (Step 14).
Counterfactual/Shapley *values* in the console output remain in model (standardised) units.

**Base-graphics title-error note:** `new_title` is written straight to `biplot_plot$Title` (no validation, matching `plot_biplotEZ()`), and biplotEZ forwards that field to `graphics::title(main = ...)` at plot-flush time. A length-1 character string (or `NA` to keep the default) is the intended contract. `title()` is tolerant: it coerces most atomic and even recursive inputs (numeric, multi-element character vectors, lists, data.frames) to a label via `as.character()`, so those do *not* error — they just produce an odd/recycled title. The error case is a value base R cannot coerce to a character vector at all — a **function/closure** or an **environment** — which raises `"cannot coerce type 'closure' to vector of type 'character'"` from *inside* the biplotEZ `plot()` call, not from the plot function's own argument handling. A length-1 character is the safe, intended input.

**What it renders (8 layers):**

| Layer | Content | Data source |
|---|---|---|
| 1 | Grey axes + variable labels | `biplot_obj` patched with `Vr_rot` and `Z_train_rot` |
| 2 | Prediction grid (blue→red) | `bl_local$Zgrid` + `bl_local$col_value` |
| 3 | Training points (optional, hidden by default) | `Z_train_rot` + confusion colours |
| 4 | Axes redrawn darker (via `par(new=TRUE)`) | Same biplotEZ object |
| 5 | Decision boundary contour lines | `bl_local$ct_local` |
| 6 | Target circle (filled, confusion colour) | `Z_target`; coloured TP/TN/FP/FN if true class known |
| 7 | CF cross (×) + arrow | `B_z` (cross); arrow from `Z_target` → `B_z` |
| 8 | Console summary | Printed to console (not plotted) |

**Key rendering detail — patching the biplotEZ object:**

The biplotEZ object stored in `bl_results_v2$biplot_obj` was built in the original (unrotated) CVA space. For the local plot, it is patched in-place:
```r
biplot_plot$Z[, proj_dims]    <- Z_train_rot   # replace training coords with rotated
biplot_plot$Lmat[, proj_dims] <- Vr_rot        # replace loadings with rotated
biplot_plot$ax.one.unit       <- (1 / diag(t(tVr_rot) %*% tVr_rot)) * t(tVr_rot)
```
This means the axes, variable label positions, and tick marks are all recalculated in the rotated coordinate frame — the biplot shows the locally rotated view, not the original CVA projection.

**Console output (Step 8 of plot):**
```
--- Local CF summary ---
  Target         : 2
  Pred prob      : 0.3120  (class 0)
  Best pair      : (1, 2)
  Distance (Z)   : 0.1342
  Boundary pred  : 0.4900
------------------------
```

---

## Step 13 — `bl_shapley()` + `print()` + `plot()`

**File:** `R/shapley.R`

**Call:**
```r
bl_shapley_values <- bl_shapley(bl_local)
```

**Raw-unit feature values (when `bl_result$scaling` is set via `bl_set_scaling()`).** The
observed and counterfactual values shown by `plot.bl_shapley()`, `print.bl_shapley()`, the
sparse table (Step 14), and `print.bl_target()` are converted back to original units at display
time by private `.scale_to_raw()` (`R/data_prepare.R`). Levels (observed, `B_x`, sparse,
counterfactual) map as `value*scale + center`; the `data_to_boundary` **delta** maps as
`value*scale` (no centre term). When scaling is present the `plot.bl_shapley()` y-axis label
switches from `observed -> change` to `observed -> counterfactual` (raw units) and the y-title
notes "raw units"; the print tables gain a "(raw units)" header. `bl_select_target()` stores
`scaling` on the `bl_target` so `print.bl_target()` can convert. The Shapley *contributions*
(`shapley_cause`) stay in prediction-impact units -- not converted. No-op when scaling absent.
See `review_section4_to_6.md` for `bl_set_scaling()` and `2 implementation_summary.txt` §4.9.

### What `bl_shapley()` does

Explains **which features drove the prediction from the observed point to the counterfactual**, using the Shapley value framework from cooperative game theory.

**The game:**
- Players: the p features
- "Start": the observed feature values `x_obs` (e.g. row 2 of test data)
- "End": the counterfactual `B_x` (the boundary point back-projected to X-space)
- The model's probability change from `pred_obs` to `B_pred` must be attributed across features

**Method selection (line 178):**
- If `p <= exact_max_vars` (default 14): use **exact** computation via `.shapley_exact_one()`
- If `p > 14`: use **approximate** via `.shapley_perm_one()` with `M = 2048` permutations

The loan dataset (after pruning) has 6 features, so exact computation is used here.

---

#### Exact method: `.shapley_exact_one()`

For each subset size `k = 0, 1, …, p-1`:

1. Enumerate all `C(p, k)` subsets S of size k.
2. For each subset S, build a "mixed" data row:
   ```
   X_S[j] = x_obs[j]   if feature j NOT in S
   X_S[j] = B_x[j]     if feature j IS in S
   ```
   So S = {loan_amnt} means: use the counterfactual value for `loan_amnt`, observed values for everything else.
3. Score `X_S` through XGBoost → `f(S)`
4. For each feature i not already in S, build `S ∪ {i}` → score → `f(S ∪ {i})`
5. Marginal contribution of feature i given coalition S: `f(S ∪ {i}) - f(S)`
6. Weight by: `k! * (p-k-1)! / p!` (the Shapley weighting factor)
7. Shapley value for feature i = sum of weighted marginal contributions across all subsets not containing i

This is the unique fair allocation of the total prediction change `f(end) - f(start)` to individual features.

#### Approximate method: `.shapley_perm_one()`

For M = 2048 random permutations:
1. Sample a random ordering of all p features
2. Build `p+1` data rows: start with all-observed, sequentially switch features to counterfactual values following the permutation order
3. Score all rows; take successive differences → marginal contributions
4. Accumulate across permutations; divide by M

---

**Contribute classification (line 206):**

Each feature is classified as `"Supports"` or `"Contradicts"` the move to the boundary:

| `pred_class` | `shapley_cause` | Label |
|---|---|---|
| 0 (predicted safe) | ≥ 0 (pushes prob up, toward boundary) | `"Supports"` |
| 0 (predicted safe) | < 0 (pushes prob down, away from boundary) | `"Contradicts"` |
| 1 (predicted default) | < 0 (pushes prob down, toward boundary) | `"Supports"` |
| 1 (predicted default) | ≥ 0 (pushes prob up, away from boundary) | `"Contradicts"` |

**Row ordering:** For class-0 target: sorted ascending by `shapley_cause` (most negative = biggest supporter at top). For class-1: sorted descending.

**`varnames_p` label** (used as Y-axis in plot):
```
"loan_amnt: 15000 -> -3250.5"
```
Format: `"feature: observed_value -> change_to_reach_boundary"`

---

**Returned object — `bl_shapley_values` (S3 class `"bl_shapley"`):**

| Field | Type | Content |
|---|---|---|
| `shapley_df` | data.frame | One row per feature, sorted by contribution. Columns: `varnames`, `pred_data` (observed value), `data_to_boundary` (change needed), `shapley_cause` (Shapley value), `pred_use` (target's predicted class), `Contribute`, `varnames_p` (formatted label) |
| `pred_prob` | numeric | Target's predicted probability |
| `pred_class` | integer | Target's predicted class |
| `pred_boundary` | numeric | Model probability at counterfactual |
| `row_id` | integer | Row index (or `NA` for external) |
| `bl_local_result` | `bl_local_result` | Reference to `bl_local` |

---

**`print(bl_shapley_values)` console output:**
```
-- bl_shapley --
  Row ID         : 2
  Pred prob      : 0.3120
  Pred class     : 0
  Boundary pred  : 0.4900

  Shapley table:
  varnames  pred_data  data_to_boundary  shapley_cause  Contribute
  loan_amnt  15000.000         -3250.500         0.1234    Supports
  loan_int_rate  11.500            0.020         0.0891    Supports
  ...
```

---

**`plot(bl_shapley_values)` — Shapley bar chart:**

Horizontal bar chart (ggplot2):
- **Y-axis:** `varnames_p` — one bar per feature, labelled with `"feature: observed -> change"`; sorted by Shapley magnitude
- **X-axis:** `shapley_cause` — the Shapley value (how much this feature contributed to the prediction move)
- **Fill colour:** `"Supports"` = dark blue; `"Contradicts"` = grey
- **Value label:** rounded Shapley value displayed inside each bar in white
- **Title:** `"Shapley Contribution Plot"`
- **Subtitle:** `"ID: 2 | Class: 0 | Pred: 0.312 | Boundary pred: 0.490"`
- **X-axis label:** `"[-> change to reach boundary]: Impact on prediction"`

**Reading the chart:** Each bar shows how much that feature's change (from observed to counterfactual) accounts for the total probability shift from `pred_prob` to `pred_boundary`. The Shapley values sum to approximately `pred_boundary - pred_prob`.

---

## Step 14 — `bl_find_sparse_cf()` + `print()` + `plot()`

**File:** `R/shapley.R`

**Call:**
```r
bl_sparse <- bl_find_sparse_cf(bl_shapley_values, round_to = NULL)
```

### What `bl_find_sparse_cf()` does

Builds the **minimal counterfactual**: keep only the features whose Shapley value `"Supports"` the prediction flip; revert all other features to their observed values. This answers "what is the fewest changes needed to flip the prediction?"

**Algorithm (lines 364–385):**

1. Start with `x_sparse = x_obs` (all features at observed values).
2. For each feature:
   - If `Contribute == "Supports"`: set `x_sparse[nm] = B_x[nm]` (use counterfactual value)
   - If `Contribute == "Contradicts"` or `"Unknown"`: leave at `x_obs[nm]`
3. **Fixed-constraint override (lines 377–385):** If any feature was marked `"fixed"` in `set_filters`, it always reverts to its observed value in the sparse CF — regardless of its Shapley classification. (The `"fixed"` constraint means the feature is not actionable at all; it was included in the boundary search only to inform the direction, not to prescribe a change.)
4. **Optional rounding** (`round_to`): if not `NULL`, round each supporting feature's new value to the nearest multiple of `round_to`. E.g. `round_to = 500` would round `loan_amnt` to the nearest 500.
5. **Score** `x_sparse` through XGBoost → `pred_sparse`.
6. **Validate:** `solution_valid = (pred_class_sparse != pred_class_target)` — does the sparse CF actually flip the prediction?

**Returned object — `bl_sparse` (S3 class `"bl_sparse_result"`):**

| Field | Type | Content |
|---|---|---|
| `x_sparse` | data.frame (1 row × p) | The sparse counterfactual in X-space |
| `pred_sparse` | numeric scalar | Model probability at sparse counterfactual |
| `solution_valid` | logical | `TRUE` if sparse CF flips predicted class |
| `shapley_df` | data.frame | Annotated Shapley table with extra columns `used_in_sparse` (logical) and `x_sparse_val` (the value used in sparse CF) |
| `bl_shapley` | `bl_shapley` | Reference to `bl_shapley_values` |

---

**`print(bl_sparse)` console output:**
```
-- bl_sparse_result --
  Solution valid : TRUE

  Predictions:
    Observed     : 0.3120  (class 0)
    Full CF      : 0.4900
    Sparse CF    : 0.5230

  Variable summary:
  variable       x_obs    B_x  x_sparse  used_in_sparse
  loan_amnt   15000.0  11749.5  11749.5            TRUE
  loan_int_rate  11.5     11.5     11.5           FALSE
  ...
```

The `used_in_sparse` column shows which features actually changed. Features with `"fixed"` constraint will show `x_obs == x_sparse` regardless.

---

**`plot(bl_sparse)` — sparse CF biplot overlay:**

1. **Calls `plot.bl_local_result(x$bl_shapley$bl_local_result, print_summary = FALSE, ...)`** — renders the full local biplot (all 7 layers: grid, axes, target circle, full CF cross + arrow). The `print_summary = FALSE` suppresses the local CF console summary. Because every `plot.bl_local_result()` argument is forwarded through `...`, `new_title` passes straight through: `plot(bl_sparse, new_title = "...")` retitles the underlying local biplot (see Step 12 for the base-graphics title-error caveat).

2. **Projects `x_sparse` to rotated Z-space:**
   ```
   x_sparse_st = (x_sparse - X_center) / sv
   z_sparse    = x_sparse_st %*% Vr_rot       # uses the ROTATED loading from best pair
   ```

3. **Plots a second cross** at `z_sparse`:
   - Colour: `"green3"` if `solution_valid == TRUE` (flip succeeded), `"yellow2"` if not
   - `pch = 4L`, `cex = 1.4`, `lwd = 2` (slightly larger than the full CF cross)

4. **Adds a legend** (top-right corner):
   - `"Full CF"` (grey cross, from `plot.bl_local_result`)
   - `"Sparse CF (valid)"` or `"Sparse CF (invalid)"` (green or yellow cross)

5. **Prints sparse CF summary** to console:
   ```
   --- Sparse CF summary ---
     Target         : 2
     Observed pred  : 0.3120  (class 0)
     Full CF pred   : 0.4900
     Sparse CF pred : 0.5230
     Solution valid : TRUE
   -------------------------
   ```

---

## Step 15 (optional) — Unconstrained local search

```r
bl_local_free  <- bl_find_local_cf(bl_results_v2, tgt)
bl_shap_free   <- bl_shapley(bl_local_free)
bl_sparse_free <- bl_find_sparse_cf(bl_shap_free, round_to = NULL)
```

Identical pipeline to Steps 11–14, with `set_filters = NULL`. The unconstrained search finds the geometrically nearest boundary point without any actionability restrictions. Comparing `bl_local` (constrained) vs `bl_local_free` (unconstrained) shows the cost of the constraints: how much farther the boundary is when `loan_amnt` must decrease and `loan_int_rate` must stay fixed.

---

## Complete Object Flow: Steps 10–14

```
bl_results_v2  [bl_result]          ← Phase 1 anchor (reduced 6-feature model)
test_pts_v2    [bl_points]          ← test data in Z-space
     │
     │  bl_predict(bl_results_v2)
     ▼
pred_summary   [data.frame]
  ├── row, pred_prob, pred_class
  ├── true_class, confusion (TP/TN/FP/FN)
  └── all feature columns
     │
     │  bl_select_target(bl_results_v2, target = 2)
     │    → extract test_data[2, var_names]
     │    → project to Z-space
     │    → score through XGBoost
     ▼
tgt  [bl_target]
  ├── x_obs    (1 × p feature values)
  ├── z_obs    (1 × 2, Z-space position)
  ├── pred_prob, pred_class
  └── row_id = 2
     │
     │  set_filters(tgt, person_age="fixed", loan_int_rate="increase", credit_score="increase")
     ▼
flt  [bl_filters]
  ├── constraints = list(person_age="fixed", loan_int_rate="increase", credit_score="increase")
  └── bl_target = tgt
     │
     │  bl_find_local_cf(bl_results_v2, tgt, set_filters=flt, max_pairs=10)
     │
     │  [for each of up to 10 eigenvector pairs]:
     │    .bl_rotate() → Vr_rot, tVr_rot (SVD-based local rotation)
     │    build m×m rotated grid → Xgrid (back-projected)
     │    score through XGBoost → grid_prob
     │    extract contour lines → ct_local
     │    filter: opposing class, train_ranges, set_filters, consistency
     │    .nearest_idx_block() → nearest boundary vertex
     │    update best if dist_z improves
     ▼
bl_local  [bl_local_result]
  ├── B_z       (1 × 2, best boundary in rotated Z)
  ├── B_x       (1 × p, counterfactual in X-space)
  ├── B_pred    (≈ 0.5)
  ├── dist_z    (Z-space distance)
  ├── Z_target, best_pair, Vr_rot, tVr_rot
  ├── Z_train_rot, Zgrid, grid_prob, col_value, ct_local
  ├── all_distances, solution_found
  ├── bl_target = tgt
  ├── bl_result = bl_results_v2
  └── set_filters = flt
     │
     ├── plot(bl_local)
     │     └─→ patched local biplot (rotated axes + grid + target circle + CF cross)
     │
     │  bl_shapley(bl_local)
     │    → start = x_obs, end = B_x
     │    → exact Shapley (p=6 ≤ 14): enumerate all 2^6 = 64 subsets
     │    → classify each feature as Supports/Contradicts
     │    → sort by contribution magnitude
     ▼
bl_shapley_values  [bl_shapley]
  ├── shapley_df   (6 rows: varnames, pred_data, data_to_boundary,
  │                 shapley_cause, Contribute, varnames_p)
  ├── pred_prob, pred_class, pred_boundary
  ├── row_id = 2
  └── bl_local_result = bl_local
     │
     ├── print(bl_shapley_values) → Shapley table to console
     ├── plot(bl_shapley_values)  → horizontal bar chart coloured by Contribute
     │
     │  bl_find_sparse_cf(bl_shapley_values, round_to=NULL)
     │    → x_sparse = x_obs (start from observed)
     │    → for "Supports" features: x_sparse[nm] = B_x[nm]
     │    → for "fixed" features: always revert to x_obs[nm]
     │    → score x_sparse → pred_sparse
     │    → solution_valid = (pred_class(x_sparse) != pred_class(x_obs))
     ▼
bl_sparse  [bl_sparse_result]
  ├── x_sparse       (1 × p sparse counterfactual)
  ├── pred_sparse    (model prob at sparse CF)
  ├── solution_valid (TRUE/FALSE)
  ├── shapley_df     (+ used_in_sparse + x_sparse_val columns)
  └── bl_shapley = bl_shapley_values
     │
     ├── print(bl_sparse) → variable-level table: x_obs, B_x, x_sparse, used
     └── plot(bl_sparse)  → local biplot + second cross (green=valid, yellow=invalid)
                             + console: Observed/Full CF/Sparse CF predictions
```

---

## Key Methodological Points

1. **Local rotation vs global boundary search.** `bl_find_boundary()` (Step 7) searches for the nearest boundary across the full test set using a fixed CVA projection. `bl_find_local_cf()` (Step 12) rotates the projection specifically around the target observation and tries multiple eigenvector pairs — this local alignment finds a closer, more actionable boundary than the global search typically can.

2. **The convex hull is NOT applied in Phase 3.** Because each eigenvector pair rotation invalidates the original grid polygon, `bl_find_local_cf()` uses `train_ranges` (per-feature min/max) and `set_filters` actionability constraints instead of the Z-space hull. The hull was the right constraint in global space; in locally-rotated space it would be incorrect.

3. **Shapley values attribute the full CF movement.** The Shapley values sum to approximately `pred_boundary - pred_prob` — the total prediction change from observation to boundary. Each Shapley value is that feature's fair share of causing the model to change its mind.

4. **Sparse CF is a greedy subset selection.** It keeps only features that `"Support"` the flip. This is not guaranteed to be the minimal subset (optimal sparse CF is NP-hard), but it is a principled, Shapley-ordered approximation. If `solution_valid = FALSE`, the remaining "supporting" features alone are not sufficient to flip the prediction — the "contradicting" features' reversal to observed values pulled the probability back below 0.5.

5. **`"fixed"` constraint has two effects:** During `bl_find_local_cf()`, it limits counterfactual search to values within ±0.5 of observed. During `bl_find_sparse_cf()`, it always reverts to observed — regardless of Shapley classification. This models features the analyst knows cannot or should not change (e.g. interest rate set by lender).

6. **The chain is fully self-contained.** Each object carries references to its parents: `bl_sparse$bl_shapley$bl_local_result$bl_result`. Every downstream method can re-extract whatever it needs (projection matrices, model, feature names) without requiring the user to pass additional arguments.
