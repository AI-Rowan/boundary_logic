# Code Walkthrough: Steps 9–15 of `03_loan_status_Boundary_Logic.R`

## Context

This document continues from `review_section7_to_8.md`. At this point `bl_results_v2` is the Phase 1 anchor built on the reduced feature set (after variable pruning; see `review_section7_to_8.md`, "Variable pruning"), and `test_pts_v2` is the projected test-set `bl_points` object. The script's Phase 3 runs Steps 9–16: single-observation local interpretation — selecting one target, finding its local counterfactual, attributing feature contributions, and building the sparsest possible flip. This walkthrough covers Steps 9–15; Step 16 (an external applicant) repeats Steps 10–14 for a data-frame target.

Note: Phase 3 uses `bl_results_v2`, the 4-feature reduced model (`loan_amnt`, `loan_int_rate`, `loan_percent_income`, `credit_score`), not `bl_results` (the 9-feature model). Its features are z-scored on the training split, and `bl_set_scaling()` records the transform, so printed feature values and biplot axes appear in raw loan units while every stored value is in z-score (model) units.

---

## Step 9 — `bl_predict()` + `bl_project_points()` + `plot()` (or `plot_biplotEZ()`) with target

### 9a — `bl_predict()` → `pred_summary`

**File:** `R/project_points.R`

**Call:**
```r
pred_summary <- bl_predict(bl_results_v2)
```

**Parameters:**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `data` | `NULL` | default | Rows to score; `NULL` uses `bl_results_v2$test_data` |

Full argument list: `?bl_predict`.

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
| `loan_amnt` … `credit_score` | All feature columns (in `var_names` order; z-score model units for `bl_results_v2`) |

**Purpose:** The analyst reads this table to identify a target row — typically a False Negative (predicted 0, true 1: the model clears an applicant who will actually default).

**Where `bl_predict()` is called from:**

No internal callers — `bl_predict()` is purely a user-facing convenience wrapper. Nothing in the package consumes its output; it exists solely for analyst inspection at the "pick a target" step.

External callers — the canonical workflow pattern (per README and `1 Foundation intro documents.txt`) places `bl_predict()` immediately before `bl_select_target()`:

```
... → bl_assemble() → bl_predict() → bl_select_target() → set_filters() → bl_find_local_cf() → ...
```

This pattern is used in the development harness (`development/03_loan_status_Boundary_Logic.R`, Step 9) and in the loan-default vignette. The analyst reads the returned `pred_summary` data frame, picks an interesting row (often a False Negative), and passes that row index to `bl_select_target()`.

`bl_project_points()` vs `bl_predict()` — same projection + scoring under the hood; the only difference is the return shape. Use `bl_project_points()` when you want a `bl_points` object for plotting; use `bl_predict()` when you want a tidy data frame for table-style inspection. (`bl_predict()` calls `bl_project_points()` internally; see `bl_predict()` in `R/project_points.R`.)

---

### 9b — `plot()` with `target_point` and biplot label customisation

```r
tdp <- 13

test_points  <- bl_project_points(bl_results_v2$test_data[tdp:10, ], bl_results_v2)
target_value <- bl_results_v2$test_data[tdp, ]

plot(
  bl_results_v2,
  points       = test_points,          # rows 13, 12, 11, 10 of test_data
  plot_points  = FALSE,                # skip the data-points layer
  target_point = target_value,         # bl_results_v2$test_data[tdp, ]
  target_label = tdp,
  label_dir         = "Paral",
  label_offset_var  = c("loan_amnt", "loan_int_rate",
                        "loan_percent_income", "credit_score"),
  label_offset_dist = c(0, 0.5, 0, 0),
  ticks_v = 3
)
```

`tdp = 13` means row 13 of `bl_results_v2$test_data`; Step 10 selects it as the target.

**Parameters:**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `points` | `NULL` | `test_points` | Rows 13 down to 10 of `test_data`, projected; not drawn, because `plot_points = FALSE` |
| `plot_points` | `TRUE` | `FALSE` | Skips the data-points layer, whatever `points` holds |
| `target_point` | `NULL` | `target_value` | Feature values of the observation to highlight (`test_data[13, ]`), drawn as a yellow circle |
| `target_label` | `NULL` | `13` | Text inside the circle |
| `label_dir` | `"Paral"` | `"Paral"` | Border-adaptive label direction (the default, set explicitly) |
| `label_offset_var` | `0L` | the four feature names | Labels to move outward; `person_age` is not in the reduced model, so it is not named |
| `label_offset_dist` | `1.5` | `c(0, 0.5, 0, 0)` | Paired with `label_offset_var` by position: only `loan_int_rate` moves, by 0.5 |
| `ticks_v` | `1L` | `3` | Three tick marks on every axis |

Full argument list: `?plot_biplotEZ`. For the full pairing rule, see `review_section4_to_6.md` Section 6.4, "Pairing rule for offsets and ticks".

`test_points` projects rows 13 down to 10 (`tdp:10`) of the test data. Because `plot_points = FALSE` skips the data-points layer, they are not drawn; `plot_biplotEZ()` still prints their prediction summary (4 points, TP / TN / FP / FN = 2 / 2 / 0 / 0 in the harness run). `plot_points = FALSE` skips the layer whatever `points` resolves to, including the `points = NULL` default (which would otherwise auto-project and draw the full training set). The `label_*` arguments are the same cosmetic biplot-label controls documented in `review_section4_to_6.md` Section 6.4, passed here via `plot()`'s `...` argument to `plot_biplotEZ()`.

The biplot draws the grid, axes and contours (the data-points layer is skipped), plus:

**Layer 6 — Target circle:**
- `target_point = target_value` (i.e. `bl_results_v2$test_data[tdp, ]`): data frame row; `plot_biplotEZ()` calls `as.numeric(target_point[var_names])` internally, so no `unlist()` is needed
- Inside `plot_biplotEZ()`, this is centered (`- X_center`) and projected to Z-space: `target_z = target_st %*% V[, proj_dims]`
- Rendered as a large filled **yellow** circle (`pch = 21`, `cex = 1.8`, `bg = "yellow"`). On the global biplot the target is a highlight only. Its confusion-group colour is shown on the local biplot (Step 12), from `tgt$class_obs`.
- `target_label = 13` is displayed inside the circle

---

## Step 10 — `bl_select_target()` → `tgt`, then `set_filters()` → `flt`

### 10a — `bl_select_target()` → `tgt`

**File:** `R/local_cf.R`

**Call:**
```r
tgt <- bl_select_target(bl_results_v2, target = tdp)   # tdp = 13
```

**Parameters:**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `target` | required | `tdp` (13) | A row index into `data`, or a one-row data frame of feature values (Step 16) |
| `data` | `NULL` | default | The frame an integer `target` indexes; `NULL` uses `test_data`. `x_obs` and `class_obs` both come from this frame |

Full argument list: `?bl_select_target`.

**What it does:**

1. Extracts `x_obs = test_data[13, var_names]` — a single-row data frame of feature values (z-score model units).
2. **Projects to Z-space** using the same formula as everywhere else:
   ```
   X_centered = x_obs - X_center          # subtract training means
   Z_obs      = X_centered %*% V[, proj_dims]   # project to 2D
   ```
   (No `X_sd` scaling — CVA, `standardise = FALSE`)
3. **Scores through XGBoost** via `.pred_function()`: `pred_prob` scalar, `pred_class` integer (0 or 1).
4. **Resolves the true class** `class_obs` from the same row that supplies `x_obs`, here `test_data[13, ]`.
   - **Which column:** it reads the `class` column when present, otherwise the first non-feature column.
   - **Why at selection time:** the true class is a property of the observation, not of a row number. Resolving it here ties it to the row that was actually chosen.
   - **Any data frame:** `data =` can be `test_data` (the default), `train_data` or any other frame. Whichever frame `row_id` indexes, `class_obs` comes from that same frame.
   - **Unknown class:** when the source has no class column, or the label is not numeric, `class_obs` is `NA`.

**Selecting from another data frame:**
```r
tgt_train <- bl_select_target(bl_results_v2, target = 5L, data = bl_results_v2$train_data)
```
`x_obs` and `class_obs` both come from `train_data[5, ]`. `row_id = 5` indexes `train_data`.

**Alternative call — external data frame (Step 16):**
```r
new_applicant <- bl_apply_scaling(new_applicant, bl_results_v2)
tgt_ext <- bl_select_target(bl_results_v2, target = new_applicant)
```

**Convert to model units first — `bl_apply_scaling()`** (`R/data_prepare.R`).
- **Why:** a data-frame target must be in the units the model was fitted on. `bl_set_scaling()` is display-only, and `bl_results_v2`'s model works in z-scores (Steps 4b–6b), while the script writes the application in raw loan units (credit score 420, loan 12,000).
- **What it does:** each feature column becomes `(raw - center) / scale`, with the training `center`/`scale` recorded on `bl_results_v2$scaling`; other columns are left alone.

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `data` | required | `new_applicant` | The new observation(s), every feature in raw units |
| `bl_result` | required | `bl_results_v2` | Supplies the recorded scaling and `var_names`; a `bl_data` or `bl_filter_result` also works |

Full argument list: `?bl_apply_scaling`.

The converted applicant (model units, from the harness run) is `loan_amnt` 0.368, `loan_int_rate` 0.620, `loan_percent_income` −1.203 and `credit_score` −5.399.

When `target` is a data frame instead of an integer, the fields are filled as follows:
- **`x_obs`:** taken directly from the frame.
- **`row_id`:** `NA_integer_`.
- **`class_obs`:** read from the frame's own class column if it has one. The script's `new_applicant` has none, so `class_obs = NA` (no true label available).

**The target-range warning.** For a data-frame target, `bl_select_target()` warns when any feature lies outside `bl_result$train_ranges`. If a scaling is recorded and converting the target would bring more features into range, the warning adds that the target looks like raw units and suggests `bl_apply_scaling()`. Here the converted applicant still warns, without that hint, because a credit score of 420 is below the training range (489 to 754; −3.73 in z-scores):
```
Warning message:
The target lies outside the training range for: credit_score. The local search looks for counterfactuals within that
range.
```
In the harness log this warning appears at the end of the run, because R defers top-level warnings.

**Step 16 result** (from the harness run):
- **The applicant:** predicted 0.296, class 0.
- **The search:** it selects pair (1, 2), with Distance (Z) 0.0135 and Mahalanobis distance 1.7026.
- **The counterfactual (prediction 0.519):**
  - raises the credit score from 420 to 489.44, the edge of the training range;
  - raises the interest rate from 13.50 to 13.99;
  - lowers the loan from 12,000 to 11,642.32;
  - lowers the loan-to-income ratio from 0.0484 to 0.0437.
- **The sparse counterfactual:** changing the interest rate alone crosses the boundary (0.52).

**Returned object — `tgt` (S3 class `"bl_target"`):**

| Field | Type | Content |
|---|---|---|
| `x_obs` | data.frame (1 row × p cols) | Feature values of the selected observation (model units) |
| `z_obs` | matrix (1 × 2) | Z-space coordinates, columns `"x"` and `"y"` |
| `pred_prob` | numeric scalar | Predicted default probability |
| `pred_class` | integer (0 or 1) | Predicted class at cutoff 0.5 |
| `row_id` | integer | Row index into `data` (`test_data` by default; 13 here); `NA_integer_` for a data-frame target |
| `class_obs` | integer (0, 1 or `NA`) | True class of the same row (`test_data[13, "class"]` here); `NA_integer_` when the source has no class column |
| `scaling` | list or `NULL` | Display-only raw-unit scaling copied from `bl_results_v2$scaling` (see `review_section4_to_6.md`) |

**`print(tgt)` console output** (the script does not call it directly; `print(bl_local)` in Step 11 prints this block; from the harness run):
```
-- bl_target --
  Row ID         : 13
  Predicted prob : 0.1110
  Predicted class: 0
  Feature values (raw units) :
          loan_amnt       loan_int_rate loan_percent_income        credit_score 
           12000.00               12.69              0.0500              634.00 
```

The values print in raw loan units because `bl_results_v2` carries a scaling; `tgt$x_obs` itself stays in z-score units. They print in fixed notation with 2 decimals. Inside `print(bl_local)`, as here, the block takes the comparison table's decimals per feature, so `loan_percent_income` prints as 0.0500; `print(tgt)` on its own prints 0.05.

---

### 10b — `set_filters()` → `flt`

**File:** `R/local_cf.R`

**Call:**
```r
flt <- set_filters(
  tgt,
  loan_amnt     = "fixed",      # make not actionable
  loan_int_rate = "increase",   # set by the lender, not the applicant
  credit_score  = "increase"
)
```

The constraints name only the four retained features, because `set_filters()` stops on a variable name the target does not have. The script keeps one further constraint commented out, `loan_percent_income = "fixed"`.

**What it does:**

1. Validates `tgt` is a `"bl_target"` object.
2. Validates all named arguments are feature names in `tgt$x_obs`.
3. Validates each constraint value is one of `"decrease"`, `"increase"`, `"fixed"`, or a `c(min, max)` numeric pair.
4. Stores everything.

**Constraint semantics (applied in X-space after back-projecting contour points):**

| Value | Applied as | Effect on search |
|---|---|---|
| `"decrease"` | `B_x[nm] <= x_obs[nm]` | Counterfactual value must be at most the observed value |
| `"increase"` | `B_x[nm] >= x_obs[nm]` | Counterfactual value must be at least the observed value |
| `"fixed"` | `abs(B_x[nm] - x_obs[nm]) <= 0.5` | Counterfactual value within ±0.5 of observed; in sparse CF, always reverts to observed |
| `c(min, max)` | `B_x[nm] >= min & B_x[nm] <= max` | Counterfactual must lie within absolute range |

No constraint = unconstrained (all contour points accepted for that feature).

The comparisons run on the stored values, which are in model units. For `bl_results_v2` those are z-scores, so the `"fixed"` tolerance of ±0.5 means ±0.5 training SD of that feature, and a `c(min, max)` range must be given in z-score units. For `loan_amnt`, whose training SD is 6,002, the tolerance is ±3,001 around the observed 12,000. A `"fixed"` feature can therefore still move in the search. `bl_find_local_cf()` then tries its observed value and keeps that point when it still crosses the boundary (Step 11, Stage H), and the sparse CF (Step 14) always holds it exactly.

**Returned object — `flt` (S3 class `"bl_filters"`):**

| Field | Type | Content |
|---|---|---|
| `constraints` | named list | Each element is a constraint value (`"decrease"`, `"increase"`, `"fixed"`, or `c(min, max)`) |
| `bl_target` | `bl_target` object | The `tgt` passed in — carried through so `bl_find_sparse_cf()` can access it |

**`print(flt)` console output** (from the harness run):
```
-- bl_filters --
  loan_amnt            : fixed
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

**Parameters:**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `bl_target` | required | `tgt` | The observation to explain (Step 10a) |
| `set_filters` | `NULL` | `flt` | Actionability constraints (Step 10b); `NULL` searches without constraints (Step 15) |
| `max_pairs` | `21L` | default | Eigenvector pairs to try, capped at p(p−1)/2: 6 for the 4-feature model |
| `m` | `200L` | default | Grid resolution per pair |
| `distance` | `"mahalanobis"` | default | Cross-pair selector: squared Mahalanobis distance in X-space, or Z-space distance under `"euclidean"` |
| `verbose` | `TRUE` | default | Prints progress messages |

Full argument list: `?bl_find_local_cf`.

**Why this function exists:** The global `bl_find_boundary()` uses the original biplot projection (fixed `proj_dims = c(1, 2)`). For a single target, we can do better: for each eigenvector pair, rotate the biplot so the target observation is aligned with that pair's plane, then search across the pairs to find whichever rotation brings the boundary closest. This local rotation maximises the chance of finding a valid counterfactual.

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

**Generate pairs:** `.generate_pairs(p, max_pairs)` produces up to `max_pairs` eigenvector index pairs in canonical order: (1,2), (1,3), (2,3), (1,4), (2,4), (3,4), (1,5), … With p = 4 there are only 6 pairs, (1,2) to (3,4), so all 6 are tried. Each pair defines a 2D projection plane from the full p-dimensional loading matrix.

**For each pair `(i, j)`, the loop does:**

#### Stage A — SVD rotation via `.bl_rotate()`

**File:** `R/local_cf.R` (private helper)

The rotation aligns the biplot so the target observation is aligned with the plane of the chosen eigenvector pair. This is done via SVD:

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

**Why:** By constructing Y from `[-target, 0, target]` and finding the rotation that maps the pair-space projection of Y onto its full-space projection, the target observation is carried into the plotted plane of the rotated biplot. This places the target at a known position and makes boundary search more reliable. The rotation fixes the target's direction only; the plane's second direction is not pinned down by the construction (see `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md` §4.1).

**Key invariant:** `YVr_padded` is non-zero only in its first two columns, so the SVD rotation `A` always concentrates the target's information in columns 1 and 2 of `Vrho = V %*% t(A)` -- regardless of which `proj_pair` was used to build `Vr`. The return therefore always selects `Vrho[, c(1, 2)]`, not `Vrho[, proj_pair]`; the `proj_pair` columns of `Vrho` are near zero for any pair other than `c(1, 2)`.

#### Stage B — Build m×m grid in rotated Z-space

```
Z_target    = x_target_st %*% Vr_rot       # 1 × 2
Z_train_rot = X_train_st %*% Vr_rot        # n × 2

z_range = range(Z_train_rot)
z_pad   = diff(z_range) * 0.10             # 10% padding
xseq    = seq(z_range[1]-pad, z_range[2]+pad, length.out = m)
Zgrid   = expand.grid(xseq, xseq)          # m² × 2
```

The 10% padding is needed because decision boundary contours frequently run near the *edge* of the training data cloud — that is where the class transition occurs. Without padding the grid would be flush with the outermost training observations, and contour lines that fall at or beyond the data range would be clipped, leaving no valid contour vertex for that pair. Adding 10% of the total range on each side gives boundary segments room to exist slightly outside the training envelope. The global grid built in `bl_build_grid()` takes its range differently, from the biplot's own plotting window (`.bl_get_plot_bounds()`; see `review_section4_to_6.md` Section 6.2).

#### Stage C — Back-project grid and score

```
Xgrid = Zgrid %*% tVr_rot + X_center      # m² × p in original feature space
grid_prob = .pred_function(model, model_type, Xgrid)   # "XGB" for bl_results_v2; chunked, 50 000 rows at a time
```

#### Stage D — Extract contour lines

```
ct_local = contourLines(xseq, yseq, matrix(grid_prob, ncol=m),
                         levels = c(cutoff - b_margin, cutoff + b_margin))
```
`b_margin = 0.01`, read from `bl_results_v2$b_margin` (set in its `bl_build_result()` call; `bl_find_local_cf()` has no `b_margin` argument), so contours at 0.49 and 0.51.

#### Stage E — Filter contour segments (four successive filters)

For each contour segment in `ct_local`:

1. **Opposing-class filter:** Keep only segments where `seg_class != target_pred_class`. A class-0 target (predicted safe) needs class-1 contours, where it would flip to "default"; a class-1 target needs class-0 contours.

2. **Back-project + actionability:** Back-project the segment to X-space, then apply the `set_filters` constraints via `.apply_actionability()` (when supplied). These are usually the most selective filter, so they run first, and a segment with no surviving vertex is skipped at once:
   - `"decrease"`: `B_x[nm] <= x_obs[nm]`
   - `"increase"`: `B_x[nm] >= x_obs[nm]`
   - `"fixed"`: `abs(B_x[nm] - x_obs[nm]) <= 0.5`
   - `c(min, max)`: range check

3. **Feasibility:** apply the `train_ranges` filter to the vertices that survived step 2.

4. **Model consistency:** Re-score the surviving back-projected points; keep only rows where prediction agrees with contour's probability level.

If all segments are eliminated, the pair is skipped.

#### Stage F — Find nearest surviving point

All surviving contour rows are combined into `all_Mi`. The nearest to `Z_target` is found via `.nearest_idx_block(Z_target, all_Mi, block=1)`. It returns the exact nearest vertex (the first on an exact tie) and uses no random numbers, so repeating the call, or calling other searches first, gives the same counterfactual (`documentation/nearest-vertex-tie-break.md`). Distance = `sqrt(sum((Z_target - B_z_local)^2))`.

#### Stage G — Update best result

For every pair (whether or not it wins), the back-projection from `B_z_local` to `B_x` is computed and the **squared Mahalanobis distance in X-space** is calculated:

```
v        = B_x - x_obs                              # 1 x p in X-space
dist_mah = as.numeric(t(v) %*% metric_inv %*% v)    # scalar squared Mahalanobis
```

where `metric_inv = W^{-1}` is the inverse of the pooled within-class covariance matrix stored on `bl_result$metric_inv`. The Cholesky decomposition is used for the inversion -- see `2 implementation_summary.txt` Section 4.2.1 for the maths and singularity safeguards.

The pair selector then depends on the `distance` argument of `bl_find_local_cf()`:

- `distance = "mahalanobis"` (default): pair selected by minimum `dist_mah`. This is the PhD's `pc.distz <- v %*% Wmat_inv %*% t(v)` criterion.
- `distance = "euclidean"`: pair selected by minimum `dist_z`.

**Under CVA, as here, the two are the same criterion.**

- **The identity:** biplotEZ's CVA loadings satisfy `V'WV = I/(n-g)`, and every rotated basis inherits it, so for the search's in-plane point `dist_mah = (n - g) * dist_z^2` for every pair.
- **On `bl_results_v2`:** the ratio `all_distances_mahalanobis / all_distances^2` is 15992 (n = 15996, g = 4) for every pair. That holds in the Step 11 search, Step 15's unconstrained search and Step 16's external search. The two options pick the same pair.

They differ only under PCA, where `dist_z` is standardised Euclidean distance and `W` is the binary-class covariance. The derivation is in `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md` §3.4.

If the chosen selector for this pair beats the running `best_selector`, save everything for this pair as the new best result. The loop continues through all configured pairs, keeping the global minimum. Both `dist_z` and `dist_mah` are reported for the winning pair, plus per-pair tables `all_distances` and `all_distances_mahalanobis`.

#### Stage H — Restore fixed features (after the loop)

**Why a fixed feature can move.** The search accepts a `"fixed"` feature within ±0.5 model units of its observed value (`.fixed_tol`, Stage E). A vertex of the 2D boundary contour back-projects to a point whose features are all set by its two coordinates, so a vertex rarely holds a feature at exactly its observed value. An exact equality would usually leave no candidate.

**The restore step.** Once the best pair is chosen, the private `.restore_fixed()` checks every `"fixed"` feature the search moved:
1. It sets those features back to their observed values and re-scores the point.
2. The restored point is used only if it still crosses the boundary (its predicted class differs from the target's) and lies within `train_ranges`. The other constraints still hold, because only the fixed features change.
3. If used, the restored point replaces `B_x` and `B_pred`. `B_z` becomes its projection, `x_st %*% Vr_rot`, the same projection `bl_select_target()` uses. The point satisfies the fixed constraint exactly but lies off the 2D boundary, so it need not be the nearest valid point.
4. `fixed_check` records the outcome either way, including the search's own point (`search_point`). `dist_z`, `dist_mahalanobis` and the per-pair tables stay the search's values, because they describe the pair choice.

In the harness run the search moves `loan_amnt` from 12,000 to 12,398.5 (inside its ±3,001 tolerance). Setting it back to 12,000 still gives 0.547, because the model is insensitive to that change here, so the restored point is the counterfactual.

---

**Returned object — `bl_local` (S3 class `"bl_local_result"`):**

| Field | Type | Content |
|---|---|---|
| `B_z` | matrix (1 × 2) | The counterfactual in rotated Z-space: the nearest boundary point, or the restored point's projection (Stage H) |
| `B_x` | data.frame (1 row × p) | The counterfactual in X-space: the back-projected boundary point, or the restored point |
| `B_pred` | numeric scalar | Model probability at the counterfactual (≈ 0.5 on the boundary) |
| `dist_z` | numeric scalar | Z-space distance from target to boundary, for the winning pair |
| `dist_mahalanobis` | numeric scalar | Squared Mahalanobis distance in X-space, for the winning pair |
| `distance` | character | The selector used (`"mahalanobis"` here) |
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
| `all_distances_mahalanobis` | named numeric | Squared Mahalanobis distance per pair tried (NA if no solution) |
| `solution_found` | logical | `TRUE` if at least one valid boundary was found |
| `blocking_constraint` | character or NULL | Reason message if no solution found |
| `bl_target` | `bl_target` | The `tgt` object |
| `bl_counterfactual` | `bl_counterfactual` or NULL | The boundary point, packaged like `bl_target` (see above) |
| `bl_result` | `bl_result` | Reference to `bl_results_v2` |
| `set_filters` | `bl_filters` or NULL | The `flt` object (needed by `bl_find_sparse_cf()`) |
| `fixed_check` | list or NULL | Stage H outcome: `vars` (fixed features the search moved), `applied`, `pred_restored`, `in_range`, `search_point` (the search's `B_x`, `B_z`, `B_pred`). `NULL` when no fixed feature moved |

**`print(bl_local)` console output** (from the harness run):
```
-- bl_local_result --
  Solution found : TRUE
  Best pair      : (3, 4)
  Selector       : mahalanobis
  Distance (Z)   : 0.0046
  Mahalanobis    : 0.5838 (sqrt of squared d_M^2 = 0.3408)
  Boundary pred  : 0.5470
  Fixed restored : loan_amnt

  Target:
-- bl_target --
  Row ID         : 13
  Predicted prob : 0.1110
  Predicted class: 0
  Feature values (raw units) :
          loan_amnt       loan_int_rate loan_percent_income        credit_score 
           12000.00               12.69              0.0500              634.00 

  Counterfactual:
-- bl_counterfactual --
  Boundary prob  : 0.5470
  Predicted class: 1
  Feature values (raw units) :
          loan_amnt       loan_int_rate loan_percent_income        credit_score 
           12000.00               14.15              0.0484              634.22 

  Per-pair distances:
  pair      dist_z dist_mahalanobis
 (1,2) 0.028309256       12.8162127
 (1,3) 0.027807716       12.3661193
 (2,3) 0.028326778       12.8320826
 (1,4) 0.015557768        3.8707700
 (2,4) 0.027625664       12.2047315
 (3,4) 0.004616263        0.3407877

  Observed vs counterfactual (raw units):
                    observed counterfactual  change constraint
loan_amnt           12000.00       12000.00    0.00      fixed
loan_int_rate          12.69          14.15    1.46   increase
loan_percent_income   0.0500         0.0484 -0.0016       none
credit_score          634.00         634.22    0.22   increase

  Note: loan_amnt is "fixed", but the boundary search moved it to 12398.5032
        (+398.5032). The search lets a "fixed" feature move by up to 0.5
        model units (+/- 3001.1777 here), because a point on the 2D decision
        boundary rarely holds a feature at exactly its observed value.
        Restoring loan_amnt to its observed value still crosses the boundary
        (prediction 0.5470), so the restored point is the selected
        counterfactual.
```

Before printing, `bl_find_local_cf()` also reports the restore as a message: "Fixed feature(s) loan_amnt restored to observed value(s); the restored point still crosses the boundary (pred 0.547) and is used as the counterfactual."

Pair (3, 4) has the smallest squared Mahalanobis distance (0.3408), so it wins. Its boundary point moves the loan amount from 12,000 to 12,398.5, inside the ±3,001 `"fixed"` tolerance (Step 10b). Setting the loan amount back to 12,000 still gives 0.547, so that restored point is the counterfactual (Stage H). It raises the interest rate from 12.69 to 14.15 and the credit score from 634 to 634.22, which satisfies both `"increase"` constraints, and keeps the loan amount at 12,000. The model's probability rises from 0.111 to 0.547. The distances and the per-pair table describe the search's boundary point, which decided the pair. The search is deterministic, so repeating the call gives the same values (`documentation/nearest-vertex-tie-break.md`).

**The comparison table at the end of `print(bl_local)`.** One row per feature, with:
- **Columns:** the observed value, the counterfactual, the change (counterfactual minus observed) and the Step 10 constraint (`"none"` when unconstrained; a range prints as `[lo, hi]`).
- **Units:** raw, because `bl_results_v2` carries a scaling. Levels convert as `value * scale + center` and the change as `value * scale`, by `.scale_to_raw()`, so the counterfactual column matches the Counterfactual block above. Without a scaling the header reads "(model units)".
- **Format:** fixed notation, never scientific, with 2 decimals. A row whose change would print as 0.00 gets more decimals, up to 4: here `loan_percent_income`, whose change is −0.0016. The Target and Counterfactual blocks above use the same decimals per feature. `print(bl_shapley_values)` and `print(bl_sparse)` use the same rule for their Shapley values. The sparse path's predictions keep 4 decimals only within 0.005 of the cutoff.
- **Note:** printed only when the search moved a `"fixed"` feature. It explains the tolerance and says whether the restored point was used. When it was not, because the restored point does not cross the boundary or leaves `train_ranges`, the counterfactual keeps the moved value and the note gives the reason.
- **No table:** when no solution exists, `print()` shows the Reason line instead, and Steps 13 and 14 stop, because `bl_shapley()` needs a solution.

The private helper `.cf_compare_table()` builds this table and the Step 14 sparse table, so the two share one layout. The counterfactual is stored in `bl_local$B_x` (z-scores) and `bl_local$bl_counterfactual` (prints in raw units).

---

## Step 12 — `plot(bl_local)` → `plot.bl_local_result()`

**File:** `R/local_cf.R`

**Call:**
```r
plot(bl_local,
     label_dir         = "Paral",
     label_offset_var  = c("loan_amnt", "loan_int_rate",
                           "loan_percent_income", "credit_score"),
     label_offset_dist = c(0, 0, 1, 0),
     ticks_var         = c("loan_amnt", "loan_int_rate",
                           "loan_percent_income", "credit_score"),
     ticks_n           = c(4, 5, 5, 5))
```

**Parameters:**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `plot_points` | `FALSE` | default | Training points are hidden, so the plot shows the target and its counterfactual. `plot_biplotEZ()` defaults to `TRUE` |
| `label_dir` | `"Paral"` | `"Paral"` | Border-adaptive label direction (the default, set explicitly) |
| `label_offset_var` | `0L` | the four feature names | Labels to move outward, by name or index |
| `label_offset_dist` | `1.5` | `c(0, 0, 1, 0)` | Paired by position: only `loan_percent_income` moves, by 1 |
| `ticks_var` | `0L` | the four feature names | Axes with their own tick count |
| `ticks_n` | `5L` | `c(4, 5, 5, 5)` | Paired by position: `loan_amnt` gets 4 ticks, the others 5 |
| `no_grid`, `no_contour` | `FALSE` | default | `TRUE` hides the prediction grid or the contour lines |
| `show_arrows` | `TRUE` | default | Draws the arrow from the target to its counterfactual |
| `print_summary` | `TRUE` | default | Prints the console summary (layer 8) |

Full argument list: `?plot.bl_local_result`. The pairing rule is in `review_section4_to_6.md` Section 6.4.

**No `rotate_deg` on the local plots.** `rotate_deg` is only for viewing the global biplot (`plot_biplotEZ()`, `plot(bl_surr)`). Whenever the rotation function `.bl_rotate()` is used -- here, the local biplot for an individual target, drawn by `plot(bl_local)` and `plot(bl_sparse)` -- the plot is not turned further: `rotate_deg` is in effect always 0. Those two plot methods have no `rotate_deg` argument and stop with an error, before any device opens, if one is passed. The local biplot's orientation is fixed by the optimal rotation, which is part of the method; `rotate_deg` is never used to plot a rotated individual point.

`new_title` (default `NA`) overrides the auto-generated `"Local biplot -- target N [pair (i,j) | p=..., class ...]"` title; `NA` keeps the auto-generated title.

**Raw-unit axis relabel:** if `bl_result$scaling` is set (via `bl_set_scaling()` — see
`review_section4_to_6.md`), the rotated local biplot's axes are relabelled into original feature
units. After the SVD rotation has patched `Z`/`Lmat`/`ax.one.unit`,
`plot.bl_local_result()` calls `.bl_rescale_biplot_axes(biplot_plot, bl_result$scaling, var_names)`.
The relabel transforms only the biplotEZ `$means`/`$sd`, which are **rotation-invariant**, so it
composes cleanly with the local rotation and leaves the plotted geometry untouched. The same
applies to `plot.bl_surrogate()` and, by pass-through, `plot.bl_sparse_result()` (Step 14).
In the console output, distances (`dist_z`, `dist_mahalanobis`) and Shapley contributions stay in
model units, while feature values (target, counterfactual, sparse CF) print in raw units.

**Base-graphics title-error note:** `new_title` is written straight to `biplot_plot$Title` (no validation, matching `plot_biplotEZ()`), and biplotEZ forwards that field to `graphics::title(main = ...)` at plot-flush time. A length-1 character string (or `NA` to keep the default) is the intended contract. `title()` is tolerant: it coerces most atomic and even recursive inputs (numeric, multi-element character vectors, lists, data.frames) to a label via `as.character()`, so those do *not* error — they just produce an odd/recycled title. The error case is a value base R cannot coerce to a character vector at all — a **function/closure** or an **environment** — which raises `"cannot coerce type 'closure' to vector of type 'character'"` from *inside* the biplotEZ `plot()` call, not from the plot function's own argument handling. A length-1 character is the safe, intended input.

**What it renders (8 layers):**

| Layer | Content | Data source |
|---|---|---|
| 1 | Grey axes + variable labels | `biplot_obj` patched with `Vr_rot` and `Z_train_rot` |
| 2 | Prediction grid (blue→red) | `bl_local$Zgrid` + `bl_local$col_value` |
| 3 | Training points (optional, hidden by default) | `Z_train_rot` + confusion colours |
| 4 | Axes redrawn darker (via `par(new=TRUE)`) | Same biplotEZ object |
| 5 | Decision boundary contour lines | `bl_local$ct_local` |
| 6 | Target circle (filled, confusion colour) | `Z_target`; colour from `bl_target$class_obs` against `bl_target$pred_class` (see below) |
| 7 | CF cross (×) + arrow | `B_z` (cross); arrow from `Z_target` → `B_z` |
| 8 | Console summary | Printed to console (not plotted) |

**Target colour (layer 6):** the colour shows whether the model is right about the target. It compares the target's true class `bl_target$class_obs`, resolved by `bl_select_target()` from the target's own row, with its predicted class at the cutoff:

| `class_obs` | `pred_class` | Group | Colour |
|---|---|---|---|
| 1 | 1 | TP | red |
| 0 | 0 | TN | blue |
| 0 | 1 | FP | purple |
| 1 | 0 | FN | orange |
| `NA` | 1 / 0 | unknown | red / blue (prediction only) |

These are the same colours the training points use in layer 3.
- **Unknown class:** an external target with no class column is coloured by its prediction alone.
- **No `class_obs` field:** a `bl_target` without the field is treated the same way.
- **The label:** the target's text label is `row_id`, or nothing for a data-frame target.

**Key rendering detail — patching the biplotEZ object:**

The biplotEZ object stored in `bl_results_v2$biplot_obj` was built in the original (unrotated) CVA space. For the local plot, it is patched in-place:
```r
biplot_plot$Z[, proj_dims]    <- Z_train_rot   # replace training coords with rotated
biplot_plot$Lmat[, proj_dims] <- Vr_rot        # replace loadings with rotated
biplot_plot$ax.one.unit       <- (1 / diag(t(tVr_rot) %*% tVr_rot)) * t(tVr_rot)
```
This means the axes, variable label positions, and tick marks are all recalculated in the rotated coordinate frame — the biplot shows the locally rotated view, not the original CVA projection.

**Console output (Step 8 of plot; from the harness run):**
```
--- Local CF summary ---
  Target         : 13
  Pred prob      : 0.1110  (class 0)
  Best pair      : (3, 4)
  Selector       : mahalanobis
  Distance (Z)   : 0.0046
  Mahalanobis    : 0.5838 (sqrt)
  Boundary pred  : 0.5470
------------------------
```

---

## Step 13 — `bl_shapley()` + `print()` + `plot()`

**File:** `R/shapley.R`

**Call:**
```r
bl_shapley_values <- bl_shapley(bl_local)
```

**Parameters:**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `exact_max_vars` | `14L` | default | Exact Shapley values when p ≤ this value, a permutation estimate above it. p = 4 here, so exact |
| `approx_perm` | `2048L` | default | Permutations for the estimate; unused here |
| `seed` | `NULL` | default | Seeds the permutation estimate (it calls `set.seed()`, which resets the session RNG); the exact method ignores it |

Full argument list: `?bl_shapley`.

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
- "Start": the observed feature values `x_obs` (row 13 of the test data here)
- "End": the counterfactual `B_x` (the boundary point back-projected to X-space)
- The model's probability change from `pred_obs` to `B_pred` must be attributed across features

**Method selection** (the `if (p <= exact_max_vars)` branch of `bl_shapley()`):
- If `p <= exact_max_vars` (default 14): use **exact** computation via `.shapley_exact_one()`
- If `p > 14`: use **approximate** via `.shapley_perm_one()` with `M = 2048` permutations

The reduced loan model has 4 features, so exact computation over all 2^4 = 16 subsets is used here, and the `seed` argument has no effect (it applies only to the permutation path).

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

**Contribute classification** (the `df$Contribute <- dplyr::case_when(...)` step of `bl_shapley()`):

Each feature is classified as `"Supports"`, `"Contradicts"` or `"No effect"` on the move to the boundary. The first row is checked first:

| `pred_class` | `shapley_cause` | Label |
|---|---|---|
| either | 0 (`abs(.) < 1e-8`: no effect on the prediction) | `"No effect"` |
| 0 (predicted safe) | > 0 (pushes prob up, toward boundary) | `"Supports"` |
| 0 (predicted safe) | < 0 (pushes prob down, away from boundary) | `"Contradicts"` |
| 1 (predicted default) | < 0 (pushes prob down, toward boundary) | `"Supports"` |
| 1 (predicted default) | > 0 (pushes prob up, away from boundary) | `"Contradicts"` |

A Shapley value of exactly 0 arises when a feature's change never moves the prediction, in any coalition. Here that happens because the change crosses no XGBoost split, and predictions are floor-rounded to 3 d.p. The PhD research script counted 0 as `"Supports"`. The `"No effect"` class keeps such a feature out of the sparse counterfactual (Step 14). See `documentation/sparse-counterfactual.md` §1.

**Row ordering:** rows are sorted by signed `shapley_cause`: ascending for a class-0 target, descending for a class-1 target. The plot draws the first row at the bottom, so for either class the strongest supporter of the move to the boundary sits at the top.

**`varnames_p` label** (used as Y-axis in plot), stored in model units:
```
"feature: observed_value -> change_to_reach_boundary"
```
The observed value is rounded to 3 d.p. and the change to 2 d.p. For `bl_results_v2` these are z-score values. When `bl_result$scaling` is set (as here), `plot.bl_shapley()` replaces the label at plot time with raw units, `"feature: observed -> counterfactual"`, for example `"loan_int_rate: 12.69 -> 14.15"`. The stored `varnames_p` is not changed.

---

**Returned object — `bl_shapley_values` (S3 class `"bl_shapley"`):**

| Field | Type | Content |
|---|---|---|
| `shapley_df` | data.frame | One row per feature, sorted by signed Shapley value. Columns: `varnames`, `pred_data` (observed value), `data_to_boundary` (change needed), `shapley_cause` (Shapley value), `pred_use` (target's predicted class), `Contribute`, `varnames_p` (formatted label) |
| `pred_prob` | numeric | Target's predicted probability |
| `pred_class` | integer | Target's predicted class |
| `pred_boundary` | numeric | Model probability at counterfactual |
| `row_id` | integer | Row index (or `NA` for external) |
| `bl_local_result` | `bl_local_result` | Reference to `bl_local` |

---

**`print(bl_shapley_values)` console output** (from the harness run):
```
-- bl_shapley --
  Row ID         : 13
  Observed pred  : 0.1110  (class 0)
  Full CF pred   : 0.5470
  Shapley sum    : 0.4360

  Observed vs counterfactual, with Shapley values (raw units):
                    observed counterfactual  change constraint shapley
loan_amnt           12000.00       12000.00    0.00      fixed    0.00
loan_int_rate          12.69          14.15    1.46   increase    0.29
loan_percent_income   0.0500         0.0484 -0.0016       none    0.14
credit_score          634.00         634.22    0.22   increase    0.00
                    contribute
loan_amnt            No effect
loan_int_rate         Supports
loan_percent_income   Supports
credit_score         No effect
```

The layout matches the tables of `print(bl_local)` (Step 11) and `print(bl_sparse)` (Step 14), with two extra columns:
- **Shared builder:** the private `.cf_compare_table()` builds all three tables. Its rows follow `var_names`; the plot keeps its Shapley order.
- **Units:** feature values and the change print in raw units (`.scale_to_raw()`, levels and a delta). `shapley` is in probability units and is not converted.
- **`Shapley sum`:** the four values add up to 0.436, the move from 0.111 to 0.547.
- **This run:** at the 80-column width of the harness log, R prints `contribute` as a second block; a wider console shows it on one line. The loan amount does not change (restored in Step 11), and the credit score's 0.22-point rise does not move the prediction. Both are therefore `"No effect"`; the interest rate and the loan-to-income ratio are `"Supports"`.

---

**`plot(bl_shapley_values)` — Shapley bar chart:**

Horizontal bar chart (ggplot2):
- **Y-axis:** one bar per feature, sorted by signed Shapley value. With a scaling set (as here) the labels read `"feature: observed -> counterfactual"` in raw units and the axis title is `"Variable (observed -> counterfactual, raw units)"`; without one they are `varnames_p`, `"feature: observed -> change"` in model units
- **X-axis:** `shapley_cause` — the Shapley value (how much this feature contributed to the prediction move)
- **Fill colour:** `"Supports"` = dark blue; `"Contradicts"` = light grey (`grey70`); `"No effect"` = mid grey (`grey45`)
- **Value label:** Shapley value rounded to 3 d.p., displayed inside each bar in white
- **Title:** `"Shapley Contribution Plot"`
- **Subtitle:** `"ID: 13 | Class: 0 | Pred: 0.111 | Boundary pred: 0.547"`
- **X-axis label:** `"[-> change to reach boundary]: Impact on prediction"`

**Reading the chart:** Each bar shows how much that feature's change (from observed to counterfactual) accounts for the total probability shift from `pred_prob` to `pred_boundary`. The Shapley values sum to approximately `pred_boundary - pred_prob`.

---

## Step 14 — `bl_find_sparse_cf()` + `print()` + `plot()`

**File:** `R/shapley.R`

**Call:**
```r
bl_sparse <- bl_find_sparse_cf(bl_shapley_values)
```

**Parameters:**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `round_to` | `NULL` | default | Rounds each changed feature's value to the nearest multiple of `round_to`. It acts on model units, which are z-scores for `bl_results_v2`, so the script leaves values unrounded (item 4 below) |

Full argument list: `?bl_find_sparse_cf`.

### What `bl_find_sparse_cf()` does

Builds a **sparse counterfactual**: it changes as few features as it needs to cross the boundary and leaves the rest at their observed values. It answers "which changes are enough to flip the prediction?"

**Algorithm** (inside `bl_find_sparse_cf()`):

1. **Candidates:** the `"Supports"` features that are not `"fixed"` in `set_filters`, ordered by absolute Shapley value, largest first. `"Contradicts"`, `"No effect"`, `"Unknown"` and fixed features are never changed. A fixed feature is not actionable; it took part in the boundary search only to inform the direction.
2. **Greedy path:** start from `x_sparse = x_obs`. For each candidate in turn:
   - set `x_sparse[nm] = B_x[nm]`;
   - re-score the point with the model;
   - record the step in `sparse_path`;
   - stop as soon as the predicted class differs from the target's.

   The candidates not reached stay at their observed values.
3. **Why re-score and not sum Shapley values:** a Shapley value is a feature's average effect over all coalitions, not its effect when it is added first. Re-scoring after each step uses the model's actual prediction, and it costs at most one model call per candidate.
4. **Optional rounding** (`round_to`): if not `NULL`, each changed feature's value is rounded to the nearest multiple of `round_to`. The rounding acts on the stored model-unit value. On raw data, `round_to = 500` rounds `loan_amnt` to the nearest 500. For `bl_results_v2`, whose values are z-scores, it rounds the z-score instead, so the script leaves values unrounded. See `documentation/sparse-counterfactual.md` §3, which is open for future investigation.
5. **Validate:** `solution_valid = (pred_class_sparse != pred_class_target)`. If no step crosses, every candidate has been changed and `solution_valid` is `FALSE`.

The selection follows the Shapley order, so it is not guaranteed to be the smallest possible subset.

**Returned object — `bl_sparse` (S3 class `"bl_sparse_result"`):**

| Field | Type | Content |
|---|---|---|
| `x_sparse` | data.frame (1 row × p) | The sparse counterfactual in X-space |
| `pred_sparse` | numeric scalar | Model probability at sparse counterfactual |
| `solution_valid` | logical | `TRUE` if sparse CF flips predicted class |
| `shapley_df` | data.frame | Annotated Shapley table with extra columns `used_in_sparse` (logical: `TRUE` when the sparse CF changes the feature, i.e. it is on the greedy path) and `x_sparse_val` (the value used in sparse CF) |
| `sparse_path` | data.frame | One row per greedy step: `step`, `variable`, `shapley`, `pred` (model probability after the step), `crosses`. Zero rows when there is no candidate |
| `bl_shapley` | `bl_shapley` | Reference to `bl_shapley_values` |

---

**`print(bl_sparse)` console output** (from the harness run):
```
-- bl_sparse_result --
  Solution valid : TRUE
  Observed pred  : 0.1110  (class 0)
  Full CF pred   : 0.5470
  Sparse CF pred : 0.5470

  Observed vs sparse counterfactual (raw units):
                    observed counterfactual   sparse  change constraint  used
loan_amnt           12000.00       12000.00 12000.00    0.00      fixed FALSE
loan_int_rate          12.69          14.15    14.15    1.46   increase  TRUE
loan_percent_income   0.0500         0.0484   0.0484 -0.0016       none  TRUE
credit_score          634.00         634.22   634.00    0.00   increase FALSE
  (change = sparse - observed)

  Path (Supports features, largest Shapley value first):
 step            variable shapley pred crosses
    1       loan_int_rate    0.29 0.36   FALSE
    2 loan_percent_income    0.14 0.55    TRUE
```

The table has the same layout as the one at the end of `print(bl_local)` (Step 11):
- **Columns:** `counterfactual` is the full counterfactual `B_x`, `sparse` is `x_sparse`, and `change` is the move the sparse CF asks for (sparse minus observed).
- **`used`:** the stored `used_in_sparse`, which is `TRUE` when the sparse CF changes the feature (it is on the greedy path).
- **Row order:** follows `var_names`, as in Step 11.
- **The path:** printed below the table, from `bl_sparse$sparse_path`.

Here the two `"Supports"` features are tried in Shapley order. Raising the interest rate alone moves the prediction from 0.111 to 0.360, which is not enough. Adding the lower loan-to-income ratio reaches 0.547 and crosses the boundary, so both are needed. The loan amount (fixed, restored in Step 11) and the credit score (`"No effect"`) stay at their observed values. The sparse CF flips the prediction with two changes instead of four.

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
   - `z_sparse` is projected with `Vr_rot`, the same frame as the rest of the local plot. `plot(bl_sparse)` accepts no `rotate_deg` (see the note after the `plot(bl_local)` Parameters table).

4. **Adds a legend** (top-right corner):
   - `"Full CF"` (grey cross, from `plot.bl_local_result`)
   - `"Sparse CF (valid)"` or `"Sparse CF (invalid)"` (green or yellow cross)

5. **Prints sparse CF summary** to console:
   ```
   --- Sparse CF summary ---
     Target         : 13
     Observed pred  : 0.1110  (class 0)
     Full CF pred   : 0.5470
     Sparse CF pred : 0.5470
     Solution valid : TRUE
   -------------------------
   ```

---

## Step 15 (optional) — Unconstrained local search

```r
bl_local_free  <- bl_find_local_cf(bl_results_v2, tgt)
print(bl_local_free)
bl_shap_free   <- bl_shapley(bl_local_free)
bl_sparse_free <- bl_find_sparse_cf(bl_shap_free)
print(bl_sparse_free)
```

(The script also plots `bl_local_free`, `bl_shap_free` and `bl_sparse_free`.)

Identical pipeline to Steps 11–14, with `set_filters = NULL`. The unconstrained search finds the geometrically nearest boundary point without any actionability restrictions. Comparing `bl_local` (constrained) vs `bl_local_free` (unconstrained) shows the cost of the constraints: how much farther the boundary is when `loan_amnt` is fixed and `loan_int_rate` and `credit_score` may only increase. In the harness run the unconstrained search also picks pair (3, 4), at a slightly smaller distance (√d²_M = 0.5554 against 0.5838), with boundary probability 0.542. Its counterfactual raises the loan amount to 12,544.2 and lowers the credit score to 633.43, two moves the constraints rule out.

**`print(bl_sparse_free)` console output** (from the harness run):
```
-- bl_sparse_result --
  Solution valid : TRUE
  Observed pred  : 0.1110  (class 0)
  Full CF pred   : 0.5420
  Sparse CF pred : 0.5420

  Observed vs sparse counterfactual (raw units):
                    observed counterfactual   sparse change constraint  used
loan_amnt           12000.00       12544.23 12000.00   0.00       none FALSE
loan_int_rate          12.69          14.11    14.11   1.42       none  TRUE
loan_percent_income     0.05           0.04     0.04  -0.01       none  TRUE
credit_score          634.00         633.43   634.00   0.00       none FALSE
  (change = sparse - observed)

  Path (Supports features, largest Shapley value first):
 step            variable shapley pred crosses
    1       loan_int_rate    0.29 0.36   FALSE
    2 loan_percent_income    0.14 0.54    TRUE
```

Without constraints the sparse counterfactual makes the same two changes as the constrained one (Step 14). It raises the interest rate (to 14.11) and lowers the loan-to-income ratio (to 0.0399). The unconstrained moves in the loan amount and the credit score do not move the prediction, so both are `"No effect"` and stay at their observed values. The interest rate alone again reaches only 0.360, so both changes are needed.

---

## Complete Object Flow: Steps 10–14

```
bl_results_v2  [bl_result]          ← Phase 1 anchor (reduced 4-feature model)
test_pts_v2    [bl_points]          ← test data in Z-space
     │
     │  bl_predict(bl_results_v2)
     ▼
pred_summary   [data.frame]
  ├── row, pred_prob, pred_class
  ├── true_class, confusion (TP/TN/FP/FN)
  └── all feature columns
     │
     │  bl_select_target(bl_results_v2, target = 13)
     │    → extract test_data[13, var_names]
     │    → project to Z-space
     │    → score through XGBoost
     │    → read true class from test_data[13, "class"]
     ▼
tgt  [bl_target]
  ├── x_obs    (1 × p feature values, model units)
  ├── z_obs    (1 × 2, Z-space position)
  ├── pred_prob, pred_class   (0.111, 0)
  ├── row_id = 13
  ├── class_obs (true class of that row; NA if unknown)
  └── scaling  (display-only z-score scaling)
     │
     │  set_filters(tgt, loan_amnt="fixed", loan_int_rate="increase", credit_score="increase")
     ▼
flt  [bl_filters]
  ├── constraints = list(loan_amnt="fixed", loan_int_rate="increase", credit_score="increase")
  └── bl_target = tgt
     │
     │  bl_find_local_cf(bl_results_v2, tgt, set_filters=flt)   # max_pairs = 21
     │
     │  [for each of the 6 eigenvector pairs (p = 4)]:
     │    .bl_rotate() → Vr_rot, tVr_rot (SVD-based local rotation)
     │    build m×m rotated grid → Xgrid (back-projected)
     │    score through XGBoost → grid_prob
     │    extract contour lines → ct_local
     │    filter: opposing class, set_filters, train_ranges, consistency
     │    .nearest_idx_block() → nearest boundary vertex
     │    update best if dist_mah improves (the default selector)
     │  [after the loop]: .restore_fixed() → loan_amnt back to 12,000, still 0.547 → used
     ▼
bl_local  [bl_local_result]
  ├── B_z       (1 × 2, counterfactual in rotated Z; restored point's projection)
  ├── B_x       (1 × p, counterfactual in X-space; loan_amnt restored)
  ├── B_pred    (0.547)
  ├── fixed_check (vars = loan_amnt, applied = TRUE, search_point)
  ├── dist_z, dist_mahalanobis, distance
  ├── Z_target, best_pair = (3, 4), Vr_rot, tVr_rot
  ├── Z_train_rot, Zgrid, grid_prob, col_value, ct_local
  ├── all_distances, all_distances_mahalanobis, solution_found
  ├── bl_target = tgt
  ├── bl_counterfactual  (the boundary point, printed in raw units)
  ├── bl_result = bl_results_v2
  └── set_filters = flt
     │
     ├── print(bl_local)
     │     └─→ summary, target, counterfactual, per-pair distances; then the observed /
     │         counterfactual / change / constraint table and the fixed-restore note
     │
     ├── plot(bl_local)
     │     └─→ patched local biplot (rotated axes + grid + target circle + CF cross)
     │
     │  bl_shapley(bl_local)
     │    → start = x_obs, end = B_x
     │    → exact Shapley (p=4 ≤ 14): enumerate all 2^4 = 16 subsets
     │    → classify each feature as Supports / Contradicts / No effect
     │    → sort by signed Shapley value
     ▼
bl_shapley_values  [bl_shapley]
  ├── shapley_df   (4 rows: varnames, pred_data, data_to_boundary,
  │                 shapley_cause, pred_use, Contribute, varnames_p)
  ├── pred_prob, pred_class, pred_boundary
  ├── row_id = 13
  └── bl_local_result = bl_local
     │
     ├── print(bl_shapley_values) → observed / counterfactual / change / constraint / shapley / contribute
     ├── plot(bl_shapley_values)  → horizontal bar chart coloured by Contribute
     │
     │  bl_find_sparse_cf(bl_shapley_values)
     │    → candidates: Supports, not fixed, largest |Shapley| first
     │      (loan_int_rate, loan_percent_income)
     │    → x_sparse = x_obs; change one candidate at a time and re-score
     │    → loan_int_rate: 0.360 (no cross); + loan_percent_income: 0.547 (crosses) → stop
     │    → solution_valid = (pred_class(x_sparse) != pred_class(x_obs))
     ▼
bl_sparse  [bl_sparse_result]
  ├── x_sparse       (1 × p sparse counterfactual)
  ├── pred_sparse    (model prob at sparse CF)
  ├── solution_valid (TRUE/FALSE)
  ├── shapley_df     (+ used_in_sparse + x_sparse_val columns)
  ├── sparse_path    (one row per greedy step)
  └── bl_shapley = bl_shapley_values
     │
     ├── print(bl_sparse) → predictions; observed / counterfactual / sparse / change / constraint / used
     └── plot(bl_sparse)  → local biplot + second cross (green=valid, yellow=invalid)
                             + console: Observed/Full CF/Sparse CF predictions
```

---

## Key Methodological Points

1. **Local rotation vs global boundary search.** `bl_find_boundary()` (Step 7) searches for the nearest boundary across the full test set using a fixed CVA projection. `bl_find_local_cf()` (Step 11) rotates the projection specifically around the target observation and tries multiple eigenvector pairs — this local alignment finds a closer, more actionable boundary than the global search typically can.

2. **The convex hull is NOT applied in Phase 3.** Because each eigenvector pair rotation invalidates the original grid polygon, `bl_find_local_cf()` uses `train_ranges` (per-feature min/max) and `set_filters` actionability constraints instead of the Z-space hull. The hull was the right constraint in global space; in locally-rotated space it would be incorrect.

3. **Shapley values attribute the full CF movement.** The Shapley values sum to approximately `pred_boundary - pred_prob` — the total prediction change from observation to boundary. Each Shapley value is that feature's fair share of causing the model to change its mind.

4. **Sparse CF is a greedy subset selection.**
   - **The rule:** it changes the `"Supports"` features one at a time, largest Shapley value first, re-scores the model after each, and stops at the first point that crosses the boundary. `"Contradicts"`, `"No effect"` and fixed features never change.
   - **Not guaranteed minimal:** optimal sparse CF is NP-hard, but this is a principled, Shapley-ordered approximation.
   - **If `solution_valid = FALSE`:** even all the supporting features together do not flip the prediction. Reverting the contradicting and no-effect features to their observed values pulled the probability back below 0.5.

5. **`"fixed"` constraint has three effects.** It models features the analyst knows cannot or should not change (for example a person's age).
   - **Search:** `bl_find_local_cf()` limits the search to values within ±0.5 of observed, in model units (±0.5 SD for the z-scored `bl_results_v2`).
   - **Restore (Stage H):** it then sets each fixed feature the search moved back to its observed value, and uses that point when it still crosses the boundary.
   - **Sparse CF:** `bl_find_sparse_cf()` always reverts it to observed, regardless of Shapley classification.

   The script's Step 10 fixes `loan_amnt`. The search moves it by 398.5, inside the ±3,001 tolerance, and the restore puts it back to 12,000 at the same prediction, 0.547.

6. **The chain is fully self-contained.** Each object carries references to its parents: `bl_sparse$bl_shapley$bl_local_result$bl_result`. Every downstream method can re-extract whatever it needs (projection matrices, model, feature names) without requiring the user to pass additional arguments.
