# Code Walkthrough: Steps 3–6 of `03_loan_status_Boundary_Logic.R`

## Context

This document traces the full data, object, and value flow through Steps 3–6 of the loan-default script. It assumes `loan_encoded` is already in memory: a data.frame of ~45,000 rows × 9 numeric columns (all categoricals have been integer-encoded; `loan_intent`, `person_education`, `person_home_ownership` already dropped in Step 1).

---

## Step 3 — Domain filter (base R, no package functions)

```r
loan_filtered <- loan_encoded[loan_encoded$previous_loan_defaults_on_file == 0, ]

feature_cols <- setdiff(
  names(loan_filtered),
  c("loan_status", "previous_loan_defaults_on_file")
)
```

**What happens:**

| Operation | Input | Output |
|---|---|---|
| Row filter | `loan_encoded` (~45 k rows) | `loan_filtered` — only applicants with no prior default (~27 k rows) |
| `setdiff()` | All column names minus the outcome + the now-constant default column | `feature_cols` — character vector of predictor names |

**Why `previous_loan_defaults_on_file` is dropped from features:**
After filtering to only `== 0`, that column is constant — zero variance. Feeding a constant feature into PCA/CVA or an ML model is numerically meaningless.

**Values at end of Step 3:**
- `loan_filtered`: data.frame, ~27 k rows × 9 cols (includes `loan_status` and `previous_loan_defaults_on_file`)
- `feature_cols`: character vector, e.g. `c("person_age", "person_gender", "person_income", "person_emp_exp", "person_home_ownership_encoded", "loan_amnt", "loan_int_rate", "loan_percent_income", "cb_person_cred_hist_length", "credit_score")`  
  *(exact names depend on what survived Step 1 pruning)*

---

## Steps 4–6 — One combined code block

The block runs 7 operations in sequence. Each produces one named object used by the next.

---

### Alternative entry point: `bl_wrap_data()`

If your data is already split (e.g. from a cross-validation framework or external pipeline),
you can bypass `bl_prepare_data()` and supply pre-split frames directly. The commented-out
block at the top of the Steps 4-6 section in script 03 shows the full pattern:

```r
# Rename outcome column to "class" and split manually
data_clean <- loan_filtered[, c(feature_cols, "loan_status")]
names(data_clean)[names(data_clean) == "loan_status"] <- "class"
set.seed(121L)
n          <- nrow(data_clean)
train_idx  <- sample(n, size = floor(0.8 * n), replace = FALSE)
train_data <- data_clean[ train_idx, , drop = FALSE]
test_data  <- data_clean[-train_idx, , drop = FALSE]
rownames(train_data) <- NULL; rownames(test_data) <- NULL
bl_dat <- bl_wrap_data(
  train_data   = train_data,
  test_data    = test_data,
  var_names    = feature_cols,
  target_class = NULL       # loan_status already encoded as 0/1
)
# bl_wrap_data() does not filter outliers; call bl_filter_outliers() if needed:
# bl_dat <- bl_filter_outliers(bl_dat, hull_fraction = 0.9)
```

`bl_wrap_data()` returns `"bl_data"` (not `"bl_filter_result"`), but both classes are accepted
by all downstream functions. All Steps 5-6 and Phases 2-3 are identical regardless of which
entry point was used.

---

### 4 — `bl_prepare_data()` → `bl_dat`

**File:** `R/data_prepare.R`

**Call:**
```r
bl_dat <- bl_prepare_data(
  data           = loan_filtered,
  class_col      = "loan_status",
  feature_cols   = feature_cols,
  train_fraction = 0.8,
  seed           = 121L,
  hull_fraction  = 0.9
)
```

**What it does, step by step:**

1. **Validates** inputs (all via `stop_if_*` helpers in `utils.R`):
   - `data` is a data.frame
   - `class_col` exists in data
   - `train_fraction` is strictly in (0, 1)

2. **Selects feature columns** — uses the supplied `feature_cols` vector; validates all columns exist.

3. **Binarises the outcome**: `loan_status` is already 0/1, so it passes through as-is and is renamed to `"class"` in the output.

4. **Assembles clean data**: a data.frame with only `feature_cols` + `"class"`.

5. **Train/test split** with `set.seed(121L)`:
   - Randomly samples `floor(0.8 × n)` row indices → `train_data`
   - Remaining rows → `test_data`
   - Row names are cleared from both.

6. **Calls `bl_filter_outliers()` internally** on the interim `"bl_data"` object:
   - Standardises training features, runs lightweight 2D PCA, builds convex hull.
   - Removes training rows outside the hull; test set passed through unchanged.
   - Prints one-line summary (if `verbose = TRUE`).

**Returned object — `bl_dat` (S3 class `"bl_filter_result"`):**

| Field | Type | Content |
|---|---|---|
| `train_data` | data.frame | ~21 600 rows × (p features + `"class"`) |
| `test_data` | data.frame | ~5 400 rows × same columns |
| `var_names` | character vector | Predictor column names (copy of `feature_cols`) |
| `num_vars` | integer | Length of `var_names` |
| `target_class` | NULL | NULL because `loan_status` was pre-coded 0/1 |

---

### 4b (integrated)

Outlier filtering is now integrated into  via . See Step 4 above for the full mechanics.

---

### 5 — Model fitting → `bl_mod`

**Files:** `R/model_fit.R` (dispatcher) + `R/model_utils.R` (internal fitter) + `R/predict_utils.R` (prediction)

`bl_fit_model()` supports four parsnip/tidymodels types directly: GLM, SVM, NNET, RForrest. For XGBoost (used in this script) the model must be fitted externally and registered via `bl_wrap_model()`. Both paths return the same `"bl_model"` S3 object.

---

#### 5a — Direct fit: SVM via `bl_fit_model()`

**Call (Step 5a in script, quick baseline):**
```r
bl_mod_svm <- bl_fit_model(
  train_data = bl_dat$train_data,
  var_names  = bl_dat$var_names,
  model_type = "SVM"
)
```

**What it does, step by step:**

1. **Validates** inputs; checks `model_type` is one of GLM/SVM/NNET/RForrest.

2. **Calls `.fit_model()`** (`R/model_utils.R`):
   - Looks up `.default_model_params[["SVM"]]` — empty list (SVM uses parsnip defaults).
   - Builds parsnip model spec:
     ```r
     parsnip::svm_rbf() |>
       parsnip::set_engine("kernlab") |>
       parsnip::set_mode("classification")
     ```
   - Wraps in a workflow:
     ```r
     workflows::workflow() |>
       workflows::add_formula(class ~ var1 + var2 + ...) |>
       workflows::add_model(model_spec)
     ```
   - Calls `parsnip::fit(wf, data = train_tm)` where `train_tm$class` is a factor.
   - Returns `list(model = <workflows::workflow>, model_type = "SVM")`.

3. **Scores training data** via `.pred_function()` (`R/predict_utils.R`):
   - Workflow path: `predict(model, new_data=train_data, type="prob")$.pred_1`
   - Floor-rounds to 3 d.p.: `floor(prob * 1000) / 1000`
   - Binarises at cutoff 0.5: `pred_class = as.numeric(pred_prob >= 0.5)`

4. **Computes training metrics:**
   - `accuracy = mean(pred_class == actual)`
   - `gini = 2 * AUC - 1` via `calc_gini()` (`utils.R`)

**`bl_mod_svm` fields:** `model` = `workflows::workflow`; `model_type` = `"SVM"`.

| Use `bl_fit_model()` when... | Use `bl_wrap_model()` when... |
|---|---|
| GLM, SVM, NNET, RForrest with default or simple params | XGB, GBM, GAM, LDA, or any other type |
| Quick exploration / baseline | Custom nrounds, watchlists, regularisation |
| | Model already fitted externally |

---

#### 5b — Wrap path: XGB via `bl_wrap_model()` with explicit `predict_fn`

**Call (Step 5b in script, used for full analysis):**
```r
xgb_data <- xgboost::xgb.DMatrix(
  data  = as.matrix(bl_dat$train_data[, bl_dat$var_names]),
  label = bl_dat$train_data$class
)
xgb_fit <- xgboost::xgb.train(
  params  = list(objective     = "binary:logistic",
                 eval_metric   = "logloss",
                 max_depth     = 3,
                 learning_rate = 0.1),
  data    = xgb_data,
  nrounds = 200,
  verbose = 0
)
bl_mod <- bl_wrap_model(
  model      = xgb_fit,
  model_type = "custom",
  var_names  = bl_dat$var_names,
  predict_fn = function(m, new_data) {
    mat <- xgboost::xgb.DMatrix(as.matrix(new_data))
    as.numeric(predict(m, newdata = mat))
  },
  train_data = bl_dat$train_data
)
```

**What it does, step by step:**

1. `bl_wrap_model()` validates inputs. Because `model_type = "custom"`, it bundles the raw `xgb_fit` object and the user-supplied `predict_fn` together:
   ```r
   model_use <- list(model = xgb_fit, predict_fn = predict_fn)
   ```

2. Calls `.pred_function(model_use, "custom", train_data[, var_names])`:
   - Dispatches to the `custom` branch: `predict_fn(model_use$model, new_data)`
   - The lambda calls `xgb.DMatrix(as.matrix(new_data))` then `predict(m, newdata = mat)`
   - Floor-rounds to 3 d.p.

3. Computes `accuracy` and `gini` from training data (same as `bl_fit_model`).

> **Alternative (direct XGB path):** Instead of `model_type = "custom"` with a `predict_fn`, you can pass the raw booster with feature metadata and let the built-in XGB dispatch handle prediction:
> ```r
> bl_mod <- bl_wrap_model(
>   model      = list(model = xgb_fit, features = bl_dat$var_names),
>   model_type = "XGB",
>   var_names  = bl_dat$var_names,
>   train_data = bl_dat$train_data
> )
> ```
> The `custom` + `predict_fn` approach is shown here because it makes the prediction contract explicit and generalises to any model type not natively supported.

**`bl_mod` fields:**

| Field | Type | Content |
|---|---|---|
| `model` | list | `list(model = <xgb.Booster>, predict_fn = <function>)` |
| `model_type` | character | `"custom"` |
| `var_names` | character vector | Predictor names (copy of `bl_dat$var_names`) |
| `cutoff` | numeric | `0.5` |
| `accuracy` | numeric | Training accuracy, e.g. `0.91` |
| `gini` | numeric | Training Gini coefficient, e.g. `0.85` |

---

### 6 — `bl_build_result()` → `bl_results`

**File:** `R/result.R` (orchestrator) → calls `bl_build_projection()`, `bl_build_grid()`, `bl_assemble()`

**Call:**
```r
bl_results <- bl_build_result(
  bl_data  = bl_dat,
  bl_model = bl_mod,
  method   = "CVA",
  title    = "Loan default (prior defaulters) — XGB, CVA biplot",
  b_margin = 0.01
)
```

This is the **Phase 1 anchor function**. It orchestrates three sub-steps and returns the central `bl_result` object that everything in Phases 2 and 3 consumes.

---

#### Sub-step 6a — `bl_build_projection()` → `bl_proj`

**File:** `R/projection.R`

**What it does:**

1. Extracts feature matrix `X = train_data[, var_names]`.
2. Computes `X_center = colMeans(X)` and `X_sd = colSDs(X)`.
3. **CVA forces `standardise = FALSE`** internally — CVA works with class-structured covariance, not z-scores.
4. **Computes confusion labels** (because CVA needs class groupings):
   - Calls `.pred_function()` on training data with `bl_mod`
   - Labels each row: TP / TN / FP / FN as a 4-level factor
5. Builds the **biplotEZ object**:
   ```r
   biplotEZ::biplot(X, scaled = FALSE, Title = title) |>
     biplotEZ::CVA(classes = confusion_labels, e.vects = c(1L, 2L))
   ```
   Internally, biplotEZ solves the generalised eigenvalue problem for the between/within class covariance to find the loading matrix `V` that maximally separates the confusion classes.
6. Extracts `V = bp$Lmat` (p × p loading matrix) and computes `tV = solve(V)` (the inverse, for back-projection).
7. Calls `get_variable_ranges()` (`R/feasibility_utils.R`) to record `c(min, max)` per feature — used later to filter test data to plausible ranges.

**`bl_proj` fields (not returned to user, passed to next sub-step):**
`V`, `tV`, `X_center`, `X_sd`, `method`, `standardise`, `proj_dims`, `biplot_obj`, `train_ranges`

---

#### Sub-step 6b — `bl_build_grid()` → `bl_grid`

**File:** `R/biplot_grid.R`

**What it does:**

1. Extracts the 2-column submatrices: `Vr = V[, 1:2]` (p×2) and `tVr = tV[1:2, ]` (2×p).

2. **Projects training data to Z-space:**
   ```
   Z_train = X_centered %*% Vr    (n × 2)
   ```
   (No scaling because CVA uses `standardise = FALSE`)

3. **Gets plot bounds** by rendering the biplot off-screen, reading `par("usr")` — gives the square `[min_val, max_val]` of the visible plot area.

4. **Builds the m×m prediction grid** in Z-space:
   ```
   xseq = seq(min_val, max_val, length.out = m)   # m = 200 default
   Zgrid = expand.grid(xseq, xseq)                # 40 000 × 2
   ```

5. **Back-projects grid to X-space:**
   ```
   Xgrid = Zgrid %*% tVr              # 40 000 × p
   Xgrid = Xgrid + X_center           # add back training means
   ```
   This is the key operation that lets us ask: "if a point is at Z-space coordinates (a, b), what feature values does it correspond to?"

6. **Computes convex hull** of `Z_train` using `aplpack::plothulls` (at `outlie=1`, so all training points included — the `hull_fraction=0.9` used in `bl_filter_outliers` was already applied to the training data).

7. **Scores all 40 000 grid points** through the XGBoost model in 50 000-row chunks:
   ```
   grid_prob = .pred_function(bl_mod$model, "XGB", Xgrid)
   ```
   Each grid probability is floor-rounded to 3 d.p.

8. **Assigns colours** for the grid: a 101-colour ramp (blue → white → red) indexed by `floor(grid_prob * 100) + 1`.

9. **Extracts contour lines** at `cutoff ± b_margin`:
   - `b_margin = 0.01` (set directly via the `b_margin` parameter)
   - `ct`: contours at `[0.49, 0.51]` — used for boundary search in Phase 2
   - `ct_surrogate`: hull-clipped contours at `[0.49, 0.50, 0.51]` — used only by `bl_surrogate()`

**`bl_grid` fields (not returned to user directly, assembled into `bl_results`):**
`Zgrid`, `Xgrid`, `grid_prob`, `col_value`, `min_val`, `max_val`, `polygon`, `hull_fraction`, `ct`, `ct_surrogate`, `xseq`, `yseq`, `b_margin`

---

#### Sub-step 6c — `bl_assemble()` → `bl_results`

**File:** `R/result.R`

Combines all artifacts into the final `bl_result` object.

**Returned object — `bl_results` (S3 class `c("bl_result", "list")`):**

| Category | Field | Content |
|---|---|---|
| **Data** | `train_data` | Filtered training data (from `bl_dat`) |
| | `test_data` | Test data filtered to training variable ranges |
| | `var_names` | Predictor names |
| | `num_vars` | Number of predictors |
| **Model** | `model` | Fitted XGBoost workflow |
| | `model_type` | `"XGB"` |
| | `cutoff` | `0.5` |
| **Projection** | `V` | p×p CVA loading matrix |
| | `tV` | Inverse of V (for back-projection) |
| | `X_center` | Training column means |
| | `X_sd` | Training column SDs |
| | `method` | `"CVA"` |
| | `standardise` | `FALSE` (CVA always) |
| | `proj_dims` | `c(1L, 2L)` |
| | `biplot_obj` | biplotEZ S3 object (used for axis rendering) |
| **Feasibility** | `train_ranges` | Named list of `c(min, max)` per feature |
| | `polygon` | Convex hull of Z-space training points (`sp::SpatialPolygons`) |
| | `hull_fraction` | `1` (grid-level hull fraction) |
| **Grid** | `biplot_grid` | Full `bl_grid` list (Zgrid, Xgrid, grid_prob, ct, ct_surrogate, ...) |
| **Performance** | `accuracy` | Training accuracy |
| | `gini` | Training Gini |
| | `b_margin` | `0.01` |
| **Metadata** | `call` | The `bl_build_result()` call expression |
| | `created_at` | Timestamp |

---

### 7 — `plot_biplotEZ(bl_results)` — renders the training biplot

**File:** `R/plot_biplot.R`

**Call:**
```r
plot_biplotEZ(
  bl_results,
  label_dir         = "Hor",
  label_offset_var  = 0L,
  label_offset_dist = 0.5
)
```

**What it renders (in layer order):**

| Layer | Source data | Visual element |
|---|---|---|
| 1 | `bl_results$biplot_obj` + `V[, proj_dims]` | Grey coordinate axes + variable labels |
| 2 | `biplot_grid$Zgrid` + `col_value` | 40 000 coloured squares (blue→red probability surface) |
| 3 | `points$Z` + `points$pred_col` | Training obs as coloured dots: TP=red, TN=blue, FP=purple, FN=orange |
| 4 | `biplot_obj` again | Darker axes re-drawn on top of the grid and points |
| 5 | `biplot_grid$ct` | Black contour lines marking the decision boundary (probability ≈ 0.5) |

**`points` here:** Because no `points=` argument is passed, `plot_biplotEZ()` automatically calls `bl_project_points(bl_results$train_data, bl_results)` internally to project the training data.

**Returns:** `bl_results` invisibly (for pipe-compatibility). Side effect: renders plot to active device.

---

### 8 — `bl_project_points()` → `test_pts`

**File:** `R/project_points.R`

**Call:**
```r
test_pts <- bl_project_points(bl_results$test_data, bl_results)
```

**What it does:**

1. Extracts from `bl_results`: `V`, `X_center`, `X_sd`, `standardise` (FALSE for CVA), `cutoff`, `proj_dims`, `polygon`.

2. **Projects test data to Z-space** using the **same** loading matrix and centering as training:
   ```
   X_centered = X_test - X_center          # subtract training means
   # (no scaling — CVA, standardise = FALSE)
   Z = X_centered %*% V[, proj_dims]       # n_test × 2
   ```
   This is critical: using training means ensures test points are placed in the *same coordinate frame* as training points.

3. **Checks polygon membership** — tests whether each test point falls inside the training convex hull using `sp::over()`.

4. **Scores through the model:**
   ```
   pred_prob  = .pred_function(bl_results$model, "XGB", X_test[, var_names])
   pred_class = as.numeric(pred_prob >= 0.5)
   ```

5. **Assigns confusion colours**: TP=red, TN=blue, FP=purple, FN=orange (test set has a `"class"` column, so true labels are known).

**Returned object — `test_pts` (S3 class `"bl_points"`):**

| Field | Type | Content |
|---|---|---|
| `Z` | numeric matrix (n_test × 2) | Z-space coordinates, columns `"x"` and `"y"` |
| `pred_prob` | numeric vector | Predicted probability of default per test obs |
| `pred_class` | numeric vector (0/1) | Predicted class at cutoff 0.5 |
| `pred_col` | character vector | Point colour per obs (TP/TN/FP/FN) |
| `class` | numeric vector (0/1) | True loan status labels |
| `inside_polygon` | logical vector | TRUE if obs falls inside training hull |

**Where `bl_project_points()` is called from:**

Internal callers (in package `R/`):

| Caller | Location | Trigger |
|---|---|---|
| `bl_predict()` | `R/project_points.R` | Always — `bl_predict()` is a tabular wrapper that calls `bl_project_points()` and reshapes the result into a data frame |
| `plot_biplotEZ()` | `R/plot_biplot.R` | When the user does not pass `points = ...` — auto-projects `bl_result$train_data` so `plot(bl_results)` works without explicit projection |
| `bl_pick_point()` | `R/pick_point.R` | Once at startup, before the interactive click loop matches click coordinates against the projected `Z` |

External (user-facing) callers — the canonical pattern is to overlay non-training data on the biplot:

```r
test_pts <- bl_project_points(bl_results$test_data, bl_results)
plot_biplotEZ(bl_results, points = test_pts)
```

Used this way in all loan scripts (03–06), the iris and Pima scripts, both vignettes, and the README. Variants: `filter_to_polygon = TRUE` (only one explicit use, in script 03 line 185) drops out-of-hull observations before plotting; `filter_to_train_ranges = TRUE` (no current scripted use) drops X-space-extrapolated rows. Without an explicit `bl_project_points()` call, `plot_biplotEZ()` only ever shows training data because of the auto-project on line ~145 of `plot_biplot.R`.

---

### 9 — `plot_biplotEZ(bl_results, points = test_pts)` — overlays test data

**Same function, different arguments.** Renders identical layers 1–5 as above, but uses `test_pts$Z` and `test_pts$pred_col` for the data points (layer 3) instead of projecting training data automatically. This shows where the held-out test observations fall on the decision surface.

---

## Complete Object Flow Summary

```
loan_encoded (data.frame, ~45 k rows × 9 cols)
     │
     │  Step 3: row filter, setdiff()
     ▼
loan_filtered (data.frame, ~27 k rows × 9 cols)
feature_cols  (character vector, p predictor names)
     │
     │  bl_prepare_data(hull_fraction = 0.9)
     │  → split → bl_filter_outliers() → standardise → SVD → 2D hull → point-in-polygon
     ▼
bl_dat  [bl_filter_result]
  ├── train_data  (filtered, ~19-20 k rows)
  ├── test_data   (unchanged, 5 400 rows)
  └── var_names, polygon, hull_fraction, n_retained, n_removed
     │
     │  bl_wrap_model(model_type = "custom", predict_fn = ...)
     │  → external xgboost fit → accuracy + Gini on train
     ▼
bl_mod  [bl_model]
  ├── model       (list(model = xgb.Booster, predict_fn = fn))
  ├── model_type  "custom"
  ├── var_names
  ├── cutoff      0.5
  ├── accuracy    (e.g. 0.90)
  └── gini        (e.g. 0.84)
     │
     │  bl_build_result(method = "CVA", b_margin = 0.01)
     │  → bl_build_projection() → V, tV, biplot_obj
     │  → bl_build_grid()       → 200×200 grid, scored, contours
     │  → bl_assemble()         → combines all
     ▼
bl_results  [bl_result]              ← THE CENTRAL ANCHOR OBJECT
  ├── train_data, test_data, var_names, num_vars
  ├── model, model_type, cutoff
  ├── V, tV, X_center, X_sd, method="CVA", standardise=FALSE
  ├── biplot_obj       (biplotEZ — for axis rendering)
  ├── train_ranges     (min/max per feature)
  ├── polygon          (hull in CVA Z-space)
  ├── biplot_grid      (Zgrid, Xgrid, grid_prob, ct, ct_surrogate)
  ├── accuracy, gini, b_margin=0.01
  └── call, created_at
     │
     │  plot_biplotEZ(bl_results)   → renders training biplot
     │
     │  bl_project_points(bl_results$test_data, bl_results)
     │  → project using same V, X_center → Z_test
     │  → score through XGBoost → pred_prob, pred_col
     ▼
test_pts  [bl_points]
  ├── Z              (n_test × 2, CVA Z-space coordinates)
  ├── pred_prob      (predicted default probability)
  ├── pred_class     (0/1 at cutoff 0.5)
  ├── pred_col       (TP/TN/FP/FN colour)
  ├── class          (true loan_status)
  └── inside_polygon (TRUE/FALSE per obs)
     │
     │  plot_biplotEZ(bl_results, points = test_pts)
     └─→ renders same biplot with test points overlaid
```

---

## Key Design Principles to Note

1. **One anchor object:** `bl_results` carries everything. All Phase 2 and Phase 3 functions take `bl_results` as their first argument and extract what they need.

2. **Two-way street between X and Z:** `V` projects forward (X → Z) and `tV = solve(V)` projects backward (Z → X). The entire boundary and counterfactual machinery depends on this invertibility.

3. **Consistent coordinate frame:** `X_center` and `X_sd` from training are applied to every subsequent projection (test data, grid, target points) — so all Z-space coordinates are directly comparable.

4. **CVA vs PCA:** CVA (`method = "CVA"`) maximises separation between the four confusion categories (TP/TN/FP/FN), producing a biplot where the decision boundary is most clearly visible. PCA would instead maximise total variance, which may not align with the class boundary.

5. **`b_margin` controls the contour band half-width**, not the predictions themselves. `b_margin = 0.01` means the boundary search looks for contour lines at probability 0.49 and 0.51 (a 0.02-wide band total). All predicted values are always stored at 3 d.p. regardless.
