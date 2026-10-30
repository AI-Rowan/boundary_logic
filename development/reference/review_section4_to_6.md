# Code Walkthrough: Steps 3–6 of `03_loan_status_Boundary_Logic.R`

## Context

This document traces the full data, object, and value flow through Steps 3–6 of the loan-default script. It assumes `loan_encoded` is already in memory: a data.frame of 45,000 rows × 11 numeric columns (the raw file has 14 columns; Step 1 integer-encodes the categoricals and drops `loan_intent`, `person_education` and `person_home_ownership`).

Section numbers follow the script's own step headings. Step 6 covers `bl_build_result()` (6.1–6.3 are its internal sub-steps) and the three plotting and projection calls the script runs at the end of its Steps 4–6 block (6.4–6.6).

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
| Row filter | `loan_encoded` (45,000 rows) | `loan_filtered` — only applicants with no prior default (22,142 rows) |
| `setdiff()` | All column names minus the outcome + the now-constant default column | `feature_cols` — character vector of predictor names |

**Why `previous_loan_defaults_on_file` is dropped from features:**
After filtering to only `== 0`, that column is constant — zero variance. Feeding a constant feature into PCA/CVA or an ML model is numerically meaningless.

**Values at end of Step 3:**
- `loan_filtered`: data.frame, 22,142 rows × 11 cols (includes `loan_status` and `previous_loan_defaults_on_file`)
- `feature_cols`: character vector of the 9 predictors, `c("person_age", "person_gender", "person_income", "person_emp_exp", "loan_amnt", "loan_int_rate", "loan_percent_income", "cb_person_cred_hist_length", "credit_score")`

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

### Optional: raw-unit biplot axes via `bl_set_scaling()`

**File:** `R/data_prepare.R` (setter), `R/plot_biplot.R` (`.bl_rescale_biplot_axes()`),
`R/result.R` / `R/projection.R` (propagation).

When features are standardised *before* model fitting (common for scale-sensitive models),
the data the package holds is in standardised units, so biplotEZ labels every axis in those
units. In the loan script this is the v2 (reduced-model) flow: the `feature_cols_v2` columns
(chosen in the variable-pruning step; see `review_section7_to_8.md`) are standardised with `scale()`, the XGB model and projection are built on the standardised
data, and the transform is recorded so the axes display in raw loan units:

```r
# Split first on raw data (the hull filter self-standardises), then fit the
# scaling on the TRAINING rows only and apply those stats to the test split --
# computing the scaling on the full data would leak test information.
bl_dat_v2  <- bl_prepare_data(loan_filtered, "loan_status",
                              feature_cols = feature_cols_v2, ...)
train_scl  <- scale(bl_dat_v2$train_data[, feature_cols_v2])
scl_center <- attr(train_scl, "scaled:center")
scl_scale  <- attr(train_scl, "scaled:scale")
bl_dat_v2$train_data[, feature_cols_v2] <- train_scl
bl_dat_v2$test_data[,  feature_cols_v2] <-
  scale(bl_dat_v2$test_data[, feature_cols_v2],
        center = scl_center, scale = scl_scale)
bl_dat_v2  <- bl_set_scaling(bl_dat_v2, center = scl_center,
                             scale = scl_scale, method = "z-score")
```

`bl_set_scaling(x, center, scale, method)` accepts a `bl_data` or `bl_filter_result`, validates
and reorders the per-feature vectors to `var_names` (via private `.align_scaling_vec()`), and
stores `x$scaling = list(center, scale, method)`. It rejects zero/`NA`/non-finite scales and
names that do not cover every feature. The scaling passes through `bl_filter_outliers()` and is
copied by `bl_assemble()` into `bl_result$scaling`, including the model-free `bl_result` that
`bl_build_result()` returns for the exploratory biplot.

**It is display-only:** the model, grid, and `predict_fn` still operate on the standardised
data the package holds; only the plot methods consume `scaling`. At plot time each biplot
method calls `.bl_rescale_biplot_axes(biplot_obj, scaling, var_names)`, which transforms only
the biplotEZ `$means`/`$sd` (`means_new = means*scale + center`, `sd_new = sd*scale`; `$scaled`/
`$center` forced TRUE) and leaves `$Z`/`$Lmat`/`$ax.one.unit` untouched -- so axis labels read
in raw units while plotted positions are byte-identical. No-op when `scaling` is `NULL`. Note
this is distinct from `bl_result$X_center`/`X_sd` (the internal PCA standardisation) and the
`standardise` flag. See `2 implementation_summary.txt` §4.4.1 for the derivation, and the
`plot.bl_*` walkthroughs in `review_section9_to_15.md` for the rotated local/sparse cases.

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

**Parameters:**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `feature_cols` | `NULL` | `feature_cols` (9 names) | Predictor columns; `NULL` uses every column except `class_col` |
| `target_class` | `NULL` | default | The `class_col` value coded as `1`; `NULL` requires the column to be 0/1 already, as `loan_status` is |
| `train_fraction` | `0.8` | `0.8` | Share of rows sampled for training; `1` keeps every row, and `bl_assemble()` then sets `test_data` to `train_data` |
| `seed` | `121L` | `121L` | Fixes the train/test split, so every downstream count is reproducible |
| `hull_fraction` | `0.9` | `0.9` | Share of training rows the outlier hull keeps (15,970 of 17,713 here); `1` removes nothing |
| `verbose` | `TRUE` | default | Prints the one-line hull-filter summary |

Full argument list: `?bl_prepare_data`.

**What it does, step by step:**

1. **Validates** inputs (all via `stop_if_*` helpers in `utils.R`):
   - `data` is a data.frame
   - `class_col` exists in data
   - `train_fraction` is in (0, 1]; `1` keeps every row for training, and `bl_assemble()` then sets `test_data` to `train_data` (the script's exploratory Step 2 uses this)

2. **Selects feature columns** — uses the supplied `feature_cols` vector; validates all columns exist.

3. **Binarises the outcome**: `loan_status` is already 0/1, so it passes through as-is and is renamed to `"class"` in the output.

4. **Assembles clean data**: a data.frame with only `feature_cols` + `"class"`.

5. **Train/test split** with `set.seed(121L)`:
   - Randomly samples `floor(0.8 × n)` row indices → `train_data`
   - Remaining rows → `test_data`
   - Row names are cleared from both.

6. **Calls `bl_filter_outliers()` internally** (`R/outlier_filter.R`) on the interim `"bl_data"` object:
   - **Hull plane:** it standardises the training features and runs a lightweight 2D PCA. The outlier hull is drawn in the plane where the training cloud is widest.
   - **Building the hull:** `.build_hull_polygon()` (`R/hull_utils.R`) builds it.
     - `aplpack::plothulls(fraction = hull_fraction)` peels convex hulls from the outside in, and returns the hull that encloses roughly `hull_fraction` of the training points (here 0.9).
     - Training rows outside that hull are removed as outliers, so the biplot and the later boundary search are fitted to the bulk of the data.
     - The test set passes through unchanged.
     - At `hull_fraction = 1` the hull is still built but encloses every point, so nothing is removed.
   - **Graphics device:** `plothulls()` computes the hull by drawing it, so `.build_hull_polygon()` draws on a temporary off-screen `png()`. It closes only that device and then re-activates whichever device the caller had active. The user's open plots stay where they are, however many devices are open.
   - **Summary:** it prints a one-line summary when `verbose = TRUE`.

**Returned object — `bl_dat` (S3 class `"bl_filter_result"`):**

| Field | Type | Content |
|---|---|---|
| `train_data` | data.frame | 15,970 rows × (9 features + `"class"`): the 17,713-row training split less the 1,743 rows outside the 0.9 hull |
| `test_data` | data.frame | 4,429 rows × same columns (not filtered) |
| `var_names` | character vector | Predictor column names (copy of `feature_cols`) |
| `num_vars` | integer | Length of `var_names` (9) |
| `target_class` | NULL | NULL because `loan_status` was pre-coded 0/1 |
| `scaling` | NULL | Display-only raw-unit scaling; `NULL` until `bl_set_scaling()` is used |
| `polygon` | `SpatialPolygons` | The outlier-filter hull, in the lightweight standardised 2-PC plane |
| `hull_fraction` | numeric | `0.9` |
| `n_retained`, `n_removed` | integer | `15970`, `1743` |

---

### 5 — Model fitting → `bl_mod`

**Files:** `R/model_fit.R` (dispatcher) + `R/model_utils.R` (internal fitter) + `R/predict_utils.R` (prediction)

`bl_fit_model()` supports four parsnip/tidymodels types directly: GLM, SVM, NNET, RForest. For XGBoost (used in this script) the model must be fitted externally and registered via `bl_wrap_model()`. Both paths return the same `"bl_model"` S3 object.

---

#### 5a — Direct fit: SVM via `bl_fit_model()`

**Call (Step 5a in the script, a quick baseline; commented out there):**
```r
bl_mod_svm <- bl_fit_model(
  train_data = bl_dat$train_data,
  var_names  = bl_dat$var_names,
  model_type = "SVM"
)
```

**Parameters:**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `model_type` | `"GLM"` | `"SVM"` | One of GLM, SVM, NNET, RForest; selects the parsnip model spec in `.fit_model()` |
| `model_params` | `list()` | default | Hyperparameter overrides with parsnip names (e.g. NNET `hidden_units`); unknown names warn |
| `cutoff` | `0.5` | default | Threshold for the training accuracy; copied to `bl_result$cutoff`. Only `0.5` is valid under the current method |

Full argument list: `?bl_fit_model`.

**What it does, step by step:**

1. **Validates** inputs; checks `model_type` is one of GLM/SVM/NNET/RForest (the old spelling `"RForrest"` is accepted with a deprecation warning).

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
| GLM, SVM, NNET, RForest with default or simple params | XGB, GBM, GAM, LDA, or any other type |
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

**Parameters:**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `model_type` | required | `"custom"` | `"custom"` predicts through `predict_fn`; a built-in type such as `"XGB"` uses its own branch in `.pred_function()` |
| `predict_fn` | `NULL` | the `xgb.DMatrix` function | `function(model, new_data)` returning probabilities; required for `"custom"`, ignored otherwise |
| `train_data` | `NULL` | `bl_dat$train_data` | Rows scored for `accuracy` and `gini`; `NULL` leaves both `NA` |
| `cutoff` | `0.5` | default | Threshold for the training accuracy; copied to `bl_result$cutoff`. Only `0.5` is valid under the current method |

Full argument list: `?bl_wrap_model`.

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
> The `custom` + `predict_fn` approach is shown here because it makes the prediction contract explicit and generalises to any model type not natively supported. The script's reduced-model refit (Steps 4b–6b) uses this direct `"XGB"` path. The booster must come from `xgboost::xgb.train()`; models from `xgboost::xgboost()` are not supported by the XGB branch.

The script's Step 5b comments also document two ways to substitute a previously-saved booster instead of fitting `xgb_fit` fresh: **Option A** loads an XGBoost binary file via `xgboost::xgb.load()` (saved earlier with `xgboost::xgb.save()`), and **Option B** loads an R object via `readRDS()` (saved earlier with `saveRDS()`). Either option's result can replace the demo `xgb.train()` call; the rest of the `bl_wrap_model()` call (`predict_fn`, `var_names`, `train_data`) is unchanged.

**`bl_mod` fields:**

| Field | Type | Content |
|---|---|---|
| `model` | list | `list(model = <xgb.Booster>, predict_fn = <function>)` |
| `model_type` | character | `"custom"` |
| `var_names` | character vector | Predictor names (copy of `bl_dat$var_names`) |
| `cutoff` | numeric | `0.5` |
| `accuracy` | numeric | Training accuracy (`0.8280` in the harness run) |
| `gini` | numeric | Training Gini coefficient (`0.8093` in the harness run) |

---

### 6 — `bl_build_result()` → `bl_results`

**File:** `R/result.R` (orchestrator) → calls `bl_build_projection()`, `bl_build_grid()`, `bl_assemble()`

**Call:**
```r
bl_results <- bl_build_result(
  bl_data  = bl_dat,
  bl_model = bl_mod,
  method   = "CVA",
  title    = "Loan default (prior defaulters) -- XGB, CVA biplot",
  b_margin = 0.001
)
```

**Parameters:**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `method` | `"PCA"` | `"CVA"` | CVA separates the four confusion groups (TP/TN/FP/FN); PCA maximises total variance |
| `title` | `""` | `"Loan default ... CVA biplot"` | Biplot title; `plot(..., new_title = )` overrides it at plot time |
| `b_margin` | `0.001` | `0.001` | Half-width of the boundary band: contours at `0.5 ± b_margin` (6.2, item 9) |
| `m` | `200L` | default | Grid resolution: m × m = 40,000 scored points |
| `proj_dims` | `c(1L, 2L)` | default | Eigenvectors that span the biplot plane |
| `standardise` | `TRUE` | default | Ignored here: CVA forces `FALSE` |
| `cva_classes` | `NULL` | default | CVA groups; `NULL` with a model uses the TP/TN/FP/FN labels |
| `outlie` | `1` | default | Fraction of training points inside the grid hull; `1` encloses all of them |
| `calc_hull` | `TRUE` | default | Trims the stored grid to the hull (6.2, item 10); visual only |

`bl_build_result()` has no `cutoff` argument. The grid draws its contours around `bl_build_grid()`'s default cutoff, `0.5`, while `bl_results$cutoff` is copied from `bl_mod$cutoff`. The two agree here because `bl_mod` uses `0.5`.

Full argument list: `?bl_build_result`.

This is the **Phase 1 anchor function**. It orchestrates three sub-steps and returns the central `bl_result` object that everything in Phases 2 and 3 consumes.

---

#### 6.1 — `bl_build_projection()` → `bl_proj`

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
7. Calls `get_variable_ranges()` (`R/feasibility_utils.R`) to record `c(min, max)` per feature. These `train_ranges` are the feasibility check in both counterfactual searches, and the optional `filter_to_train_ranges` filter in `bl_project_points()`.
8. Computes the within-class metric `W` from the same four confusion groups, and its inverse by Cholesky (`.compute_metric_inverse()`); these are the Mahalanobis metric of Phases 2 and 3 (see `2 implementation_summary.txt` §4.2.1).

**`bl_proj` fields (not returned to user, passed to next sub-step):**
`V`, `tV`, `X_center`, `X_sd`, `method`, `standardise`, `proj_dims`, `biplot_obj`, `cva_classes` (the TP/TN/FP/FN factor), `point_col`, `train_ranges`, `metric`, `metric_inv`, `metric_type`, `condition_number`

---

#### 6.2 — `bl_build_grid()` → `bl_grid`

**File:** `R/biplot_grid.R`

**What it does:**

1. Extracts the 2-column submatrices: `Vr = V[, 1:2]` (p×2) and `tVr = tV[1:2, ]` (2×p).

2. **Projects training data to Z-space:**
   ```
   Z_train = X_centered %*% Vr    (n × 2)
   ```
   (No scaling because CVA uses `standardise = FALSE`)

3. **Gets plot bounds** with `.bl_get_plot_bounds()`.
   - **Why it draws:** biplotEZ fixes its axis extents only when it draws. The helper therefore renders the biplot axes on a temporary off-screen `png()` and reads `par("usr")`.
   - **The bounds:** the smaller and larger of the x and y limits give the square `[min_val, max_val]`. This square sets the grid's range, so the prediction grid covers exactly the area the biplot shows.
   - **Graphics device:** the helper closes only its temporary device and then re-activates the caller's device, so the user's open plots are untouched.
   - **Fallback:** if the render fails, it uses the range of `Z_train`, with a warning.

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

6. **Computes the convex hull** of `Z_train` with `.build_hull_polygon()`.
   - **No trimming here:** it runs `aplpack::plothulls()` at `outlie = 1`, so every training point is inside. Trimming already happened in Step 4 through `hull_fraction`.
   - **Where it is used:** this hull is built in the final biplot plane. It is stored as `bl_grid$polygon` and then `bl_results$polygon`, and Phase 2 clips the boundary contours to it. A global counterfactual therefore stays in the region the training data supports.
   - **Graphics device:** it draws on a temporary off-screen device, handled as in item 3.

7. **Scores all 40 000 grid points** through the XGBoost model in 50 000-row chunks:
   ```
   grid_prob = .pred_function(bl_mod$model, bl_mod$model_type, Xgrid)   # model_type "custom"
   ```
   The `"custom"` branch calls the stored `predict_fn`. Each grid probability is floor-rounded to 3 d.p., and the grid probabilities are clamped to `[0, 1]`.

8. **Assigns colours** for the grid: a 101-colour ramp (`deepskyblue` → `white` → `lightsalmon`) indexed by `floor(grid_prob * 100) + 1`.

9. **Extracts contour lines** at `cutoff ± b_margin`:
   - `b_margin = 0.001` (from the `bl_build_result()` call)
   - `ct`: contours at `0.499` and `0.501` — used for boundary search in Phase 2. The two levels bracket the 0.5 boundary as a thin band.
   - `ct_surrogate`: the same two levels, extracted after every grid cell outside the hull is set to `cutoff`, so each contour closes inside the hull. The `0.5` level is extracted too and then dropped. Used only by `bl_surrogate()`.

10. **Trims the stored grid to the hull** (`calc_hull = TRUE`): `Zgrid`, `Xgrid`, `grid_prob` and `col_value` keep only the points inside the hull, so the prediction surface is drawn in the region the training data covers. The contours above are extracted before this trim.

**`bl_grid` fields (not returned to user directly, assembled into `bl_results`):**
`Zgrid`, `Xgrid`, `grid_prob`, `col_value`, `min_val`, `max_val`, `polygon`, `hull_fraction`, `ct`, `ct_surrogate`, `xseq`, `yseq`, `b_margin`

---

#### 6.3 — `bl_assemble()` → `bl_results`

**File:** `R/result.R`

Combines all artifacts into the final `bl_result` object.

**Returned object — `bl_results` (S3 class `c("bl_result", "list")`):**

| Category | Field | Content |
|---|---|---|
| **Data** | `train_data` | Filtered training data (from `bl_dat`; 15,970 rows) |
| | `test_data` | Test data, passed through unfiltered (4,429 rows) |
| | `var_names` | Predictor names |
| | `num_vars` | Number of predictors (9) |
| **Model** | `model` | `list(model = <xgb.Booster>, predict_fn = <function>)` |
| | `model_type` | `"custom"` |
| | `cutoff` | `0.5` |
| | `b_margin` | `0.001` |
| **Projection** | `V` | p×p CVA loading matrix |
| | `tV` | Inverse of V (for back-projection) |
| | `X_center` | Training column means |
| | `X_sd` | Training column SDs |
| | `method` | `"CVA"` |
| | `standardise` | `FALSE` (CVA always) |
| | `proj_dims` | `c(1L, 2L)` |
| | `biplot_obj` | biplotEZ S3 object (used for axis rendering) |
| **Metric** | `metric` | Within-class covariance `W` (p×p) of the four confusion groups |
| | `metric_inv` | `W^{-1}` by Cholesky; the Mahalanobis metric of Phases 2 and 3 |
| | `metric_type` | `"W_cva"` |
| **Feasibility** | `train_ranges` | Named list of `c(min, max)` per feature |
| | `polygon` | Convex hull of Z-space training points (`sp::SpatialPolygons`) |
| | `hull_fraction` | `1` (the grid hull's `outlie`; printed as "Grid hull") |
| **Display** | `scaling` | `NULL` here (no `bl_set_scaling()` in the full-model flow) |
| **Grid** | `biplot_grid` | Full `bl_grid` list (Zgrid, Xgrid, grid_prob, ct, ct_surrogate, ...) |
| **Performance** | `accuracy` | Training accuracy (`0.8280`) |
| | `gini` | Training Gini (`0.8093`) |
| **Metadata** | `call` | The `bl_build_result()` call expression |
| | `created_at` | Timestamp |

---

### 6.4 — `plot_biplotEZ(bl_results)` / `plot(bl_results)` — renders the training biplot

**File:** `R/plot_biplot.R`

`plot(bl_results)` dispatches to `plot.bl_result()`, a one-line S3 wrapper that forwards all arguments to `plot_biplotEZ()`.

**Call:**
```r
plot(
  bl_results,
  label_dir         = "Hor",  # "Hor" = horizontal, "Orthog" = orthogonal to axis
  label_offset_var  = 0L,     # variable index/indices to shift, e.g. c(1L, 3L)
  label_offset_dist = 0.5     # outward distance per shifted label
)
```

**Parameters:**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `points` | `NULL` | default | `NULL` projects and draws `bl_results$train_data`; a `bl_points` object draws other rows (6.6) |
| `plot_points` | `TRUE` | default | `FALSE` skips the data-points layer |
| `label_dir` | `"Paral"` | `"Hor"` | `"Paral"` is border-adaptive (left/right labels vertical, top/bottom horizontal); `"Hor"` is all horizontal; `"Orthog"` is perpendicular to each axis |
| `label_offset_var` | `0L` | `0L` | Variables whose labels move outward, by name (`"loan_amnt"`) or index; `0L` moves none |
| `label_offset_dist` | `1.5` | `0.5` | Outward distance, in margin lines, for each listed variable; no effect here, as no variable is listed |
| `ticks_v` | `1L` | default | Tick marks on each axis not listed in `ticks_var` |
| `ticks_var`, `ticks_n` | `0L`, `5L` | default | Axes with their own tick count, and those counts; `0L` overrides none |
| `rotate_deg` | `0` | default | Rotates the whole plot clockwise; grid, points, contours and axes turn together |
| `new_title` | `NA` | default | Replaces the title stored by `bl_build_result()` |

Full argument list: `?plot_biplotEZ`.

**Pairing rule for offsets and ticks.** `label_offset_dist` pairs with `label_offset_var` position by position, and `ticks_n` pairs with `ticks_var` in the same way. The private helpers `.make_label_line_vec()` and `.make_ticks_vec()` (`R/plot_biplot.R`) build these into one value per axis. They are shared by `plot_biplotEZ()`, `plot.bl_local_result()` and `plot.bl_surrogate()`.
- **A single value** applies to every listed variable.
- **A name the model does not have,** such as a variable pruned from a reduced model, is dropped with a warning, together with its value. The remaining values stay with their own variables.
- **A value vector of any other length** (neither 1 nor the length of the variable vector) is recycled with a warning.
- **Unlisted axes** keep the defaults: label line `0.1`, and `ticks_v` ticks.

**What it renders (in layer order):**

| Layer | Source data | Visual element |
|---|---|---|
| 1 | `bl_results$biplot_obj` + `V[, proj_dims]` | Grey coordinate axes + variable labels |
| 2 | `biplot_grid$Zgrid` + `col_value` | Coloured squares for the grid points inside the hull (the 200×200 grid trimmed by `calc_hull`; blue→red probability surface) |
| 3 | `points$Z` + `points$pred_col` | Training obs as coloured dots: TP=red, TN=blue, FP=purple, FN=orange |
| 4 | `biplot_obj` again | Darker axes re-drawn on top of the grid and points |
| 5 | `biplot_grid$ct` | Black contour lines marking the decision boundary (probability ≈ 0.5) |

#### biplotEZ visual features: what `plot_biplotEZ()` uses, and what's available for later

`plot_biplotEZ()` deliberately uses only a thin slice of biplotEZ's rendering
API — enough to draw the coordinate-axis skeleton (`axes()`) on a blank canvas
(`samples(opacity = 0)`), then layers the package's own grid, points, contours,
and target marker on top with `graphics::points()`/`lines()`/`text()`. This is
intentional: `boundarylogic` needs full control over confusion-category colours
(TP/TN/FP/FN), the prediction-probability surface, decision-boundary contours,
and synchronised rotation — none of which map cleanly onto biplotEZ's
group-aesthetic-driven sample styling.

**Currently used:** `biplot()`, `PCA()`/`CVA()` (projection — `R/projection.R`),
`samples()` (canvas only, `opacity = 0`), `axes()` (full axis styling, exposed
via `label_dir`, `label_cex`, `tick_label_cex`, `ticks_v`/`ticks_var`/`ticks_n`,
`label_offset_var`/`label_offset_dist`, `which`, `X_names`).

**Available in biplotEZ but not currently used — candidates for future
additions:**

| biplotEZ function | What it would add | Notes |
|---|---|---|
| `means()` | Class/group mean markers + labels | Would need confusion-category-aware grouping to fit the package's colour scheme |
| `alpha.bags()` | Alpha-bag contours -- region containing alpha% of a class's points | Density-based summary region per class |
| `ellipses()` | Kappa concentration ellipses per class | Parametric alternative to alpha-bags |
| `density1D()` / `density2D()` | Kernel density curves / heatmaps over the biplot | Visual alternative/complement to the existing probability-coloured grid |
| `legend.type()` | In-plot legend for samples/means/bags/ellipses/regions | Package currently prints a text summary to the console instead |
| `newsamples()` / `newaxes()` | Supplementary points/axes styled distinctly (e.g. orange) from the primary biplot | Could visually distinguish genuinely external data from train/test overlays |
| `interpolate()` | Projects a new dataset onto an *existing* fitted biplot object | Different mechanism from `bl_project_points()`, which re-derives coordinates from `bl_result$V`/`X_center` directly |
| `classification()` / `prediction()` | biplotEZ's own classification-region and sample/mean prediction overlays | Alternative to this package's custom prediction-grid + contour approach |
| `regress()` | Regression biplot (axes fitted via linear regression or B-splines) | An alternative *construction*, not an addition to the current PCA/CVA biplot |
| `rotate()` / `reflect()` / `translate_axes()` | Generic geometric transforms | `boundarylogic` already implements its own (`.apply_biplot_rotation()`) because rotation must keep the prediction grid, contours, and target point in sync -- biplotEZ's generic versions don't know about these package-specific elements |
| `fit.measures()` | Quality-of-fit diagnostics (overall quality, axis/sample/class predictivity) | Not visual, but could help decide which additions above are worth prioritising |
| `CA()`, `CATPCA()`, `PCO()`, `AoD()`, `CLPs()`, `CLRs()` | Alternative biplot constructions for categorical / distance-matrix data | Not applicable to this package's continuous-feature binary-classification setting |

None of the above are wired into `plot_biplotEZ()`, `plot.bl_local_result()`,
or `plot.bl_surrogate()` today. Adding any of them would be new, scoped feature
work (see CLAUDE.md Section 9, "Deferred / Future Work").

**`points` here:** Because no `points=` argument is passed, `plot_biplotEZ()` automatically calls `bl_project_points(bl_results$train_data, bl_results)` internally to project the training data.

**Returns:** `bl_results` invisibly (for pipe-compatibility). Side effect: renders plot to active device, and prints a prediction summary (in the harness run: 15,970 points, accuracy 82.8%, TP / TN / FP / FN = 5472 / 7751 / 1080 / 1667).

---

### 6.5 — `bl_project_points()` → `test_pts`

**File:** `R/project_points.R`

**Call:**
```r
test_pts <- bl_project_points(bl_results$test_data, bl_results, filter_to_polygon = TRUE)
```

**Parameters:**

| Parameter | Default | This call | Effect |
|---|---|---|---|
| `filter_to_polygon` | `FALSE` | `TRUE` | Drops rows that project outside the training hull (`bl_results$polygon`), so the overlay shows only the region the training data supports: 9 of 4,429 here, leaving 4,420 |
| `filter_to_train_ranges` | `FALSE` | default | Drops rows with any feature outside `train_ranges`, before projecting |

Full argument list: `?bl_project_points`.

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
   pred_prob  = .pred_function(bl_results$model, bl_results$model_type, X_test[, var_names])   # "custom"
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
plot(bl_results, points = test_pts)   # or plot_biplotEZ(...)
```

This pattern is used in the development harness (this step, and Steps 4b–6b for the reduced model), the loan-default vignette, and the README. Variants: `filter_to_polygon = TRUE` (the harness here, and the vignette) drops out-of-hull observations before plotting; `filter_to_train_ranges = TRUE` (no current scripted use) drops X-space-extrapolated rows. Without an explicit `bl_project_points()` call, `plot_biplotEZ()` only ever shows training data, because it auto-projects `bl_result$train_data` when `points` is `NULL`.

---

### 6.6 — `plot_biplotEZ(bl_results, points = test_pts)` / `plot(bl_results, points = test_pts)` — overlays test data

**Same function, different arguments.** `plot(bl_results, points = ...)` works via the `plot.bl_result()` S3 wrapper, as in 6.4. Renders identical layers 1–5 as above, but uses `test_pts$Z` and `test_pts$pred_col` for the data points (layer 3) instead of projecting training data automatically. This shows where the held-out test observations fall on the decision surface. Its prediction summary in the harness run: 4,420 points, accuracy 82.2%, TP / TN / FP / FN = 1494 / 2138 / 345 / 443.

---

## Complete Object Flow Summary

```
loan_encoded (data.frame, 45,000 rows × 11 cols)
     │
     │  Step 3: row filter, setdiff()
     ▼
loan_filtered (data.frame, 22,142 rows × 11 cols)
feature_cols  (character vector, 9 predictor names)
     │
     │  bl_prepare_data(hull_fraction = 0.9)
     │  → split → bl_filter_outliers() → standardise → SVD → 2D hull → point-in-polygon
     ▼
bl_dat  [bl_filter_result]
  ├── train_data  (filtered, 15,970 rows)
  ├── test_data   (unchanged, 4,429 rows)
  └── var_names, num_vars, target_class, scaling, polygon, hull_fraction, n_retained, n_removed
     │
     │  bl_wrap_model(model_type = "custom", predict_fn = ...)
     │  → external xgboost fit → accuracy + Gini on train
     ▼
bl_mod  [bl_model]
  ├── model       (list(model = xgb.Booster, predict_fn = fn))
  ├── model_type  "custom"
  ├── var_names
  ├── cutoff      0.5
  ├── accuracy    0.8280
  └── gini        0.8093
     │
     │  bl_build_result(method = "CVA", b_margin = 0.001)
     │  → bl_build_projection() → V, tV, biplot_obj, train_ranges, W metric
     │  → bl_build_grid()       → 200×200 grid, scored, contours, hull-trimmed
     │  → bl_assemble()         → combines all
     ▼
bl_results  [bl_result]              ← THE CENTRAL ANCHOR OBJECT
  ├── train_data, test_data, var_names, num_vars
  ├── model, model_type="custom", cutoff, b_margin=0.001
  ├── V, tV, X_center, X_sd, method="CVA", standardise=FALSE, proj_dims
  ├── biplot_obj       (biplotEZ — for axis rendering)
  ├── metric, metric_inv, metric_type  (within-class W, Mahalanobis metric)
  ├── train_ranges     (min/max per feature)
  ├── polygon, hull_fraction  (grid hull in CVA Z-space, outlie = 1)
  ├── scaling          (NULL here)
  ├── biplot_grid      (Zgrid, Xgrid, grid_prob, ct, ct_surrogate)
  ├── accuracy, gini
  └── call, created_at
     │
     │  plot(bl_results)   → renders training biplot   # or plot_biplotEZ(...)
     │
     │  bl_project_points(bl_results$test_data, bl_results, filter_to_polygon = TRUE)
     │  → project using same V, X_center → Z_test; drop out-of-hull rows
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
     │  plot(bl_results, points = test_pts)   # or plot_biplotEZ(...)
     └─→ renders same biplot with test points overlaid
```

---

## Key Design Principles to Note

1. **One anchor object:** `bl_results` carries everything. All Phase 2 and Phase 3 functions take `bl_results` as their first argument and extract what they need.

2. **Two-way street between X and Z:** `V` projects forward (X → Z) and `tV = solve(V)` projects backward (Z → X). The entire boundary and counterfactual machinery depends on this invertibility.

3. **Consistent coordinate frame:** `X_center` and `X_sd` from training are applied to every subsequent projection (test data, grid, target points) — so all Z-space coordinates are directly comparable.

4. **CVA vs PCA:** CVA (`method = "CVA"`) maximises separation between the four confusion categories (TP/TN/FP/FN), producing a biplot where the decision boundary is most clearly visible. PCA would instead maximise total variance, which may not align with the class boundary.

5. **`b_margin` controls the contour band half-width**, not the predictions themselves. `b_margin = 0.001` (the package default, used for `bl_results`) means the boundary search looks for contour lines at probability 0.499 and 0.501 (a 0.002-wide band total). The reduced-model `bl_results_v2` uses `b_margin = 0.01` (0.49 and 0.51). All predicted values are always stored at 3 d.p. regardless.
