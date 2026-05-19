---
name: reference-data-prep-functions
description: "How bl_prepare_data, bl_wrap_data, and bl_filter_outliers work, why they are separate, and how to combine them — including the train_fraction=1 vs bl_wrap_data comparison"
metadata: 
  node_type: memory
  type: reference
  originSessionId: 5d1f2f11-cc4c-460f-b0fd-25216426084d
---

# Data preparation functions: bl_prepare_data, bl_wrap_data, bl_filter_outliers

Source files: `R/data_prepare.R`, `R/outlier_filter.R`

---

## bl_prepare_data()

**Purpose:** All-in-one entry point for raw data. Handles feature selection, binary
class encoding, and a seeded train/test split in a single call.

```r
bl_prepare_data(data, class_col, target_class = NULL,
                feature_cols = NULL, train_fraction = 0.8, seed = 121L)
```

- `class_col` — name of the target column (any name; renamed to `"class"` internally)
- `target_class` — if the target is multiclass, the value to map to 1 (all others → 0);
  if `NULL`, `class_col` must already be numeric 0/1
- `train_fraction = 1` — retains all rows in `train_data`, produces an empty `test_data`;
  used for exploratory biplots (script 03 Step 2 pattern)
- **Always shuffles rows** via `set.seed(seed); sample(n, floor(train_fraction * n))`
  even when `train_fraction = 1`

Returns: `"bl_data"` object.

---

## bl_wrap_data()

**Purpose:** Bypass entry point for users who have already split their data (e.g. from
a CV framework, external pipeline, or pre-defined holdout). Validates and packages
pre-split frames without re-splitting or re-encoding.

```r
bl_wrap_data(train_data, test_data = NULL, var_names = NULL, target_class = NULL)
```

- `train_data` must already contain a column named **`"class"`** with numeric 0/1 values
- `target_class` is stored informally only — no conversion is performed
- **Preserves original row order** — no shuffle
- If `test_data = NULL`, creates an empty frame with matching columns

Returns: `"bl_data"` object — structurally identical to `bl_prepare_data()` output.

### bl_prepare_data(train_fraction=1) vs bl_wrap_data()

| | `bl_prepare_data(train_fraction=1)` | `bl_wrap_data(test_data=NULL)` |
|---|---|---|
| Row order | **Shuffled** (sample runs) | **Preserved** |
| test_data | Empty (0 rows) | Empty (0 rows) |
| Class encoding | Handles multiclass via `target_class` | Must already be 0/1 |
| Return class | `"bl_data"` | `"bl_data"` |
| Downstream compatible | Yes | Yes |

For biplots and model fitting, the shuffle makes no practical difference.
For time-ordered or index-matched data, prefer `bl_wrap_data()`.

---

## bl_filter_outliers()

**Purpose:** Separate modelling decision — trims extreme training observations using
a convex hull in standardised 2D PCA-proxy space.

```r
bl_filter_outliers(bl_data, hull_fraction = 0.9, verbose = TRUE)
```

- Accepts a `"bl_data"` object from **either** `bl_prepare_data()` or `bl_wrap_data()`
- `hull_fraction = 0.9` removes the outermost ~10 % of training points
- `hull_fraction = 1` retains all training points but still builds the polygon
- **Test data is passed through unchanged** — only training rows are filtered
- Returns `"bl_filter_result"` (adds `polygon`, `hull_fraction`, `n_retained`, `n_removed`)

The polygon produced here flows into `bl_build_grid()` to clip prediction contours
to the observed-data region.

---

## Why they are separate

1. **Iterative outlier tuning** — users commonly call `bl_filter_outliers()` multiple
   times with different `hull_fraction` values to inspect retention without re-splitting
   the data. Merging would force a re-split each time.

2. **`bl_wrap_data()` compatibility** — users bringing pre-split data via `bl_wrap_data()`
   can still apply outlier filtering afterwards. Merging filtering into `bl_prepare_data()`
   would break this path.

3. **Separation of concerns** — `bl_prepare_data()` / `bl_wrap_data()` are about data
   structure; `bl_filter_outliers()` is a modelling decision. Inspecting the split before
   deciding to filter is intentional.

4. **Different return classes** — `bl_filter_result` carries extra fields (`polygon`,
   counts) that are meaningful outputs in their own right.

---

## How to combine them

### Standard path (raw data, auto split)
```r
bl_dat  <- bl_prepare_data(data, class_col = "loan_status",
                           feature_cols = feature_cols,
                           train_fraction = 0.8, seed = 121L)
bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 0.9)
# pass bl_filt to bl_fit_model() and bl_build_result()
```

### Pre-split path (bring your own split)
```r
# Manual steps replacing bl_prepare_data:
data_clean <- raw_data[, c(feature_cols, "target_col")]
names(data_clean)[names(data_clean) == "target_col"] <- "class"
set.seed(121L)
train_idx  <- sample(nrow(data_clean), floor(0.8 * nrow(data_clean)))
train_data <- data_clean[ train_idx, ]; rownames(train_data) <- NULL
test_data  <- data_clean[-train_idx, ]; rownames(test_data)  <- NULL

bl_dat  <- bl_wrap_data(train_data, test_data, var_names = feature_cols)
bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 0.9)
# pass bl_filt to bl_fit_model() and bl_build_result()
```

### Exploratory biplot (no model, all data)
```r
bl_dat_exp  <- bl_prepare_data(data, class_col = "loan_status", train_fraction = 1)
bl_filt_exp <- bl_filter_outliers(bl_dat_exp, hull_fraction = 1)
bl_proj     <- bl_build_result(bl_data = bl_filt_exp, method = "PCA", title = "...")
plot_biplotEZ(bl_proj)
```

See `scripts/04_loan_wrap_data_demo.R` for a full working example of the pre-split path.
