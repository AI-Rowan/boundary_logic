---
name: reference-data-prep-functions
description: "bl_prepare_data (integrated split + filter, returns bl_filter_result) vs bl_wrap_data (external path, no filter, returns bl_data) — signatures, return fields, when to use each"
metadata:
  node_type: memory
  type: reference
  originSessionId: 5d1f2f11-cc4c-460f-b0fd-25216426084d
---

# Data preparation: bl_prepare_data vs bl_wrap_data

Source files: `R/data_prepare.R`, `R/outlier_filter.R`

There are two entry points depending on where your data comes from. Pick one; they are
not meant to be chained.

---

## bl_prepare_data() — in-package path

Use when starting from a raw data frame inside the package workflow. Does everything in
one call: feature selection, binary class encoding, train/test split, and convex hull
outlier filtering. Returns a `"bl_filter_result"` ready for `bl_fit_model()` and
`bl_build_result()`.

```r
bl_prepare_data(data, class_col, target_class = NULL,
                feature_cols = NULL, train_fraction = 0.8, seed = 121L,
                hull_fraction = 0.9, verbose = TRUE)
```

Key parameters:
- `class_col` — target column name; renamed to `"class"` internally
- `target_class` — value to map to 1 (all others -> 0); if `NULL`, column must already be 0/1
- `hull_fraction = 0.9` — fraction of training points retained by the convex hull filter;
  `1` keeps all points (polygon still built); default removes outermost ~10%
- `train_fraction = 1` — all rows go to `train_data`, `test_data` is empty (exploratory use)
- **Always shuffles** rows regardless of `train_fraction`

How it works internally (final block of `R/data_prepare.R`):
```r
bl_data_raw <- structure(
  list(train_data = train_data, test_data = test_data,
       var_names = feature_cols, num_vars = length(feature_cols),
       target_class = target_class),
  class = "bl_data"
)
bl_filter_outliers(bl_data_raw, hull_fraction = hull_fraction, verbose = verbose)
```
The interim `"bl_data"` object is never exposed to the caller. `bl_filter_outliers()` is
an internal implementation detail here, not a step the caller takes.

Returns `"bl_filter_result"` with 9 fields:
`train_data`, `test_data`, `var_names`, `num_vars`, `target_class`,
`polygon`, `hull_fraction`, `n_retained`, `n_removed`.

Result variable: `bl_dat`.

---

## bl_wrap_data() — external-data path

Use when data is already split and class-encoded outside the package (e.g. from a CV
framework, a pre-defined holdout, or an externally trained model workflow). Packages the
frames without re-splitting, re-encoding, or filtering.

```r
bl_wrap_data(train_data, test_data = NULL, var_names = NULL, target_class = NULL)
```

- `train_data` must already have a column named `"class"` with numeric 0/1 values
- `target_class` stored for reference only — no conversion performed
- **No shuffle** — original row order preserved
- **No outlier filtering** — by design; downstream steps use `bl_dat` directly
- If `test_data = NULL`, an empty frame with matching columns is created

Returns `"bl_data"` with 5 fields:
`train_data`, `test_data`, `var_names`, `num_vars`, `target_class`.

Result variable: `bl_dat`.

---

## Comparison

| | `bl_prepare_data()` | `bl_wrap_data()` |
|---|---|---|
| Row order | Shuffled | Preserved |
| Class encoding | Handles multiclass via `target_class` | Must already be 0/1 |
| Outlier filtering | Integrated via `hull_fraction` | None |
| Return class | `"bl_filter_result"` (9 fields) | `"bl_data"` (5 fields) |
| Downstream compatible | Yes | Yes |

Both return types are accepted by `bl_assemble()` and `bl_build_result()`.

---

## bl_filter_outliers() — internal / power-user only

`bl_filter_outliers()` is what `bl_prepare_data()` calls internally. Callers of
`bl_prepare_data()` never interact with it directly.

The only reason to call it explicitly is after `bl_wrap_data()`, when you want to
iterate on different `hull_fraction` values without re-splitting the data.

```r
bl_filter_outliers(bl_data, hull_fraction = 0.9, verbose = TRUE)
```

- Input must be a `"bl_data"` object (from `bl_wrap_data()`)
- Standardises features, projects to first 2 PCs, builds convex hull via
  `aplpack::plothulls()`, removes points outside the hull
- Test data passes through unchanged
- Returns `"bl_filter_result"` (same 9 fields as `bl_prepare_data()`)

Not shown in scripts, vignettes, or examples.

---

## Usage patterns

### Standard workflow
```r
bl_dat <- bl_prepare_data(data, class_col = "loan_status",
                           feature_cols = feature_cols,
                           train_fraction = 0.8, seed = 121L,
                           hull_fraction = 0.9)
# bl_dat is "bl_filter_result" — pass directly to bl_fit_model() / bl_build_result()
```

### External-data workflow
```r
bl_dat <- bl_wrap_data(train_data, test_data, var_names = feature_cols)
# bl_dat is "bl_data" — pass directly to bl_fit_model() / bl_build_result()
# No filtering step; bl_filter_outliers() is not called
```

### Exploratory biplot (all data, no model)
```r
bl_dat <- bl_prepare_data(data, class_col = "loan_status",
                           train_fraction = 1, hull_fraction = 1)
bl_proj <- bl_build_result(bl_data = bl_dat, method = "PCA", title = "...")
plot_biplotEZ(bl_proj)
```

### Power user — iterate hull fractions on pre-split data
```r
# Only pattern where bl_filter_outliers() is called explicitly
bl_dat_raw <- bl_wrap_data(train_data, test_data, var_names = feature_cols)
bl_dat_90  <- bl_filter_outliers(bl_dat_raw, hull_fraction = 0.90)
bl_dat_95  <- bl_filter_outliers(bl_dat_raw, hull_fraction = 0.95)
```
