---
name: reference-data-prep-functions
description: "bl_prepare_data (integrated split + filter, returns bl_filter_result) vs bl_wrap_data (external path, no filter, returns bl_data) — signatures, return fields, when to use each; bl_set_scaling and bl_apply_scaling"
metadata:
  node_type: memory
  type: reference
  originSessionId: 5d1f2f11-cc4c-460f-b0fd-25216426084d
  modified: 2026-11-15T23:13:22.933Z
---

# Data preparation: bl_prepare_data vs bl_wrap_data

Source files: `R/data_prepare.R`, `R/outlier_filter.R`. This file has an identical copy in
`.claude/reference/` (CLAUDE.md §2, memory-to-reference sync rule).

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
                hull_fraction = 1, verbose = TRUE)
```

Key parameters:
- `class_col` — target column name; renamed to `"class"` internally
- `target_class` — value to map to 1 (all others -> 0); if `NULL`, column must already be 0/1
- `hull_fraction = 1` — fraction of training points retained by the convex hull filter. The
  default `1` keeps all points (the polygon is still built; no outlier removal); `0.9` removes
  the outermost ~10%, as the loan harness and vignette set it
- `train_fraction = 1` — all rows go to `train_data`, `test_data` is empty (exploratory use)
- **Always shuffles** rows regardless of `train_fraction`; the split depends only on the rows and
  the seed, so the same rows with fewer `feature_cols` give the same split

How it works internally (final block of `bl_prepare_data()`): it builds an interim `"bl_data"`
object (`train_data`, `test_data`, `var_names`, `num_vars`, `target_class`) and passes it to
`bl_filter_outliers(bl_data_raw, hull_fraction = hull_fraction, verbose = verbose)`. The interim
object is never exposed to the caller. `bl_filter_outliers()` is an internal implementation
detail here, not a step the caller takes.

Returns `"bl_filter_result"` with 10 fields:
`train_data`, `test_data`, `var_names`, `num_vars`, `target_class`, `scaling` (NULL unless
`bl_set_scaling()` was applied first), `polygon`, `hull_fraction`, `n_retained`, `n_removed`.

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
`train_data`, `test_data`, `var_names`, `num_vars`, `target_class` (plus `scaling` once
`bl_set_scaling()` is applied).

Result variable: `bl_dat`.

---

## Comparison

| | `bl_prepare_data()` | `bl_wrap_data()` |
|---|---|---|
| Row order | Shuffled | Preserved |
| Class encoding | Handles multiclass via `target_class` | Must already be 0/1 |
| Outlier filtering | Integrated via `hull_fraction` | None |
| Return class | `"bl_filter_result"` (10 fields) | `"bl_data"` (5 fields) |
| Downstream compatible | Yes | Yes |

Both return types are accepted by `bl_assemble()` and `bl_build_result()`.

---

## bl_set_scaling() — record raw-unit display scaling (optional)

Use when features were standardised **before** model fitting, so the biplot axes would
otherwise read in standardised units. Attach it as the final data-prep step (after
`bl_prepare_data()` / `bl_wrap_data()`), passing the per-feature transform
`standardised = (raw - center) / scale`.

```r
bl_set_scaling(x, center, scale, method = "z-score")
# x: a "bl_data" or "bl_filter_result"; center/scale: named or var_names-ordered numeric
bl_dat <- bl_set_scaling(bl_dat, center = attr(z, "scaled:center"),
                                 scale  = attr(z, "scaled:scale"))
```

- Stores `x$scaling = list(center, scale, method)` (vectors reordered to `var_names`).
- Validated/reordered by private `.align_scaling_vec()`; rejects zero/NA scales and names
  that don't cover every feature.
- Propagates through `bl_filter_outliers()` and `bl_assemble()` into `bl_result$scaling`
  (and onto `bl_projection$scaling` when the projection is built through `bl_build_result()`).
- **Display-only:** never changes the projection, grid, or what the model/`predict_fn` sees.
  At plot time each biplot method calls private `.bl_rescale_biplot_axes()` (in
  `R/plot_biplot.R`), which transforms only the biplotEZ `$means`/`$sd` so axis labels read
  in raw units while geometry (`$Z`/`$Lmat`/`$ax.one.unit`) is byte-identical.
- Distinct from `X_center`/`X_sd` (internal PCA standardisation) and the `standardise` flag.
- One affine `method` for all features; per-variable transform families are a future
  development (`documentation/future-improvements-review.md`).
- Observed/counterfactual/sparse **feature values** also display in raw units (the Shapley plot
  labels `observed -> counterfactual`, `print.bl_shapley()`, `print.bl_sparse_result()`,
  `print.bl_target()`, `print.bl_counterfactual()`, `print.bl_local_result()`), via private
  `.scale_to_raw()` (level = `v*scale+center`, delta = `v*scale`). Shapley *contributions* stay
  in prediction-impact units. See `documentation/scaling_axis_relabel_note.md` and
  `2 implementation_summary.txt` §4.4.1 and §4.9.
- In scripts, fit the scaling on the **training** rows only and apply those stats to the test
  split (split first via `bl_prepare_data()`, then `scale()` train, then apply to test) to
  avoid test-set leakage. See the harness's v2 block (`development/03_loan_status_Boundary_Logic.R`).

## bl_apply_scaling() — new raw-unit data into model units

```r
bl_apply_scaling(data, bl_result)   # bl_result: a bl_result, bl_data or bl_filter_result
```

- Returns `data` with each feature converted as `(raw - center) / scale`, the input-side
  counterpart of `.scale_to_raw()`. Every function that takes new data expects model units.
- `bl_select_target()` warns when a data-frame target is outside `train_ranges`, and suggests
  `bl_apply_scaling()` when converting would bring it into range. The harness's Step 16
  converts its external applicant this way.

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
- Returns `"bl_filter_result"` (the same 10 fields as `bl_prepare_data()`)

Not shown in scripts, vignettes, or examples (CLAUDE.md §8).

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
