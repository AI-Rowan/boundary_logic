> **IMPLEMENTED 2026-06-16.** Archived copy of the approved plan from
> `.claude/plans/currently-the-biplot-does-greedy-kurzweil.md`. Delivered: private
> `.scale_to_raw()` (level vs delta), raw observed->counterfactual values across
> `plot.bl_shapley()` / `print.bl_shapley()` / `print.bl_sparse_result()` /
> `print.bl_target()` (scaling stored on `bl_target`), and the script 03
> training-only scaling fix. Shapley contributions left in prediction-impact
> units. Tests 127 PASS / 0 FAIL. The *why* lives in
> `documentation/scaling_axis_relabel_note.md`; the *how* in
> `2 implementation_summary.txt` §4.9.

# Plan: Raw-unit feature values in Shapley / sparse / target displays

Date created: 2026-06-16
Date last edited: 2026-06-16

## Context

The earlier change (`bl_set_scaling()` + `.bl_rescale_biplot_axes()`) relabels biplot
**axes** into raw units but deliberately left **printed/plotted feature values** in model
(standardised) units -- recorded as deferred work in CLAUDE.md S9. The Shapley contribution
plot still shows each variable's observed value and change in standardised units (its
per-variable y-axis label, built as `"<var>: <observed> -> <change>"`), which is unreadable
when the user standardised features before fitting.

This change converts the **current (observed) and counterfactual feature values** back to
raw units everywhere they are displayed, gated on `bl_result$scaling`. It is display-only:
the stored `shapley_df` / `x_sparse` / `x_obs` stay in model units (the model and geometry
are untouched), exactly mirroring the axis-relabel philosophy. Shapley *contributions*
(`shapley_cause`) are prediction-impact units, **not** feature units, and are left unchanged.

Per the user's decisions:
- Shapley plot label becomes **observed -> counterfactual** (both levels, raw), e.g.
  `income: 45000 -> 57000`, replacing the old `observed -> change`.
- All related print tables also convert: `print.bl_shapley`, `print.bl_sparse_result`,
  `print.bl_target` (and therefore `print.bl_local_result`, which delegates to it).

## Key arithmetic (must be exact)

For per-feature scaling `std = (raw - center) / scale`:
- **Level** value (observed, full CF `B_x`, sparse CF, counterfactual end):
  `raw = value * scale + center`.
- **Delta** value (`data_to_boundary = end - start`): `raw = delta * scale` (NO centre term).

The Shapley plot/print currently store `pred_data` (level = observed) and `data_to_boundary`
(delta). The counterfactual level is `pred_data + data_to_boundary` (model units) -> convert
as a level.

## Files to read first
- `R/shapley.R` -- `bl_shapley()` (builds `shapley_df`, `varnames_p` at lines ~206-228),
  `plot.bl_shapley()` (~260), `print.bl_shapley()` (~303), `bl_find_sparse_cf()`,
  `print.bl_sparse_result()` (~512).
- `R/local_cf.R` -- `bl_select_target()` (`bl_target` construction ~184), `print.bl_target()`
  (~1031), `print.bl_local_result()` (delegates to `print.bl_target`).
- `R/data_prepare.R` -- `.align_scaling_vec()` (the existing private scaling helper to sit beside).

## Design

### New private helper: `.scale_to_raw()` -- `R/data_prepare.R`
Beside `.align_scaling_vec()` (both private scaling utilities):
```r
.scale_to_raw(values, scaling, var_names, kind = c("level", "delta"))
```
- `values`: numeric aligned to `var_names`.
- No-op (returns `values`) when `scaling` is `NULL`.
- Aligns `scaling$center`/`scale` to `var_names`; if any feature is missing, returns `values`
  unchanged with a warning (same defensive pattern as `.bl_rescale_biplot_axes()`).
- `kind = "level"` -> `values * scale + center`; `kind = "delta"` -> `values * scale`.

### Apply at display time (storage stays in model units)

| Location | File | Change |
|---|---|---|
| `plot.bl_shapley()` | `R/shapley.R` | Read `scaling <- x$bl_local_result$bl_result$scaling`. When non-NULL: `obs_raw <- .scale_to_raw(df$pred_data, ..., "level")`, `cf_raw <- .scale_to_raw(df$pred_data + df$data_to_boundary, ..., "level")`, rebuild the y label as `"<var>: <obs_raw> -> <cf_raw>"` (ordered factor preserving the existing sort), and switch the y-axis title to note "observed -> counterfactual (raw units)". When NULL, keep the current `varnames_p` and title unchanged. |
| `print.bl_shapley()` | `R/shapley.R` | When scaling present, convert the printed `pred_data` (level) and `data_to_boundary` (delta) columns; add a one-line "(raw units)" note. `shapley_cause`/`Contribute` unchanged. |
| `print.bl_sparse_result()` | `R/shapley.R` | When scaling present, convert `x_obs`, `B_x`, `x_sparse` (all levels) in the variable-summary table. |
| `print.bl_target()` | `R/local_cf.R` | When scaling present, print the feature values in raw units. Requires the scaling on the object (next item). |
| `bl_select_target()` | `R/local_cf.R` | Add `scaling = bl_result$scaling` to the `bl_target` structure so `print.bl_target()` can reach it (additive field; removes no metadata). |

Notes:
- Gate the Shapley plot's *semantic* change (observed->counterfactual) on scaling presence:
  with no scaling the label keeps the existing observed->change form, so non-scaled users and
  existing output are unaffected.
- `plot.bl_sparse_result()`'s console summary shows only probabilities -- no feature values,
  no change needed.

## Error cases
- `scaling = NULL` (older objects / no scaling): every path is a no-op; output identical to today.
- Scaling not covering all features: `.scale_to_raw()` warns once and leaves values in model
  units (defensive, matches `.bl_rescale_biplot_axes()`).
- Rounding preserved: `round(., 3)` for plot labels, `round(., 4)` for print tables (as today).

## Test cases (`tests/testthat/test-set_scaling.R`, extend)
- `.scale_to_raw()`: level vs delta arithmetic; NULL no-op; named/ordered alignment; missing
  feature -> warning + unchanged.
- Integration on iris with a GLM model + `bl_set_scaling()`:
  build `bl_find_local_cf()` -> `bl_shapley()` -> `bl_find_sparse_cf()`;
  - assert the `plot.bl_shapley()` ggplot's y labels contain raw-scale observed/CF numbers
    (inspect the built plot's data / label levels), not standardised ones;
  - `capture.output(print(sparse))` and `print(bl_shapley)` contain raw-scale values;
  - with `scaling = NULL`, outputs match the pre-change form.

## Documentation to update
- `2 implementation_summary.txt` Section 4.9 (Shapley module) -- document the display-time
  raw-unit conversion of observed/counterfactual/sparse feature values and `.scale_to_raw()`;
  note contributions stay in prediction-impact units.
- `.claude/reference/review_section9_to_15.md` -- Steps 13/14: note the raw-unit value
  conversion in the Shapley plot label and the print tables.
- `documentation/scaling_axis_relabel_note.md` -- update "Scope and limitations": value
  displays (observed/CF/sparse) now convert; only Shapley contributions remain in model units.
- `CLAUDE.md` S9 -- update the existing "Raw-unit reporting of counterfactual/target values"
  deferred bullet to mark it implemented (feature-value displays now convert; Shapley
  contributions intentionally not converted).
- Memory sync: `~/.claude/.../memory/reference_data_prep_functions.md` `bl_set_scaling`
  section currently says "printed counterfactual/Shapley values stay in model units
  (deferred)" -- update to reflect that observed/CF feature values now convert.
- Archive this plan to `.claude/reference/` with an IMPLEMENTED banner once executed.

## Additional: fix script 03 to scale on the training data only

**Problem.** The current v2 block (`scripts/03_loan_status_Boundary_Logic.R` ~331-350)
standardises with `scale(loan_std[, feature_cols_v2])` on the **full** `loan_filtered`
(train + test) *before* `bl_prepare_data()` splits it. The centre/scale therefore leak
test-set information -- the standardisation should be derived from the training rows only and
then *applied* to the test rows.

**Fix.** Split first on raw data (the internal hull filter self-standardises, so filtering on
raw vs standardised is equivalent), derive the transform from the (filtered) **training**
features, apply the same training centre/scale to both train and test in place, then record it:

```r
# Split + outlier-filter on raw data first
bl_dat_v2 <- bl_prepare_data(
  data           = loan_filtered,
  class_col      = "loan_status",
  feature_cols   = feature_cols_v2,
  train_fraction = 0.8, seed = 121L, hull_fraction = 0.9
)

# Scaling derived from TRAINING rows only (no test leakage)
train_scl  <- scale(bl_dat_v2$train_data[, feature_cols_v2])
scl_center <- attr(train_scl, "scaled:center")
scl_scale  <- attr(train_scl, "scaled:scale")

# Apply the training transform to both splits (test uses training stats)
bl_dat_v2$train_data[, feature_cols_v2] <- train_scl
bl_dat_v2$test_data[,  feature_cols_v2] <-
  scale(bl_dat_v2$test_data[, feature_cols_v2],
        center = scl_center, scale = scl_scale)

bl_dat_v2 <- bl_set_scaling(bl_dat_v2, center = scl_center,
                            scale = scl_scale, method = "z-score")
```

Update the block comment to state the scaling is fit on the training split and applied to the
test split. This is a script-only change (the package already accepts whatever centre/scale the
user records); no package code changes for this part. Mention training-only scaling as the
recommended pattern in `.claude/reference/review_section4_to_6.md` where `bl_set_scaling()` is
documented.

## Verification
- `devtools::document()` clean; `devtools::test()` >= current 113 PASS / 0 FAIL with new tests.
- Manual: write `verify_shapley_scaling.R` (run via `Rscript verify_shapley_scaling.R`, not
  `-e`, per the biplotEZ Windows rule) building the iris local-CF -> shapley pipeline with a
  known `(center, scale)`; print the Shapley table and sparse summary and eyeball that
  observed/CF values match raw ranges. Delete the scratch file after use (diff against the
  session-start `git status` snapshot first).
- Report PASS/FAIL and stop; no commit/push without explicit instruction.
