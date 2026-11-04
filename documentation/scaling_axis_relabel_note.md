# Technical note: raw-unit biplot axes for pre-standardised features

## Motivation

A user who standardises features before fitting a model (common for SVM, neural
nets, penalised logistic regression, or any pipeline that scales once up front)
hands `boundarylogic` data that is already in standardised units. The PCA/CVA
projection, the prediction grid, and -- visibly -- the biplot axis tick labels
are then all expressed in standardised units (e.g. an age axis reading
`-1.5, 0, 1.5` rather than `25, 35, 45`). This makes reading values off the
biplot, the whole point of the method, meaningless.

The goal is to display axis ticks in the original (raw) units **without** changing
the projection, the grid, or what the model sees -- the geometry must stay in the
standardised space the model was trained on. This is therefore a pure display
relabel, not a re-projection.

## API

`bl_set_scaling(x, center, scale, method = "z-score")` records the per-feature
affine transform

    standardised = (raw - center) / scale            =>    raw = standardised * scale + center

on a `bl_data` / `bl_filter_result` object as `x$scaling = list(center, scale, method)`.
It propagates through `bl_filter_outliers()` and `bl_assemble()` into
`bl_result$scaling` (and onto `bl_projection$scaling`). It is display-only metadata;
a custom model wrapped with `bl_wrap_model(predict_fn = ...)` keeps receiving the
standardised data it was trained on, because the grid back-projects to the data's
own (standardised) units.

## Why a `$means`/`$sd` transform suffices (and is exact)

biplotEZ computes axis tick **labels** in `.calibrate.axis()` (`plot2D.R`) as an
affine function of the biplot object's `$means` and `$sd`:

    nice ticks  :  pretty(range(Xhat))                  # Xhat reconstructed from means/sd/scaled/center
    position    :  axis.vals = (label - means) / sd     # along ax.one.unit
    label drawn :  means + sd * axis.vals

The 2D geometry -- `$Z`, `$Lmat`, `$ax.one.unit`, `$e.vects` -- is independent of
`$means`/`$sd`, and a data point's displacement along an axis (`axis.vals`) is
fixed by geometry alone. For a pre-standardisation with per-feature `(center c,
scale s)`, set

    means_new = means * s + c
    sd_new    = sd * s
    $scaled   = TRUE ,  $center = TRUE          # geometry fields untouched

Then for any tick:

    label_new = means_new + sd_new * axis.vals
              = (means + sd * axis.vals) * s + c
              = label_old * s + c
              = raw

and `Xhat` reconstructs onto the raw range so `pretty()` places ticks on nice raw
numbers. The plotted positions are byte-identical (verified by `all.equal()` on
`$Z`/`$Lmat`/`$ax.one.unit`).

### Uniformity across PCA and CVA

- **PCA `standardise = TRUE`**: `sd` = column SDs of the (standardised) input;
  the transform is the textbook case above.
- **PCA `standardise = FALSE` and CVA** (`$scaled = FALSE`, so `sd = 1`): then
  `sd_new = s`. This is exactly correct: `ax.one.unit` was built for one unit of
  the internal (standardised) variable, and multiplying the label map by `s`
  supplies the missing `d(raw)/d(std) = s` factor. Forcing `$scaled = TRUE` makes
  the `Xhat` reconstruction land on the raw range so the pretty ticks are nice raw
  numbers, again without moving any plotted point.

The implementation (`.bl_rescale_biplot_axes()` in `R/plot_biplot.R`) is a no-op
when `scaling` is `NULL` (back-compatible with older `bl_result` objects) and
warns, leaving axes in standardised units, if the scaling does not cover every
biplot feature.

## Scope and limitations

- Applied to every biplot plot method: `plot_biplotEZ()`/`plot.bl_result()`,
  `plot.bl_projection()`, `plot.bl_local_result()` (the SVD-rotated local biplot --
  `$means`/`$sd` are rotation-invariant, so the relabel composes with the rotation),
  `plot.bl_surrogate()`, and `plot.bl_sparse_result()` (via pass-through).
- **Feature-value displays now also convert (2026-06-16).** The observed and
  counterfactual feature values in `plot.bl_shapley()`, `print.bl_shapley()`,
  `print.bl_sparse_result()`, and `print.bl_target()` are converted back to raw
  units at display time via the private helper `.scale_to_raw(values, scaling,
  var_names, kind)` (`R/data_prepare.R`). A **level** (observed, full CF, sparse,
  counterfactual) converts as `value * scale + center`; a **delta** (the Shapley
  `data_to_boundary` change) converts as `value * scale` only. The Shapley plot
  label reads `observed -> counterfactual` in raw units. Conversion is
  display-only -- stored objects stay in model units, mirroring the axis relabel.
  Shapley *contributions* (`shapley_cause`) are in prediction-impact units, not
  feature units, and are deliberately left unconverted. The rules for how these
  values print are in the next section.
- **One affine method for all features.** `method` is a single label; each feature
  has its own `(center, scale)` but shares the transform family. Per-variable
  transform families (z-score / min-max / log / none mixed) are future work -- and
  non-affine inverses (e.g. log) would need biplotEZ's spline-calibrated axis path
  (`PCOaxes = "splines"`) since a linear biplot axis cannot render non-linear ticks
  exactly. See `CLAUDE.md` Section 9.

## Display conventions

This section is the one place for the rules that govern how feature values and
Shapley values are displayed. The roxygen, the implementation summary, the
walkthroughs and `CLAUDE.md` point here rather than restating them.

**1. Units.** A displayed feature value is in raw units when `bl_result$scaling`
is set, and in model units otherwise. The conversion is `.scale_to_raw()`
(`R/data_prepare.R`):

| quantity | `kind` | conversion |
|---|---|---|
| a level: an observed, counterfactual, sparse or target value | `"level"` | `value * scale + center` |
| a delta: a change such as counterfactual minus observed | `"delta"` | `value * scale` |

Stored objects always stay in model units. A scaling that does not cover every
feature leaves the values in model units, with a warning. Shapley values are
prediction units and are never converted.

**2. Decimals.** Numbers print in fixed notation (never scientific, never
`-0.00`) with 2 decimals. A feature whose change would print as 0.00 gets more
decimals, 2 significant digits and at most 4, so that a small move such as a
loan-to-income change of -0.0016 stays visible. The helpers are in
`R/local_cf.R`:

- `.fmt_digits(x)` returns the decimals for each value: 2, or more for a
  non-zero value that rounds to 0 at 2;
- `.fmt_fixed(x, d)` formats with those decimals;
- `.cf_compare_table()` sets one number of decimals per feature row, from that
  row's changes (counterfactual minus observed, in display units), and
  `.format_cf_table()` applies it to every column of the row.

The stored tables stay unrounded.

**3. Where the rules apply.**

- **The observed-versus-counterfactual tables** of `print()` for
  `bl_local_result`, `bl_shapley` and `bl_sparse_result` (one builder,
  `.cf_compare_table()`), the Target and Counterfactual blocks, which take the
  table's row decimals, and the sparse path.
- **The Shapley plot labels** (`plot.bl_shapley()`) read
  `"variable: observed -> counterfactual"`, both values levels, for example
  `"Glucose: 179.00 -> 170.63"`. Each label takes its decimals from its change,
  as a table row does (`.shapley_labels()`, `R/shapley.R`). The same wording is
  used with and without a scaling; the y-axis title adds ", raw units" when the
  values are converted. `bl_shapley()` stores the labels in model units
  (`shapley_df$varnames_p`), and `plot.bl_shapley()` rebuilds them in raw units
  when a scaling is set.
- **The Shapley values**, on the plot's bars and in the `shapley` column of
  `print.bl_shapley()`, use `.fmt_fixed(s, .fmt_digits(s))`: 2 decimals, more
  for a small non-zero value.
- **Predictions** are a separate quantity: `.pred_function()` floors them to
  3 decimals (`floor(x * 1000) / 1000`), and the headers of the print methods
  show them to 4. The sparse path's `pred` column takes its decimals from
  `pred - cutoff`, so a prediction just across the cutoff stays visible.

## See also

- `2 implementation_summary.txt` Section 4.4.1 -- function-level mechanics.
- `.claude/reference/scaling_axis_relabel_plan.md` -- the approved implementation plan.
- `development/reference/review_section4_to_6.md` -- `bl_set_scaling()` in the loan workflow.
