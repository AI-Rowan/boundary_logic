> **IMPLEMENTED 2026-06-16.** Archived copy of the approved plan from
> `.claude/plans/biplot-raw-scale-axes.md`. Delivered: `bl_set_scaling()`,
> `.bl_rescale_biplot_axes()`, the `scaling` field threaded through
> `bl_filter_outliers()` / `bl_assemble()` / `bl_build_result()`, raw-unit axes
> wired into all biplot plot methods, script 03 standardisation step, tests
> (113 PASS / 0 FAIL), and docs. The *secondary* raw-unit value reporting was
> deliberately **deferred** (see CLAUDE.md §9). The *why* lives in
> `documentation/scaling_axis_relabel_note.md`; the *how* in
> `2 implementation_summary.txt` §4.4.1.

# Plan: Raw-scale biplot axes via stored feature scaling

Date created: 2026-06-16
Date last edited: 2026-06-16

## Context

When a user standardises their features **before** fitting a model (a common step to
improve fit for scale-sensitive models — logistic/SVM/NNET, and routinely done in
pipelines), the data frame handed to `boundarylogic` already holds standardised values.
Every downstream artifact — the PCA/CVA projection, the prediction grid, and crucially the
biplot **axis tick labels** — is therefore expressed in standardised units (e.g.
`person_age` reads `-0.5, 0, 0.5` instead of `25, 35, 45`). This makes reading values off
the biplot meaningless.

The fix: let the user record the scaling transform `(center, scale)` as a data-prep step,
and use it **only at plot time** to relabel biplot axis ticks back into the original
(raw) units. The biplot geometry, the projection, the grid, and the model all stay in
standardised space — nothing about the computation changes. This is a pure
display/relabel feature.

**How it works with a custom modelling step (answer to the open question):**
The stored scaling is *display-only metadata*. A custom model wrapped via
`bl_wrap_model(predict_fn = ...)` continues to receive data in **standardised** units —
exactly what it was trained on — because `bl_build_grid()` back-projects to the data's own
units (standardised) and passes those to `predict_fn`. The scaling never touches model
fitting or grid scoring; it is consumed solely by the plot methods (and, secondarily, by
the print methods that report counterfactual values). So there is zero interaction with,
or constraint on, the custom model path.

### Why this is feasible (biplotEZ mechanism — verified by source reading)

Both biplotEZ's draw path (`plotting.R:188-215`) and `axes_coordinates()`
(`calibrate_axes.R`) compute axis tick **labels** via `.calibrate.axis()`
(`plot2D.R:305-351`):

```
std.ax.tick.label <- pretty(range(Xhat[, j]), n = tick)   # nice values in Xhat units
interval          <- (std.ax.tick.label - means[j]) / sd[j]
label printed     <- interval * sd[j] + means[j]          # affine: means + sd * position
```

`Xhat` is reconstructed from `x$means`, `x$sd`, `x$scaled`, `x$center`; the **geometry**
(`x$Z`, `x$Lmat`, `x$ax.one.unit`, `x$e.vects`) is independent of those fields. A data
point's displacement along an axis (`interval`/`axis.vals`) is fixed by geometry alone.
Therefore an affine relabel of every axis is achieved purely by transforming `$means` and
`$sd` (with `$scaled`/`$center` forced `TRUE`), leaving geometry untouched.

Let the user pre-standardisation be `std = (raw - center) / scale`, i.e.
`raw = std * scale + center`. With `m0 = x$means`, `s0 = x$sd`, set:

```
means_new = m0 * scale + center
sd_new    = s0 * scale
x$scaled  = TRUE ; x$center = TRUE      # geometry fields untouched
```

Then `label_new = means_new + sd_new * axis.vals = (m0 + s0*axis.vals)*scale + center
              = label_old * scale + center = raw`, and `Xhat` reconstructs onto the raw
range so `pretty()` lands ticks on nice raw numbers. This holds **uniformly** for
PCA `standardise=TRUE` (`s0` = column SDs), PCA `standardise=FALSE` and CVA (`s0 = 1`),
because for the latter `sd_new = scale` exactly supplies the `d(raw)/d(std)` factor that
`ax.one.unit` (built for unit `std`) otherwise lacks. **Implementation must still confirm
this with the round-trip verification test below** (belt-and-suspenders against biplotEZ
version drift).

## Design

A new optional field `scaling` on the data-prep object, threaded into `bl_result`, plus a
private rescaling helper invoked by every biplot plot method.

### New field: `scaling`
`list(center = <named numeric>, scale = <named numeric>, method = <character>)` or `NULL`.
Stored on `bl_data` / `bl_filter_result`, copied into `bl_result` by `bl_assemble()`.
Distinct from `X_center`/`X_sd` (the package-internal PCA standardisation of whatever data
was supplied) and from the `standardise` flag — these are not conflated.

### New exported function: `bl_set_scaling()` — `R/data_prepare.R`
```r
bl_set_scaling(x, center, scale, method = "z-score")
```
- `x`: a `"bl_data"` or `"bl_filter_result"` object (same acceptance set as `bl_assemble()`).
- `center`, `scale`: named numeric vectors; names must cover `x$var_names`. Reordered to
  `x$var_names` internally so column order always matches the biplot.
- Returns `x` with `x$scaling` populated; class unchanged.
- Convenience: if `center`/`scale` carry the attributes produced by base `scale()`
  (`scaled:center` / `scaled:scale`) the user can pass `attr(z, "scaled:center")` etc.

### New private helper: `.bl_rescale_biplot_axes()` — `R/plot_biplot.R`
```r
.bl_rescale_biplot_axes(biplot_obj, scaling, var_names)
```
- No-op (returns `biplot_obj` unchanged) when `scaling` is `NULL`.
- Aligns `scaling$center`/`scaling$scale` to `var_names`, applies the `means_new`/`sd_new`
  transform above, sets `$scaled <- TRUE`, `$center <- TRUE`. Touches no geometry field.
- Lives beside the other private plot helpers (`.make_label_line_vec()`,
  `.make_ticks_vec()`, `.apply_biplot_rotation()`).

### Wiring into all biplot plot methods
In each method, fetch the biplot object as today, then **after** any rotation
(`.apply_biplot_rotation()` — `$means`/`$sd` are unaffected by rotation) and **after** the
`new_title` override, call `.bl_rescale_biplot_axes(obj, scaling, var_names)` before the
first `plot()`/overlay. Source the scaling from the reachable `bl_result`:

| Method | File | Scaling source |
|---|---|---|
| `plot_biplotEZ()` / `plot.bl_result()` | `R/plot_biplot.R` | `bl_result$scaling` |
| `plot.bl_projection()` | `R/projection.R` | `x$scaling` (add field to `bl_projection`; see note) |
| `plot.bl_local_result()` | `R/local_cf.R` | `x$bl_result$scaling` |
| `plot.bl_surrogate()` | `R/surrogate.R` | `x$bl_result$scaling` |
| `plot.bl_sparse_result()` | `R/shapley.R` (forwards to local) | inherited |

Note on `plot.bl_projection()`: the exploratory biplot (`bl_build_result(bl_model = NULL)`)
returns a `bl_projection`, which is built before any `bl_data` is in scope inside
`bl_build_projection()`. Simplest: have `bl_assemble()`/`bl_build_result()` copy `scaling`
onto the projection too, OR accept that the no-model exploratory biplot shows standardised
axes (acceptable, since scaling is most meaningful once a model exists). **Decision: copy
`scaling` onto `bl_projection` via a small assignment in `bl_build_result()` after the
projection is built**, so the exploratory `plot()` also relabels.

### Secondary (recommended for coherence): raw-unit reporting of values
Axis labels in raw units but printed counterfactual/target tables in standardised units
would be confusing. Add raw-unit display to the value-reporting paths, gated on `scaling`:
- `plot_biplotEZ()` Step 6 target-point label, and `print.bl_local_result()` /
  `print.bl_sparse_result()` / the Shapley value table — inverse-transform displayed
  feature values via `value * scale + center`.
This is a clearly separable step; it can be deferred without affecting the core axis fix,
but is included here so the feature is internally consistent.

## Files to read first
- `R/data_prepare.R` (add setter; `bl_data`/print)
- `R/result.R` (`bl_assemble()` field list; `bl_build_result()` projection copy)
- `R/projection.R` (`bl_build_projection`, `plot.bl_projection`, `print.bl_projection`)
- `R/plot_biplot.R` (private helpers + `plot_biplotEZ()` render order)
- `R/local_cf.R`, `R/surrogate.R`, `R/shapley.R` (sub-object plot methods; all carry `bl_result`)
- `R/outlier_filter.R` (confirm `bl_filter_outliers()` copies through any new `bl_data` field — `scaling` must survive filtering)

## Implementation steps
1. `bl_set_scaling()` + validation (`stop_if_*` helpers, `call. = FALSE`) and roxygen
   (`@export`). Accept `bl_data` and `bl_filter_result`.
2. Ensure `bl_filter_outliers()` preserves `$scaling` when it rebuilds the object (so the
   `bl_prepare_data` path, which filters internally, keeps a scaling set afterward).
3. `.bl_rescale_biplot_axes()` private helper in `R/plot_biplot.R`.
4. Add `scaling` to the `bl_assemble()` return list (read from `bl_data$scaling`); document
   in the Value section. Copy `scaling` onto `bl_projection` in `bl_build_result()`.
5. Wire the helper into the five plot methods per the table.
6. (Secondary) raw-unit value display in target-label + print/Shapley paths.
7. Update `print.bl_data` / `print.bl_result` to note when a scaling is attached.
8. `devtools::document()` + `devtools::test()`; then the verification script.

## Function design notes
- `.bl_rescale_biplot_axes()` must match `scaling` vectors to the biplot's variable order
  (`var_names`), not to `names(biplot_obj$means)` (which may be `V1..Vp`).
- Never mutate the caller's object in place beyond the local copy used for rendering;
  return a modified copy (consistent with `.apply_biplot_rotation()`).

## Error cases
- `center`/`scale` not numeric, differing lengths, or names not covering `var_names` -> stop.
- `scale` containing `0` (or `NA`) -> stop (division/inversion would be invalid).
- `bl_set_scaling()` on a wrong class -> stop with the `bl_assemble()`-style message.
- Older `bl_result` with no `scaling` field -> helper no-ops; all plots behave as today
  (back-compatible; no warning needed since absence is the normal default).

## Test cases (`tests/testthat/`)
- `bl_set_scaling()` stores correctly, reorders to `var_names`, rejects bad inputs.
- Round-trip: build a `bl_result` on standardised iris with known `(center, scale)`;
  assert `axes_coordinates()` (or `.calibrate.axis` label column) on the rescaled biplot
  equals `raw = std*scale + center` at matching positions, for **PCA standardise=TRUE,
  PCA standardise=FALSE, and CVA**. Assert `$Z`/`$Lmat`/`$ax.one.unit` are byte-identical
  pre/post rescale (geometry invariance).
- `scaling = NULL` path identical to current behaviour (snapshot of one plot's
  `axes_coordinates()`).

## Script 03 changes (`scripts/03_loan_status_Boundary_Logic.R`)
Add a standardisation step to the existing reduced-model (v2) flow:
1. Before fitting `xgb_fit_v2`, standardise the `feature_cols_v2` columns of the training
   features (`z <- scale(...)`), capturing `center`/`scale`.
2. Fit XGB on the standardised features (note in a comment: XGB is scale-invariant
   per-feature, so this is illustrative — the payoff is the raw-unit biplot axes, which the
   projection/CVA geometry *does* depend on).
3. `bl_dat_v2 <- bl_set_scaling(bl_dat_v2, center = ..., scale = ..., method = "z-score")`.
4. Show `plot(bl_results_v2)` now renders axes in raw loan units; keep one commented
   before/after note.

## Documentation to update (mandatory — part of this change, not optional follow-up)

Per the CLAUDE.md S2 triplet rule, this is a methodologically significant change, so the
documentation updates below ship **with** the code, not afterward:

- **`2 implementation_summary.txt`** (REQUIRED) — add: (a) `bl_set_scaling()` to the
  function inventory; (b) the new `scaling` field to the `bl_data` / `bl_filter_result` /
  `bl_result` / `bl_projection` S3 field lists, explicitly distinguished from `X_center` /
  `X_sd` / `standardise`; (c) a subsection describing the means/sd affine-relabel mechanism
  and its geometry-invariance, cross-referencing the biplotEZ `.calibrate.axis()` path.
  This is the canonical "how the code works" record — do not skip it.
- **`.claude/reference/review_section4_to_6.md`** (REQUIRED) — covers `data_prepare.R`,
  `projection.R`, `result.R`, `plot_biplot.R`: document `bl_set_scaling()`, the `scaling`
  field flow (set on data-prep object -> `bl_assemble()` -> `bl_result` / `bl_projection`),
  and `.bl_rescale_biplot_axes()` invoked in the render path. (Update-trigger rule: these
  source files are changing.)
- **`.claude/reference/review_section9_to_15.md`** (REQUIRED) — note the raw-axis relabel
  in the local / surrogate / sparse plot methods (these source files change too) and the
  secondary raw-unit value reporting in the print/Shapley paths.
- `documentation/` — short technical note on the affine relabel derivation and why
  geometry is invariant (the *why*, completing the triplet alongside the plan archive).
- `CLAUDE.md` — three edits:
  1. R/ file map entry for `bl_set_scaling()` and the new `.bl_rescale_biplot_axes()`
     private helper.
  2. An "Always Do" (S5) note that every biplot plot method must call
     `.bl_rescale_biplot_axes()` when `bl_result$scaling` is present (mirrors the existing
     `new_title` convention bullet).
  3. A new "Deferred / Future Work" (S9) entry — see the next section.
- Memory sync: if any updated reference doc has a memory mirror in
  `~/.claude/projects/.../memory/`, apply the same change there (CLAUDE.md S2 sync rule).

## Future work to record in CLAUDE.md Section 9: per-variable standardisation methods

The current change supports a **per-variable** `(center, scale)` pair, but a **single
affine method** applied uniformly across all features (`method` is one scalar label, e.g.
`"z-score"`). Add a Section 9 deferred-work bullet capturing the natural next step:

> **Variable-specific standardisation methods.** Allow a *different* transform family per
> feature (e.g. z-score for some, min-max for others, log / Box-Cox for skewed features,
> and "none" for already-interpretable features) rather than one affine method for all.
> This needs `scaling` generalised from two numeric vectors + one method label to a
> per-variable transform spec carrying each feature's forward and inverse function (so the
> biplot relabel can invert non-affine transforms). Note: non-affine inverses (e.g. log)
> will produce **non-linearly spaced** axis ticks; the biplotEZ axis is linear, so the
> relabel would need either spline-calibrated axes (biplotEZ `PCOaxes = "splines"` path) or
> a documented limitation that only affine transforms get exact raw-unit ticks. Keep the
> current uniform-affine `bl_set_scaling()` as the simple default.

(Recorded here so it is not implemented now without a fresh instruction.)

## Verification
- Write `verify_scaling.R` (invoke via `Rscript verify_scaling.R`, not `-e`, per the
  biplotEZ Windows-segfault rule). It: builds a standardised-iris `bl_result` for PCA and
  CVA, attaches a known scaling, and prints `axes_coordinates()` label ranges before vs
  after rescale to confirm raw-unit ticks and unchanged geometry. Delete after use (diff
  against the session-start `git status` snapshot before any cleanup).
- `devtools::document()` clean; `devtools::test()` >= current 77 PASS / 0 FAIL with the new
  tests added. Report PASS/FAIL and stop (no commit/push without explicit instruction).
