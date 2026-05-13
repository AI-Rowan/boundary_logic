# CLAUDE.md — boundarylogic

## 1. Project Identity

`boundarylogic` is an R package that operationalises a PhD methodology for post-hoc interpretability of binary classification models. It projects high-dimensional data into a 2D biplot (PCA or CVA), approximates the model's decision boundary in that space, and finds counterfactuals — the minimum change to an observation that flips the prediction.

**Three-phase workflow:**
- Phase 1: data prep → model → biplot projection → prediction grid → `bl_result` (the central anchor object)
- Phase 2 (global): nearest boundary point per observation, surrogate model, robustness
- Phase 3 (local): single-observation biplot rotation → counterfactual → Shapley attribution → sparse counterfactual

The `bl_result` object is what everything downstream consumes.

## 2. Reference Documents

Before modifying any core function, read the relevant documents below.

### Generic references (function-level, not script-specific)

- `1 Foundation intro documents.txt` (project root) — mathematical intent, design principles, and coding standards. Takes precedence on all methodological questions.
- `2 implementation_summary.txt` (project root) — function inventory, S3 class fields, and design decisions. **Update when any function signature, return object, or S3 field changes.**

### Script-specific references (tied to `scripts/03_loan_status_Boundary_Logic.R`)

Detailed walkthroughs tracing every function call, parameter, and returned object field for the loan-default workflow. Read the relevant document before modifying the covered source files. **Update the document when the covered source files change.**

| Document | Steps | Source files covered |
|---|---|---|
| `.claude/reference/review_section4_to_6.md` | Steps 3–6 (Phase 1) | `data_prepare.R`, `outlier_filter.R`, `model_fit.R`, `model_utils.R`, `result.R`, `biplot_grid.R`, `projection.R`, `plot_biplot.R`, `project_points.R` |
| `.claude/reference/review_section7_to_8.md` | Steps 7–8 (Phase 2) | `boundary.R`, `boundary_plot.R`, `hull_utils.R`, `feasibility_utils.R` |
| `.claude/reference/review_section11_to_17.md` | Steps 11–17 (Phase 3) | `local_cf.R`, `shapley.R`, `project_points.R` |

Script: `scripts/03_loan_status_Boundary_Logic.R`

## 3. R/ File Map

| File | Responsibility |
|---|---|
| `utils.R` | Input validators (`stop_if_*`), `calc_gini()` |
| `hull_utils.R` | Convex hull construction and polygon clipping |
| `feasibility_utils.R` | Training range extraction and row-wise filtering |
| `model_utils.R` | Internal model fitting dispatcher (`.fit_model()`) |
| `predict_utils.R` | Unified prediction interface (`.pred_function()`) |
| `data_prepare.R` | `bl_prepare_data()`, `bl_wrap_data()` |
| `outlier_filter.R` | `bl_filter_outliers()` |
| `model_fit.R` | `bl_fit_model()`, `bl_wrap_model()` |
| `projection.R` | `bl_build_projection()` — PCA / CVA loading matrix |
| `biplot_grid.R` | `bl_build_grid()` — m×m prediction grid and contours |
| `result.R` | `bl_assemble()`, `bl_build_result()` — Phase 1 anchor |
| `plot_biplot.R` | `plot_biplotEZ()` — main biplot renderer |
| `project_points.R` | `bl_project_points()`, `bl_predict()` |
| `boundary.R` | `bl_find_boundary()` — global counterfactual search |
| `boundary_plot.R` | `plot.bl_boundary()`, `bl_robustness()` |
| `surrogate.R` | `bl_surrogate()` |
| `local_cf.R` | `bl_select_target()`, `set_filters()`, `bl_find_local_cf()` |
| `shapley.R` | `bl_shapley()`, `bl_find_sparse_cf()` |
| `pick_point.R` | `bl_pick_point()` — interactive biplot point picker |

## 4. Critical Constraints — Never Do

- **Never simplify the maths to make code shorter.** Replace mathematical logic only when outputs are proven equivalent.
- **Never alter default scaling, centring, or projection conventions** without explicit instruction.
- **Never apply the Z-space hull polygon inside `bl_find_local_cf()`.** The SVD rotation invalidates the original polygon coordinates. Only `train_ranges` and `set_filters` constrain the local search.
- **Never pass `rounding` into `.pred_function()`.** It controls only the contour band width (`b_margin`). Predictions are always floor-rounded to 3 d.p. inside `.pred_function()` regardless.
- **Never entangle plotting code with computation code.** `plot_biplotEZ()` renders; it does not compute boundary points or projections.
- **Never remove metadata fields** from S3 objects (`V`, `tV`, `X_center`, `X_sd`, `train_ranges`, etc.). These are required for inverse projection and downstream analysis.
- **Never change the `cutoff` default from 0.5.** Other values are accepted but are methodologically invalid under the current method.
- **Never change the XGB model format.** The `bl_model$model` field for XGB must be `list(model = <xgb.Booster>, features = <character vector>)`.
- **Never expose `outlie`, `calc_hull`, or `bl_robustness()` as sequential steps in examples or vignettes.** These are internal or redundant.

## 5. Always Do

- Add `roxygen2` documentation (`#' @param`, `#' @return`, `#' @export`) to every exported function.
- Validate inputs with the `stop_if_*()` helpers in `utils.R`, always with `call. = FALSE`.
- Name the assembled pipeline result `bl_results` in scripts (not `result` or `bl_result`).
- Name the Shapley object `bl_shapley_values` in scripts (not `bl_shap`, to avoid confusion with the Python SHAP library).
- Run `devtools::test()` after any change to a core function.
- Confirm changes against both the iris and Pima datasets when modifying core functions.
- Before modifying a core function: read its roxygen header, grep for callers, check what `bl_result` fields flow from it, then make a minimal change without refactoring surrounding code.

## 6. Key Architectural Facts

- `rounding` controls only the contour band width (`b_margin = 1 / 10^rounding`). Predictions are always stored at 3 d.p. via `floor(x * 1000) / 1000` in `.pred_function()`. These are separate concerns.
- The hull polygon in `bl_result` always comes from `bl_grid` (final Z-space), never from `bl_filter_result` (raw PCA space).
- CVA forces `standardise = FALSE` internally inside `bl_build_projection()`, regardless of the user-supplied value.
- In `bl_find_boundary()`, `train_ranges` are applied post-selection (after the nearest boundary point is identified), not as a pre-filter on contour segments.
- `bl_build_result()` returns a `bl_projection` object (not `bl_result`) when `bl_model = NULL`.
- `.poly_clip()` exists in both `hull_utils.R` and `boundary.R`. The authoritative copy is in `hull_utils.R`. Do not add a third copy.

## 7. Adding a New Model Type

Touches three files — do all three or do none:

1. Add the type string to `valid_types` in `bl_fit_model()` and `bl_wrap_model()` (`R/model_fit.R`)
2. Add a fitting branch in `.fit_model()` (`R/model_utils.R`)
3. Add a prediction branch in `.pred_function()` (`R/predict_utils.R`)
4. Test with the iris two-class binary dataset
5. Add the new type to the supported model types table in the `bl_fit_model()` roxygen header

## 8. What Not to Show in Examples or Vignettes

- `outlie` parameter in `bl_build_grid()` — internal visual tuning only
- `calc_hull` parameter in `bl_build_grid()` — visual-only
- `bl_robustness()` called after `plot(bl_bnd)` — redundant; `plot(bl_bnd)` already prints the robustness summary to the console

## 9. Deferred / Future Work

Do not implement these without a new instruction:

- Mathematical method vignette (derives biplot projection and counterfactual geometry from first principles)
- Mahalanobis distance-to-boundary alternative (see `memory/future_mahalanobis_distance.md`)
- Unit tests for Phase 3 functions (`bl_select_target`, `set_filters`, `bl_find_local_cf`, `bl_shapley`, `bl_find_sparse_cf`)
- Unit tests for `bl_build_result()` and `bl_assemble()`
- CRAN submission preparation
- Per-variable label direction in `plot_biplotEZ()` — currently `label_dir` accepts only a single scalar (`"Hor"` or `"Orthog"`) applied to all labels, as biplotEZ::axes() does not support per-variable direction. A future enhancement would add a `label_dir_var` vector parameter with a second-pass redraw for specific variables.
- GAM fitting via `bl_fit_model()` — the parsnip/workflows two-formula workaround is broken in some configurations. Until fixed, use `bl_wrap_model()` with `mgcv::gam()` directly (see `scripts/00_pima_Boundary_Logic.R` Step 3).
- Zero-length arrow warning in `plot_biplotEZ()` boundary overlay — the warning "zero-length arrow is of indeterminate angle and so skipped" is generated by R's graphics device at deferred flush time (not at call time), making it resistant to `suppressWarnings()` and `withCallingHandlers()`. Three attempted fixes all failed: (1) pre-filter zero-length rows before `graphics::arrows()`, (2) `suppressWarnings()` on biplotEZ `plot()` calls, (3) `withCallingHandlers()` in `plot.bl_boundary()`. The warning fires in the graphics device layer outside the R call stack. **Recommended alternative to investigate:** change the `show_arrows` default to `FALSE` in `plot_biplotEZ()` and remove arrow drawing from the boundary overlay entirely — at large observation counts (1000+) the crosses at boundary points suffice and the arrows add visual noise.

## 10. Branches and Repo

- `main` — stable, public-facing
- `method_developments` — active development

**Always start every session on `method_developments`.** Run `git checkout method_developments` at the beginning of any session before making changes.

Develop on `method_developments`, merge to `main` when stable. Both branches push to `origin` (GitHub: `AI-Rowan/boundary_logic`).
