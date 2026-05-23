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
- `documentation/mahalanobis_technical_note.md` — theoretical motivation, mathematical derivation, and trade-offs for the Mahalanobis-aware distance measures (`bl_find_local_cf()` cross-pair selector and `bl_robustness()` / `plot.bl_boundary()` per-feature denominator). Read this before changing how `metric`/`metric_inv` are computed, propagated, or consumed.

### Implementation plan archive

- `.claude/reference/mahalanobis_implementation_plan.md` — the approved plan that drove the Mahalanobis migration (2026-05-20). Records *what* was done; the technical note above records *why*. Reference document; not used during code execution.

### Documentation folder convention

The repo has three places for narrative documentation; each has a specific role.

| Location | Purpose | When to add a new file |
|---|---|---|
| Project root (numbered: `1 ...`, `2 ...`) | Foundational design docs that take precedence on methodological questions. | Almost never — reserved for the foundation/implementation-summary canon. |
| `documentation/` | Methodological technical notes, the PhD thesis PDF, workflow HTML exports, presentations. | When you implement a methodologically significant change, write a technical note here explaining *why* (use plain, un-numbered filenames). |
| `.claude/reference/` | Code walkthrough docs (per source-file or per-script) and archived implementation plans (after approval + execution). Tracked by git. | When a plan from `.claude/plans/` has been implemented, copy it here with an IMPLEMENTED banner. |

For methodologically significant changes, the standard triplet is:
1. **Technical note** in `documentation/` — explains *why* (theory, derivation, trade-offs).
2. **Plan archive** in `.claude/reference/` — explains *what was done* (the approved plan, verbatim).
3. **Implementation summary update** to `2 implementation_summary.txt` — explains *how the code works* (function-level mechanics).

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
| `pick_point.R` | `bl_pick_point()` — interactive biplot point picker; optionally draws boundary counterfactual arrow per click via `bl_boundary =` |

## 4. Critical Constraints — Never Do

- **Never simplify the maths to make code shorter.** Replace mathematical logic only when outputs are proven equivalent.
- **Never alter default scaling, centring, or projection conventions** without explicit instruction.
- **Never apply the Z-space hull polygon inside `bl_find_local_cf()`.** The SVD rotation invalidates the original polygon coordinates. Only `train_ranges` and `set_filters` constrain the local search.
- **Never entangle plotting code with computation code.** `plot_biplotEZ()` renders; it does not compute boundary points or projections.
- **Never remove metadata fields** from S3 objects (`V`, `tV`, `X_center`, `X_sd`, `train_ranges`, etc.). These are required for inverse projection and downstream analysis.
- **Never change the `cutoff` default from 0.5.** Other values are accepted but are methodologically invalid under the current method.
- **Never change the XGB model format.** When calling `bl_wrap_model(model_type = "XGB")`, the `model` argument must be `list(model = <xgb.Booster>, features = <character vector>)`. XGB is no longer supported by `bl_fit_model()`.
- **Never expose `outlie`, `calc_hull`, or `bl_robustness()` as sequential steps in examples or vignettes.** These are internal or redundant.
- **Never add a `boundary =` parameter back to `plot_biplotEZ()`.** Counterfactual arrows belong exclusively in `bl_pick_point(bl_result, bl_boundary = bl_bnd)` — drawing n arrows at plot-flush time caused unsuppressable zero-length arrow warnings.
- **Never use literal non-ASCII characters in R source files.** Use ASCII equivalents (`--` for em dash, `x` for ×, `>=` for ≥) or `\uXXXX` escape sequences in code strings. The Edit tool cannot reliably write `\uXXXX` escapes — use a Python one-liner via Bash for bulk replacement.

## 5. Always Do

- **After implementing any plan: stop after `devtools::document()` and `devtools::test()` pass.** Report the PASS/FAIL count and wait. Do NOT proceed to `git add`, `git commit`, or `git push` unless the user explicitly asks — even in auto-accept mode.
- Add `roxygen2` documentation (`#' @param`, `#' @return`, `#' @export`) to every exported function.
- Validate inputs with the `stop_if_*()` helpers in `utils.R`, always with `call. = FALSE`.
- Name the assembled pipeline result `bl_results` in scripts (not `result` or `bl_result`).
- Name the Shapley object `bl_shapley_values` in scripts (not `bl_shap`, to avoid confusion with the Python SHAP library).
- Run `devtools::test()` after any change to a core function.
- Confirm changes against both the iris and Pima datasets when modifying core functions.
- Before modifying a core function: read its roxygen header, grep for callers, check what `bl_result` fields flow from it, then make a minimal change without refactoring surrounding code.
- Every new base-R plotting function (one that calls `graphics::*` directly) must open a graphics device if none is active: `if (grDevices::dev.cur() == 1L) grDevices::dev.new()` at the top of the function body, after input validation. Add `@importFrom grDevices dev.cur dev.new` to its roxygen block.
- **Distance-measure functions** use the `distance = c("mahalanobis", "euclidean")` parameter convention with `"mahalanobis"` as the default. The Mahalanobis path reads `bl_result$metric_inv` (full `W^{-1}`) or `bl_result$metric` (raw `W`, e.g. `sqrt(diag(W))` for per-feature standardisation); the Euclidean path preserves pre-Mahalanobis behaviour for reproducibility. Always include the `is.null(bl_result$metric_inv)` fallback-with-warning to handle older `bl_result` objects.
- **Matrix inversion of within-class scatter or any SPD covariance** must use `chol2inv(chol(M))` (not `solve(M)`) with a `tryCatch` ridge fallback `chol2inv(chol(M + lambda * I))` (`lambda = 1e-6 * mean(diag(M))`). Rationale: numerical stability on SPD matrices, automatic singularity detection via `chol()`, ~2x speed. See `documentation/mahalanobis_technical_note.md` Section 3.2.
- **Verification scripts that touch biplotEZ** must be written to a `.R` file and invoked via `Rscript verify.R`, not `Rscript -e '...'`. The `-e` mode segfaults inside `biplotEZ::PCA()` / `biplotEZ::CVA()` on Windows. Tests under `devtools::test()` are unaffected.

## 6. Key Architectural Facts

- `b_margin` controls the decision boundary contour band half-width directly (default `0.001`, valid range `(0, 0.5)`). It is stored on `bl_grid$b_margin` and propagates to `bl_result$b_margin`. Predictions are independently stored at 3 d.p. via `floor(x * 1000) / 1000` in `.pred_function()`. These are separate concerns.
- The hull polygon in `bl_result` always comes from `bl_grid` (final Z-space), never from `bl_filter_result` (raw PCA space).
- CVA forces `standardise = FALSE` internally inside `bl_build_projection()`, regardless of the user-supplied value.
- In `bl_find_boundary()`, `train_ranges` are applied post-selection (after the nearest boundary point is identified), not as a pre-filter on contour segments.
- In `bl_find_local_cf()`, `set_filters` (actionability) runs before `train_ranges` (feasibility) within each candidate eigenvector-pair loop. Actionability constraints are typically more selective, so running them first enables per-segment early exit (`if (!any(keep)) next`). Filter 4 (model re-score) is not batched -- it stays per-segment to avoid O(n) vertex accumulation across pairs.
- `bl_build_result()` always returns a `bl_result` object (class `c("bl_result", "list")`), even when `bl_model = NULL`. In that case the model and grid fields are NULL but the object is still a `bl_result`, not a `bl_projection`.
- `.poly_clip()` exists in both `hull_utils.R` and `boundary.R`. The authoritative copy is in `hull_utils.R`. Do not add a third copy.
- `bl_build_projection()` computes and stores the within-class metric matrix `W` (`metric`), its inverse via Cholesky (`metric_inv`), and the `metric_type` tag. These propagate through `bl_assemble()` into `bl_result`. The metric is the same `W` used by CVA's eigenproblem -- distances measured with it are consistent with the biplot axes. See `2 implementation_summary.txt` Section 4.2.1.
- `bl_find_local_cf()` selects the best eigenvector pair via squared Mahalanobis distance in X-space by default (`distance = "mahalanobis"`). The pre-existing Z-space Euclidean behaviour is available as `distance = "euclidean"`. Within-pair selection (nearest valid contour vertex) is unchanged -- still 2D Euclidean.
- `plot.bl_boundary()` and `bl_robustness()` use `sqrt(diag(W))` (within-class SD per feature) as the per-variable standardisation denominator by default (`distance = "mahalanobis"`). This correctly amplifies features that are good class separators. The legacy `X_sd` (total SD) denominator is available as `distance = "euclidean"`. A cross-correlation diagnostic is printed when Mahalanobis is used; see `2 implementation_summary.txt` Section 4.6.1 for why the full Shapley decomposition is theoretically better but not used in Phase 2.

## 7. Adding a New Model Type

Two distinct paths depending on where the type should be supported.

**`bl_wrap_model()` only** (externally fitted model, complex training, or non-parsnip engine):
1. Add the type string to `valid_types` in `bl_wrap_model()` (`R/model_fit.R`)
2. Add a prediction branch in `.pred_function()` (`R/predict_utils.R`)
3. Test with the iris two-class binary dataset

**`bl_fit_model()` + `bl_wrap_model()`** (simple parsnip-compatible type; current four are GLM, SVM, NNET, RForrest):
1. Add the type string to `valid_types` in both `bl_fit_model()` and `bl_wrap_model()` (`R/model_fit.R`)
2. Add a fitting branch in `.fit_model()` (`R/model_utils.R`)
3. Add a prediction branch in `.pred_function()` (`R/predict_utils.R`)
4. Test with the iris two-class binary dataset
5. Add the new type to the supported model types table in the `bl_fit_model()` roxygen header

## 8. What Not to Show in Examples or Vignettes

- `outlie` parameter in `bl_build_grid()` — internal visual tuning only
- `calc_hull` parameter in `bl_build_grid()` — visual-only
- `bl_robustness()` called after `plot(bl_bnd)` — redundant; `plot(bl_bnd)` already prints the robustness summary to the console
- `boundary =` in `plot_biplotEZ()` — this parameter no longer exists; the correct pattern is `bl_pick_point(bl_results, bl_boundary = bl_bnd)` in an interactive session (shown with `eval=FALSE` in vignettes)

## 9. Deferred / Future Work

Do not implement these without a new instruction:

- Mathematical method vignette (derives biplot projection and counterfactual geometry from first principles)
- Shapley attribution of full Mahalanobis distance in Phase 2 (`plot.bl_boundary()` / `bl_robustness()`). Currently uses diagonal-only Mahalanobis (`sqrt(diag(W))` per feature) as a deliberate `O(n)` trade-off. Full Shapley would be `O(n * 2^p)`. See `2 implementation_summary.txt` Section 4.6.1.
- Unit tests for Phase 3 functions (`bl_select_target`, `set_filters`, `bl_find_local_cf`, `bl_shapley`, `bl_find_sparse_cf`)
- Unit tests for `bl_build_result()` and `bl_assemble()`
- CRAN submission preparation
- Per-variable label direction in `plot_biplotEZ()` — currently `label_dir` accepts only a single scalar (`"Hor"` or `"Orthog"`) applied to all labels, as biplotEZ::axes() does not support per-variable direction. A future enhancement would add a `label_dir_var` vector parameter with a second-pass redraw for specific variables.
- ~~GAM fitting via `bl_fit_model()`~~ — **REMOVED.** The parsnip/workflows two-formula workaround was broken and unmaintainable; GAM support was removed from `bl_fit_model()` entirely. Use `bl_wrap_model()` with `mgcv::gam()` directly (see `scripts/00_pima_Boundary_Logic.R` Step 3 for the pattern).
- LDA multi-class biplots (k > 2 classes) — LDA naturally produces one discriminant axis per class boundary and could support multi-class interpretation. The current package is restricted to binary 0/1 outcomes. A future enhancement could add an LDA path in `bl_wrap_model()` (using `MASS::lda()` directly) and a dedicated biplot variant where each class pair yields a separate boundary. See `2 implementation_summary.txt` for the design note.
- ~~**`.bl_rotate()` bug** — when `best_pair != c(1, 2)`, the target point appears at the biplot origin~~ — **RESOLVED.** `.bl_rotate()` now returns `Vrho[, c(1L, 2L)]` and `tVrho[c(1L, 2L), ]` unconditionally. The SVD construction guarantees target information concentrates in columns 1-2 of `Vrho` regardless of `proj_pair`.
- ~~Zero-length arrow warning in `plot_biplotEZ()` boundary overlay~~ — **RESOLVED.** The `boundary`, `show_arrows`, and `arrow_col` parameters have been removed from `plot_biplotEZ()`. Arrow drawing has moved to `bl_pick_point(bl_result, bl_boundary = bl_bnd)`, where it fires once per interactively selected observation. The global n-arrow overlay (which fired the warning n times at plot-flush time) is gone. The per-click code in `bl_pick_point()` guards with `all(is.finite(bz_row))` before drawing; if an observation is exactly on the boundary the warning may fire at most once per click, in interactive mode. Known limitation: if the biplot was rendered with `rotate_deg != 0`, the arrow will be misaligned (boundary coordinates are in unrotated space). Use `rotate_deg = 0` when combining `bl_pick_point()` with boundary overlay.

## 10. Branches and Repo

- `main` — stable, public-facing
- `method_developments` — active development

**Always start every session on `method_developments`.** Run `git checkout method_developments` at the beginning of any session before making changes.

Develop on `method_developments`, merge to `main` when stable. Both branches push to `origin` (GitHub: `AI-Rowan/boundary_logic`).
