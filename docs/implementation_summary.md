# Boundary Logic — Implementation Summary

**Package:** `boundarylogic`
**Branch:** `feature/phase3-local`
**Last updated:** 2026-03-17

---

## 1. Purpose and Scope

The `boundarylogic` package operationalises a PhD research methodology for post hoc interpretability of binary classification models. The core idea is that a two-dimensional biplot projection of high-dimensional data preserves enough geometric structure to:

- approximate model decision boundaries in the reduced space;
- identify the nearest boundary-crossing point (counterfactual) for any observation;
- explain which variable changes drive the prediction across the boundary; and
- support both local (single observation) and global (full dataset) interpretability analyses.

This document records: the package file structure, the function inventory, how each component maps to the instruction manual modules A–G, and the key design decisions made during implementation.

---

## 2. Package File Structure

```
boundarylogic/
├── DESCRIPTION              # Package metadata, dependencies
├── NAMESPACE                # Exports and importFrom directives (roxygen-managed)
├── R/
│   ├── utils.R              # Input validation helpers, Gini coefficient
│   ├── hull_utils.R         # Convex hull polygon construction and clipping
│   ├── feasibility_utils.R  # Training range extraction and filtering
│   ├── model_utils.R        # Internal model fitting dispatcher
│   ├── predict_utils.R      # Unified prediction interface across model types
│   ├── data_prepare.R       # Module A: Data ingestion and class conversion
│   ├── outlier_filter.R     # Module A: Iterative polygon outlier filter
│   ├── model_fit.R          # Module A: Model fitting and wrapping
│   ├── projection.R         # Module B: PCA/CVA projection matrix
│   ├── biplot_grid.R        # Module C: Prediction grid and contours
│   ├── result.R             # Phase 1 assembly object
│   ├── plot_biplot.R        # Phase 1 biplot visualisation
│   ├── project_points.R     # Projection utility + bl_predict()
│   ├── boundary.R           # Module D (global): Counterfactual search
│   ├── boundary_plot.R      # Module E: Distance plots and robustness
│   ├── surrogate.R          # Module F: Surrogate model
│   ├── local_cf.R           # Module B+D (local): Rotation, target, CF search
│   └── shapley.R            # Module G: Shapley values and sparse CF
├── scripts/
│   ├── 00_pima_GAM_CVA.R           # Full workflow: Steps 1–18
│   ├── 00_pima_phase3_local.R      # Standalone Phase 3 walkthrough
│   ├── 00_pima_external_test_comparison.R
│   ├── 00_pima_phase2_unlabeled.R
│   └── 00_inspect_biplot_CVA.R
└── docs/
    └── implementation_summary.md   # This document
```

---

## 3. S3 Class Hierarchy

Each phase produces a named S3 object. Later phases consume objects from earlier phases, forming a chain:

```
bl_data
  └─► bl_filter_result
        └─► bl_model
              └─► bl_projection
                    └─► bl_grid
                          └─► bl_result          ← Phase 1 anchor
                                ├─► bl_points     (bl_project_points / bl_predict)
                                │
                                ├── PHASE 2 (global) ─────────────────────────
                                ├─► bl_boundary   (bl_find_boundary)
                                ├─► bl_surrogate  (bl_surrogate)
                                │
                                └── PHASE 3 (local) ──────────────────────────
                                    ├─► bl_target        (bl_select_target)
                                    ├─► bl_filters       (set_filters)
                                    ├─► bl_local_result  (bl_find_boundary_local)
                                    ├─► bl_shapley       (bl_shapley)
                                    └─► bl_sparse_result (bl_sparse_cf)
```

All objects have `print()` methods. Objects with visual outputs have `plot()` methods.

---

## 4. Module-by-Module Function Inventory

### 4.1 Module A — Data Preprocessing and Model Fitting
*Instruction manual §3A*

**`R/data_prepare.R`**

| Function | Signature | Returns |
|---|---|---|
| `bl_prepare_data()` | `(data, class_col, target_class=NULL, feature_cols=NULL, train_fraction=0.8, seed=121L)` | `bl_data` |
| `bl_wrap_data()` | `(train_data, test_data=NULL, var_names=NULL, target_class=NULL)` | `bl_data` |

`bl_prepare_data()` handles the full ingestion path: converts multiclass factors to binary 0/1 indicators, stratified train/test split, stores `var_names`, `class_col`, and `target_class` metadata. `bl_wrap_data()` accepts pre-split data frames directly.

**`R/outlier_filter.R`**

| Function | Signature | Returns |
|---|---|---|
| `bl_filter_outliers()` | `(bl_data, hull_fraction=0.9, verbose=TRUE)` | `bl_filter_result` |

Implements the iterative convex hull polygon filter. At each iteration: removes points outside the hull, refits the model, reports Gini coefficient. The fraction of the minimum bounding box retained by the hull is controlled by `hull_fraction`. The retained training data from the final iteration is stored for downstream use.

**`R/model_fit.R`**

| Function | Signature | Returns |
|---|---|---|
| `bl_fit_model()` | `(train_data, var_names, model_type="GLM", cutoff=0.5, rounding=2L, model_params=list())` | `bl_model` |
| `bl_wrap_model()` | `(model, model_type, var_names, predict_fn=NULL, train_data=NULL, cutoff=0.5, rounding=2L)` | `bl_model` |

`bl_fit_model()` supports: GLM, GAM, GBM, SVM, LDA, NNET, decision tree (via parsnip/tidymodels or direct calls). `bl_wrap_model()` accepts any external fitted model with a custom `predict_fn`, enabling GAMs (mgcv), XGBoost, or any user-defined classifier. Both store `model_type`, `cutoff`, `rounding`, and `var_names` for consistent downstream use.

---

### 4.2 Module B — Projection and Biplot Construction
*Instruction manual §3B*

**`R/projection.R`**

| Function | Signature | Returns |
|---|---|---|
| `bl_build_projection()` | `(train_data, var_names, method="PCA", standardise=TRUE, proj_dims=c(1L,2L), bl_model=NULL, cva_classes=NULL, title="")` | `bl_projection` |

Constructs the loading matrix `V` (p × p) and its inverse `tV` (p × p) for either PCA or CVA. For PCA, standardisation is applied by default. For CVA, standardisation is bypassed (biplotEZ handles the within-class normalisation internally). The projected training coordinates `Z_train` are stored alongside `X_center`, `X_sd`, `standardise`, `proj_dims`, and the biplotEZ S3 object (`biplot_obj`) used for subsequent plotting. The `biplot_obj` field is the key link between the mathematical projection and the biplotEZ rendering pipeline.

**Design note (Q9):** CVA with only 2 classes is technically a degenerate case but is supported via biplotEZ. This is documented in the function help.

**`R/local_cf.R` — rotation (Phase 3 extension of Module B)**

The private function `.bl_rotate(x_target_st, V, proj_pair)` adapts the `Biplot_rotation()` logic from the research script `1.3 Optimal Rotation.R`. Given a standardised target point and an eigenvector pair, it computes an SVD-based orthogonal rotation matrix `A` such that the rotated loading matrix `Vrho = V %*% t(A)` aligns the target with the chosen projection plane. Returns `Vr_rot` (p × 2) and `tVr_rot` (2 × p) for the selected pair.

---

### 4.3 Module C — Grid Construction and Boundary Approximation
*Instruction manual §3C*

**`R/biplot_grid.R`**

| Function | Signature | Returns |
|---|---|---|
| `bl_build_grid()` | `(train_data, bl_projection, bl_model, m=200L, cutoff=0.5, b_margin=NULL, rounding=2L, polygon=NULL, outlie=1, calc_hull=TRUE)` | `bl_grid` |

Builds an `m × m` regular grid over the Z-space, back-projects each grid point to X-space using `tV`, scores it through the model, and extracts contour lines at `cutoff ± b_margin` using `grDevices::contourLines`. Stores:

- `ct` — full contour list (for plotting and global CF search)
- `ct_surrogate` — hull-clipped contours (for surrogate model only)
- `col_value` — per-grid-point colour vector (blue→white→red by probability)
- `polygon` — the Z-space convex hull as a `SpatialPolygons` object
- `hull_fraction` — the `outlie` value used

The hull polygon (`calc_hull = TRUE` by default) is visual-only and does not affect CF search. The `bl_surrogate()` function reads `ct_surrogate` directly; user code never needs to distinguish the two.

---

### 4.4 Phase 1 Assembly and Plotting

**`R/result.R`**

| Function | Signature | Returns |
|---|---|---|
| `bl_assemble()` | `(bl_data, bl_filter_result=NULL, bl_model, bl_projection, bl_grid)` | `bl_result` |

The central Phase 1 anchor object. Consolidates all metadata into a single reproducible object: training/test data, model, projection matrices (`V`, `tV`), centring/scaling (`X_center`, `X_sd`), grid and contours, polygon, `train_ranges`, `var_names`, `proj_dims`, `cutoff`, and `rounding`. The `biplot_obj` from `bl_projection` is carried forward for plotting.

**`R/plot_biplot.R`**

| Function | Signature | Returns |
|---|---|---|
| `plot_biplotEZ()` | `(bl_result, points=NULL, boundary=NULL, show_arrows=TRUE, arrow_col, target_point=NULL, target_label=NULL, no_grid, no_points, no_contour, new_title, cex_z, label_dir, tick_label_cex, ticks_v, which, X_names, label_offset_var, label_offset_dist, rotate_deg, contour_col, contour_lwd, contour_lty)` | `bl_result` (invisible) |

Eight-step rendering pipeline using biplotEZ:

1. **Base axes** — biplotEZ renders the axis skeleton (invisible samples, `opacity=0`)
2. **Prediction grid** — coloured squares per grid point
3. **Data points** — confusion-coloured (TP=red, TN=blue, FP=purple, FN=orange) or prediction-coloured for unlabelled data
4. **Axes redrawn on top** — so variable axes are visible over the grid
5. **Contour lines** — decision boundary approximation
6. **Target** — optional highlighted observation (yellow circle)
7. **Boundary overlay** — optional `bl_boundary` crosses and arrows
8. **Prediction summary** — printed to console (accuracy if labels known, class split if not)

The `rotate_deg` parameter patches the biplotEZ object's `Z`, `Lmat`, and `ax.one.unit` fields to rotate all rendered elements without rerunning the projection pipeline.

**`R/project_points.R`**

| Function | Signature | Returns |
|---|---|---|
| `bl_project_points()` | `(data, bl_result, filter_to_polygon=FALSE, filter_to_train_ranges=FALSE)` | `bl_points` |
| `bl_predict()` | `(bl_result, data=NULL)` | `data.frame` |

`bl_project_points()` projects any data frame (train, test, external) into Z-space and scores it through the model. Returns a `bl_points` object suitable for passing to `plot_biplotEZ()` via the `points` parameter. `bl_predict()` wraps `bl_project_points()` and returns a tidy data frame with columns `row`, `pred_prob`, `pred_class`, `true_class` (if known), `confusion` (TP/TN/FP/FN if known), and all feature columns — intended for inspection and target selection in interactive scripts.

---

### 4.5 Module D (Global) — Global Counterfactual Search
*Instruction manual §3D*

**`R/boundary.R`**

| Function | Signature | Returns |
|---|---|---|
| `bl_find_boundary()` | `(bl_result, data=NULL, tdp=NULL)` | `bl_boundary` |

Searches for the nearest decision boundary point for each observation in `data`. For each observation, iterates over all contour segments in `bl_grid$ct`, clips segments to the training hull polygon, applies `train_ranges` feasibility constraints, verifies model consistency at the boundary, then finds the nearest contour point using a block-wise nearest-neighbour search (`.nearest_idx_block()`). Stores `Z_obs`, `B_z`, `B_x`, `B_pred`, and `dist_z` per observation.

**Design note:** `train_ranges` are always applied in global search — this is not optional. The hull polygon constraint is also always applied. Both correspond to the "rows2" filter in the original research code.

---

### 4.6 Module E — Distance and Robustness Analysis
*Instruction manual §3E*

**`R/boundary_plot.R`**

| Function | Signature | Returns |
|---|---|---|
| `plot.bl_boundary()` | `(x, type=c("jitter","boxplot"), ...)` | list (invisible) |
| `bl_robustness()` | `(bl_boundary)` | list |

`plot.bl_boundary()` visualises the per-variable distance from each observation to its boundary counterfactual, using ggplot2 jitter or boxplot. Points are coloured by confusion type (TP/TN/FP/FN). `bl_robustness()` returns a per-variable distance summary and a scalar total robustness score.

---

### 4.7 Module F — Surrogate Model
*Instruction manual §3F*

**`R/surrogate.R`**

| Function | Signature | Returns |
|---|---|---|
| `bl_surrogate()` | `(bl_result, data=NULL)` | `bl_surrogate` |

Assigns each observation to the nearest biplot contour region (using `ct_surrogate`, the hull-clipped contours from `bl_grid`). Observations outside the training polygon receive `surrogate_pred = NA`. Reports accuracy of surrogate predictions vs the fitted model and vs true labels (when available). The plot renders the hull-clipped grid with surrogate-coloured points over the contour regions.

---

### 4.8 Module D (Local) — Local Counterfactual Search
*Instruction manual §3D, §3B (rotation)*

**`R/local_cf.R`**

| Function | Signature | Returns |
|---|---|---|
| `bl_select_target()` | `(bl_result, target, data=NULL)` | `bl_target` |
| `set_filters()` | `(bl_target, ...)` | `bl_filters` |
| `bl_find_boundary_local()` | `(bl_result, bl_target, set_filters=NULL, max_pairs=10L, m=200L, verbose=TRUE)` | `bl_local_result` |

**`bl_select_target()`** wraps a single observation (integer row index or single-row data frame) into a `bl_target` object, projecting it into Z-space and scoring it through the model. External data frames (no class column, not in train/test) are fully supported.

**`set_filters()`** defines per-variable actionability constraints for the local search:

| Constraint | Meaning during search | Meaning in sparse CF |
|---|---|---|
| `"decrease"` | CF value must be ≤ observed | Not special — Shapley decides |
| `"increase"` | CF value must be ≥ observed | Not special — Shapley decides |
| `"fixed"` | CF value within ±0.5 of observed | Always reverts to observed value |
| `c(min, max)` | CF value must lie in `[min, max]` | Not special — Shapley decides |

**`bl_find_boundary_local()`** iterates over up to `max_pairs` eigenvector pairs in canonical order — (1,2), (1,3), (2,3), (1,4), … — testing a maximum of 10 combinations (configurable). For each pair:

1. Rotates the biplot via `.bl_rotate()` to align the target with the projection plane
2. Projects the target and all training observations into the rotated Z-space
3. Builds an `m × m` grid, scores through the model (chunked for memory)
4. Extracts contour lines at `cutoff ± b_margin`
5. Back-projects contour segments to X-space
6. Applies `train_ranges` (feasibility) and `set_filters` (actionability)
7. Verifies model consistency at retained contour points
8. Finds the nearest valid boundary point

The pair producing the globally minimum Z-space distance is selected. The result stores the full grid, contours, rotation matrices (`Vr_rot`, `tVr_rot`), and the `set_filters` object (needed by `bl_sparse_cf()`).

**Important constraint difference from global search:**

| Constraint | Global (`bl_find_boundary`) | Local (`bl_find_boundary_local`) |
|---|---|---|
| Z-space hull polygon | Applied | **Not applied** — rotation invalidates it |
| `train_ranges` | Applied | Applied |
| Actionability filters | Not available | Applied via `set_filters` |

**`plot.bl_local_result()`** renders using the same biplotEZ pipeline as `plot_biplotEZ()`. The biplotEZ object from `bl_result` is patched in-place with the rotated local coordinates before rendering:

```r
biplot_plot$Z[, proj_dims]    <- Z_train_rot   # rotated training coordinates
biplot_plot$Lmat[, proj_dims] <- Vr_rot         # rotated loading matrix
biplot_plot$ax.one.unit       <- (1 / diag(t(tVr_rot) %*% tVr_rot)) * t(tVr_rot)
```

This formula for `ax.one.unit` is identical to the `rotate_deg` path in `plot_biplotEZ()`. The target point is coloured by the confusion scheme (TP/TN/FP/FN) when its true class can be determined from `bl_result$test_data`; otherwise red (predicted 1) or blue (predicted 0). Parameters match `plot_biplotEZ()` exactly: `no_grid`, `no_points`, `no_contour`, `label_dir`, `which`, `X_names`, `label_offset_var/dist`, `contour_col/lwd/lty`, `show_arrows`, `arrow_col`.

---

### 4.9 Module G — Shapley Values and Sparse Counterfactual
*Instruction manual §3G*

**`R/shapley.R`**

| Function | Signature | Returns |
|---|---|---|
| `bl_shapley()` | `(bl_local_result, exact_max_vars=14L, approx_perm=2048L, seed=1L)` | `bl_shapley` |
| `bl_sparse_cf()` | `(bl_shapley_result, round_to=NULL)` | `bl_sparse_result` |

**`bl_shapley()`** computes Shapley values for the path from the observed point to the boundary counterfactual. Two computation paths:

- **Exact** (`p ≤ exact_max_vars`, default 14): evaluates all 2^p subsets using the weighted marginal contribution formula, adapted from research script `1.2 Run_Shapley_optimized.R`.
- **Permutation** (`p > exact_max_vars`): samples `approx_perm` random permutations and averages marginal contributions.

Variables are classified as `"Supports"` or `"Contradicts"` using the following logic from research script `4.4 Local Interpretation Shapley.R`:

```r
class 0 observation:  positive Shapley → Supports,  negative → Contradicts
class 1 observation:  negative Shapley → Supports,  positive → Contradicts
```

`print.bl_shapley()` shows the per-variable Shapley table. `plot.bl_shapley()` renders a ggplot2 horizontal bar chart (darkblue = Supports, grey = Contradicts) with y-axis labels `"variable: observed → change_to_boundary"`.

**`bl_sparse_cf()`** builds a sparse version of the counterfactual:

1. Retain CF value only for `"Supports"` variables (optionally rounded by resolution `round_to`)
2. Revert `"Contradicts"` and `"Unknown"` variables to observed value
3. Revert `"fixed"`-constrained variables to observed value regardless of Shapley class
4. Score the sparse CF through the model
5. Set `solution_valid = TRUE` if the predicted class flips

The `round_to` parameter is **resolution-based** (round to nearest multiple), not decimal places:
- `round_to = 0.5` → rounds to nearest 0.5
- `round_to = 1` → rounds to nearest integer
- `round_to = NULL` (default) → no rounding

`print.bl_sparse_result()` shows three predictions side-by-side — Observed, Full CF, Sparse CF — to allow verification of the classification change at each stage. `plot.bl_sparse_result()` renders the full local biplot with both CFs overlaid: grey cross = full CF, green = sparse CF valid, yellow = sparse CF invalid.

---

## 5. Internal Utility Modules

### `R/utils.R`
Input validation functions used by all exported functions: `stop_if_not_data_frame()`, `stop_if_col_missing()`, `stop_if_not_in_range()`, `stop_if_not_positive_integer()`, `stop_if_not_scalar_numeric()`. Also `calc_gini()` (2 × AUC − 1) used for per-iteration accuracy reporting in the outlier filter.

### `R/hull_utils.R`
`.build_hull_polygon()` — constructs a convex hull `SpatialPolygons` object from the Z-space training coordinates, scaled by `outlie` to control tightness. Uses `aplpack::plothulls()` for the hull computation and `sp::SpatialPolygons` for the polygon object. `.poly_clip()` clips a matrix of Z-space points to the interior of the polygon. `.points_in_polygon()` returns a logical membership vector.

### `R/feasibility_utils.R`
`get_variable_ranges()` extracts per-variable min/max bounds from training data as `rlang` filter expressions. `get_filter_logical_vector()` evaluates these expressions row-wise against any data frame — used in both global and local CF search to enforce `train_ranges`.

### `R/model_utils.R`
`.fit_model()` dispatches to the appropriate fitting backend by `model_type`: GLM/GAM/SVM/LDA/NNET/rpart via parsnip workflows; GBM via direct `gbm::gbm()` call; XGBoost via `xgboost`. Returns a fitted model object wrapped in a tidymodels workflow where applicable.

### `R/predict_utils.R`
`.pred_function()` is the single unified prediction interface used by all boundary, surrogate, and Shapley functions. Dispatches to: tidymodels `predict()`, `gbm::predict.gbm()`, `e1071` SVM probability extraction, or a user-supplied `predict_fn`. Always returns a numeric vector in [0, 1], optionally rounded to `rounding` decimal places.

---

## 6. Design Decisions and Manual Mapping

### 6.1 Two-phase workflow (Manual §4, Architecture Q2–Q3)

The manual specifies a two-phase workflow. This is implemented as:

- **Phase 1** (`bl_prepare_data` → `bl_filter_outliers` → `bl_fit_model/bl_wrap_model` → `bl_build_projection` → `bl_build_grid` → `bl_assemble`): produces `bl_result`, the anchor object for all downstream analysis. All Phase 1 choices (model, projection, grid resolution, polygon) are fixed at assembly.
- **Phase 2 global** (`bl_find_boundary`, `bl_robustness`, `bl_surrogate`): operates on datasets, produces per-observation summaries.
- **Phase 2 local / Phase 3** (`bl_select_target`, `set_filters`, `bl_find_boundary_local`, `bl_shapley`, `bl_sparse_cf`): operates on single observations.

### 6.2 biplotEZ dependency (Architecture Q4)

biplotEZ is kept for plotting only. All mathematical operations — projection, grid scoring, contour extraction, inverse projection, Shapley computation — are performed without biplotEZ. The plotting pipeline uses biplotEZ solely to draw the axis skeleton (`samples(opacity=0)` + `axes()`), then overlays grid, points, contours, and markers manually using `graphics::points()` and `graphics::lines()`. This decoupling means the computation pipeline is not affected by biplotEZ API changes.

### 6.3 CVA with two classes (Architecture Q9)

CVA is technically defined for three or more classes, but biplotEZ permits two-class CVA. The package supports this; the limitation is documented in `bl_build_projection()` help. The Pima diabetes example (two classes) uses CVA throughout.

### 6.4 Cutoff parameter (Architecture Q10)

The `cutoff` parameter (default 0.5) is stored in `bl_result` and propagated through all downstream functions. Only `cutoff = 0.5` produces valid results under the current methodology. The parameter is retained to avoid breaking the interface in future extensions, but non-default values are clearly warned against in the documentation.

### 6.5 Hull polygon sourcing

The hull polygon in `bl_result` always comes from `bl_grid`, not from `bl_filter_result`. Reason: the filter polygon is constructed in the PCA space of the raw training data, while the grid polygon is constructed in the final Z-space (which may be CVA, or PCA with different `proj_dims`). Only the grid polygon is consistent with the rendered biplot. `bl_assemble()` enforces this by reading the polygon from `bl_grid`.

### 6.6 Local rotation invalidates the hull polygon

The per-pair SVD rotation in `bl_find_boundary_local()` produces a new Z-space for each eigenvector pair. The convex hull polygon from `bl_grid` is expressed in the original (unrotated) Z-space and is no longer geometrically meaningful after rotation. Therefore:

- Global search: hull polygon applied
- Local search: hull polygon **not applied**; only `train_ranges` (X-space bounds) and `set_filters` (actionability) constrain the search

This is the most important constraint difference between the two search modes.

### 6.7 Actionability filters (Manual §3D, Architecture Q8)

`set_filters()` exposes per-variable constraints that were deferred during initial architecture planning (Q8: "develop later"). These are Phase 3 only — they are not available in `bl_find_boundary()`. The `"fixed"` constraint allows a ±0.5 search window during boundary search (not strict equality) to avoid filtering out all contour points due to floating-point imprecision in the back-projection. In sparse CF generation, fixed variables always revert to the observed value regardless of their Shapley classification, because the intent of `"fixed"` is that the variable does not change in the final recommendation.

### 6.8 Shapley computation (Manual §3G)

The Shapley implementation adapts `1.2 Run_Shapley_optimized.R` (exact) and `4.4 Local Interpretation Shapley.R` (Supports/Contradicts classification). The exact algorithm is used when `p ≤ 14` (the default threshold), producing provably correct attribution. Beyond that, the permutation approximation is used with `M = 2048` samples by default. The Supports/Contradicts classification logic is class-direction-aware: for class-0 observations, positive Shapley values push toward class 1 (support crossing the boundary); for class-1 observations, negative values push toward class 0.

### 6.9 Custom model support (`bl_wrap_model`)

The instruction manual requires model-agnostic operation. The package achieves this via `bl_wrap_model()` + a user-supplied `predict_fn`. This was the mechanism used for the GAM example (which cannot be wrapped by parsnip without special handling): the user fits `mgcv::gam()` directly and passes `predict.gam()` as the predict function. The `predict_fn` is stored in `bl_model` and called uniformly via `.pred_function()` throughout the package.

### 6.10 Naming conventions

All exported functions use the `bl_` prefix. S3 method naming follows R conventions (`print.bl_result`, `plot.bl_boundary`). Private helpers use a `.` prefix (`.build_hull_polygon`, `.bl_rotate`). Internal utility functions (`stop_if_*`, `calc_gini`) have no prefix as they are not user-facing. This aligns with instruction manual §7.4 and §9.4.

---

## 7. Implementation Status

| Component | Status | Commit |
|---|---|---|
| Data preparation (`bl_prepare_data`, `bl_wrap_data`) | Complete | `feature/phase1-final` |
| Outlier filter (`bl_filter_outliers`) | Complete | `feature/phase1-final` |
| Model fitting (`bl_fit_model`, `bl_wrap_model`) | Complete | `feature/phase1-final` |
| Projection (`bl_build_projection`) | Complete | `feature/phase1-final` |
| Grid and contours (`bl_build_grid`) | Complete | `feature/phase1-final` |
| Assembly (`bl_assemble`) | Complete | `feature/phase1-final` |
| Biplot plotting (`plot_biplotEZ`) | Complete | `feature/phase1-final` |
| Point projection (`bl_project_points`, `bl_predict`) | Complete | `feature/phase3-local` |
| Global CF search (`bl_find_boundary`) | Complete | `feature/phase2-global` |
| Robustness (`bl_robustness`, `plot.bl_boundary`) | Complete | `feature/phase2-global` |
| Surrogate model (`bl_surrogate`) | Complete | `feature/phase2-global` |
| Local target selection (`bl_select_target`) | Complete | `feature/phase3-local` |
| Actionability filters (`set_filters`) | Complete | `feature/phase3-local` |
| Local CF search (`bl_find_boundary_local`) | Complete | `feature/phase3-local` |
| Shapley values (`bl_shapley`) | Complete | `feature/phase3-local` |
| Sparse CF (`bl_sparse_cf`) | Complete | `feature/phase3-local` |
| **Unit tests** | **Pending** | — |
| **Vignettes** | **Pending** | — |

---

## 8. Dependencies

| Package | Role | Used in |
|---|---|---|
| `biplotEZ` | Biplot rendering pipeline | `projection.R`, `plot_biplot.R`, `local_cf.R` |
| `dplyr` | `case_when()` for colour/classification logic | `project_points.R`, `boundary_plot.R`, `local_cf.R`, `shapley.R` |
| `ggplot2` | Shapley bar chart, boundary plots | `shapley.R`, `boundary_plot.R` |
| `sp` | SpatialPolygons hull, point-in-polygon | `hull_utils.R`, `project_points.R` |
| `aplpack` | Convex hull computation (`plothulls`) | `hull_utils.R` |
| `mgcv` | GAM support (direct fit path) | `model_utils.R` |
| `e1071` | SVM support | `model_utils.R`, `predict_utils.R` |
| `gbm` | GBM support | `model_utils.R` |
| `nnet` | Neural net support | `model_utils.R` |
| `rpart` | Decision tree support | `model_utils.R` |
| `xgboost` | XGBoost support | `model_utils.R` |
| `parsnip`/`workflows`/`discrim` | tidymodels workflow wrapping | `model_utils.R`, `predict_utils.R` |
| `kernlab` | SVM via parsnip | `model_utils.R` |
| `MASS` | LDA (`lda()`) | `model_utils.R` |
| `rlang` | Tidy evaluation for range filters | `feasibility_utils.R` |
| `grDevices` | `contourLines()`, colour ramps | `biplot_grid.R`, `local_cf.R` |

All dependencies are listed in `DESCRIPTION` under `Imports`. `testthat`, `knitr`, and `rmarkdown` are in `Suggests`.
