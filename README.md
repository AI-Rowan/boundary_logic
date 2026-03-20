# boundarylogic

**Biplot-based decision boundary interpretability for binary classification models**

`boundarylogic` is an R package that provides post hoc explanations for binary
classification models. It projects high-dimensional data into a 2D PCA or CVA
biplot space, approximates the model's decision boundary in that reduced space,
and identifies counterfactual observations — the nearest points that would flip
the predicted class.

## Key ideas

- The decision boundary of any classifier (GLM, GAM, GBM, SVM, XGB, …) is
  approximated in a 2D biplot plane, making it visually interpretable regardless
  of the number of features.
- Counterfactuals are identified as the nearest boundary points in biplot space,
  then inverse-projected back to the original feature space.
- **Local** interpretations (per-observation: what would need to change?) and
  **global** interpretations (variable importance by distance to boundary,
  surrogate models) are both supported.
- Shapley values decompose the path from an observation to its counterfactual
  into per-variable contributions.
- Actionability constraints let you restrict counterfactual searches to
  clinically or practically plausible directions.
- The 2D boundary is an approximation — the package is designed to communicate
  this clearly while still providing actionable insight.

---

## Installation

```r
# Install from GitHub
devtools::install_github("AI-Rowan/boundary_logic")
```

For development use (recommended while the package is actively evolving):

```r
# Clone the repository, open boundary_logic.Rproj in RStudio, then:
devtools::load_all()   # sources all R/ files instantly, no install needed
```

---

## Three-phase workflow

The package is organised around three phases. All Phase 1 outputs are stored in
a single `bl_result` object that acts as the reproducible anchor for Phases 2
and 3.

```
Phase 1 ─ Build the pipeline
  bl_prepare_data()  →  bl_filter_outliers()  →  bl_fit_model() / bl_wrap_model()
  →  bl_build_result()  →  plot_biplotEZ()

Phase 2 ─ Global interpretation (population level)
  bl_find_boundary()  →  plot(bl_bnd)  →  bl_surrogate()

Phase 3 ─ Local interpretation (single observation)
  bl_predict()  →  bl_select_target()  →  set_filters()  →  bl_find_local_cf()
  →  bl_shapley()  →  bl_find_sparse_cf()
```

---

## Quick start (iris dataset)

```r
library(boundarylogic)

# Phase 1: prepare data, fit model, build biplot
bl_dat <- bl_prepare_data(
  data           = iris,
  class_col      = "Species",
  target_class   = "versicolor",   # 1 = versicolor, 0 = all others
  train_fraction = 0.8,
  seed           = 42L
)

bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 0.91)

bl_mod <- bl_fit_model(
  train_data = bl_filt$train_data,
  var_names  = bl_filt$var_names,
  model_type = "GLM"
)

bl_results <- bl_build_result(
  bl_data  = bl_filt,
  bl_model = bl_mod,
  method   = "PCA",
  rounding = 3L
)

plot_biplotEZ(bl_results)

# Overlay test data
test_pts <- bl_project_points(bl_results$test_data, bl_results)
plot_biplotEZ(bl_results, points = test_pts)
```

```r
# Phase 2: global counterfactual search
bl_bnd <- bl_find_boundary(bl_results)
plot(bl_bnd)               # jitter plot — per-variable distances + robustness
plot(bl_bnd, type = "box") # boxplot by confusion category

bl_sur <- bl_surrogate(bl_results)
plot(bl_sur)
```

```r
# Phase 3: local interpretation of one observation
pred_summary <- bl_predict(bl_results)
tdp <- which(pred_summary$pred_class == 1)[1]   # first predicted class-1 row
tgt <- bl_select_target(bl_results, target = tdp)

bl_local <- bl_find_local_cf(bl_results, tgt)
plot(bl_local)

bl_shapley_values <- bl_shapley(bl_local)
plot(bl_shapley_values)

bl_sparse <- bl_find_sparse_cf(bl_shapley_values, round_to = NULL)
print(bl_sparse)
plot(bl_sparse)
```

---

## Bringing your own model

Use `bl_wrap_model()` when you have fitted your own model externally — for
example, using custom hyperparameters or a model type not natively supported.

```r
# Example: custom GAM with spline terms (Pima diabetes dataset)
library(mgcv)

custom_gam <- mgcv::gam(
  class ~ s(Glucose) + s(BMI) + s(BloodPressure) + Insulin + Age,
  data   = bl_filt$train_data,
  family = binomial(link = "logit")
)

bl_mod <- bl_wrap_model(
  model      = custom_gam,
  model_type = "custom",
  var_names  = bl_filt$var_names,
  predict_fn = function(m, new_data) {
    as.numeric(mgcv::predict.gam(m, newdata = new_data, type = "response"))
  },
  train_data = bl_filt$train_data
)

bl_results <- bl_build_result(bl_filt, bl_mod, method = "CVA", rounding = 2L)
plot_biplotEZ(bl_results)
```

---

## Actionability constraints

`set_filters()` restricts the counterfactual search to clinically or practically
plausible directions:

```r
flt <- set_filters(
  tgt,
  Glucose     = "decrease",    # can only be reduced
  Age         = "fixed",       # not actionable
  Pregnancies = "fixed"        # not actionable
)

bl_local <- bl_find_local_cf(bl_results, tgt, set_filters = flt)
```

| Constraint | Meaning |
|---|---|
| `"decrease"` | Counterfactual value must be ≤ observed |
| `"increase"` | Counterfactual value must be ≥ observed |
| `"fixed"` | Search allows ±0.5; always reverts to observed in sparse CF |
| `c(min, max)` | Counterfactual must lie within this absolute range |

---

## Supported model types

| `model_type` | Algorithm | Backend |
|---|---|---|
| `"GLM"` | Logistic regression | parsnip → `glm` |
| `"GAM"` | Generalised additive model | parsnip → `mgcv` |
| `"GBM"` | Gradient boosting | `gbm` (direct) |
| `"LDA"` | Linear discriminant analysis | parsnip → `MASS` |
| `"SVM"` | Support vector machine | parsnip → `kernlab` |
| `"NNET"` | Neural network | parsnip → `nnet` |
| `"RForrest"` | Decision tree | parsnip → `rpart` |
| `"XGB"` | XGBoost | parsnip → `xgboost` |
| `"custom"` | Any model with a `predict_fn` | `bl_wrap_model()` |

---

## Projection methods

| Method | When to use |
|---|---|
| `"PCA"` | General purpose; separates by overall variance |
| `"CVA"` | Maximises between-class separation; recommended when class structure matters |

---

## Function reference

### Phase 1 — Build the pipeline

| Function | Purpose | Returns |
|---|---|---|
| `bl_prepare_data()` | Convert raw data to binary 0/1, train/test split | `bl_data` |
| `bl_wrap_data()` | Register pre-split, pre-coded data | `bl_data` |
| `bl_filter_outliers()` | Remove outliers via convex hull filter | `bl_filter_result` |
| `bl_fit_model()` | Fit a supported classifier | `bl_model` |
| `bl_wrap_model()` | Register an externally fitted model | `bl_model` |
| `bl_build_result()` | Steps 4–6 combined (preferred) | `bl_result` |
| `bl_build_projection()` | Compute PCA or CVA loading matrix | `bl_projection` |
| `bl_build_grid()` | Build prediction grid, extract boundary contours | `bl_grid` |
| `bl_assemble()` | Collect all Phase 1 artifacts | `bl_result` |
| `plot_biplotEZ()` | Plot biplot with grid, points, and boundary | (plot) |

### Phase 2 — Global interpretation

| Function | Purpose | Returns |
|---|---|---|
| `bl_project_points()` | Project any data frame into Z-space | `bl_points` |
| `bl_predict()` | Score all observations; tidy prediction data frame | `data.frame` |
| `bl_find_boundary()` | Nearest boundary point per observation | `bl_boundary` |
| `bl_robustness()` | Per-variable distance totals and robustness scalar | named vector |
| `bl_surrogate()` | Region-based spatial surrogate model | `bl_surrogate` |

### Phase 3 — Local interpretation

| Function | Purpose | Returns |
|---|---|---|
| `bl_select_target()` | Wrap one observation as a target | `bl_target` |
| `set_filters()` | Encode actionability constraints | `bl_filters` |
| `bl_find_local_cf()` | Nearest boundary via biplot rotation | `bl_local_result` |
| `bl_shapley()` | Shapley attribution of the CF path | `bl_shapley` |
| `bl_find_sparse_cf()` | Parsimonious counterfactual (Contradicts reverted) | `bl_sparse_result` |

---

## Vignettes

Two vignettes cover the complete workflow:

- **Iris dataset** (`Boundary_Logic-workflow.Rmd`) — introductory walkthrough
  using the built-in iris dataset; covers Phases 1–3 step by step.
- **Pima diabetes dataset** (`Boundary_Logic_Pima_diabetes_workflow.Rmd`) —
  clinical context with a custom GAM, CVA biplot, actionability constraints,
  and external patient analysis (Phases 1–3).

```r
vignette("Boundary_Logic-workflow",              package = "boundarylogic")
vignette("Boundary_Logic_Pima_diabetes_workflow", package = "boundarylogic")
```

---

## Important notes

- **`cutoff = 0.5` only:** non-default cutoff values are accepted but do not
  give valid boundary results in the current methodology.
- **The 2D boundary is an approximation:** `boundarylogic` is designed for
  exploratory interpretability, not exact boundary reconstruction. The 2D
  projection captures the dominant structure of the decision surface but cannot
  represent every dimension of a high-dimensional model.
- **`rounding` controls contour band width only:** the `rounding` parameter in
  `bl_build_result()` / `bl_build_grid()` determines how tightly the boundary
  contour is extracted (`rounding = 2L` → band of ±0.01; `rounding = 3L` →
  band of ±0.001). All model predictions are always stored to 3 decimal places
  regardless of this setting.
- **Local search does not apply the convex hull:** `bl_find_local_cf()` rotates
  the biplot plane, which invalidates the original hull polygon. The local
  search is constrained only by `train_ranges` and `set_filters`.

---

## Running the tests

```r
devtools::test()
```

---

## Repository structure

```
boundary_logic/
├── R/                   ← Package source
├── tests/testthat/      ← Unit tests
├── vignettes/           ← Iris and Pima diabetes workflows
├── inst/extdata/        ← Bundled datasets (pima_diabetes.csv)
├── scripts/             ← Interactive example scripts
└── docs/                ← Flow diagram, implementation summary
```

Branches:

| Branch | Purpose |
|---|---|
| `main` | Current stable code — use this for general access |
| `method_developments` | Active development branch |
| `original_PhD_code` | Preserved original PhD research scripts (read-only) |
