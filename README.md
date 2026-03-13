# boundarylogic

**Biplot-based decision boundary interpretability for classification models**

`boundarylogic` is an R package that provides post hoc explanations for binary classification models by projecting high-dimensional data into a 2D PCA or CVA biplot space, approximating the model's decision boundary in that reduced space, and identifying counterfactual observations — the nearest points that would change the predicted class.

## Key ideas

- The decision boundary of any classifier (GLM, GAM, GBM, SVM, …) is approximated in a 2D biplot plane, making it visually interpretable regardless of the number of features.
- Counterfactuals are identified as the nearest boundary points in Z-space, then inverse-projected back to the original feature space.
- Both **local** interpretations (per-observation: what would need to change?) and **global** interpretations (variable importance by distance to boundary, surrogate models) are supported.
- The 2D boundary is an approximation — the package is designed to communicate this clearly while still providing actionable insight.

## Two-phase workflow

**Phase 1 — Build and store the biplot anchor:**
1. Prepare data (`bl_prepare_data`)
2. Optionally filter outliers (`bl_filter_outliers`)
3. Fit a classification model (`bl_fit_model`)
4. Build a PCA or CVA projection (`bl_build_projection`)
5. Build the prediction grid and decision boundary contours (`bl_build_grid`)
6. Assemble into a reproducible result object (`bl_assemble`)

**Phase 2 — Interpret (not yet implemented):**
- Local: counterfactual search, Shapley attribution, biplot rotation to target
- Global: variable importance by distance, surrogate model

## Installation

```r
# Development version
devtools::install_github("your-org/boundarylogic")
```

## Quick start

```r
library(boundarylogic)

# 1. Prepare data — multiclass Iris converted to binary (versicolor vs rest)
bl_dat <- bl_prepare_data(datasets::iris,
                           class_col    = "Species",
                           target_class = "versicolor")

# 2. Fit a logistic regression model
bl_mod <- bl_fit_model(bl_dat$train_data, bl_dat$var_names, model_type = "GLM")

# 3. Build PCA projection
bl_proj <- bl_build_projection(bl_dat$train_data, bl_dat$var_names, method = "PCA")

# 4. Build prediction grid (decision boundary contours)
bl_grid <- bl_build_grid(bl_dat$train_data, bl_proj, bl_mod, m = 150L)

# 5. Assemble Phase 1 result
result <- bl_assemble(bl_dat,
                      bl_model      = bl_mod,
                      bl_projection = bl_proj,
                      bl_grid       = bl_grid)

# 6. Plot the biplot with decision boundary
plot_biplotEZ(result)

# 7. Overlay test data
plot_biplotEZ(result, points = bl_project_points(result$test_data, result))
```

## Piped workflow (from Step 3)

```r
bl_dat <- bl_prepare_data(datasets::iris,
                           class_col    = "Species",
                           target_class = "versicolor")
bl_mod <- bl_fit_model(bl_dat$train_data, bl_dat$var_names, model_type = "GLM")

result <- bl_build_projection(bl_dat$train_data, bl_dat$var_names, method = "PCA") |>
  (\(proj) list(
    proj = proj,
    grid = bl_build_grid(bl_dat$train_data, proj, bl_mod, m = 150L)
  ))() |>
  (\(x) bl_assemble(bl_dat,
                    bl_model      = bl_mod,
                    bl_projection = x$proj,
                    bl_grid       = x$grid))()

result |> plot_biplotEZ(rotate_deg = 45)
```

## Supported model types

| Code | Model |
|------|-------|
| `"GLM"` | Logistic regression (`glm`) |
| `"GAM"` | Generalised additive model (`mgcv`) |
| `"GBM"` | Gradient boosting (`gbm`) |
| `"LDA"` | Linear discriminant analysis (`MASS`) |
| `"SVM"` | Support vector machine (`kernlab`) |
| `"NNET"` | Neural network (`nnet`) |
| `"RForrest"` | Decision tree (`rpart`) |
| `"XGB"` | XGBoost (`xgboost`) |
| `"custom"` | User-supplied `predict_fn` via `bl_wrap_model()` |

## Projection methods

| Method | When to use |
|--------|-------------|
| `"PCA"` | General purpose; separates by overall variance |
| `"CVA"` | Maximises between-class separation; requires a fitted model to derive TP/TN/FP/FN class labels (4 levels) |

## Phase 1 function reference

| Function | Purpose |
|----------|---------|
| `bl_prepare_data()` | Load data, convert multiclass to binary 0/1, train/test split |
| `bl_wrap_data()` | Wrap already-split data into a `bl_data` object |
| `bl_filter_outliers()` | Iterative convex hull filter to remove outliers from training data |
| `bl_fit_model()` | Fit a supported classification model |
| `bl_wrap_model()` | Register an externally fitted model |
| `bl_build_projection()` | Compute PCA or CVA projection matrix and biplotEZ object |
| `bl_build_grid()` | Build m×m prediction grid and extract decision boundary contours |
| `bl_assemble()` | Assemble all Phase 1 artifacts into a reproducible `bl_result` |
| `bl_project_points()` | Project any dataset into the biplot Z-space |
| `plot_biplotEZ()` | Plot the biplot with grid, points, and contours |
| `plot.bl_projection()` | Quick biplot render directly from a `bl_projection` object |

## Important notes

- **cutoff = 0.5 only**: non-default cutoff values are accepted but do not give valid boundary results in the current methodology.
- **The 2D boundary is an approximation**: the package is designed for exploratory interpretability, not exact boundary reconstruction.
- **CVA with 2 classes**: biplotEZ supports CVA with 2 class levels, but 4-level confusion labels (TP/TN/FP/FN) derived from a fitted model are recommended for better separation.
