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
# Install from GitHub (build_vignettes = TRUE makes vignette() available)
devtools::install_github("AI-Rowan/boundary_logic", build_vignettes = TRUE)
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
  bl_prepare_data()  (split + outlier filter via hull_fraction)
  →  bl_fit_model() / bl_wrap_model()  →  bl_build_result()  →  plot_biplotEZ()

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
  seed           = 42L,
  hull_fraction  = 0.9             # drop the outermost 10% of training points
)

bl_mod <- bl_fit_model(
  train_data = bl_dat$train_data,
  var_names  = bl_dat$var_names,
  model_type = "GLM"
)

bl_results <- bl_build_result(
  bl_data  = bl_dat,
  bl_model = bl_mod,
  method   = "PCA",
  b_margin = 0.001                 # contour band: cutoff +/- 0.001 (the default)
)

plot_biplotEZ(bl_results)

# Overlay test data
test_pts <- bl_project_points(bl_results$test_data, bl_results)
plot_biplotEZ(bl_results, points = test_pts)
```

```r
# Phase 2: global counterfactual search
bl_bnd <- bl_find_boundary(bl_results)
plot(bl_bnd)                   # jitter plot — per-variable distances + robustness
plot(bl_bnd, type = "boxplot") # boxplot by confusion category

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

bl_sparse <- bl_find_sparse_cf(bl_shapley_values)
print(bl_sparse)
plot(bl_sparse)
```

---

## Bringing your own model

Use `bl_wrap_model()` when you have fitted your own model externally — for
example, using custom hyperparameters or a model type not natively supported.

```r
# Example: custom GAM with spline terms (Pima diabetes dataset, bundled)
library(mgcv)

pima_raw <- read.csv(system.file("extdata", "pima_diabetes.csv",
                                 package = "boundarylogic"))
# Zeros in these columns are missing values; keep complete rows only
pima <- pima_raw[with(pima_raw, Glucose != 0 & BloodPressure != 0 &
                        SkinThickness != 0 & Insulin != 0 & BMI != 0), ]

bl_dat_pima <- bl_prepare_data(
  data          = pima,
  class_col     = "Outcome",       # already 0/1
  feature_cols  = c("Pregnancies", "Glucose", "BloodPressure",
                    "Insulin", "BMI", "Age"),
  hull_fraction = 1                # no outlier trimming
)

custom_gam <- mgcv::gam(
  class ~ s(Glucose) + s(BMI) + s(BloodPressure) + Insulin + Age + Pregnancies,
  data   = bl_dat_pima$train_data,
  family = binomial(link = "logit")
)

bl_mod_pima <- bl_wrap_model(
  model      = custom_gam,
  model_type = "custom",           # or model_type = "GAM", with no predict_fn
  var_names  = bl_dat_pima$var_names,
  predict_fn = function(m, new_data) {
    as.numeric(mgcv::predict.gam(m, newdata = new_data, type = "response"))
  },
  train_data = bl_dat_pima$train_data
)

bl_results_pima <- bl_build_result(bl_dat_pima, bl_mod_pima,
                                   method = "CVA", b_margin = 0.01)
plot_biplotEZ(bl_results_pima)
```

---

## Actionability constraints

`set_filters()` restricts the counterfactual search to clinically or practically
plausible directions. Continuing the Pima example with a patient the model
predicts as diabetic:

```r
pred_pima <- bl_predict(bl_results_pima)
tgt_pima  <- bl_select_target(bl_results_pima,
                              target = which(pred_pima$pred_class == 1)[1])

flt <- set_filters(
  tgt_pima,
  Glucose     = "decrease",    # can only be reduced
  Age         = "fixed",       # not actionable
  Pregnancies = "fixed"        # not actionable
)

bl_local_pima <- bl_find_local_cf(bl_results_pima, tgt_pima, set_filters = flt)
```

Constraints are checked on the values the model sees (here raw clinical units).
If the features were standardised before fitting, give `c(min, max)` ranges on
that standardised scale, and read the `"fixed"` tolerance as ±0.5 of it. The
exception is `"integer"`: it rounds to whole numbers in raw units when a
scaling was recorded with `bl_set_scaling()`.

| Constraint | Meaning |
|---|---|
| `"decrease"` | Counterfactual value must be ≤ observed |
| `"increase"` | Counterfactual value must be ≥ observed |
| `"fixed"` | Search allows ±0.5, then restores the observed value; if the restored point does not cross the boundary or leaves the training ranges, no counterfactual is reported. Always at observed in the sparse CF |
| `c(min, max)` | Counterfactual must lie within this absolute range |
| `"integer"` | A whole number in raw units (a count or a score); combine as `c("increase", "integer")` or `list(c(min, max), "integer")`. The counterfactual is rounded after the search and re-checked |

---

## Supported model types

**Fitted by the package** with `bl_fit_model()` (parsnip / tidymodels workflows),
and also accepted by `bl_wrap_model()` as externally fitted models:

| `model_type` | Algorithm | Backend |
|---|---|---|
| `"GLM"` | Logistic regression | parsnip → `glm` |
| `"SVM"` | Support vector machine | parsnip → `kernlab` |
| `"NNET"` | Neural network | parsnip → `nnet` |
| `"RForest"` | Decision tree | parsnip → `rpart` |

**Externally fitted only**, registered with `bl_wrap_model()`:

| `model_type` | Algorithm | Model object to pass |
|---|---|---|
| `"GAM"` | Generalised additive model | an `mgcv::gam()` fit |
| `"GBM"` | Gradient boosting | a `gbm::gbm()` fit |
| `"LDA"` | Linear discriminant analysis | a `MASS::lda()` fit |
| `"XGB"` | XGBoost | `list(model = <booster from xgboost::xgb.train() or xgboost::xgboost()>, features = <feature names>)` |
| `"custom"` | Any model | the model plus a `predict_fn(model, new_data)` returning probabilities |

The older spelling `"RForrest"` is still accepted, with a deprecation warning.
An `xgboost::xgboost()` model must be fitted with a factor 0/1 response, e.g.
`y = factor(class)`; with a numeric response `xgboost()` fits a regression, and
`bl_wrap_model()` stops.

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
| `bl_prepare_data()` | Convert raw data to binary 0/1, train/test split, convex hull outlier filter | `bl_filter_result` |
| `bl_wrap_data()` | Register pre-split, pre-coded data (no filtering) | `bl_data` |
| `bl_filter_outliers()` | Convex hull outlier filter for `bl_wrap_data()` data | `bl_filter_result` |
| `bl_set_scaling()` | Record a pre-model standardisation so plots and prints show raw units (display only) | same class as input |
| `bl_fit_model()` | Fit a supported classifier | `bl_model` |
| `bl_wrap_model()` | Register an externally fitted model | `bl_model` |
| `bl_build_result()` | Projection + grid + assembly in one call (preferred) | `bl_result` |
| `bl_build_projection()` | Compute PCA or CVA loading matrix | `bl_projection` |
| `bl_build_grid()` | Build prediction grid, extract boundary contours | `bl_grid` |
| `bl_assemble()` | Collect all Phase 1 artifacts | `bl_result` |
| `plot_biplotEZ()` | Plot biplot with grid, points, and boundary contours | (plot) |

### Phase 2 — Global interpretation

| Function | Purpose | Returns |
|---|---|---|
| `bl_project_points()` | Project any data frame into Z-space | `bl_points` |
| `bl_predict()` | Score all observations; tidy prediction data frame | `data.frame` |
| `bl_find_boundary()` | Nearest boundary point per observation | `bl_boundary` |
| `bl_robustness()` | Per-variable distance totals and robustness scalar | list (`sum_of_distance`, `robustness`, `distance`) |
| `bl_surrogate()` | Region-based spatial surrogate model | `bl_surrogate` |
| `bl_pick_point()` | Interactively pick points on a drawn biplot; with `bl_boundary`, draws each pick's counterfactual arrow (pass the same `data` as the boundary search) | `data.frame` |

### Phase 3 — Local interpretation

| Function | Purpose | Returns |
|---|---|---|
| `bl_select_target()` | Wrap one observation as a target | `bl_target` |
| `set_filters()` | Encode actionability constraints | `bl_filters` |
| `bl_find_local_cf()` | Nearest boundary via biplot rotation | `bl_local_result` |
| `bl_shapley()` | Shapley attribution of the CF path | `bl_shapley` |
| `bl_find_sparse_cf()` | Parsimonious counterfactual (changes Supports features, largest Shapley first, until the boundary is crossed) | `bl_sparse_result` |

---

## Vignettes

The package vignette covers the complete workflow:

- **Loan-default dataset** (`Boundary_Logic_loan_default_workflow.Rmd`) —
  credit-risk walkthrough with an XGBoost model and CVA biplot: Phase 2
  feature pruning, actionability constraints, external-applicant analysis,
  and a side-by-side comparison of the global and local interpretation
  against SHAP (`fastshap` / `shapviz`) for the same model (Phases 1–3).

```r
vignette("Boundary_Logic_loan_default_workflow", package = "boundarylogic")
```

The vignette is available after installing with `build_vignettes = TRUE` (see
Installation), and on the package's pkgdown site.

---

## Important notes

- **`cutoff = 0.5` only:** non-default cutoff values are accepted but do not
  give valid boundary results in the current methodology.
- **The 2D boundary is an approximation:** `boundarylogic` is designed for
  exploratory interpretability, not exact boundary reconstruction. The 2D
  projection captures the dominant structure of the decision surface but cannot
  represent every dimension of a high-dimensional model.
- **`b_margin` controls contour band width only:** the `b_margin` parameter in
  `bl_build_result()` / `bl_build_grid()` is the half-width of the band around
  the cutoff in which the boundary contour is extracted (`b_margin = 0.001`,
  the default → contours at 0.499 and 0.501; `b_margin = 0.01` → 0.49 and
  0.51). All model predictions are always stored to 3 decimal places
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
├── vignettes/           ← Loan-default workflow
├── inst/extdata/        ← Bundled datasets (loan_data.csv, pima_diabetes.csv)
├── development/         ← Dev harness + code walkthroughs
├── documentation/       ← Technical notes (and the scripts behind them)
└── docs/                ← pkgdown site (generated)
```

The implementation summary is `2 implementation_summary.txt` in the repository root.

Branches:

| Branch | Purpose |
|---|---|
| `main` | Current stable code — use this for general access |
| `method_developments` | Active development branch |
| `paper_nr1_changes` | Package fixes for the first journal paper (in progress) |
| `original_PhD_code` | Preserved original PhD research scripts (read-only) |
