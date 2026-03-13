############################################################
# Example: biplot rotation — piped workflow from Step 3
#
# Steps 1-2 (data prep and model fit) produce bl_dat and
# bl_mod which are needed by multiple downstream functions,
# so they are kept as named objects.
#
# From Step 3 onwards (projection → grid → assemble → plot)
# results are piped using R native lambdas (\(x) ...) to
# carry intermediate outputs forward without naming them.
#
# Dataset: Iris (built-in), binary versicolor vs other.
############################################################

# ---- Load the package (development workflow) ---------------------------
rm(list = ls())
devtools::load_all()

# ---- Step 1: Prepare data -------------------------------------------
bl_dat <- bl_prepare_data(
  datasets::iris,
  class_col    = "Species",
  target_class = "versicolor"
)

# ---- Step 2: Fit model ----------------------------------------------
bl_mod <- bl_fit_model(
  bl_dat$train_data,
  bl_dat$var_names,
  model_type = "GLM"
)

# ---- Steps 3-5 piped: PCA projection → grid → assemble -------------
# bl_build_projection is piped into a lambda that also builds the grid,
# returning both as a list so bl_assemble can receive proj and grid.
result <- bl_build_projection(
    bl_dat$train_data,
    bl_dat$var_names,
    method = "PCA"
  ) |>
  (\(proj) list(
    proj = proj,
    grid = bl_build_grid(bl_dat$train_data, proj, bl_mod, m = 150L)
  ))() |>
  (\(x) bl_assemble(
    bl_dat,
    bl_model      = bl_mod,
    bl_projection = x$proj,
    bl_grid       = x$grid
  ))()

# ---- Plot without rotation (baseline) --------------------------------
plot_biplotEZ(result, new_title = "No rotation (0 deg)")

# ---- Plot with rotation — pipe result directly into plot -------------
result |> plot_biplotEZ(rotate_deg = 45,  new_title = "Rotated 45 deg clockwise")
result |> plot_biplotEZ(rotate_deg = 90,  new_title = "Rotated 90 deg clockwise")
result |> plot_biplotEZ(rotate_deg = 270, new_title = "Rotated 270 deg clockwise")

# ---- Rotation with target point --------------------------------------
result |> plot_biplotEZ(
  rotate_deg   = 45,
  target_point = unlist(result$test_data[1L, result$var_names]),
  target_label = 1L,
  new_title    = "Rotated 45 deg with target point"
)

# ---- Steps 3-5 piped: CVA projection --------------------------------
result_cva <- bl_build_projection(
    bl_dat$train_data,
    bl_dat$var_names,
    method   = "CVA",
    bl_model = bl_mod
  ) |>
  (\(proj) list(
    proj = proj,
    grid = bl_build_grid(bl_dat$train_data, proj, bl_mod, m = 150L)
  ))() |>
  (\(x) bl_assemble(
    bl_dat,
    bl_model      = bl_mod,
    bl_projection = x$proj,
    bl_grid       = x$grid
  ))()

result_cva |> plot_biplotEZ(new_title = "CVA — no rotation")
result_cva |> plot_biplotEZ(rotate_deg = 270, new_title = "CVA — rotated 270 deg")
