############################################################
# Example: biplot rotation with plot_biplotEZ()
#
# Demonstrates how to use the rotate_deg parameter in
# plot_biplotEZ() to rotate the entire biplot — grid,
# data points, contour lines, and variable axes — clockwise
# by a specified angle.
#
# Dataset: Iris (built-in), binary versicolor vs other.
############################################################

# ---- Load the package (development workflow) ---------------------------
rm()
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

# ---- Step 3: Build projection (PCA) ---------------------------------
bl_proj <- bl_build_projection(
  bl_dat$train_data,
  bl_dat$var_names,
  method = "PCA"
)

# ---- Step 4: Build prediction grid ----------------------------------
bl_grid <- bl_build_grid(
  bl_dat$train_data,
  bl_proj,
  bl_mod,
  m = 150L
)

# ---- Step 5: Assemble -----------------------------------------------
result <- bl_assemble(
  bl_dat,
  bl_model      = bl_mod,
  bl_projection = bl_proj,
  bl_grid       = bl_grid
)

# ---- Step 6: Plot without rotation (baseline) -----------------------
plot_biplotEZ(result, new_title = "No rotation (0 deg)")

# ---- Step 7: Plot with 45-degree clockwise rotation -----------------
# rotate_deg rotates the entire plot clockwise:
#   - Variable axis arrows and labels move with the rotation
#   - Data points, grid, and contour lines all rotate consistently
#   - The underlying model and projection are unchanged

plot_biplotEZ(result, rotate_deg = 45, new_title = "Rotated 45 deg clockwise")

# ---- Step 8: Plot with 90-degree rotation ---------------------------
plot_biplotEZ(result, rotate_deg = 90, new_title = "Rotated 90 deg clockwise")

# ---- Step 9: Rotation with a target point ---------------------------
# The target point is also rotated to the correct position.
target <- unlist(result$test_data[1L, result$var_names])

plot_biplotEZ(
  result,
  rotate_deg   = 45,
  target_point = target,
  target_label = 1L,
  new_title    = "Rotated 45 deg with target point"
)

# ---- Step 10: Rotation with CVA biplot ------------------------------
bl_proj_cva <- bl_build_projection(
  bl_dat$train_data,
  bl_dat$var_names,
  method   = "CVA",
  bl_model = bl_mod
)

bl_grid_cva <- bl_build_grid(
  bl_dat$train_data,
  bl_proj_cva,
  bl_mod,
  m = 150L
)

result_cva <- bl_assemble(
  bl_dat,
  bl_model      = bl_mod,
  bl_projection = bl_proj_cva,
  bl_grid       = bl_grid_cva
)

plot_biplotEZ(result_cva, new_title = "CVA — no rotation")
plot_biplotEZ(result_cva, rotate_deg = 270, new_title = "CVA — rotated 30 deg")

