############################################################
# Data inspection biplot — CVA, no model required
#
# Purpose: visualise training data in the CVA biplot plane
#   using the true class labels (red = class 1, blue = class 0)
#   BEFORE fitting any model.
#
# CVA note: CVA works best with >= 3 class levels. Here we
#   pass the binary class factor (2 levels) directly via
#   cva_classes. biplotEZ will produce a warning about this,
#   which is expected and safe to ignore for exploration.
#   The result may be less well-separated than a 4-label CVA
#   (TP/TN/FP/FN), but is useful for initial data inspection.
#
# Run interactively: Ctrl+Enter line by line
############################################################

# ---- Load the package (development workflow) ---------------------------
devtools::load_all()


# ---- Step 1: Prepare data ---------------------------------------------
bl_dat <- bl_prepare_data(
  data           = datasets::iris,
  class_col      = "Species",
  target_class   = "versicolor",
  train_fraction = 0.8,
  seed           = 121L
)

print(bl_dat)


# ---- Step 2: CVA projection — binary class labels ---------------------
# cva_classes is supplied directly from the training data class column.
# No model is needed. A warning about 2 levels is expected.
bl_proj <- bl_build_projection(
  train_data  = bl_dat$train_data,
  var_names   = bl_dat$var_names,
  method      = "PCA",
  proj_dims   = c(1L, 2L),
  cva_classes = as.factor(bl_dat$train_data$cldata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==ass),
  title       = "CVA biplot — class 1 (red) vs class 0 (blue)"
)

print(bl_proj)


# ---- Step 3: Plot — colour points by true class -----------------------
# Colours (red = class 1, blue = class 0) are derived automatically from
# the binary cva_classes stored inside bl_proj.
plot(bl_proj)

# Add a legend manually
legend("topright",
       legend = c("Class 1 (target)", "Class 0"),
       col    = c("red", "blue"),
       pch    = 16,
       bty    = "n")
