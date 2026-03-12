############################################################
# Phase 1 walkthrough — boundary_logic package
# Run each section interactively (Ctrl+Enter line by line)
#
# Models are fitted via tidymodels (parsnip + workflows) for all
# types except "GBM". New package dependencies required:
#   install.packages(c("parsnip", "workflows", "discrim", "kernlab"))
############################################################

# ---- Load the package (development workflow) ---------------------------
devtools::load_all()   # loads from R/ source directly; no install needed


# ---- Step 1: Prepare data ---------------------------------------------
bl_dat <- bl_prepare_data(
  data           = iris,
  class_col      = "Species",
  target_class   = "versicolor",
  train_fraction = 0.8,
  seed           = 121L
)

print(bl_dat)


# ---- Step 2 (optional): Filter outliers --------------------------------
# hull_fraction = 1 keeps all points (no filtering)
# hull_fraction = 0.9 trims the outermost 10%

##A This step includes filtering the data for 
# train_ranges <- get_variable_ranges(train_filtered[, var_names, drop = FALSE])
# test_rows    <- get_filter_logical_vector(test_data[, var_names, drop = FALSE],
#                                           train_ranges)

bl_filt <- bl_filter_outliers(bl_dat, hull_fraction = 1)

print(bl_filt)


# ---- Step 3: Fit model -------------------------------------------------
# All types except "GBM" are fitted via parsnip + workflows.
# Hyperparameter overrides use parsnip parameter names (except GBM):
#   NNET: hidden_units, penalty, epochs
#   XGB:  trees, tree_depth, learn_rate, sample_size
#   RForrest: min_n
bl_mod <- bl_fit_model(
  train_data = bl_filt$train_data,
  var_names  = bl_filt$var_names,
  model_type = "GLM"
)

print(bl_mod)

## Steps 1-3 are background setup steps.
## Any fitted model and prepared data can be substituted here.
## The interpretability workflow proper starts at Step 4.


# ---- Step 4: Build projection ------------------------------------------
# Choose ONE of the two options below (PCA or CVA).

# Option A — PCA (default; works for any data)
bl_proj <- bl_build_projection(
  train_data = bl_filt$train_data,
  var_names  = bl_filt$var_names,
  method     = "PCA"
)

# Option B — CVA (maximises class separation; requires bl_model)
# Uncomment to use CVA instead of PCA.
# Note: set calc_ct_in_hull = TRUE in Step 5 when using CVA.
bl_proj <- bl_build_projection(
  train_data = bl_filt$train_data,
  var_names  = bl_filt$var_names,
  method     = "CVA",
  proj_dims  = c(1L, 2L),
  bl_model   = bl_mod        # derives TP/TN/FP/FN labels automatically
)

print(bl_proj)


# ---- Step 5: Build prediction grid ------------------------------------

bl_grid <- bl_build_grid(
  train_data      = bl_filt$train_data,
  bl_projection   = bl_proj,
  bl_model        = bl_mod,
  m               = 200L,
  calc_ct_in_hull = FALSE   # set TRUE if using CVA
)

print(bl_grid)


# ---- Step 6: Assemble result object -----------------------------------
result <- bl_assemble(
  bl_data          = bl_dat,
  bl_filter_result = bl_filt,
  bl_model         = bl_mod,
  bl_projection    = bl_proj,
  bl_grid          = bl_grid
)

print(result)
summary(result)


# ---- Step 7: Plot biplot -----------------------------------------------
plot_biplotEZ(
  result,
  label_offset_var  = c(1L, 4L),
  label_offset_dist = c(0.6,1.2)
  #rotate_deg = 90
)
