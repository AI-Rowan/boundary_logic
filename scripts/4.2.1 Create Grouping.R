############################################################
# Script: 4.2.1 Create Grouping 
# Description:
#   Build the biplot grid for surrogate grouping and prepare inputs
#   for assigning points to the nearest contour (ct). 
#
# Notes:
#   * Uses training-data V (Vmat/tVmat), not any rotated V.
#   * Key difference for grouping: calc_ct_in_hull = TRUE and no_polygon = TRUE.
############################################################

# --- Rounding for contour grid -------------------------------------------
rounding_ct <- rounding

# --- Build grouping grid (biplot_grouping) ----------------------------------
#   Based on training data biplot_data
biplot_grouping <- biplot_plane(
  X = biplot_data_surrogate,
  biplot_plot = biplot_matrix,
  model_use, var_names, numvars = num_vars, cutoff = cut_off,
  model_select = model_select,
  # Standardisation handling (FALSE automatically if CVA_method == TRUE)
  standardise = standardise_use,
  Xcenter = X_center, Xsd = X_sd,
  CVA_method = CVA_method,
  proj = pc.values,
  # Use starting V from training (never the rotated values)
  V = Vmat, tV = tVmat,
  # Ignore any pre-existing polygon; create a new one for grouping
  # no_polygon must be TRUE for grouping; calc_hull must be TRUE
  no_polygon = TRUE,
  calc_hull = TRUE,
  # NB: outer ct hull value used to group by ct value
  calc_ct_in_hull = TRUE,   # Set TRUE for biplot grouping, together with no_polygon = TRUE
  # outlie = 0.9,           # (optional) keep as-is, commented in original
  m = m_grid,
  b_margin = 1/(10^rounding_ct), rounding = rounding_ct,
  return_values = TRUE
)

# --- Prepare inputs for grouping by contours ------------------------------
{
  # 1) Data to classify based on the nearest ct-line
  #    (X.st not needed here; we take the 2D Z.st from biplot_grouping)
  Z.st <- biplot_grouping$Z.st

  # 2) Contours and polygon from training grid
  ct <- biplot_grouping$ct
  data_polygon <- biplot_grouping$polygon  # select CT within the polygon

  # 3) Get grouping based on contours (produces boundary_set_nr, boundary_used)
  source("4.2.2 Biplot Grouping CT.R")

  # --- Optional: derive surrogate class from contour levels ---------------
  #     boundary_class aligns each retained boundary to class 0/1 via cutoff,
  #     then indexes by boundary_set_nr for point-level labels.
  boundary_class <- c()
  for (i in 1:length(boundary_used)) boundary_class[i] <- ct[[boundary_used[i]]]$level
  boundary_class <- ifelse(boundary_class > cut_off, 1, 0)
  boundary_class <- boundary_class[boundary_set_nr]

  # (Optional) Alternative subclass grouping ideas kept from original:
  # original_class <- unique(train_data$Original_class)
  # nr_original_class <- length(original_class)
  # new_class <- original_class
  # for (i in 1:nr_original_class) {
  #   j <- train_data$Original_class == original_class[i]
  #   new_class[j] <- boundary_set_nr[j] + original_class[i] * nr_original_class
  # }
}

#################################################################
#### Plotting of grouping results based on boundary_set_nr   ####
#################################################################
# (Plotting is done in the surrogate script; this file only prepares inputs.)
