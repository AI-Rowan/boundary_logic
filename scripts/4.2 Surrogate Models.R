############################################################
# Section 4.2: Surrogate models 
# Description:
#   Surrogate classification via biplot contour grouping.
############################################################

##################################
####### option 1  ################
##################################

# Select nearest boundary set number, then determine that contour value's class
# (Assumes objects from prior steps exist in the environment.)
# - Expects that 4.2.1 Create Grouping has constructed biplot_grouping with fields:
#   Z.st, ct, polygon, pred_use, class, etc.


# !!!! NB !!!! set the biplot_grid to search
biplot_grid_search <- biplot_grid_Test # Set the biplot to search
biplot_data_surrogate <- biplot_dataT # Set the data points to search
biplot_matrix <- biplot_input$biplot_out # Training matrix

source("4.2.1 Create Grouping.R")

{
  # Contour value and class/colour based on cutoff
  ct_value <- z_boundary_type[boundary_set_nr]
  B_boundary_class <- ifelse(ct_value < cut_off, 0, 1)       # classify by score
  B_boundary_col   <- ifelse(ct_value < cut_off, "blue", "red")  # colour by classification
  # pch_z <- as.matrix(if_else(biplot_data$class == 0, 16, 17))  # shape by actual class (optional)

  # Surrogate model accuracy:
  # Compare surrogate class vs. biplot_grouping predictions
  ConfusionMatrix(B_boundary_class, biplot_grouping$pred_use)
  Accuracy(B_boundary_class, biplot_grouping$pred_use)

  # biplot_grouping may filter rows; compare with its retained data/class as well
  ConfusionMatrix(B_boundary_class, biplot_grouping$class)
  Accuracy(B_boundary_class, biplot_grouping$class)
}


# Update the colours on the grouping grid based on surrogate class
biplot_grouping$pred_col <- B_boundary_col

# Plot (kept as in original; tweak aesthetic params upstream if needed)
plot_biplotEZ(
  grid = biplot_grouping,   
  biplot_plot = biplot_matrix,# biplot function sed
  # k = c(1, 4),
  # label_dist = c(0.5),
  
  cex_z = 0.65,
  no_points = FALSE,
  tick.label.cex = 0.9,
  label_dir = "Paral",
  ticks_v = 5,
  new_title = "CVA Biplot GBM Surrogate Model"

)

##################################
####### option 2  ################
##################################

# Allocate class based on prediction score of closest grid point:
# - Compute distance from each Z.st point to each grid point
# - Select gridscore of the closest grid point
# - If gridscore > cutoff -> class 1, else class 0
# (Implementation not changed here; see original notes.)

