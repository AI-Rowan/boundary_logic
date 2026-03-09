############################################################
# Section 4.1: Global Interpretations – Robust & VIP
# Description:
#   Plot global variable importance-like diagnostics based on
#   distances to the decision boundary on the biplot plane.
#   Guards against undefined vars
#
# Notes:
#   * Assumes objects from the main driver exist (e.g., Vmat, tVmat,
#     biplot_grid, train_ranges, model_use, model_select, etc.).

############################################################

####################################
############  Plot change  #########
####################################

## Use unrotated biplots with no additional filters
## The biplot polygon should come from TRAIN to avoid implausible CFs.
## biplot_plane already provides a polygon; use no_polygon = FALSE and calc_hull = TRUE.

# Filters: keep structure but allow easy extension
set_filters <- list(TRUE)

# Main training biplot: no rotation, select grid to search
{

  ## biplot:
  pc.values <- c(1, 2)
  #biplot_grid_search <- biplot_grid # Set the biplot to search
  biplot_grid_search <- biplot_grid_Test # Set the biplot to search
  
  ## CFs:
  #biplot_CF_search <- biplot_search
  biplot_CF_search <- biplot_search_T
}

# Run Boundary Search on the existing (training) grid
# {
#   biplot_CF_search <- Biplot_boundary_search(
#     X = biplot_grid_search, # No require biplot_data, only the results of the biplot_grid
#     model_select = model_select,
#     standardise = standardise_use,   # FALSE automatically if CVA_method == TRUE
#     Xcenter = X_center,
#     Xsd = X_sd,
#     CVA_method = CVA_method,
#     proj = pc.values,
#     V = V,
#     tV = tV,
#     set_filters = set_filters,
#     train_ranges = train_ranges
#   )
# }

########################################
## VIP-style plot (distance to boundary)
########################################
tdp <- 1:nrow(biplot_dataT) #
#tdp <- sample(nrow(biplot_dataT),100000,replace = F)

{
  # Color semantics for points (predicted class)
  class_data <- biplot_grid_search$pred_use[tdp]
  # Alternative options:
  # class_data <- biplot_data[ , "pred_col"]
  # class_data <- biplot_grid$class


  # Distances are measured on the biplot plane, not in original X-space
  # Vector from each selected point to its nearest boundary point in Z-space
  Z_st_subset <- as.matrix(biplot_grid_search$Z.st[tdp,]) #data points
  Bz_subset   <- as.matrix(biplot_CF_search$B_boundary_z_value[tdp,]) # counterfactuals
  #vec_to_boundary_Z <- Bz_subset - Z_st_subset
  vec_to_boundary_Z <- Z_st_subset - Bz_subset # show which side of the boundary you are: move from B to point, so if negative then you are below

  # Map to X-space using the inverse loading rows for the plotted components
  # (2×p block of tV)
  tVr <- tV[pc.values,]
  vec_to_boundary <- vec_to_boundary_Z %*% tVr
  colnames(vec_to_boundary) <- var_names
  vec_to_boundary <- as.matrix(vec_to_boundary)

  # If PCA was standardised, convert back by multiplying by X_sd
  if (isTRUE(standardise_use)) {
    vec_to_boundary <- sweep(vec_to_boundary, 2, X_sd, "*")
  }

  # Standardise per-variable for comparability (like z-scores by sd)
  vec_to_boundary_sd <- sweep(vec_to_boundary, 2, X_sd, "/")

  # Combine with class labels for plotting
  boundary_vec_class <- data.frame(
    class = class_data,
    vec_to_boundary_sd,
    check.names = FALSE
  )

  # Aggregate absolute distances per variable across all selected points
  sum_of_distance <- colSums(abs(vec_to_boundary_sd), na.rm = TRUE)

  # Rename columns to include totals (rounded)
  distance_names <- paste0(colnames(boundary_vec_class)[-1],
                           " : ", round(sum_of_distance, 0))
  colnames(boundary_vec_class)[-1] <- distance_names

  # Order variables by total absolute distance
  ord <- order(sum_of_distance)
  boundary_vec_class <- boundary_vec_class[, c(1, ord + 1L)]

  # Convert class for ggplot colouring (keep 0/1 semantics)
  boundary_vec_class[["class"]] <- factor(boundary_vec_class[["class"]] + 1,
                                           labels = c("0", "1"))

  # Long format for ggplot
  dfp1 <- melt(boundary_vec_class,
                         id.vars = c("class"),
                         variable.name = "Variable")
  dfp1$value <- round(as.numeric(dfp1$value), 2)

  # Plot
  distance_plot <- ggplot(data = dfp1, aes(x = value)) +
    geom_jitter(aes(y = Variable, colour = class), size = 1.5) +
    labs(
      x = "Distance to Boundary – Standardised",
      y = "Variable : sum of distance to boundary",
      title = "Variable Distance to Boundary"
    ) +
    scale_color_manual(values = c("0" = "deepskyblue", "1" = "#F8766D")) +
    geom_vline(xintercept = 0) +
    labs(color = "Prediction Class") +
    theme_light() +
    theme(legend.position = "bottom")
}

print(sum_of_distance)
print(distance_plot)


################################################
#### Robustness (total distance) ################
################################################

# Total absolute distance across variables (rough global robustness)
sum(sum_of_distance, na.rm = TRUE)

# Compare to model-based VIP
library(vip)
vip(model_use)

####################################
###########  END  ##################
####################################
