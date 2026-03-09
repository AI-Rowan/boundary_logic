############################################################
# Main section: Biplot Search of local decision boundaries
# Description:
#   Main driver for loading trainig data, setting global flags, fitting/using a
#   selected model, set tthe g and V functions, constructing the biplot (PCA/CVA) and
#   generating the prediction grid and contours, projecting points

#
# Notes:
#   * Requires the helper scripts and RData files referenced via source()/load().
#   * Process is very iterative - repeat until varaibles, filtered data and MLM
#       is selected to generate the final biplot projection function used in XAI
############################################################

############################################################
# Step 0: Load functions and set model parameters
############################################################
rm(list = ls())

{
source("scripts/1 Functions for biplot.R")
  # replace with two new functions
source("scripts/1.1 biplot_plane_optimized.R")
source("scripts/1.2 Run_Shapley_optimized.R")
  
# load functions
source("scripts/1.3 Optimal Rotation.R")   # Provides Biplot_rotation
source("scripts/1.4 Biplot Boundary Search_optimized.R")  # Provides Biplot_boundary_search
source("scripts/1.5 Data point Picker.R") # select a point on the biplot - default 2 points

}

{
  # Global switches
  first_standardise <- FALSE         # If FALSE: modeling uses standardized data in helper scripts
  no_polygon <- FALSE                # If TRUE: do not limit to convex hull polygon
  train_polygon <- NA
  no_pca <- FALSE                    # If TRUE: return rotation = identity (original space)
  cut_off <- 0.5                     # Classification threshold

  # Biplot method and standardization
  CVA_method <- T
  standardise_use <- TRUE
  if (CVA_method == TRUE) standardise_use <- FALSE

  pc.values <- c(1, 2)
  rounding <- 2
  m_grid <- 200                      # Grid density (reduce for heavy models)

  # LDB filter defaults
  set_filters <- list(TRUE)

}

# Select the MLM type used
# to name the plot in graphs later automatically: Set default GLM
{
  glm_use <- F
 ##>> Remove this option: lda <- F
  lda <- F
  svm_u <- F
  gam_used <- F
  gbm_used <- F
  rpart_used <- F
  nnet_used <- F
  xgb_used <- T
  
}

# Whether to re-run with polygon filters later
Polygon_filter <- F

############################################################
## Step 1 and 2: Load data and Fit model f(x) (or load pre-fit model)

############################################################


# load / create the "data" data.frame - including the test data, saved as "model_data"
# Open "2 Data File.R" and select the correct data: create "data" data.frame
# Always need to pre-calc this part: 
## source("2 Data File.R")

# prepare the data file for use by MLM and Biplots
source("scripts/2.1 Data Setup.R")          # Updates data object with var_names, num_vars, etc.

# fit the MLM model , or load existing
{
  # include the prediction funtion 
  source("scripts/2.2 Model_use fitting v2.R") # run model, Requires "train_data" data.frame
  
  if(Polygon_filter == F) model_results <- model_fitting(train_data = train_data)  
  if(Polygon_filter == T) model_results <- model_fitting(train_data = biplot_data[,1:(num_vars+1)])  


  # Allocate results, or load exising model
  model_select <- model_results$model_select
  model_use <- model_results$model_use
  summary(model_use)
  
}



############################################################
# Step 3.1: Select train/test data and set feasibility filters

## Start visual investigation and data region selection, without changing the MLM
## different data will create different biplot: 
## >> so decide how to use the filters to remove outliers
############################################################

Polygon_filter

{
  # find the ranges
  if (Polygon_filter == FALSE) train_ranges <- get_variable_ranges(train_data[, 1:num_vars])
  if (Polygon_filter == TRUE) train_ranges <- get_variable_ranges(biplot_data[, 1:num_vars])

  # filter by ranges: Row-selection logical vector: a
  if (Polygon_filter == FALSE) rows_to_exclude <- get_filter_logical_vector(train_data[, 1:num_vars], train_ranges)
  if (Polygon_filter == TRUE) rows_to_exclude <- get_filter_logical_vector(biplot_data[, 1:num_vars], train_ranges)
}

# set start data
{
  if (Polygon_filter == FALSE) start_data <- train_data[rows_to_exclude, ]
  if (Polygon_filter == TRUE) start_data <- biplot_data[, c(1:(num_vars + 1))]
}

#inspect the data
str(start_data)

############################################################
# Step 3.2: Create data input for biplot: Determine point colours (classification results)

# Option 1: Re-run with polygon filters
#   ##> Can Remove data once the initial polygon data filter is available
# round 1
# Polygon_filter <- F
# train_polygon <- NA

############################################################


# round 2
 Polygon_filter <- F
## train_polygon

biplot_data <- set_model_data(
  model_use,
  data_use = start_data,
  num_vars = num_vars,
  cutoff = cut_off,
  rounding = 2,
  target_class = tdp_class,
  model_select = model_select,
  
  ##> REmove data once the initial polygon data filter is available
  polygon = NA,#train_polygon,
  proj = pc.values,
  CVA_method = CVA_method,
  standardise = standardise_use
)$data_use

#rm(start_data)

## set colours to target data values for PCA
# biplot_data$pred_col <- if_else(biplot_data$class == 1, "red", "blue")

############################################################
# Step 4: Get projection matrix V (PCA/CVA) and standardisation function g

## Target class determined by biplot_data[, "pred_col"]
## if CVA biplot, put prediction accuracy is 100%, then only two class available, instead of 3 needed for 2-d CVA biplot
############################################################

{
  biplot_input <- biplot_input_calc(
    X = biplot_data[, c(1:num_vars)],
    standardise = standardise_use,
    CVA_method = CVA_method,
    main_heading = paste0(ifelse(CVA_method == TRUE, "CVA", "PCA"),
                          " Biplot based ", model_select),
    CVA_classes = as.factor(biplot_data[, "pred_col"])   # Allow 4 classes for CVA biplot
  )

  V <- Vmat <- biplot_input$biplot_out$Lmat
  tV <- tVmat <- solve(Vmat)

  Wmat <- biplot_input$biplot_out$Wmat
  biplot_ax.one.unit <- biplot_input$biplot_out$ax.one.unit

  X_center <- biplot_input$X_center
  X_sd <- biplot_input$X_sd
  
}

############################################################
# Steps 5–9 (Section 1): Build initial biplot grid & contours
############################################################
rounding <- 2
outlie <- 1
outlie <- ifelse(Polygon_filter == TRUE, 1, 0.9)

biplot_grid <- biplot_plane(
  X = biplot_data,
  biplot_plot = biplot_input$biplot_out,
  model_use,
  var_names,
  numvars = num_vars,
  model_select = model_select,
  cutoff = cut_off,
  standardise = standardise_use,
  Xcenter = X_center,
  Xsd = X_sd,
  CVA_method = CVA_method,
  proj = pc.values,
  V = V,
  tV = tV,
  polygon =NA,# final_polygon,           # Could also use biplot_train$polygon
  no_polygon = FALSE,
  calc_hull = F,
  calc_ct_in_hull = F,
  outlie = outlie,
  m = 200,#m_grid,
  b_margin = 1 / (10^rounding),
  rounding = rounding
)

# Keep polygon for later searches if needed
if (Polygon_filter == FALSE) train_polygon <- biplot_grid$polygon
if (Polygon_filter == TRUE) final_polygon <- biplot_grid$polygon

############################################################
# Step 9 (Part 2): Plot grid & points

## Now visually investigate the plot to explain the model
## can decide to run again with the data filtered to create zoomed in biplot
## >> return to 3.1
## or select different variables and data
## >> return to "2.Data File", rerun the MLM and create new data, and rerun this file

############################################################
# var_names

plot_biplotEZ(
  grid = biplot_grid,
  biplot_plot = biplot_input$biplot_out,
  # Label offsets (tune as needed)
  
  #Loans PCA
   # k = c(1,2,3,4,5,6,7,8),
   # label_dist = c(0,0.5, 1.5,0,0, 0.5,1,1.5),
  
  #Loans initial CVA
  k = c(1,2,3,4,5,6,7,8),
  label_dist = c(0,0.5, 0.5,1.5,0, 0.5,2,1.5),
  
    #Loans updated CVA
  # k = c(3,5,1,4,8,9),
  # label_dist = c(0, 0.5,2,1.5,1,0.5),
  
  tick.label.cex = 0.75,

  # points and features
  cex_z = 0.7,
  no_grid = F,
  no_points = FALSE,
  no_contour = F,
  label_dir = "Paral",
  ticks_v = 3,
  new_title = ""#CVA biplot of training data"
)

############################################################
# Steps 10–11: Run CF search (unconstrained or constrained)

## Can plot individual data points
## or select data points for further manual investigation
## once CFs are selected, can proceed to 4 for model explainion steps

############################################################

biplot_search <- Biplot_boundary_search(
  X = biplot_grid,
  model_select = model_select,
  #tdp = NA,#,tdp,            # Use NA to search all points
  standardise = standardise_use,
  Xcenter = X_center,
  Xsd = X_sd,
  CVA_method = CVA_method,
  proj = pc.values,
  V = V,
  tV = tV,
  set_filters = set_filters,
  train_ranges = NA#,train_ranges
)

# Plot TDP or CF 
ztdp <- point_inter_proj(
  #X = biplot_data[, 1:num_vars],
  X = biplot_search$B_boundary_x_value,
  V = V,
  #tdp = 1,
  plot = TRUE,
  col_x = "green",
  #col_x = biplot_data$pred_col,
  cex_x = 1.00,
  standardise = standardise_use,
  CVA_method = CVA_method
)

# Click to pick points (press ESC to stop early)
picked <- pick_biplot_points(
  X = biplot_data[, 1:num_vars, drop = FALSE],
  proj = pc.values,
  standardise = standardise_use,
  CVA_method = CVA_method,
  V = V,
  n = 2,                 # pick two points
  method = "identify"    # or "locator"
)

picked$idx      # row numbers in your_data
picked$picked   # data.frame with row, z_x, z_y, and the original data row
