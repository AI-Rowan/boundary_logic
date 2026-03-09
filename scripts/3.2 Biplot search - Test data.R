
############################################################
# Section 3.2 : Biplot plot and search with Test data
# Description:
#   Continues from file 3.1, with the test data now generated
#   uses the biplot projection functions generated with the Training data
#   Have the option in step 4 to retaining the biplot on the test data
#   Have the option in step 4.1 to rotate the biplot for the tdp selected
#   Generates biplot_dataT and biplot_grid_Test, which can be used later
#   Identifies individual counterfactual

# Notes:
#   * Requires the helper scripts and RData files referenced via source()/load().
############################################################

############################################################
# Restart Step 3: Use TEST data
############################################################

  model_data <- test_data #[sample(nrow(test_data),10000,replace = F),] 
  head(model_data)
  
{
  # Ensure TEST data is within TRAIN ranges
  rows_to_exclude <- get_filter_logical_vector(model_data, train_ranges)
  model_data <- model_data[rows_to_exclude, ]

  # Create biplot-ready test data, limited by training polygon if provided
  biplot_dataT <- set_model_data(
    model_use,
    data_use = model_data,
    num_vars = num_vars,
    cutoff = cut_off,
    rounding = 2,
    target_class = tdp_class,
    model_select = model_select,
    polygon = NA,#final_polygon,          # If provided; NA by default, can also use final_polygon at end of training process
    proj = pc.values,
    CVA_method = CVA_method,
    standardise = standardise_use
  )$data_use
}

############################################################
# Step 4 (optional): Recompute V from TEST (if desired)
############################################################
# If you want a new V/tV for TEST, re-run biplot_input_calc here.
# This section is commented out to preserve the original flow.
# {
#   biplot_inputT <- biplot_input_calc(
#     X = biplot_dataT[, c(1:num_vars)],
#     standardise = standardise_use,
#     CVA_method = CVA_method,
#     #main_heading = "Create Biplot using projection matrix V2",
#     CVA_classes = as.factor(biplot_dataT[, "pred_col"])
#   )
# }

############################################################
# Step 4 (part 2): Select rotation and compute Vrho
############################################################

# allocate results to biplot    
biplot_input_search <- biplot_input
# biplot_input_search <- biplot_inputT

# allocate the correct biplot values to the underlying V and g
# {  
#   V  <- Vmat <- biplot_input_search$biplot_out$Lmat
#   tV <- tVmat <- solve(Vmat)
#   pc.values <- c(1, 2)
#   Vr <- V[,pc.values]
#   tVr <- tV[pc.values,]
#   
#   Wmat <- biplot_input_search$biplot_out$Wmat
#   biplot_ax.one.unit <- biplot_input_search$biplot_out$ax.one.unit
#   
#   X_center <- biplot_input_search$X_center
#   X_sd <- biplot_input_search$X_sd
#   
# }  

# Choose target data point index for rotation (unconstrained solution)
  #Iris: 30, PIMA: 12,  loans: tdp <- c(3, 13, 26, 2)[1] 
tdp <- 3#,nrow(biplot_dataT)
biplot_dataT[tdp, ]
#biplot_dataT %>% filter(pred_confusion == 1)

pc.values <- c(1, 2)
Rotate <- F
{
  if (Rotate == TRUE) {
    V_rotate <- Biplot_rotation(
      X = biplot_dataT,
      tdp = tdp,
      standardise = standardise_use,
      Xcenter = X_center,
      Xsd = X_sd,
      CVA_method = CVA_method,
      proj = pc.values,
      V = Vmat  # from selected biplot input data
    )

    # update the biplot for the new projection matrix
    
      V <- V_rotate$Vrho
      tV <- V_rotate$tVrho
      Vr <- V[,pc.values]
      tVr <- tV[pc.values,]


      biplot_input_search$biplot_out$ax.one.unit <- 1/(diag(t(tVr) %*% tVr)) * t(tVr)
      biplot_input_search$biplot_out$Lmat <- V
      biplot_input_search$biplot_out$Mr <- Vr
      biplot_input_search$biplot_out$Linv
      biplot_input_search$biplot_out$Z <- biplot_input_search$biplot_out$X %*% Vr
      
    }
  
   if (Rotate == FALSE) {
     # allocate results to biplot    
     biplot_input_search <- biplot_input
     
     {  
       V  <- Vmat <- biplot_input_search$biplot_out$Lmat
       tV <- tVmat <- solve(Vmat)
       Vr <- V[,pc.values]
       tVr <- tV[pc.values,]
       
     #  biplot_input_search$ax.one.unit <- 1/(diag(t(tVr) %*% tVr)) * t(tVr)
     #  biplot_input_search$biplot_out$Lmat <- V
    #   biplot_input_search$biplot_out$Z <- biplot_input_search$biplot_out$X %*% Vr
     }  

  }
}

############################################################
# Step 5-9 (rotated): Build grid & contours for rotated V
############################################################

biplot_grid_Test <- biplot_plane(
  X = biplot_dataT,
  biplot_plot = biplot_input_search$biplot_out,
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
  polygon = NA,#train_polygon,     # Limit test to training area if desired: not applicable after a rotation
  no_polygon = F,           # Override polygon (use full area)
  calc_hull = F,           # Irrelevant when no_polygon = TRUE: Need to be TRUE for Biplot CF search
  calc_ct_in_hull = FALSE,
  m = 250,
  rounding = 2,
  b_margin = 1 / (10^2)
)

############################################################
# Step 9 (plot rotated grid & points)
############################################################
var_names

plot_biplotEZ(
  grid = biplot_grid_Test,
  biplot_plot = biplot_input_search$biplot_out,
  
#  k = c(1,2,3,4,5,6,7,8),
 # label_dist = c(0,0.5, 0.5,1.5,0, 0.5,2,1.5),
  
  k = c(1,2,3,4,5,6,7,8,9),
  label_dist = c(0,0.5, 1,1.5,2.5, 0.5,2,1.5,2),
  
  tick.label.cex = 0.85,
  
  #tdp = tdp,
  cex_z = 0.65,

  # 
  no_grid = F,
  no_contour = F,
  no_points = F,
  label_dir = "Paral",
  ticks_v = 3,
  new_title = ""
  #new_title = "Rotated Component biplot"
)

# Plot Ttdp# Plot TDP in rotated space
ztdp <- point_inter_proj(
  X = biplot_dataT[, 1:num_vars],
  V = V,
  tdp = tdp,
  plot = TRUE,
  col_x = "green",
  #col_x = biplot_dataT$pred_col,
  cex_x = 1.25,
  standardise = standardise_use,
  CVA_method = CVA_method
)

############################################################
# Steps 10–11: Run CF search (unconstrained or constrained)
############################################################
set_filters <- list(TRUE)

# Example constrained filters (commented):
# set_filters <- list(
#   expr(Age >= 33  & Age <= 61),
#   expr(Pregnancies >= 4 & Pregnancies <= 21)
# )

biplot_search_T <- Biplot_boundary_search(
  X = biplot_grid_Test, # contains the biplot polygon and determines if remove empty CTs
  model_select = model_select,
 # tdp = tdp,            # Use NA to search all points
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

# Plot the CF solution
plot_search <- point_inter_proj(
  X = biplot_search_T$B_boundary_x_value,
  V = V,
  tdp = tdp,
  plot = TRUE,
  col_x = "yellow",
  cex_x = 0.65,
  standardise = standardise_use,
  CVA_method = CVA_method
)

## connect tdp to cf
# segments(ztdp[,1], ztdp[,2], plot_search[,1], plot_search[,2], col = "grey70",lwd = 2)

# Click to pick points (press ESC to stop early)
picked <- pick_biplot_points(
  X = biplot_dataT[, 1:num_vars, drop = FALSE],
  proj = pc.values,
  standardise = standardise_use,
  CVA_method = CVA_method,
  V = V,
  n = 2,                 # pick two points
  method = "identify"    # or "locator"
)

picked$idx      # row numbers in your_data
picked$picked   # data.frame with row, z_x, z_y, and the original data row

# Inspect results
biplot_dataT[tdp, ]
biplot_search$B_boundary_x_value[1,]

############################################################
# End of Main Section
############################################################


