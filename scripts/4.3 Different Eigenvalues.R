############################################################
# Section 4.3: Different eigenvalues (Cleaned, Same Logic)
# Description:
#   Explore counterfactuals across all 2D eigenvector (PC) combinations.
#   For each PC pair, rotate onto the plane, generate a grid, run the
#   boundary search for the target data point, and record distances.
#   Select the nearest pair according to Mahalanobis-like distance.
#
# Notes:
#   * Logic unchanged from the original script, only formatting,
#     comments, and safety guards added for clarity.
#   * Expects objects from previous steps: Vmat, tVmat, model_use,
#     model_select, biplot_data, Wmat, etc.
############################################################

###############################
# 0.1) Setup: target data point #
###############################

#tdp <- c(3, 13, 26, 2)[4]   # choose a test data point index
#tdp <- 1                  # PIMA

biplot_dataT[tdp, ]          # inspect chosen row


# Test + training data (filter by training convex hull)
rotate_data <- rbind(biplot_dataT[tdp, ], biplot_dataT)

################
# 0.2) Filters #
################

set_filters <- list(TRUE)   # default

# Examples retained from original:
# set_filters <- list(expr(Sepal.Width > 3))

# set_filters <- list(expr(Age >= 33 & Age <= 61),
#                     expr(Pregnancies >= 4 & Pregnancies <= 21))

# set_filters <- list(expr(loan_amnt >= 33000))

############################
# 1) Prep structures, defs #
############################

{
  # All 2D PC combinations
  n_contrib <- t(combn(num_vars, 2))
  n_pcs <- nrow(n_contrib)

  # Storage
  pc.cf    <- as.data.frame(matrix(0, ncol = num_vars, nrow = n_pcs, dimnames = list(NULL, var_names)))
  pc.cfz   <- matrix(0, ncol = 2, nrow = n_pcs)
  pc.z     <- matrix(0, ncol = 2, nrow = n_pcs)
  pc.distz <- matrix(0, ncol = 1, nrow = n_pcs)
  pc.distz_eucl <- matrix(0, ncol = 1, nrow = n_pcs)

  rounding <- 2
  Rotate <- TRUE

  Wmat_inv <- solve(Wmat)
}

#############################################
# 2) Iterate over all 2D eigenvector pairs  #
#############################################

for (i in seq_len(n_pcs)) {
  cat("PC pair index:", i, "\n")
  pc.vects <- n_contrib[i, ]

  tryCatch({
    #############################
    # 2.1) Rotation to new plane
    #############################
    V_rotate <- Biplot_rotation(
      X = rotate_data,
      tdp = 1,
      standardise = standardise_use,
      Xcenter = X_center, Xsd = X_sd,
      CVA_method = CVA_method,
      proj = pc.vects,
      V = Vmat
    )

    V <- V_rotate$Vrho
    tV <- V_rotate$tVrho
    
    # Use base indices for 2D block projection
    pc.vects <- c(1, 2)
    # update the biplot for the new projection matrix
    {
      Vr <- V[,pc.vects , drop = FALSE]
      tVr <- tV[pc.vects, , drop = FALSE]
      
      biplot_input_search$biplot_out$ax.one.unit <- 1/(diag(t(tVr) %*% tVr)) * t(tVr)
      biplot_input_search$biplot_out$Lmat <- V
      biplot_input_search$biplot_out$Z <- as.matrix(rotate_data[,1:num_vars]) %*% Vr
      
    #   biplot_input$biplot_out$ax.one.unit <- 1/(diag(t(tVr) %*% tVr)) * t(tVr)
    #   biplot_input$biplot_out$Lmat<- V
    # # biplot_input$biplot_out$Vr <- Vr
    #   biplot_input$biplot_out$Z <- biplot_input$biplot_out$X %*% biplot_input$biplot_out$Vr
    }
    
    #######################################
    # 2.2) Build biplot grid for rotation #
    #######################################
    biplot_grid_rotate <- biplot_plane(
      X = rotate_data,
      biplot_plot = biplot_input_search$biplot_out,
      model_use, var_names, numvars = num_vars,
      model_select = model_select,
      cutoff = cut_off,

      standardise = standardise_use,
      Xcenter = X_center, Xsd = X_sd,
      CVA_method = CVA_method,

      proj = pc.vects,
      V = V, tV = tV,

      polygon = NA,
      no_polygon = FALSE, # F = create a polygon , with includes the training data based on rotated data
      calc_hull = TRUE, # T = only select boundaries within the polygon
      calc_ct_in_hull = FALSE,
      outlie = 1,

      m = m_grid,
      b_margin = 1 / (10^rounding), rounding = rounding
    )

    ####################################
    # 2.3) Boundary search in Z-space  #
    ####################################
    biplot_search_ev <- Biplot_boundary_search(
      biplot_grid_rotate,
      
      model_select = model_select,
      tdp = 1,

      standardise = standardise_use,
      Xcenter = X_center, Xsd = X_sd,
      CVA_method = CVA_method,

      proj = pc.vects,
      V = V, tV = tV,

      set_filters = set_filters,
      train_ranges = train_ranges
    )

    ############################################################
    # 2.4) Save CF results for this pair (X-space and Z-space) #
    ############################################################
    pc.cf[i, ]  <- as.numeric(biplot_search_ev$B_boundary_x_value)
    pc.cfz[i, ] <- as.numeric(biplot_search_ev$B_boundary_z_value)
    pc.z[i, ]   <- as.numeric(biplot_grid_rotate$Z.st[1, ])

    # Distance vector in Z-space
    vec_to_boundary_Z <- as.matrix(biplot_search_ev$B_boundary_z_value - biplot_grid_rotate$Z.st[1, , drop = FALSE])
    vec_to_boundary   <- vec_to_boundary_Z %*% tVr

    # Mahalanobis-like distance
    pc.distz[i] <- (as.matrix(vec_to_boundary, nrow = 1)) %*% Wmat_inv %*% t(as.matrix(vec_to_boundary, nrow = 1))
    # pc.distz_eucl[i] <- apply((vec_to_boundary)^2, 1, "sum")
  }, error = function(e) {
    message("Skipping index ", i, " due to error: ", e$message)
  })
}

#########################
# 3) Collate and select #
#########################

{
  pc.distz_row <- cbind(seq_len(n_pcs), pc.distz)
  pc.distz_row <- pc.distz_row[pc.distz_row[, 2] != 0, , drop = FALSE]

  nearest_pc <- ifelse(nrow(pc.distz_row) == 1,
                       pc.distz_row[which.min(pc.distz_row[, 2])],
                       pc.distz_row[which.min(pc.distz_row[, 2]), 1])

  boundary <- pc.cf[nearest_pc, ]
  pc.rotate <- n_contrib[nearest_pc, ]
  print(pc.rotate)
}

##########################
# 4) Apply final rotation #
##########################

{
  V_rotate <- Biplot_rotation(
    X = rotate_data,
    tdp = 1,
    standardise = standardise_use,
    Xcenter = X_center, Xsd = X_sd,
    CVA_method = CVA_method,
    proj = pc.rotate,
    V = Vmat
  )

  V <- V_rotate$Vrho
  tV <- V_rotate$tVrho
  
}

###############################
# 5) Final biplot grid & plot #
###############################

# Use base indices for 2D block projection
pc.values <- c(1, 2)

    # update the biplot for the new projection matrix
    {
      Vr <- V[,pc.vects , drop = FALSE]
      tVr <- tV[pc.vects, , drop = FALSE]
      
      biplot_input_search$biplot_out$ax.one.unit <- 1/(diag(t(tVr) %*% tVr)) * t(tVr)
      biplot_input_search$biplot_out$Lmat<- V
    #  biplot_input$biplot_out$Vr <- Vr
      biplot_input_search$biplot_out$Z <- biplot_input_search$biplot_out$X %*% Vr
    }


biplot_grid_final <- biplot_plane(
  X = rotate_data,
  biplot_plot = biplot_input_search$biplot_out,
  model_use, var_names, numvars = num_vars,
  model_select = model_select,
  cutoff = cut_off,

  standardise = standardise_use,
  Xcenter = X_center, Xsd = X_sd,
  CVA_method = CVA_method,

  proj = pc.values,
  V = V, tV = tV,

  polygon = NA,
  no_polygon = FALSE,
  calc_hull = TRUE,
  calc_ct_in_hull = FALSE,
  outlie = 1,

  m = m_grid,
  b_margin = 1 / (10^rounding), rounding = rounding
)

var_names

plot_biplotEZ(
  grid = biplot_grid_final,
  biplot_plot = biplot_input_search$biplot_out,
  k = c(1,2,3,4,5,6,7,8,9),
  label_dist = c(0,0.5,0.5,0,0.5,0.5,0.5,0.5,0.5),
  tdp = 1,
  cex_z = 1.5,
  tick.label.cex = 0.7,
  label_dir = "Paral",
  ticks_v = 4,
  new_title = "Rotated 2 and 4 Component biplot"
)

point_inter_proj(
  X = rotate_data[, 1:num_vars],
  V = V,
  tdp = 1,
  plot = TRUE,
  proj = pc.values,
  col_x = "green",
  cex_x = 1.5,
  standardise = standardise_use,
  CVA_method = CVA_method
)

#########################################
# 6) Final boundary search on final grid #
#########################################

biplot_search <- Biplot_boundary_search(
  biplot_grid_final,
  model_select = model_select,
  tdp = 1,

  standardise = standardise_use,
  Xcenter = X_center, Xsd = X_sd,
  CVA_method = CVA_method,

  proj = pc.values,
  V = V, tV = tV,

  set_filters = set_filters,
  train_ranges = train_ranges
)

point_inter_proj(
  X = biplot_search$B_boundary_x_value,
  V = V,
  plot = TRUE,
  proj = pc.values,
  col_x = "darkgreen",
  cex_x = 1.5,
  standardise = standardise_use,
  CVA_method = CVA_method
)

# Inspect
rotate_data[1, 1:num_vars]
biplot_search$B_boundary_x_value
pc.cf[nearest_pc, ]

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

