############################################################
# Section 1: Biplot-based Model Utilities
# Description:
#   Utilities for model prediction wrappers, data preparation for
#   biplot visualizations (PCA/CVA), grid coloring and contouring on the
#   biplot plane, plotting helpers, indicator matrices, and Shapley-value
#   approximations for local feature contributions. 
#
# Notes:
#   * Function parameters are documented using roxygen-style comments.
############################################################

#############################
# 1) Libraries
#############################
suppressPackageStartupMessages(
{
  library(MASS)
  library(tidyverse)
  library(pdist)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(grid)

  library(e1071)
  library(gbm)
  library(mgcv)
  library(rpart)
  library(MLmetrics)

  ## Shapley
  library(dplyr)
  library(combinat)

  # for coordinates in biplot LDB
  library(sp)
  library(devtools)

  library(biplotEZ)
})


###########################################################
# 4) Data Preparation for Model Prediction & Evaluation
###########################################################
#' Prepare dataset for prediction, evaluation, and optional polygon filtering
#'
#' Adds predictions, confusion labels, and plotting attributes.
#'
#' @param model_use Model object
#' @param data_use Data frame to process
#' @param num_vars Number of feature columns (assumed to be the first num_vars)
#' @param cutoff Classification cutoff (default 0.5)
#' @param model_select Character indicating model family (see pred_function)
#' @param target_class Target class index for LDA 1-vs-all reduction (default 2)
#' @param rounding Rounding for prediction values (default 2)
#' @param proj Projection indices (for biplot; default pc.values)
#' @param CVA_method Logical; whether CVA is used (affects scaling)
#' @param standardise Logical; standardize features for PCA (default standardise_use)
#' @param polygon Optional SpatialPolygons object to filter points inside
#' @param Xcenter Numeric vector of feature means
#' @param Xsd Numeric vector of feature sds
#' @param V Projection/loadings matrix
#'
#' @return list with element data_use (augmented data.frame)
set_model_data <- function(model_use, data_use, num_vars = num_vars,
                           cutoff = 0.5,
                           model_select = model_select,
                           target_class = 2,
                           rounding = 2,
                           proj = pc.values,
                           CVA_method = CVA_method,
                           standardise = standardise_use, # set to FALSE for CVA plots
                           polygon = NA, # train_polygon
                           Xcenter = X_center,
                           Xsd = X_sd,
                           V = Vmat) {

  # Remove data outside polygon, if provided
  if(is(polygon, "SpatialPolygons"))  {
    sv <- Xsd
    if(CVA_method == TRUE) standardise <- FALSE # set to FALSE if CVA is TRUE
    if(standardise  == FALSE) sv <- FALSE
    X.st <- scale(data_use[, 1:num_vars], center = Xcenter, scale = sv)
    Z_values <- X.st %*% V[, proj]
    colnames(Z_values) <- c("x", "y")

    # Select Z data points inside hull polygon
    Z_data_points <- as.data.frame(Z_values)
    coordinates(Z_data_points) <- ~x + y
    Z_data_points$inside_polygon <- !is.na(over(Z_data_points, polygon))
    data_use <- data_use[Z_data_points$inside_polygon == TRUE, ]
  }

  # Get predictions
  data_use$pred_value <- pred_function(model_use = model_use,
                                       model_select = model_select,
                                       rounding = rounding,
                                       new_data = data_use)

  # Convert LDA multiclass to 0/1 for final comparison; retain original class
  if(model_select  == "LDA")
    data_use$class <- as.numeric(ifelse(data_use$class == target_class, 1, 0))

  # Round prediction
  data_use <- data_use %>% mutate(pred_value = round(pred_value, 6))

  # Predicted class based on cut-off value
  data_use$pred_use <- as.numeric(ifelse(data_use$pred_value < cutoff, 0, 1))

  # Correct and incorrect prediction
  data_use <- data_use %>% mutate(pred_correct = if_else(pred_use == class, 1, 0))

  # Store original order
  data_use <- data_use %>% mutate(Original_order = row_number())

  # Colour allocation for confusion categories (for plotting)
  data_use <- data_use %>% mutate(
    pred_confusion = case_when(
      data_use[, "class"] == data_use[, "pred_use"] & data_use[, "class"] == 1 ~ 1,   # TP
      data_use[, "class"] == data_use[, "pred_use"] & data_use[, "class"] == 0 ~ 2,   # TN
      data_use[, "class"] != data_use[, "pred_use"] & data_use[, "class"] == 0 ~ 3,   # FP
      TRUE ~ 4                                                                           # FN
    )
  )

  data_use <- data_use %>% mutate(
    pred_col = case_when(
      data_use[, "class"] == data_use[, "pred_use"] & data_use[, "class"] == 1 ~ "red",     # TP
      data_use[, "class"] == data_use[, "pred_use"] & data_use[, "class"] == 0 ~ "blue",    # TN
      data_use[, "class"] != data_use[, "pred_use"] & data_use[, "class"] == 0 ~ "purple",  # FP
      TRUE ~ "orange"                                                                           # FN
    )
  )

  # Ensure different classes have different symbols. Class 1 = triangle
  data_use <- data_use %>% mutate(pch_x = if_else(data_use$class == 0, 16, 17))

  # Print quick metrics
  print("Confusion matrix: ")
  print(ConfusionMatrix(data_use$pred_use, data_use$class))
  print("Mean prediction: ")
  print(1 - mean(data_use$class)) # for comparison against accuracy: lift provided
  print("")
  print("Accuracy: ")
  print(Accuracy(data_use$pred_use, data_use$class))

  return(list(data_use = data_use))
}


#############################################
# 5) Biplot Input Calculation (PCA / CVA)
#############################################
#' Compute PCA/CVA biplot objects and centering/scaling
#'
#' @param X Feature matrix/data.frame (default train_data[, 1:num_vars])
#' @param standardise Logical; standardize features for PCA
#' @param main_heading Title for biplot
#' @param proj Integer vector of eigenvector indices to retain
#' @param CVA_method Logical; use CVA (if FALSE, PCA)
#' @param CVA_classes Factor of class labels for CVA
#'
#' @return list with X, X_center, X_sd, and biplot_out (biplotEZ object)
biplot_input_calc <- function(X = train_data[, 1:num_vars],
                              standardise = standardise_use, # for PCA
                              main_heading = "Test",
                              proj = 1:2,
                              CVA_method = FALSE,
                              CVA_classes = train_data$Original_class) { # default when CVA FALSE

  X_center <- apply(X, 2, "mean")
  X_sd <- apply(X, 2, "sd")

  # PCA
  if (CVA_method == FALSE) {
    create_bp <- biplot(X, scaled = standardise, Title = main_heading) |>
      PCA(e.vects = proj)
  }

  # CVA
  if (CVA_method == TRUE) {
    create_bp <- biplot(X, scaled = FALSE, Title = main_heading) |>
      CVA(e.vects = proj, classes = CVA_classes, show.class.means = FALSE)
  }

  return(list(
    #X = X,
    X_center = X_center,
    X_sd = X_sd,
    biplot_out = create_bp
  ))
}


########################################################
# 6) Biplot Plane Coloring, Grid, and Contours
########################################################

## >> Replace with "1.1 biplot_plane_optimized.R"

#' library(aplpack)
#' library(sp) # spatial data on the 2-D surface for polygon calculation
#' 
#' #' Create color-coded prediction grid on the biplot plane
#' #'
#' #' @param X Data with at least the first num_vars feature columns and plotting columns
#' #' @param model Model object used for predictions
#' #' @param varnames Character vector of feature names
#' #' @param numvars Number of features
#' #' @param model_select Character indicating model family
#' #' @param cutoff Decision threshold for contours
#' #' @param standardise Logical; standardization flag for back-projection
#' #' @param proj Indices of projection vectors
#' #' @param V Matrix of loadings (columns are eigenvectors)
#' #' @param tV Transpose of V
#' #' @param CVA_method Logical; whether CVA is used
#' #' @param Xcenter Numeric vector of means
#' #' @param Xsd Numeric vector of sds
#' #' @param m Integer grid density (per side)
#' #' @param polygon Optional SpatialPolygons to restrict grid area
#' #' @param no_polygon Logical; if TRUE, use full grid as polygon
#' #' @param calc_hull Logical; if TRUE, remove grid points outside polygon
#' #' @param calc_ct_in_hull Logical; if TRUE, compute contours with hull edge at cutoff
#' #' @param outlie Fraction for convex hull trimming (aplpack::plothulls)
#' #' @param b_margin Band around cutoff for contours
#' #' @param rounding Floor-based rounding for prediction scores
#' #' @param remove_filter Logical; retain for compatibility (not used)
#' #' @param return_values Logical; return list of grid and metadata
#' #'
#' #' @return list with Zgrid, Xgrid, grid.prob, Z.st, pred_col, pred_use, class,
#' #'   col.value, min_val, max_val, polygon, ct (if return_values=TRUE)
#' biplot_plane <- function(X = biplot_data,
#'                          biplot_plot = biplot_input$biplot_out,
#'                          model = model_use, varnames, numvars,
#'                          model_select = model_select,
#'                          cutoff = 0.5,
#'                          standardise = FALSE, # default is FALSE to include in CVA method
#'                          proj = c(1, 2),
#'                          V = svdx$v, # data from X train to calibrate the biplot
#'                          tV = t(svdx$v),
#'                          CVA_method = FALSE,
#'                          Xcenter = X_center, Xsd = X_sd,
#'                          m = 300,
#'                          polygon = NA,
#'                          no_polygon = FALSE, # X.Poly not currently used, density not specified
#'                          calc_hull = FALSE,
#'                          calc_ct_in_hull = FALSE,
#'                          outlie = 1,
#'                          b_margin = 0.01, # contour lines and margin around cutoff
#'                          rounding = 2, # must be consistent with the b-margin value
#'                          remove_filter = FALSE,
#'                          return_values = TRUE) {
#' 
#'   Vr <- V[, proj]
#'   tVr <- tV[proj, ]
#' 
#'   {
#'     sv <- Xsd
#'     if(CVA_method == TRUE) standardise <- FALSE # set to FALSE if CVA is TRUE
#'     if(standardise  == FALSE) sv <- FALSE
#'     X.st <- scale(X[, 1:num_vars], center = Xcenter, scale = sv)
#'   }
#' 
#'   # In case of blank, create a plot to ensure new plot steps work correctly
#'   if(dev.cur() == 1)  plot(1, 1, type = "n", asp = TRUE)
#' 
#'   ## Plot to get the grid area
#'   dev.new()
#' 
#'   {
#'     # First plot to get plot dimensions of x and y axes
#'     Z.st <- (X.st %*% Vr)
#'     colnames(Z.st) <- c("x", "y")
#' 
#'     if(is(polygon, "SpatialPolygons")) plot(polygon, asp = 1)
#'     ## NB: create the biplot
#'     # First plot to get plot dimensions of x and y axes: use direct, as will be used later
#'     if(!is(polygon, "SpatialPolygons")) plot(Z.st, asp = 1, ylab = "Y[2]", xlab = "Y[1]")
#'       # biplot_plot |>
#'       #   samples(opacity = 0, which = NULL) |>
#'       #   axes(col = "grey")  |> plot()
#' 
#'     # Get the plot dimensions
#'     dev.set(dev.cur())
#' 
#'     minx <- par("usr")[1]
#'     miny <- par("usr")[3]
#'     maxx <- par("usr")[2]
#'     maxy <- par("usr")[4]
#' 
#'     dev.off()
#'   }
#' 
#'   ## Get the grid points
#'   {
#'     # min/max values of x and y dimensions, allow square plot later
#'     min_val <- (min(minx, miny))
#'     max_val <- (max(maxx, maxy))
#' 
#'     # Grid points in plot area
#'     yseq <- xseq <- seq(from = min_val, to = max_val, length.out = m)
#'     Zgrid <- as.matrix(expand.grid(xseq, yseq))
#'     colnames(Zgrid) <- c("x", "y")
#' 
#'     # Convert back to X dimension space
#'     Xgrid <- as.data.frame(Zgrid %*% tVr)
#' 
#'     # Return the unstandardised Xgrid
#'     if(standardise == TRUE) Xgrid <- sweep(Xgrid, 2, Xsd, "*")
#'     Xgrid <- sweep(Xgrid, 2, Xcenter, "+")
#'     Xgrid <- as.data.frame(Xgrid)
#' 
#'     # Set feature names
#'     colnames(Xgrid) <- varnames
#'   }
#' 
#'   ## Create convex hull (polygon) for the prediction area
#'   # use existing polygon provided unless none provided or no_polygon == TRUE
#'   {
#'     if(!is(polygon, "SpatialPolygons")) {
#'       # create polygon if no existing polygon was provided
#'       dev.new()
#'       polygon_points <- plothulls(
#'         x = Z.st[, 1], y = Z.st[, 2], asp = 1,
#'         xlim = c(min_val, max_val), ylim = c(min_val, max_val),
#'         fraction = outlie, n.hull = 1, main, add = FALSE, col.hull = 2,
#'         lty.hull = 1, lwd.hull = 2, density = 1
#'       )
#'       dev.off()
#'       polygon <- SpatialPolygons(list(Polygons(list(Polygon(polygon_points)), ID = "polygon")))
#'     }
#' 
#'     if(no_polygon == TRUE){
#'       # Use entire grid if not limiting to polygon: overwrite the existing polygon
#'       dev.new()
#'       polygon_points <- plothulls(
#'         x = Zgrid[, 1], y = Zgrid[, 2], asp = 1,
#'         xlim = c(min_val, max_val), ylim = c(min_val, max_val),
#'         fraction = 1, n.hull = 1, main, add = FALSE, col.hull = 2,
#'         lty.hull = 1, lwd.hull = 2, density = 1
#'       )
#'       dev.off()
#'       polygon <- SpatialPolygons(list(Polygons(list(Polygon(polygon_points)), ID = "polygon")))
#'     }
#' 
#'     ## use polygon to set the convex hull filter region: only apply if convex hull is true
#'     # Select grid points inside hull polygon
#'     G_data_points <- as.data.frame(Zgrid)
#'     coordinates(G_data_points) <- ~x + y
#'     G_data_points$inside_polygon <- !is.na(over(G_data_points, polygon))
#' 
#'     # Remove data points outside the polygon
#'     Z_data_points <- as.data.frame(Z.st)
#'     coordinates(Z_data_points) <- ~x + y
#'     Z_data_points$inside_polygon <- !is.na(over(Z_data_points, polygon))
#'   }
#' 
#'   # Score each grid point using the selected model
#'   {
#'     grid.prob <- pred_function(model_use = model_use,
#'                                model_select = model_select,
#'                                rounding = rounding,
#'                                new_data = Xgrid)
#'   }
#' 
#'   # Color palette for grid values
#'   {
#'     col.vec <- colorRampPalette(c("deepskyblue", "white", "lightsalmon"))(100 + 1)
#'     col.value <- col.vec[floor(grid.prob * 100) + 1]
#'   }
#' 
#'   ###############
#'   ## Contours  ##
#'   ###############
#'   {
#'     grid_contour <- matrix(grid.prob, ncol = m, byrow = FALSE)
#' 
#'     if(calc_ct_in_hull == TRUE) {
#'       # Compute contours in hull; make edge = cutoff to avoid double grouping
#'       grid_contour[G_data_points$inside_polygon == FALSE] <- cutoff
#'       ct <- contourLines(xseq, yseq, grid_contour, levels = c(cutoff - b_margin, cutoff, cutoff + b_margin))
#'       ct <- ct[sapply(ct, function(cl) cl$level != cut_off)]
#'     }
#' 
#'     if(calc_ct_in_hull == FALSE)
#'       ct <- contourLines(xseq, yseq, grid_contour, levels = c(cutoff - b_margin, cutoff + b_margin))
#'   }
#' 
#'   ################################################
#'   # Remove grid points outside the polygon
#'   ################################################
#'   
#'   ## use convex hull to remove point, based on polygon: No polygon implies no convex hull filter even if true
#'   
#'   pred_col <- X$pred_col
#'   pred_use <- X$pred_use
#'   class <- X$class
#' 
#'   if(calc_hull == TRUE) {
#'     grid.prob <-  grid.prob[G_data_points$inside_polygon == TRUE]
#'     col.value <-  col.value[G_data_points$inside_polygon == TRUE]
#'     Zgrid <- Zgrid[G_data_points$inside_polygon == TRUE, ]
#'     Xgrid <- Xgrid[G_data_points$inside_polygon == TRUE, ]
#' 
#'     Z.st <- Z.st[Z_data_points$inside_polygon == TRUE, ]
#'     pred_col <- pred_col[Z_data_points$inside_polygon == TRUE]
#'     pred_use <- pred_use[Z_data_points$inside_polygon == TRUE]
#'     class <- class[Z_data_points$inside_polygon == TRUE]
#'   }
#' 
#'   ## base plot
#'   dev.off()
#' 
#'   # Retained for compatibility
#'   if(remove_filter == TRUE) {
#'     col.value <- col.value[rows_to_exclude]
#'     grid.prob <- grid.prob[rows_to_exclude]
#'     Xgrid <- Xgrid[rows_to_exclude, ]
#'     Zgrid <- Zgrid[rows_to_exclude, ]
#'   }
#' 
#'   if(return_values == TRUE)
#'     return(list(
#'       Zgrid = Zgrid, Xgrid = Xgrid, grid.prob = grid.prob, Z.st = Z.st,
#'       pred_col = pred_col, pred_use = pred_use, class = class,
#'       col.value = col.value,
#'       min_val = min_val, max_val = max_val, # for 3d plots
#'       polygon = polygon, ct = ct
#'     ))
#' }


###############################################
# 7) Plot Helpers for biplotEZ Grid & Points
###############################################
#' Plot biplot grid, samples, axes, and optional contours using biplotEZ
#'
#' @param grid Output from biplot_plane
#' @param k Integer index of variable to offset label
#' @param label_dist Numeric label offset distance
#' @param tdp Index/indices of points to plot (default all)
#' @param V Loadings matrix (passed through)
#' @param no_grid Logical; hide colored grid
#' @param no_points Logical; hide sample points
#' @param no_contour Logical; hide contours
#' @param new_title Optional new biplot title
#' @param ticks_v Number of ticks on axes
#' @param cex_z Point size for samples
#' @param label_dir "Hor" or "Rad" label direction
#' @param tick.label.cex Axis tick label size
#' @param which Variable indices to label
#' @param X.names Custom variable names for axes
#'
#' @return Invisibly returns NULL (plot function)
plot_biplotEZ <- function(grid = biplot_grid,
                          biplot_plot = biplot_input$biplot_out,
                          k = 0,
                          label_dist = 0.5,
                          tdp = NA,
                          no_grid = FALSE,
                          no_points = FALSE,
                          no_contour = FALSE,
                          new_title = NA,
                          ticks_v = 1,
                          cex_z = 0.5,
                          label_dir = "Hor",
                          tick.label.cex = 0.6,
                          which = 1:num_vars,
                          X.names = var_names) {

  ## Set data from the biplot function
  Z.st <- grid$Z.st
  pred_col <- grid$pred_col
  Zgrid <- grid$Zgrid
  grid.value <- grid$col.value
  ct <- grid$ct

  # Update title if provided
  if(!is.na(new_title))
    biplot_plot$Title <- new_title

  # Move position of labels
  label.line.vec <- rep(0, num_vars)
  label.line.vec[k] <- label_dist

  # Base plot
  biplot_plot |>
    samples(opacity = 0, which = NULL) |>
    axes(col = "grey",
         label.dir = label_dir,
         which = which,
         X.names = X.names,
         tick.label.cex = tick.label.cex,
         ticks = ticks_v,
         label.line = label.line.vec)  |> plot()

  ## Add the grid/dense map
  if(no_grid == FALSE) {
    points(Zgrid, type = "p", col = grid.value, pch = 15, cex = 0.5)
  }

  # Add points
  if(is.na(tdp[1])) tdp <- 1:nrow(Z.st)

  if(no_points == FALSE)
    points(x = Z.st[tdp, 1],
           y = Z.st[tdp, 2],
           type = "p",
           col = pred_col[tdp],
           pch = 16, cex = cex_z)

  # Add axes on top
  biplot_plot |>
    samples(opacity = 0, which = NULL) |>
    axes(col = "grey22",
         label.dir = label_dir,
         which = which,
         X.names = X.names,
         tick.label.cex = tick.label.cex,
         ticks = ticks_v,
         label.line = label.line.vec)  |> plot(add = TRUE)

  # Add contour lines
  if(no_contour == FALSE) {
    for(i in 1:length(ct)) {
      lines(ct[[i]]$x, ct[[i]]$y)
    }
  }

  invisible(NULL)
}


#########################################################
# 8) Project Points into Biplot Coordinates (PCA/CVA)
#########################################################
#' Project observations into biplot space
#'
#' @param X Feature data.frame/matrix
#' @param proj Projection indices (default pc.values)
#' @param standardise Logical; whether to standardize by Xsd
#' @param CVA_method Logical; if TRUE, force standardise = FALSE
#' @param tdp Indices to plot (if plot = TRUE)
#' @param Xcenter Numeric means
#' @param Xsd Numeric sds
#' @param V Loadings/eigenvector matrix
#' @param plot Logical; draw points (base graphics)
#' @param col_x Point color
#' @param pch_x Point pch
#' @param cex_x Point size
#'
#' @return Matrix with columns (x, y) of projected coordinates
point_inter_proj <- function(X, proj = pc.values,
                             standardise = FALSE, CVA_method = FALSE,
                             tdp = NA,
                             Xcenter = X_center, Xsd = X_sd,
                             V = V,
                             plot = FALSE,
                             col_x = 8, pch_x = 16, cex_x = 1.5) {

  # Apply the standardisation of X: class or testing data
  sv <- Xsd
  if(CVA_method == TRUE) standardise <- FALSE # set to FALSE if CVA is TRUE
  if(standardise  == FALSE) sv <- FALSE
  X.st <- scale(X, center = Xcenter, scale = sv)

  # Project to biplot coordinates
  Z_values <- X.st %*% V[, proj]
  colnames(Z_values) <- c("x", "y")

  if(is.na(tdp[1])) {
    tdp <- 1:nrow(Z_values)
  }
  
  if (plot == TRUE)
    points(x = Z_values[tdp, 1], y = Z_values[tdp, 2], cex = cex_x, col = col_x, pch = pch_x)

  Z_values <- Z_values[tdp,]
  
  return(Z_values = Z_values)
}


###############################################
# 9) Filter Helpers (Dynamic Range & AND-combine)
###############################################
library(dplyr)
library(rlang)
library(mlbench)

#' Get the range of values in each variable in a data.frame
#' Remove new data points that are outside the range of the training data
#'
#' @param data data.frame
#'
#' @return List of expressions like var >= min & var <= max
get_variable_ranges <- function(data) {
  vars <- names(data)

  filters <- map(vars, function(var) {
    rng <- range(data[[var]], na.rm = TRUE)
    expr((!!sym(var)) >= !!rng[1] & (!!sym(var)) <= !!rng[2])
  })

  filters
}


#' Evaluate a list of expressions row-wise with AND logic
#'
#' @param data data.frame to filter
#' @param filter_exprs List of expressions (e.g., from get_variable_ranges)
#'
#' @return Logical vector indicating rows that pass all filters
get_filter_logical_vector <- function(data, filter_exprs) {
  # Combine multiple expressions with AND
  row_filter_expr <- reduce(filter_exprs, function(x, y) expr((!!x) & (!!y)))
  # Evaluate the logical expression on the dataset
  eval_tidy(row_filter_expr, data = data)
}


#########################################################
# 10) Shapley Value Calculation (Exact via Permutations)
#########################################################

## >> replace with "1.2 Run_Shapley_optimized.R"

#' Compute exact Shapley values for one or all observations
#'
#' WARNING: Complexity is factorial in the number of variables.
#'
#' @param model_use Model object for prediction
#' @param num_vars Number of variables (columns of test_data used)
#' @param test_data Data frame of observations (X)
#' @param B_boundary Data frame of boundary values (same cols as test_data)
#' @param individual Logical; if TRUE, compute only for tdp
#' @param tdp Integer row index to compute for when individual=TRUE
#' @param lda,rpart_used,svm_u Flags retained for compatibility
#' @param rounding Integer decimals for floor-based rounding
#'
#' @return list with shapley_store (matrix: nrow(test_data) x num_vars)

# Run_Shapley <- function(model_use, num_vars, test_data,
#                         B_boundary,
#                         individual = FALSE, tdp = 1,
#                         lda = FALSE, rpart_used = FALSE, svm_u = FALSE,
#                         rounding = rounding) {
# 
#   n_orders <- factorial(num_vars)
#   n_contrib <- matrix(unlist(permn(num_vars)), ncol = num_vars, byrow = TRUE)
#   margin_calc <- rep(0, num_vars)
#   margin_store <- matrix(0, nrow = n_orders, ncol = num_vars)
#   shapley_store <- matrix(0, ncol = num_vars, nrow = nrow(test_data))
# 
#   run_s <- 1
#   run_e <- nrow(test_data)
# 
#   # Single data point
#   if(individual  == TRUE) {
#     run_s <- tdp
#     run_e <- tdp
#   }
# 
#   for(t in run_s:run_e) {
#     start_time = Sys.time()
#     print(t)
# 
#     start_data <- as.numeric(test_data[t, ])
#     new_data <- as.numeric(B_boundary[t, ])
# 
#     for(j in 1:n_orders){
#       s_order <- n_contrib[j, ]
#       s_contrib <- matrix(0, nrow = num_vars, ncol = num_vars)
# 
#       # Build marginal contribution mask for this permutation
#       for(i in 1:num_vars){
#         s_contrib[, s_order[i]] <- c(rep(0, i - 1), rep(1, num_vars - i + 1))
#       }
# 
#       # Add opening row of zeros for starting value
#       s_contrib <- rbind(c(rep(0, num_vars)), s_contrib)
# 
#       # Replace 0/1s with actual values
#       s_contrib <- sweep(1 - s_contrib, 2, start_data, FUN = "*") + sweep(s_contrib, 2, new_data, FUN = "*")
#       colnames(s_contrib) <- var_names
#       s_contrib <- as.data.frame(s_contrib)
# 
#       # Predict
#       s_value <- pred_function(model_use = model_use,
#                                model_select = model_select,
#                                rounding = rounding,
#                                new_data = s_contrib)
# 
#       # Marginal differences
#       for(i in 1:num_vars){
#         margin_calc[i] <-  s_value[i + 1] - s_value[i]
#       }
#       margin_store[j, s_order] <- margin_calc
#     }
# 
#     s_value[-4] - s_value[-1]  # (Retained from original script)
# 
#     # Sum marginal contributions per variable across permutations
#     shapley_store[t, ] <- apply(margin_store, 2, sum) / n_orders
# 
#     end_time = Sys.time()
#     print(end_time - start_time)
#     print("...")
#   }
# 
#   return(list(shapley_store = shapley_store))
# }
