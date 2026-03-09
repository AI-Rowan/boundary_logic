
########################################################
# 6) Biplot Plane Coloring, Grid, and Contours
########################################################
library(aplpack)
library(sp) # spatial data on the 2-D surface for polygon calculation

#' Create color-coded prediction grid on the biplot plane
#'
#' @param X Data with at least the first num_vars feature columns and plotting columns
#' @param model Model object used for predictions
#' @param varnames Character vector of feature names
#' @param numvars Number of features
#' @param model_select Character indicating model family
#' @param cutoff Decision threshold for contours
#' @param standardise Logical; standardization flag for back-projection
#' @param proj Indices of projection vectors
#' @param V Matrix of loadings (columns are eigenvectors)
#' @param tV Transpose of V
#' @param CVA_method Logical; whether CVA is used
#' @param Xcenter Numeric vector of means
#' @param Xsd Numeric vector of sds
#' @param m Integer grid density (per side)
#' @param polygon Optional SpatialPolygons to restrict grid area
#' @param no_polygon Logical; if TRUE, use full grid as polygon
#' @param calc_hull Logical; if TRUE, remove grid points outside polygon
#' @param calc_ct_in_hull Logical; if TRUE, compute contours with hull edge at cutoff
#' @param outlie Fraction for convex hull trimming (aplpack::plothulls)
#' @param b_margin Band around cutoff for contours
#' @param rounding Floor-based rounding for prediction scores
#' @param remove_filter Logical; retain for compatibility (not used)
#' @param return_values Logical; return list of grid and metadata
#'
#' @return list with Zgrid, Xgrid, grid.prob, Z.st, pred_col, pred_use, class,
#'   col.value, min_val, max_val, polygon, ct (if return_values=TRUE)
biplot_plane <- function(X = biplot_data,
                         biplot_plot = biplot_input$biplot_out,
                         model = model_use, varnames, numvars,
                         model_select = model_select,
                         cutoff = 0.5,
                         standardise = FALSE, # default is FALSE to include in CVA method
                         proj = c(1, 2),
                         V = svdx$v, # data from X train to calibrate the biplot
                         tV = t(svdx$v),
                         CVA_method = FALSE,
                         Xcenter = X_center, Xsd = X_sd,
                         m = 300,
                         polygon = NA,
                         no_polygon = FALSE, # X.Poly not currently used, density not specified
                         calc_hull = FALSE,
                         calc_ct_in_hull = FALSE,
                         outlie = 1,
                         b_margin = 0.01, # contour lines and margin around cutoff
                         rounding = 2, # must be consistent with the b-margin value
                         remove_filter = FALSE,
                         return_values = TRUE) {
  
  Vr <- V[, proj]
  tVr <- tV[proj, ]
  
  {
    sv <- Xsd
    if(CVA_method == TRUE) standardise <- FALSE # set to FALSE if CVA is TRUE
    if(standardise  == FALSE) sv <- FALSE
    X.st <- scale(X[, 1:num_vars], center = Xcenter, scale = sv)
  }

  Z.st <- (X.st %*% Vr)
  colnames(Z.st) <- c("x", "y")  

  # First plot to get plot dimensions of x and y axes
  biplot_plot |>
    samples(opacity = 0, which = NULL) |>
    axes(col = "grey", which = NULL)  |> plot()
  
  # {
  #   # First plot to get plot dimensions of x and y axes: only the min and max to get plot size correct
  #   if(!is(polygon, "SpatialPolygons")) {
  #     
  #     rx <- range(Z.st[, 1])
  #     ry <- range(Z.st[, 2])
  #     min_val <- c(min(rx[1]), min(ry[1]))
  #     max_val <- c(max(rx[2]), max(ry[2])) 
  #     
  #     plot(rbind(min_val, max_val), asp = 1, ylab = "Y[2]", xlab = "Y[1]")
  #   }
  #   
  #   # if have polygon, follwow these steps
  #   if(is(polygon, "SpatialPolygons")) {
  #     plot(polygon, asp = 1, ylab = "Y[2]", xlab = "Y[1]")
  #     # coord_points <- polygon@polygons[[1]]@Polygons[[1]]@coords #final_polygon@polygons[[1]]@Polygons[[1]]@coords
  #   
  #   }
  #   
  # }
  
  ## Get the grid points
  {
    # min/max values of x and y dimensions, allow square plot later
    # Get the plot dimensions
    minx <- par("usr")[1]
    miny <- par("usr")[3]
    maxx <- par("usr")[2]
    maxy <- par("usr")[4]
    dev.off()
    
    min_val <- (min(minx, miny))
    max_val <- (max(maxx, maxy))
    
    # Grid points in plot area
    yseq <- xseq <- seq(from = min_val, to = max_val, length.out = m)
    Zgrid <- as.matrix(expand.grid(xseq, yseq))
    colnames(Zgrid) <- c("x", "y")
    
  }

  # Convert back to X dimension space  
  {  
    Xgrid <- as.data.frame(Zgrid %*% tVr)
    # Return the unstandardised Xgrid
    if(standardise == TRUE) Xgrid <- sweep(Xgrid, 2, Xsd, "*")
    Xgrid <- sweep(Xgrid, 2, Xcenter, "+")
    Xgrid <- as.data.frame(Xgrid)
    # Set feature names
    colnames(Xgrid) <- varnames
  }
  
  ## Create convex hull (polygon) for the prediction area
  # use existing polygon provided unless none provided or no_polygon == TRUE
  {
    if(!is(polygon, "SpatialPolygons")) {
      # create polygon if no existing polygon was provided
      dev.new()
      polygon_points <- plothulls(
        x = Z.st[, 1], y = Z.st[, 2], asp = 1,
        xlim = c(min_val, max_val), ylim = c(min_val, max_val),
        fraction = outlie, n.hull = 1, main, add = FALSE, col.hull = 2,
        lty.hull = 1, lwd.hull = 2, density = 1
      )
      dev.off()
      polygon <- SpatialPolygons(list(Polygons(list(Polygon(polygon_points)), ID = "polygon")))
    }
    
    if(no_polygon == TRUE){
      # Use entire grid if not limiting to polygon: overwrite the existing polygon
      dev.new()
      polygon_points <- plothulls(
        x = Zgrid[, 1], y = Zgrid[, 2], asp = 1,
        xlim = c(min_val, max_val), ylim = c(min_val, max_val),
        fraction = 1, n.hull = 1, main, add = FALSE, col.hull = 2,
        lty.hull = 1, lwd.hull = 2, density = 1
      )
      dev.off()
      polygon <- SpatialPolygons(list(Polygons(list(Polygon(polygon_points)), ID = "polygon")))
    }
    
    ## use polygon to set the convex hull filter region: only apply if convex hull is true
    # Select grid points inside hull polygon
    G_data_points <- as.data.frame(Zgrid)
    coordinates(G_data_points) <- ~x + y
    G_data_points$inside_polygon <- !is.na(over(G_data_points, polygon))
    
    # Remove data points outside the polygon
    Z_data_points <- as.data.frame(Z.st)
    coordinates(Z_data_points) <- ~x + y
    Z_data_points$inside_polygon <- !is.na(over(Z_data_points, polygon))
    
    inside_grid <- G_data_points$inside_polygon == TRUE
    inside_data <- Z_data_points$inside_polygon == TRUE
    
  }
  
  # Score each grid point using the selected model
  # Score each grid point using the selected model
  # Optional chunking to reduce peak memory
  chunk <- 50000L
  grid.prob <- numeric(nrow(Xgrid))
  if (nrow(Xgrid) > 0) {
    for (i0 in seq(1L, nrow(Xgrid), by = chunk)) {
      j0 <- min(i0 + chunk - 1L, nrow(Xgrid))
      grid.prob[i0:j0] <- pred_function(
        model_use = model_use,
        model_select = model_select,
        rounding = rounding,
        new_data = Xgrid[i0:j0, , drop = FALSE]
      )
    }
  }
  
  # Color palette for grid values
  {
    col.vec <- colorRampPalette(c("deepskyblue", "white", "lightsalmon"))(100 + 1)
    col.value <- col.vec[floor(grid.prob * 100) + 1]
  }
  
  ###############
  ## Contours  ##
  ###############
  {
    grid_contour <- matrix(grid.prob, ncol = m, byrow = FALSE)
    
    if(calc_ct_in_hull == TRUE) {
      # Compute contours in hull; make edge = cutoff to avoid double grouping
      grid_contour[G_data_points$inside_polygon == FALSE] <- cutoff
      ct <- contourLines(xseq, yseq, grid_contour, levels = c(cutoff - b_margin, cutoff, cutoff + b_margin))
      ct <- ct[sapply(ct, function(cl) cl$level != cut_off)]
    }
    
    if(calc_ct_in_hull == FALSE)
      ct <- contourLines(xseq, yseq, grid_contour, levels = c(cutoff - b_margin, cutoff + b_margin))
  }
  
  ################################################
  # Remove grid points outside the polygon
  ################################################
  
  ## use convex hull to remove point, based on polygon: No polygon implies no convex hull filter even if true
  
  pred_col <- X$pred_col
  pred_use <- X$pred_use
  class <- X$class
  
  if(calc_hull == TRUE) {
    grid.prob <-  grid.prob[inside_grid]
    col.value <-  col.value[inside_grid]
    Zgrid <- Zgrid[inside_grid, ]
    Xgrid <- Xgrid[inside_grid, ]
    
    Z.st <- Z.st[inside_data, ]
    pred_col <- pred_col[inside_data]
    pred_use <- pred_use[inside_data]
    class <- class[inside_data]
  }
  
  # Retained for compatibility
  if(remove_filter == TRUE) {
    col.value <- col.value[rows_to_exclude]
    grid.prob <- grid.prob[rows_to_exclude]
    Xgrid <- Xgrid[rows_to_exclude, ]
    Zgrid <- Zgrid[rows_to_exclude, ]
  }
  
  if(return_values == TRUE)
    return(list(
      Zgrid = Zgrid, Xgrid = Xgrid, grid.prob = grid.prob, Z.st = Z.st,
      pred_col = pred_col, pred_use = pred_use, class = class,
      col.value = col.value,
      min_val = min_val, max_val = max_val, # for 3d plots
      polygon = polygon, ct = ct
    ))
}
