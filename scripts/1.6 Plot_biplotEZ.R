##############################################
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
