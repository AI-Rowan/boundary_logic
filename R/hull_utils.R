############################################################
# Internal convex hull / polygon helpers
# Refactored from: scripts/1.1 biplot_plane_optimized.R
# Isolates the aplpack::plothulls + sp::SpatialPolygons side-effect pattern.
############################################################

#' Build a convex hull polygon from projected 2D coordinates
#'
#' Uses `aplpack::plothulls()` to compute the hull boundary and wraps
#' the result as an `sp::SpatialPolygons` object. The graphics device
#' is opened and closed internally; no plot is visible to the user.
#'
#' @param Z       Numeric matrix with exactly 2 columns (x, y). Typically
#'   the projected training data (`Z_train`).
#' @param min_val Numeric; lower bound for the square plot region passed to
#'   `plothulls()`.
#' @param max_val Numeric; upper bound for the square plot region.
#' @param outlie  Numeric in (0, 1]; the `fraction` argument for
#'   `aplpack::plothulls()`. Values below 1 trim the most extreme points.
#'   Default `0.9`.
#'
#' @return An `sp::SpatialPolygons` object representing the convex hull.
#' @keywords internal
.build_hull_polygon <- function(Z, min_val, max_val, outlie = 0.9) {
  # Suppress the graphics output from plothulls using an off-screen device
  tmp_file <- tempfile(fileext = ".png")
  grDevices::png(tmp_file, width = 600L, height = 600L)
  # Minimise margins so plot.new() does not fail on small devices
  graphics::par(mar = c(1, 1, 1, 1))
  on.exit({
    grDevices::dev.off()
    unlink(tmp_file)
  }, add = TRUE)

  polygon_points <- aplpack::plothulls(
    x        = Z[, 1],
    y        = Z[, 2],
    asp      = 1,
    xlim     = c(min_val, max_val),
    ylim     = c(min_val, max_val),
    fraction = outlie,
    n.hull   = 1L,
    add      = FALSE,
    col.hull = 2L,
    lty.hull = 1L,
    lwd.hull = 2L,
    density  = 1L
  )

  sp::SpatialPolygons(
    list(sp::Polygons(list(sp::Polygon(polygon_points)), ID = "polygon"))
  )
}


#' Filter a matrix of (x, y) points to those inside a SpatialPolygons object
#'
#' @param Mxy  Numeric matrix with 2 columns (x, y).
#' @param poly An `sp::SpatialPolygons` object.
#'
#' @return The subset of rows in `Mxy` that lie inside `poly`. If `poly`
#'   is not a `SpatialPolygons` object, `Mxy` is returned unchanged.
#' @keywords internal
.poly_clip <- function(Mxy, poly) {
  if (is.null(Mxy) || nrow(Mxy) == 0L) return(matrix(numeric(0L), ncol = 2L))
  if (!inherits(poly, "SpatialPolygons"))  return(as.matrix(Mxy))

  coords  <- poly@polygons[[1L]]@Polygons[[1L]]@coords
  inside  <- sp::point.in.polygon(Mxy[, 1L], Mxy[, 2L],
                                  coords[, 1L], coords[, 2L])
  Mxy[inside > 0L, , drop = FALSE]
}


#' Test whether a set of points lies inside a SpatialPolygons hull
#'
#' @param Z    n x 2 matrix of (x, y) coordinates.
#' @param poly An `sp::SpatialPolygons` object.
#'
#' @return Logical vector of length n.
#' @keywords internal
.points_in_polygon <- function(Z, poly) {
  if (!inherits(poly, "SpatialPolygons")) return(rep(TRUE, nrow(Z)))
  Z_sp <- as.data.frame(Z)
  colnames(Z_sp) <- c("x", "y")
  sp::coordinates(Z_sp) <- ~x + y
  !is.na(sp::over(Z_sp, poly))
}
