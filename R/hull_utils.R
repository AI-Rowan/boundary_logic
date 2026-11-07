############################################################
# Internal convex hull / polygon helpers
# Refactored from: scripts/1.1 biplot_plane_optimized.R
# Isolates the aplpack::plothulls + sp::SpatialPolygons side-effect pattern.
############################################################

#' Build a convex hull polygon from projected 2D coordinates
#'
#' Uses `aplpack::plothulls()` to compute the hull boundary and wraps
#' the result as an `sp::SpatialPolygons` object. A temporary off-screen
#' device is opened and closed internally; no plot is visible to the user,
#' and the caller's active graphics device is restored on exit.
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
#' @importFrom aplpack plothulls
#' @importFrom sp Polygon Polygons SpatialPolygons coordinates over point.in.polygon
#' @importFrom grDevices png dev.off dev.cur dev.set
#' @importFrom graphics par
#' @return An `sp::SpatialPolygons` object representing the convex hull.
#' @keywords internal
.build_hull_polygon <- function(Z, min_val, max_val, outlie = 0.9) {
  # Suppress the graphics output from plothulls using an off-screen device.
  # Close only that device and re-activate the caller's one: a bare dev.off()
  # would leave R's *next* open device active, not the caller's.
  prev_dev <- grDevices::dev.cur()
  tmp_file <- tempfile(fileext = ".png")
  grDevices::png(tmp_file, width = 600L, height = 600L)
  hull_dev <- grDevices::dev.cur()
  # Minimise margins so plot.new() does not fail on small devices
  graphics::par(mar = c(1, 1, 1, 1))
  on.exit({
    grDevices::dev.off(hull_dev)
    # dev.cur() is 1 when no device was open; dev.set(1) would open a new one
    if (prev_dev > 1L) grDevices::dev.set(prev_dev)
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


#' Test whether points lie inside a hull given by its vertices
#'
#' Rebuilds the `SpatialPolygons` object from stored hull vertices (as
#' `.build_hull_polygon()` builds it) and applies `.points_in_polygon()`, the
#' rule the global grid uses to keep its points inside the hull. Used by
#' `plot.bl_local_result()` to colour the local grid only inside the hull of
#' its plane (`bl_local_result$hull_polygon`).
#'
#' @param Z           n x 2 matrix of (x, y) coordinates.
#' @param hull_coords Two-column matrix of hull vertices, or `NULL`.
#'
#' @return Logical vector of length n; all `TRUE` when `hull_coords` is `NULL`.
#' @keywords internal
.grid_in_hull <- function(Z, hull_coords) {
  Z <- as.matrix(Z)
  if (is.null(hull_coords)) return(rep(TRUE, nrow(Z)))
  poly <- sp::SpatialPolygons(
    list(sp::Polygons(list(sp::Polygon(as.matrix(hull_coords))), ID = "polygon"))
  )
  unname(.points_in_polygon(Z, poly))   # sp::over() names its result
}


#' Does a boundary contour enclose training data?
#'
#' The empty-region rule for decision boundaries (thesis section 3.10, step 4:
#' a boundary area that contains no data is discarded). A contour that is
#' still closed after clipping to the hull (first vertex equal to the last,
#' `.is_closed()`) is kept only if at least one training point lies inside or
#' on it; an open contour is always kept. Used by `bl_find_boundary()` (the
#' training data on the global plane) and `bl_find_local_cf()` (the training
#' data on each rotated plane), so both searches apply the same rule.
#'
#' @param Mi    Two-column matrix of contour vertices (x, y).
#' @param Z_ref Two-column matrix of the training points on the same plane.
#'
#' @return Logical scalar: `TRUE` to keep the contour.
#' @keywords internal
.encloses_data <- function(Mi, Z_ref) {
  if (!.is_closed(Mi)) return(TRUE)
  inside <- sp::point.in.polygon(Z_ref[, 1L], Z_ref[, 2L], Mi[, 1L], Mi[, 2L])
  any(inside > 0L)
}
