# .poly_clip() keeps the rows of a 2-column matrix that lie inside a
# SpatialPolygons hull. bl_find_boundary() clips the global contours with it,
# and bl_find_local_cf() the contours on each rotated plane.

unit_square <- function() {
  sp::SpatialPolygons(list(sp::Polygons(
    list(sp::Polygon(cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)))),
    ID = "polygon"
  )))
}

xy <- function(...) {
  M <- rbind(...)
  colnames(M) <- c("x", "y")
  M
}

test_that(".poly_clip keeps inside, edge and vertex points and drops outside ones", {
  M   <- xy(c(0.5, 0.5), c(2, 2), c(1, 0.5), c(0, 0), c(-0.1, 0.5))
  out <- boundarylogic:::.poly_clip(M, unit_square())
  expect_equal(out, M[c(1L, 3L, 4L), , drop = FALSE])
  expect_identical(colnames(out), c("x", "y"))
})

test_that(".poly_clip keeps a single kept row as a 1 x 2 matrix", {
  out <- boundarylogic:::.poly_clip(xy(c(0.5, 0.5), c(2, 2)), unit_square())
  expect_true(is.matrix(out))
  expect_identical(dim(out), c(1L, 2L))
})

test_that(".poly_clip returns a 0 x 2 matrix when every point is outside", {
  out <- boundarylogic:::.poly_clip(xy(c(2, 2), c(3, 3)), unit_square())
  expect_identical(dim(out), c(0L, 2L))
})

test_that(".poly_clip returns a 0 x 2 matrix for empty input", {
  expect_identical(dim(boundarylogic:::.poly_clip(NULL, unit_square())),
                   c(0L, 2L))
  expect_identical(
    dim(boundarylogic:::.poly_clip(matrix(numeric(0L), ncol = 2L),
                                   unit_square())),
    c(0L, 2L)
  )
})

test_that(".poly_clip returns the points unchanged when poly is not a polygon", {
  M <- xy(c(0.5, 0.5), c(2, 2), c(1, 0.5))
  expect_identical(boundarylogic:::.poly_clip(M, NULL), M)
})

test_that(".grid_in_hull tests points against stored hull vertices", {
  sq <- cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0))
  Z  <- rbind(c(0.5, 0.5), c(2, 2), c(0.2, 0.9), c(-0.1, 0.5))
  expect_identical(boundarylogic:::.grid_in_hull(Z, sq),
                   c(TRUE, FALSE, TRUE, FALSE))
  expect_identical(boundarylogic:::.grid_in_hull(Z, NULL), rep(TRUE, 4L))
})
