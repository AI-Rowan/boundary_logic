# The internal helpers .build_hull_polygon() and .bl_get_plot_bounds() each
# draw on a temporary png() device. They must close only that device and leave
# the caller's active device (and the list of open devices) unchanged.
# pdf(NULL) devices only, so the tests run without a screen.

hull_Z <- function() {
  set.seed(1)
  matrix(rnorm(200), ncol = 2)
}

dev_bl_dat <- bl_prepare_data(iris, class_col = "Species",
                              target_class = "versicolor", verbose = FALSE)
dev_bl_mod <- bl_fit_model(dev_bl_dat$train_data, dev_bl_dat$var_names,
                           model_type = "GLM")
dev_bl_proj <- bl_build_projection(dev_bl_dat$train_data, dev_bl_dat$var_names,
                                   method = "PCA")

open_two_devices <- function() {
  grDevices::graphics.off()
  grDevices::pdf(NULL); d1 <- grDevices::dev.cur()
  grDevices::pdf(NULL); d2 <- grDevices::dev.cur()
  c(d1, d2)
}

test_that(".build_hull_polygon keeps the caller's device with two open", {
  on.exit(grDevices::graphics.off(), add = TRUE)
  devs <- open_two_devices()
  boundarylogic:::.build_hull_polygon(hull_Z(), -4, 4, outlie = 0.9)
  expect_identical(grDevices::dev.cur(), devs[2L])
  expect_identical(unname(grDevices::dev.list()), unname(devs))
})

test_that(".bl_get_plot_bounds keeps the caller's device with two open", {
  on.exit(grDevices::graphics.off(), add = TRUE)
  Z_train <- hull_Z()
  devs <- open_two_devices()
  boundarylogic:::.bl_get_plot_bounds(dev_bl_proj$biplot_obj, Z_train)
  expect_identical(grDevices::dev.cur(), devs[2L])
  expect_identical(unname(grDevices::dev.list()), unname(devs))
})

test_that("the helpers leave no device open when none was open", {
  on.exit(grDevices::graphics.off(), add = TRUE)
  grDevices::graphics.off()
  boundarylogic:::.build_hull_polygon(hull_Z(), -4, 4, outlie = 0.9)
  expect_null(grDevices::dev.list())
  boundarylogic:::.bl_get_plot_bounds(dev_bl_proj$biplot_obj, hull_Z())
  expect_null(grDevices::dev.list())
})

test_that("the hull polygon is the same whether or not other devices are open", {
  on.exit(grDevices::graphics.off(), add = TRUE)
  grDevices::graphics.off()
  hull_none <- boundarylogic:::.build_hull_polygon(hull_Z(), -4, 4, outlie = 0.9)
  open_two_devices()
  hull_two <- boundarylogic:::.build_hull_polygon(hull_Z(), -4, 4, outlie = 0.9)
  expect_identical(sp::coordinates(hull_none@polygons[[1L]]@Polygons[[1L]]),
                   sp::coordinates(hull_two@polygons[[1L]]@Polygons[[1L]]))
})

test_that("bl_prepare_data and bl_build_result keep the caller's device", {
  on.exit(grDevices::graphics.off(), add = TRUE)
  devs <- open_two_devices()

  bl_dat <- bl_prepare_data(iris, class_col = "Species",
                            target_class = "versicolor", verbose = FALSE)
  expect_identical(grDevices::dev.cur(), devs[2L])
  expect_identical(unname(grDevices::dev.list()), unname(devs))

  bl_results <- suppressWarnings(
    bl_build_result(bl_dat, dev_bl_mod, method = "PCA", m = 30L))
  expect_s3_class(bl_results, "bl_result")
  expect_identical(grDevices::dev.cur(), devs[2L])
  expect_identical(unname(grDevices::dev.list()), unname(devs))
})
