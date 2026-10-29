# .make_label_line_vec() and .make_ticks_vec() pair each value with its own
# variable, position by position. An unknown variable name is dropped together
# with its value; a single value applies to every listed variable; any other
# length mismatch is recycled with a warning.

ph_vars <- c("a", "b", "c")

test_that("label offsets land on their own variables when all names are valid", {
  expect_no_warning(
    vec <- boundarylogic:::.make_label_line_vec(c("a", "c"), c(1, 2), 3L, ph_vars))
  expect_equal(vec, c(1, 0.1, 2))
})

test_that("an unknown label name is dropped together with its offset", {
  expect_warning(
    vec <- boundarylogic:::.make_label_line_vec(c("zz", "b", "c"), c(9, 1, 2),
                                                3L, ph_vars),
    "not found: zz")
  expect_equal(vec, c(0.1, 1, 2))
})

test_that("a single label offset applies to every listed variable", {
  expect_no_warning(
    vec <- boundarylogic:::.make_label_line_vec(c("a", "b"), 1.5, 3L, ph_vars))
  expect_equal(vec, c(1.5, 1.5, 0.1))
})

test_that("a label offset length mismatch warns and recycles", {
  expect_warning(
    vec <- boundarylogic:::.make_label_line_vec(c("a", "b"), c(1, 2, 3),
                                                3L, ph_vars),
    "values are recycled")
  expect_equal(vec, c(1, 2, 0.1))
})

test_that("the label offset defaults give no warning and no offset", {
  expect_no_warning(
    vec <- boundarylogic:::.make_label_line_vec(0L, 1.5, 3L, ph_vars))
  expect_equal(vec, rep(0.1, 3L))
})

test_that("tick counts land on their own variables when all names are valid", {
  expect_no_warning(
    vec <- boundarylogic:::.make_ticks_vec(1L, c("a", "c"), c(4L, 6L),
                                           3L, ph_vars))
  expect_equal(vec, c(4L, 1L, 6L))
})

test_that("an unknown tick name is dropped together with its count", {
  expect_warning(
    vec <- boundarylogic:::.make_ticks_vec(1L, c("zz", "b", "c"), c(9L, 4L, 6L),
                                           3L, ph_vars),
    "not found: zz")
  expect_equal(vec, c(1L, 4L, 6L))
})

test_that("a single tick count applies to every listed variable", {
  expect_no_warning(
    vec <- boundarylogic:::.make_ticks_vec(1L, c("a", "b"), 5L, 3L, ph_vars))
  expect_equal(vec, c(5L, 5L, 1L))
})

test_that("a tick count length mismatch warns and recycles", {
  expect_warning(
    vec <- boundarylogic:::.make_ticks_vec(1L, c("a", "b"), c(4L, 5L, 6L),
                                           3L, ph_vars),
    "values are recycled")
  expect_equal(vec, c(4L, 5L, 1L))
})
