# Branch-order pipeline + the downstream table functions. Build the graph once,
# then exercise the consumers on it.

make_order_map <- function() {
  skip_if_not_installed("terra")
  skel <- make_binary_spatraster()
  # baseR distance-transform fallback exists, so imager is not required here.
  branch_order_map(skel = skel, mask = skel, order = "branch_order",
                   unit = "px", return_map = TRUE)
}

test_that("branch_order_map returns the documented object", {
  res <- make_order_map()
  expect_s3_class(res, "branchOrderMap")
  expect_true(all(c("edges", "summary", "order", "unit") %in% names(res)))
  expect_s3_class(res$edges, "data.frame")
  expect_true(all(c("tip_order", "root_order", "branch_order") %in% names(res$edges)))
})

test_that("convert_root_units rescales lengths and records attributes", {
  res <- make_order_map()
  cm  <- convert_root_units(res$edges, unit = "cm", dpi = 300)
  expect_equal(attr(cm, "unit"), "cm")
  expect_equal(cm$length, res$edges$length * 2.54 / 300, tolerance = 1e-6)
})

test_that("summarize_orders and order_metrics agree on shared columns", {
  res <- make_order_map()
  so <- summarize_orders(res$edges, order_col = "branch_order")
  om <- order_metrics(res, focal = NULL)
  expect_s3_class(so, "data.frame")
  expect_s3_class(om, "data.frame")
  # same per-order totals (order_metrics adds length_fraction + col name differs)
  expect_equal(sort(so$total_length), sort(om$total_length), tolerance = 1e-6)
})

test_that("order_metrics focal split returns focal + rest", {
  res <- make_order_map()
  om <- order_metrics(res, focal = "thinnest")
  expect_true(all(om$group %in% c("focal", "rest")))
})

test_that("order_classification_map rasterises onto the template", {
  res <- make_order_map()
  cm  <- order_classification_map(res$edges, template = make_binary_spatraster(),
                                  value = "branch_order")
  expect_s4_class(cm, "SpatRaster")
})

test_that("summarize_orders handles NULL edge table without erroring", {
  expect_s3_class(summarize_orders(NULL), "data.frame")
})
