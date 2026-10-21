# Branch-order pipeline + the downstream table functions. Build the graph once,
# then exercise the consumers on it.

make_order_map <- function() {
  skip_if_not_installed("terra")
  skel <- make_binary_spatraster()
  # baseR distance-transform fallback exists, so imager is not required here.
  # `order` is given in the 1.1.0 spelling, so every test built on this fixture
  # also covers the alias path.
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

test_that("the edge table carries all four order columns and the 1.1.0 copies", {
  res <- make_order_map()
  et  <- res$edges
  expect_true(all(c("order_from_tip", "strahler_order",
                    "root_max_order_from_tip", "order_from_base") %in% names(et)))
  expect_equal(et$tip_order,    et$order_from_tip)
  expect_equal(et$root_order,   et$root_max_order_from_tip)
  expect_equal(et$branch_order, et$order_from_base)
})

test_that("a 1.1.0 order name is resolved to the column that is computed", {
  res <- make_order_map()                                  # order = "branch_order"
  expect_equal(res$order, "order_from_base")
  expect_equal(names(order_metrics(res$edges, order_col = "tip_order"))[1],
               "order_from_tip")
  expect_equal(attr(summarize_orders(res$edges, order_col = "root_order"), "order_col"),
               "root_max_order_from_tip")
})

test_that("an unknown order scheme is rejected by name", {
  skip_if_not_installed("terra")
  skel <- make_binary_spatraster()
  expect_error(branch_order_map(skel = skel, mask = skel, order = "soil_depth",
                                unit = "px", return_map = FALSE),
               "must be one of")
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
  om <- order_metrics(res)
  expect_s3_class(so, "data.frame")
  expect_s3_class(om, "data.frame")
  # same per-order totals (order_metrics adds length_fraction + col name differs)
  expect_equal(sort(so$total_length), sort(om$total_length), tolerance = 1e-6)
})

test_that("groups = \"thinnest\" returns that class and a rest row", {
  res <- make_order_map()
  om  <- order_metrics(res, groups = "thinnest")
  expect_true(all(om$group %in% c("thinnest", "rest")))
  expect_true("thinnest" %in% om$group)
  expect_true(is.list(attr(om, "groups")))
})

test_that("a named groups list keeps every ordered segment in exactly one row", {
  res  <- make_order_map()
  per  <- order_metrics(res)
  lv   <- per[[1]]
  skip_if(length(lv) < 2, "need at least two order classes")
  om   <- order_metrics(res, groups = list(fine = lv[1]))
  expect_equal(sort(om$group), c("fine", "rest"))
  expect_equal(sum(om$total_length), sum(per$total_length), tolerance = 1e-8)
  expect_equal(sum(om$n_segments), sum(per$n_segments))
  expect_equal(attr(om, "groups")$rest, lv[-1])
})

test_that("an order named in two groups is an error", {
  res <- make_order_map()
  expect_error(order_metrics(res, groups = list(a = 1:2, b = 2:3)),
               "more than one group")
})

test_that("an unnamed groups list and a stray group name are rejected", {
  res <- make_order_map()
  expect_error(order_metrics(res, groups = list(1:2)), "named list")
  expect_error(order_metrics(res, groups = list(rest = 1)), "another name")
  expect_error(order_metrics(res, groups = "fattest"), "thinnest")
})

test_that("a group that catches no segment gets no row", {
  res <- make_order_map()
  om  <- order_metrics(res, groups = list(absent = 999, present = 1))
  expect_false("absent" %in% om$group)
  expect_true("present" %in% om$group)
})

test_that("focal still works and warns", {
  res <- make_order_map()
  expect_warning(om <- order_metrics(res, focal = 1), "deprecated")
  expect_true("focal" %in% om$group)
  expect_warning(ot <- order_metrics(res, focal = "thinnest"), "deprecated")
  expect_true("thinnest" %in% ot$group)
})

test_that("cumulative length runs over the orders present, with no zero rows", {
  res <- make_order_map()
  om  <- order_metrics(res)
  expect_true(all(c("cumulative_length", "cumulative_fraction") %in% names(om)))
  expect_equal(om[[1]], sort(om[[1]]))                      # ascending
  expect_equal(om$cumulative_length, cumsum(om$total_length), tolerance = 1e-8)
  expect_equal(utils::tail(om$cumulative_fraction, 1), 1, tolerance = 1e-8)
  expect_true(all(om$total_length > 0))                     # no synthetic zeros
  expect_false(any(c("cumulative_length", "cumulative_fraction") %in%
                     names(order_metrics(res, groups = list(all = om[[1]])))))
})

test_that("order_classification_map rasterises onto the template", {
  res <- make_order_map()
  cm  <- order_classification_map(res$edges, template = make_binary_spatraster(),
                                  value = "branch_order")
  expect_s4_class(cm, "SpatRaster")
  # The layer is named after the resolved column, not the spelling passed in.
  expect_equal(names(cm), "order_from_base")
})

test_that("summarize_orders handles NULL edge table without erroring", {
  expect_s3_class(summarize_orders(NULL), "data.frame")
})
