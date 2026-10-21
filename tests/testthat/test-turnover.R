# Turnover is scored on masks whose pixel counts are known by construction, so
# the expected production and ratios are arithmetic rather than recorded output.

skip_if_not_installed("terra")

test_that("turnover_tc counts root pixels and reports production against both dates", {
  m1 <- matrix(0, 10, 10); m1[3, 2:6] <- 1      # 5 root pixels
  m2 <- matrix(0, 10, 10); m2[3, 2:9] <- 1      # 8 root pixels
  out <- turnover_tc(terra::rast(m1), terra::rast(m2),
                     method = "rootpx", select_layer = NULL)

  expect_equal(out$standingroot_t1, 5)
  expect_equal(out$standingroot_t2, 8)
  expect_equal(out$production, 3)
  # The two percentages take the same production over each date's standing root.
  # data.frame() makes the names syntactic, so the percent sign reads as a dot.
  expect_equal(out[["newroot.per_t1"]], round(3 / 5, 4))
  expect_equal(out[["newroot.per_t2"]], round(3 / 8, 4))
  expect_named(out, c("standingroot_t1", "standingroot_t2", "production",
                      "newroot.per_t1", "newroot.per_t2"))
})

test_that("turnover_tc with kimura measures length, so a 5-pixel line is 4 steps", {
  m1 <- matrix(0, 10, 10); m1[3, 2:6] <- 1      # 5 pixels in a row
  m2 <- matrix(0, 10, 10); m2[3, 2:9] <- 1      # 8 pixels in a row
  # At 254 dpi one pixel edge is 2.54/254 = 0.01 cm. A run of n pixels has
  # n - 1 steps between them, all orthogonal, so the length is (n - 1) * 0.01.
  out <- turnover_tc(terra::rast(m1), terra::rast(m2),
                     method = "kimura", select_layer = NULL, dpi = 254, unit = "cm")

  expect_equal(out$standingroot_t1, 0.04)
  expect_equal(out$standingroot_t2, 0.07)
  expect_equal(out$production, 0.03, tolerance = 1e-12)
})

test_that("include_virtualroots changes the denominator, not the pixel counts", {
  set.seed(5)
  mk <- function(p) terra::rast(matrix(sample(c(0, 255), 400, TRUE, c(1 - p, p)), 20, 20))
  img <- c(mk(0.3), mk(0.4), mk(0.5))

  a <- turnover_dpc(img, include_virtualroots = FALSE)
  b <- turnover_dpc(img, include_virtualroots = TRUE)

  # The classification is identical; only how the ratios are normalized moves.
  expect_equal(a[, c("tape", "constant", "production", "decay")],
               b[, c("tape", "constant", "production", "decay")])

  den_all <- a$production + a$constant + a$decay
  expect_equal(a$newgrowth.ratio, round(a$production / (a$production + a$constant), 4))
  expect_equal(b$newgrowth.ratio, round(b$production / den_all, 4))
  expect_equal(a$decay.ratio,     round(a$decay / (a$decay + a$constant), 4))
  expect_equal(b$decay.ratio,     round(b$decay / den_all, 4))

  # constant.ratio already divides by every class, so it does not move at all.
  expect_equal(a$constant.ratio, b$constant.ratio)
  # And with decay present, counting it in the denominator can only lower the
  # reported new-growth share.
  expect_gt(a$decay, 0)
  expect_lt(b$newgrowth.ratio, a$newgrowth.ratio)
})

test_that("turnover_dpc can return the classified layers instead of the counts", {
  set.seed(5)
  mk <- function(p) terra::rast(matrix(sample(c(0, 255), 400, TRUE, c(1 - p, p)), 20, 20))
  img <- c(mk(0.3), mk(0.4), mk(0.5))

  ims <- turnover_dpc(img, im_return = TRUE)
  expect_named(ims, c("tape", "constant", "production", "decay"))
  expect_true(all(vapply(ims, inherits, logical(1), "SpatRaster")))
  # The rasters must be the ones the counts were taken from, so each layer's
  # foreground count has to match the tabular run on the same image.
  cnt <- turnover_dpc(img, im_return = FALSE)
  expect_equal(count_pixels(ims$production), cnt$production)
  expect_equal(count_pixels(ims$decay),      cnt$decay)
  expect_equal(count_pixels(ims$constant),   cnt$constant)
})

test_that("root_turnover dispatches to the two methods and guards the inputs", {
  m1 <- matrix(0, 10, 10); m1[3, 2:6] <- 1
  m2 <- matrix(0, 10, 10); m2[3, 2:9] <- 1
  r1 <- terra::rast(m1); r2 <- terra::rast(m2)

  expect_equal(root_turnover(r1, r2, method = "tc", tc_method = "rootpx")$production, 3)
  expect_error(root_turnover(r1, method = "tc"), "requires two images")

  set.seed(5)
  mk <- function(p) terra::rast(matrix(sample(c(0, 255), 400, TRUE, c(1 - p, p)), 20, 20))
  img <- c(mk(0.3), mk(0.4), mk(0.5))
  expect_warning(root_turnover(img, img, method = "dpc"), "img2 is ignored")
})
