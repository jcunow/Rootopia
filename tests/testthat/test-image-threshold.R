# image_threshold, adaptive mode.
#
# The scene is the one adaptive thresholding exists for: bright roots on a
# background that is almost three times brighter at the top-left corner than at
# the bottom-right, so no single global level is right everywhere. Roots are
# 1 px wide and the window is 21, which satisfies "window much wider than the
# root"; the truth mask is known, so precision and recall are exact rather than
# eyeballed.
adaptive_scene <- function(nr = 80L, nc = 80L, contrast = 2, cv = 0.02) {
  set.seed(7)
  g <- outer(seq(1, 0.35, length.out = nr), seq(1, 0.6, length.out = nc))
  truth <- matrix(FALSE, nr, nc)
  truth[, c(12, 30, 55, 70)] <- TRUE
  truth[c(20, 48), ] <- TRUE
  m <- g * 0.4
  m[truth] <- m[truth] * contrast
  m <- m + stats::rnorm(nr * nc, 0, cv * 0.4)
  m[] <- pmax(0, pmin(1, m))
  list(img = terra::rast(m), truth = terra::rast(truth * 1))
}

# Scores are taken from the two rasters' own cell vectors, so the comparison
# cannot be thrown off by a matrix orientation.
pr_scores <- function(pred, truth) {
  p <- terra::values(pred) > 0
  t <- terra::values(truth) > 0
  tp <- sum(p & t); fp <- sum(p & !t); fn <- sum(!p & t)
  list(f1 = 2 * tp / (2 * tp + fp + fn), precision = tp / (tp + fp),
       recall = tp / (tp + fn), frac_foreground = mean(p))
}

test_that("adaptive thresholding recovers a known mask under a gradient", {
  skip_if_not_installed("terra")
  s <- adaptive_scene()
  out <- image_threshold(s$img, method = "adaptive", adaptive_k = 1.15,
                         window_size = 21, select_layer = 1, mask_layer = NULL,
                         binary_01 = TRUE)
  sc <- pr_scores(out, s$truth)
  expect_gt(sc$f1, 0.9)
  expect_gt(sc$precision, 0.9)
  expect_gt(sc$recall, 0.9)
})

test_that("a multiplier at or below 1 cannot separate anything", {
  skip_if_not_installed("terra")
  s <- adaptive_scene()
  truth_frac <- mean(terra::values(s$truth) > 0)

  # Below 1 every pixel clears its own local mean, so the whole image is
  # foreground and precision falls to the base rate.
  low <- suppressWarnings(
    image_threshold(s$img, method = "adaptive", adaptive_k = 0.4,
                    window_size = 21, select_layer = 1, mask_layer = NULL,
                    binary_01 = TRUE))
  sc <- pr_scores(low, s$truth)
  expect_gt(sc$frac_foreground, 0.95)
  expect_lt(sc$precision, 2 * truth_frac)

  # At exactly 1 the local mean sits in the middle of the background's spread.
  # On a scene with roots in it the roots pull their own neighborhood's mean up,
  # which rescues some of the nearby background, so the over-detection here is
  # about twofold rather than the near-half of a root-free patch (covered by
  # its own test below). Precision is what collapses either way.
  at1 <- suppressWarnings(
    image_threshold(s$img, method = "adaptive", adaptive_k = 1,
                    window_size = 21, select_layer = 1, mask_layer = NULL,
                    binary_01 = TRUE))
  sc1 <- pr_scores(at1, s$truth)
  expect_gt(sc1$frac_foreground, 1.5 * truth_frac)
  expect_lt(sc1$precision, 0.7)
})

test_that("a root-free image yields almost nothing above 1 and half of itself at 1", {
  skip_if_not_installed("terra")
  set.seed(3)
  flat <- terra::rast(matrix(0.4 + stats::rnorm(60 * 60, 0, 0.02 * 0.4), 60, 60))
  frac <- function(k) {
    o <- suppressWarnings(
      image_threshold(flat, method = "adaptive", adaptive_k = k,
                      window_size = 21, select_layer = 1, mask_layer = NULL,
                      binary_01 = TRUE))
    mean(terra::values(o) > 0)
  }
  expect_gt(frac(1.0), 0.4)
  expect_lt(frac(1.15), 0.01)
})

test_that("adaptive beats the best global level when the illumination is uneven", {
  skip_if_not_installed("terra")
  s <- adaptive_scene()
  ad <- image_threshold(s$img, method = "adaptive", adaptive_k = 1.15,
                        window_size = 21, select_layer = 1, mask_layer = NULL,
                        binary_01 = TRUE)
  best_global <- max(vapply(seq(0.05, 0.95, by = 0.025), function(t) {
    o <- image_threshold(s$img, threshold = t, method = "global",
                         select_layer = 1, mask_layer = NULL, binary_01 = TRUE)
    f <- pr_scores(o, s$truth)$f1
    if (is.finite(f)) f else 0
  }, numeric(1)))
  expect_gt(pr_scores(ad, s$truth)$f1, best_global)
})

test_that("the local mean is the separable box mean terra::focal gives", {
  skip_if_not_installed("terra")
  set.seed(5)
  r <- terra::rast(matrix(stats::runif(40 * 36), 40, 36))
  direct <- terra::focal(r, w = matrix(1, 15, 15), fun = "mean", na.rm = TRUE)
  expect_equal(terra::values(.box_mean(r, 15L)), terra::values(direct),
               tolerance = 1e-10)
})

test_that("adaptive arguments are checked and a stale threshold is reported", {
  skip_if_not_installed("terra")
  s <- adaptive_scene()
  common <- list(s$img, method = "adaptive", window_size = 21,
                 select_layer = 1, mask_layer = NULL, binary_01 = TRUE)

  expect_warning(do.call(image_threshold, c(common, list(threshold = 0.4))),
                 "adaptive_k")
  expect_warning(do.call(image_threshold, c(common, list(adaptive_k = 0.9))),
                 "at or below 1")
  expect_error(do.call(image_threshold, c(common, list(adaptive_k = -1))),
               "adaptive_k")
  expect_error(do.call(image_threshold, c(common, list(adaptive_k = c(1.1, 1.2)))),
               "adaptive_k")
  expect_error(image_threshold(s$img, method = "adaptive", window_size = 20,
                               select_layer = 1, mask_layer = NULL),
               "odd")
})

test_that("a window wider than the image is clamped, not passed to focal", {
  skip_if_not_installed("terra")
  set.seed(9)
  small <- terra::rast(matrix(stats::runif(24 * 28), 24, 28))
  expect_warning(
    out <- image_threshold(small, method = "adaptive", adaptive_k = 1.15,
                           window_size = 101, select_layer = 1,
                           mask_layer = NULL, binary_01 = TRUE),
    "exceeds the image")
  expect_s4_class(out, "SpatRaster")
  expect_true(all(terra::values(out) %in% c(0, 1)))
})

test_that("global thresholding is unchanged and both methods return 0/1", {
  skip_if_not_installed("terra")
  set.seed(1)
  r <- terra::rast(matrix(stats::runif(24 * 28), 24, 28))
  g <- image_threshold(r, threshold = 0.4, method = "global",
                       select_layer = 1, mask_layer = NULL, binary_01 = TRUE)
  expect_equal(terra::values(g),
               terra::values((r >= 0.4 * max(terra::values(r))) * 1))

  a <- image_threshold(r, method = "adaptive", adaptive_k = 1.15,
                       window_size = 21, select_layer = 1, mask_layer = NULL,
                       binary_01 = TRUE)
  for (out in list(g, a)) {
    expect_s4_class(out, "SpatRaster")
    v <- terra::values(out)
    expect_true(all(v %in% c(0, 1) | is.na(v)))
  }
})

test_that("on a multilayer image the mask layer gets its own local level", {
  skip_if_not_installed("terra")
  set.seed(4)
  r <- c(terra::rast(matrix(stats::runif(60 * 60), 60, 60)),
         terra::rast(matrix(stats::runif(60 * 60), 60, 60)),
         terra::rast(matrix(stats::runif(60 * 60), 60, 60)))
  out <- image_threshold(r, method = "adaptive", adaptive_k = 1.15,
                         window_size = 21, select_layer = 2, mask_layer = 1,
                         binary_01 = TRUE)
  expect_equal(terra::nlyr(out), 3L)

  # layer 1 is judged against layer 1's own neighborhood
  expect_equal(terra::values(out)[, 1],
               terra::values((r[[1]] >= .box_mean(r[[1]], 21L) * 1.15) * 1)[, 1])
  # every other layer carries the mask built from select_layer
  expect_equal(terra::values(out)[, 2],
               terra::values((r[[2]] >= .box_mean(r[[2]], 21L) * 1.15) * 1)[, 1])
  expect_equal(terra::values(out)[, 2], terra::values(out)[, 3])
})

test_that("an unusable method stops instead of falling through", {
  skip_if_not_installed("terra")
  r <- terra::rast(matrix(stats::runif(30 * 30), 30, 30))
  expect_error(image_threshold(r, method = "nope", select_layer = 1,
                               mask_layer = NULL), "global")
  expect_error(image_threshold(r, method = NULL, select_layer = 1,
                               mask_layer = NULL), "global")
})

test_that("the documented example recovers its three roots exactly", {
  skip_if_not_installed("terra")
  # the scene from ?image_threshold: background falling 0.40 -> 0.11 across the
  # frame, three 1 px roots at twice their own local background
  m <- outer(seq(1, 0.4, length.out = 120), seq(1, 0.7, length.out = 120)) * 0.4
  expect_equal(round(c(max(m), min(m)), 2), c(0.40, 0.11))
  m[, c(20, 50, 90)] <- m[, c(20, 50, 90)] * 2

  bw <- image_threshold(terra::rast(m), method = "adaptive", adaptive_k = 1.15,
                        window_size = 31, binary_01 = TRUE)
  truth <- matrix(FALSE, 120, 120)
  truth[, c(20, 50, 90)] <- TRUE
  sc <- pr_scores(bw, terra::rast(truth * 1))
  expect_equal(sc$f1, 1)

  # and that the gradient is genuinely beyond a single global level
  best_global <- max(vapply(seq(0.05, 0.95, by = 0.025), function(t) {
    f <- pr_scores(image_threshold(terra::rast(m), threshold = t,
                                   method = "global", binary_01 = TRUE),
                   terra::rast(truth * 1))$f1
    if (is.finite(f)) f else 0
  }, numeric(1)))
  expect_lt(best_global, 0.9)
})


test_that("deblur separates the two layers at a fraction of select_layer's max", {
  skip_if_not_installed("terra")
  # Both output layers are worked out by hand from the documented rule: the
  # mask layer becomes (mask - select >= max(select) * threshold) * max(select),
  # and the select layer becomes (select >= max(select) * threshold) * max(select).
  l_mask   <- matrix(c(200, 200, 100, 0), 2, 2, byrow = TRUE)
  l_select <- matrix(c(100,   0,   0, 0), 2, 2, byrow = TRUE)
  img <- c(terra::rast(l_mask), terra::rast(l_select))

  out <- image_threshold(img, threshold = 0.5, select_layer = 2,
                         mask_layer = 1, deblur = TRUE)

  mx  <- max(l_select)                       # 100
  thr <- mx * 0.5                            # 50
  expect_equal(as.vector(terra::values(out[[1]])),
               as.vector(t((l_mask - l_select >= thr) * mx)))
  expect_equal(as.vector(terra::values(out[[2]])),
               as.vector(t((l_select >= thr) * mx)))
})

test_that("deblur needs two different layers to subtract", {
  skip_if_not_installed("terra")
  one <- terra::rast(matrix(1, 2, 2))
  expect_error(image_threshold(one, deblur = TRUE), "multi-layer")
  two <- c(terra::rast(matrix(1, 2, 2)), terra::rast(matrix(1, 2, 2)))
  expect_error(image_threshold(two, select_layer = 1, mask_layer = 1, deblur = TRUE),
               "cannot be the same")
})

# invert = TRUE is the dark-roots-on-a-bright-tray case: the same gradient
# background, with the roots at half their own local level instead of twice it.
# The local cut is the mean divided by adaptive_k, so a pixel has to be that
# factor darker than its surroundings, which is the mirror of what the upright
# scene asks of a bright root.
dark_adaptive_scene <- function(nr = 80L, nc = 80L, contrast = 2, cv = 0.02) {
  set.seed(7)
  g <- outer(seq(1, 0.35, length.out = nr), seq(1, 0.6, length.out = nc))
  truth <- matrix(FALSE, nr, nc)
  truth[, c(12, 30, 55, 70)] <- TRUE
  truth[c(20, 48), ] <- TRUE
  m <- g * 0.8
  m[truth] <- m[truth] / contrast
  m <- m + stats::rnorm(nr * nc, 0, cv * 0.4)
  m[] <- pmax(0, pmin(1, m))
  list(img = terra::rast(m), truth = terra::rast(truth * 1))
}

test_that("invert recovers a known mask of dark roots under a gradient", {
  skip_if_not_installed("terra")
  s <- dark_adaptive_scene()

  dn <- image_threshold(s$img, method = "adaptive", adaptive_k = 1.15,
                        window_size = 21, select_layer = NULL,
                        mask_layer = NULL, binary_01 = TRUE, invert = TRUE)
  expect_gt(pr_scores(dn, s$truth)$f1, 0.8)

  # Without the flag the same scene returns the background, not the roots.
  up <- image_threshold(s$img, method = "adaptive", adaptive_k = 1.15,
                        window_size = 21, select_layer = NULL,
                        mask_layer = NULL, binary_01 = TRUE)
  expect_lt(pr_scores(up, s$truth)$f1, pr_scores(dn, s$truth)$f1)
})

test_that("invert reverses which side of a global cut is kept", {
  skip_if_not_installed("terra")
  r <- terra::rast(matrix(c(0, 0.25, 0.75, 1), 2, 2))

  hi <- image_threshold(r, threshold = 0.5, method = "global",
                        select_layer = NULL, mask_layer = NULL,
                        binary_01 = TRUE)
  lo <- image_threshold(r, threshold = 0.5, method = "global",
                        select_layer = NULL, mask_layer = NULL,
                        binary_01 = TRUE, invert = TRUE)

  expect_equal(sort(terra::values(hi, mat = FALSE)), c(0, 0, 1, 1))
  expect_equal(sort(terra::values(lo, mat = FALSE)), c(0, 0, 1, 1))
  expect_true(all(terra::values(hi, mat = FALSE) +
                  terra::values(lo, mat = FALSE) == 1))
})

test_that("invert is validated and is ignored by deblur", {
  skip_if_not_installed("terra")
  r <- terra::rast(matrix(c(0, 0.5, 0.5, 1), 2, 2))
  expect_error(image_threshold(r, select_layer = NULL, mask_layer = NULL,
                               invert = NA), "TRUE or FALSE")
  two <- c(r, r * 0.5)
  expect_warning(image_threshold(two, select_layer = 2, mask_layer = 1,
                                 deblur = TRUE, invert = TRUE),
                 "not used when deblur")
})
