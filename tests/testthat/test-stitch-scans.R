# Stitching is scored against tiles cut from a known image, so the expected
# offset is the one imposed when cutting and the expected mosaic is the image
# the tiles came from. Nothing here is read back off the function under test.

# Noise makes a sharp correlation peak and has no periodicity to alias onto,
# which is what the estimator is entitled to assume about a real scan.
stitch_fixture <- function(h = 48L, w = 240L, seed = 11L) {
  set.seed(seed)
  array(runif(h * w * 3, 0, 255), dim = c(h, w, 3))
}

# Tile 1 is columns 1:w1; tile 2 starts `overlap` columns before tile 1 ends,
# so the two share exactly `overlap` columns and together cover the image.
cut_tiles <- function(base, w1, overlap) {
  list(base[, 1:w1, , drop = FALSE],
       base[, (w1 - overlap + 1):ncol(base), , drop = FALSE])
}

test_that("align_phase_correlation reports the overlap it was given", {
  base <- stitch_fixture()
  h <- dim(base)[1]
  # dx = edge_width - overlap is the documented contract, so it is the
  # assertion. edge_width is a search area: any band at least as wide as the
  # overlap must find it.
  for (overlap in c(40L, 60L, 90L)) {
    tl <- cut_tiles(base, 180L, overlap)
    for (edge_width in c(overlap, overlap + 20L, as.integer(1.8 * overlap))) {
      a <- align_phase_correlation(tl[[1]], tl[[2]], edge_width = edge_width,
                                   vertical_region = h, vertical_offset = 0)
      # The band is clamped to the narrower frame, so the contract holds for
      # the band actually used rather than the value passed in.
      ew_used <- min(edge_width, ncol(tl[[1]]), ncol(tl[[2]]))
      expect_equal(as.numeric(a[["dx"]]), ew_used - overlap,
                   info = sprintf("overlap=%d edge_width=%d (band %d)",
                                  overlap, edge_width, ew_used))
      expect_equal(as.numeric(a[["dy"]]), 0)
    }
  }
})

test_that("a band equal to the overlap gives a perfect correlation peak", {
  base <- stitch_fixture()
  tl   <- cut_tiles(base, 180L, 60L)
  # The two bands are then the same pixels, so the peak carries all the mass.
  a <- align_phase_correlation(tl[[1]], tl[[2]], edge_width = 60L,
                               vertical_region = dim(base)[1], vertical_offset = 0)
  expect_equal(as.numeric(a[["peak"]]), 1, tolerance = 1e-6)
})

test_that("frames that do not overlap score far below ones that do", {
  base <- stitch_fixture()
  h    <- dim(base)[1]
  tl   <- cut_tiles(base, 180L, 60L)
  disjoint <- list(base[, 1:60, , drop = FALSE], base[, 121:180, , drop = FALSE])
  # The peak share depends on the band size, so this is a comparison between
  # two fixtures of the same band size rather than an absolute threshold.
  true_peak <- align_phase_correlation(tl[[1]], tl[[2]], edge_width = 60L,
                                       vertical_region = h, vertical_offset = 0)[["peak"]]
  none_peak <- align_phase_correlation(disjoint[[1]], disjoint[[2]], edge_width = 60L,
                                       vertical_region = h, vertical_offset = 0)[["peak"]]
  expect_gt(as.numeric(true_peak), 10 * as.numeric(none_peak))
})

test_that("stitch_image_pair puts the image back together", {
  base <- stitch_fixture()
  h    <- dim(base)[1]
  for (overlap in c(40L, 60L, 90L)) {
    tl <- cut_tiles(base, 180L, overlap)
    for (edge_width in c(overlap, overlap + 20L)) {
      m <- stitch_image_pair(tl[[1]], tl[[2]], edge_width = edge_width,
                             vertical_region = h, vertical_offset = 0,
                             blend = "overlay_first")
      # overlay_first keeps tile 1 across the overlap band, so a correct
      # placement reproduces the original exactly, not approximately.
      expect_equal(dim(m), dim(base),
                   info = sprintf("overlap=%d edge_width=%d", overlap, edge_width))
      expect_equal(max(abs(m - base)), 0)
      expect_equal(as.numeric(attr(m, "offset")[["dx"]]), edge_width - overlap)
    }
  }
})

test_that("a search area many times the overlap still finds it", {
  # The point of calling edge_width a search area: when the overlap is unknown
  # you give a generous upper bound, and the estimate must not degrade.
  base <- stitch_fixture(h = 40L, w = 420L, seed = 5L)
  h    <- dim(base)[1]
  w1   <- 200L
  for (overlap in c(12L, 25L, 50L)) {
    tl <- cut_tiles(base, w1, overlap)
    for (search in c(overlap, 4L * overlap, 8L * overlap, 200L)) {
      if (search > min(w1, ncol(tl[[2]]))) next
      a <- align_phase_correlation(tl[[1]], tl[[2]], edge_width = search,
                                   vertical_region = h, vertical_offset = 0)
      expect_equal(as.numeric(a[["dx"]]), search - overlap,
                   info = sprintf("overlap=%d search=%d (%.0fx)",
                                  overlap, search, search / overlap))
      m <- stitch_image_pair(tl[[1]], tl[[2]], edge_width = search,
                             vertical_region = h, vertical_offset = 0,
                             blend = "overlay_first")
      expect_equal(dim(m), dim(base))
      expect_equal(max(abs(m - base)), 0)
    }
  }
})

test_that("a search area narrower than the overlap cannot find it", {
  # The one precondition: the band has to contain the overlap. An overlap wider
  # than the band lies outside the pixels being compared, so it is reported as a
  # smaller one. The test pins the precondition, not the wrong number it gives.
  base <- stitch_fixture(h = 40L, w = 420L, seed = 5L)
  h    <- dim(base)[1]
  tl   <- cut_tiles(base, 200L, 120L)
  a <- align_phase_correlation(tl[[1]], tl[[2]], edge_width = 60L,
                               vertical_region = h, vertical_offset = 0)
  expect_lt(60L - as.numeric(a[["dx"]]), 120L)
  m <- suppressWarnings(stitch_image_pair(tl[[1]], tl[[2]], edge_width = 60L,
                                          vertical_region = h, vertical_offset = 0,
                                          blend = "overlay_first"))
  expect_false(dim(m)[2] == dim(base)[2] && max(abs(m - base)) == 0)
})

test_that("stitch_image_sequence reassembles three tiles", {
  base <- stitch_fixture(w = 300L)
  h    <- dim(base)[1]
  ov   <- 50L
  frames <- list(base[,   1:140, , drop = FALSE],
                 base[,  91:240, , drop = FALSE],
                 base[, 191:300, , drop = FALSE])
  # One generous search area for every step, which is how a sequence is run
  # when the overlap is only known approximately.
  res <- stitch_image_sequence(frames, edge_width = 2L * ov, vertical_region = h,
                               vertical_offset = 0, blend = "overlay_first",
                               return_offsets = TRUE)
  expect_equal(dim(res$mosaic), dim(base))
  expect_equal(max(abs(res$mosaic - base)), 0)
  expect_equal(nrow(res$offsets), 2L)
  expect_true(all(res$offsets$overlap == ov))
})

test_that("direction = 'vertical' stitches along rows", {
  # The vertical path transposes, stitches horizontally and transposes back,
  # so the same construction rotated must give the same reconstruction.
  base <- stitch_fixture(h = 240L, w = 48L, seed = 3L)
  ov   <- 60L
  t1 <- base[1:180, , , drop = FALSE]
  t2 <- base[(180 - ov + 1):240, , , drop = FALSE]
  m <- stitch_image_pair(t1, t2, edge_width = ov, direction = "vertical",
                         vertical_region = dim(base)[2], vertical_offset = 0,
                         blend = "overlay_first")
  expect_equal(dim(m), dim(base))
  expect_equal(max(abs(m - base)), 0)
})

test_that("the overlap bracket is what the peak is read from", {
  base <- stitch_fixture()
  h    <- dim(base)[1]
  ov   <- 60L
  tl   <- cut_tiles(base, 180L, ov)
  band <- 120L

  # A bracket that contains the truth leaves the answer alone.
  inside <- align_phase_correlation(tl[[1]], tl[[2]], edge_width = band,
                                    vertical_region = h, vertical_offset = 0,
                                    min_overlap = 40L, max_overlap = 90L)
  expect_equal(band - as.numeric(inside[["dx"]]), ov)

  # A bracket above the truth puts the answer inside the bracket, not on the
  # truth: the columns holding the real match are not searched at all.
  above <- align_phase_correlation(tl[[1]], tl[[2]], edge_width = band,
                                   vertical_region = h, vertical_offset = 0,
                                   min_overlap = 80L, max_overlap = band)
  ov_above <- band - as.numeric(above[["dx"]])
  expect_gte(ov_above, 80L)
  expect_lte(ov_above, band)

  # Same from the other side.
  below <- align_phase_correlation(tl[[1]], tl[[2]], edge_width = band,
                                   vertical_region = h, vertical_offset = 0,
                                   max_overlap = 40L)
  ov_below <- band - as.numeric(below[["dx"]])
  expect_gte(ov_below, 1L)
  expect_lte(ov_below, 40L)
})

test_that("an overlap bracket wider than the band is refused or clamped", {
  base <- stitch_fixture()
  h    <- dim(base)[1]
  tl   <- cut_tiles(base, 180L, 60L)
  # A floor the band cannot reach has no admissible column, so it errors
  # rather than quietly searching the whole band.
  expect_error(align_phase_correlation(tl[[1]], tl[[2]], edge_width = 80L,
                                       vertical_region = h, vertical_offset = 0,
                                       min_overlap = 200L),
               "exceeds the search band")
  # A ceiling above the band is simply the whole band.
  wide <- align_phase_correlation(tl[[1]], tl[[2]], edge_width = 80L,
                                  vertical_region = h, vertical_offset = 0,
                                  max_overlap = 10000L)
  open <- align_phase_correlation(tl[[1]], tl[[2]], edge_width = 80L,
                                  vertical_region = h, vertical_offset = 0)
  expect_equal(as.numeric(wide[["dx"]]), as.numeric(open[["dx"]]))
  expect_error(align_phase_correlation(tl[[1]], tl[[2]], edge_width = 80L,
                                       vertical_region = h, vertical_offset = 0,
                                       min_overlap = 50L, max_overlap = 20L),
               "max_overlap")
  # A ceiling below one pixel admits no column at all.
  expect_error(align_phase_correlation(tl[[1]], tl[[2]], edge_width = 80L,
                                       vertical_region = h, vertical_offset = 0,
                                       max_overlap = 0),
               "at least 1 pixel")
})

test_that("stitch_image_pair passes the bracket down to the alignment", {
  base <- stitch_fixture()
  h    <- dim(base)[1]
  ov   <- 60L
  tl   <- cut_tiles(base, 180L, ov)
  # A bracket around the true overlap must still reconstruct the image, which
  # is the only evidence that the argument reaches the estimator intact.
  m <- stitch_image_pair(tl[[1]], tl[[2]], edge_width = 120L,
                         vertical_region = h, vertical_offset = 0,
                         blend = "overlay_first",
                         min_overlap = 40L, max_overlap = 90L)
  expect_equal(dim(m), dim(base))
  expect_equal(max(abs(m - base)), 0)
  expect_equal(as.numeric(attr(m, "offset")[["dx"]]), 120L - ov)
})


test_that("stitch_order_frames sorts on the key, numerically where it can", {
  # Unpadded frame numbers are exactly where the default name sort goes wrong:
  # "_10" sorts before "_2" as text.
  p <- c("scan_T067_2.tif", "scan_T067_10.tif", "scan_T067_1.tif")

  expect_equal(stitch_order_frames(p),
               c("scan_T067_1.tif", "scan_T067_10.tif", "scan_T067_2.tif"))

  expect_equal(stitch_order_frames(p, order_by = "\\d+(?=\\.)"),
               c("scan_T067_1.tif", "scan_T067_2.tif", "scan_T067_10.tif"))

  expect_equal(stitch_order_frames(p, order_by = "\\d+(?=\\.)", decreasing = TRUE),
               c("scan_T067_10.tif", "scan_T067_2.tif", "scan_T067_1.tif"))
})

test_that("stitch_order_frames falls back to text for non-numeric keys", {
  p <- c("tube_c.tif", "tube_a.tif", "tube_b.tif")
  expect_equal(stitch_order_frames(p, order_by = "[a-z](?=\\.)"),
               c("tube_a.tif", "tube_b.tif", "tube_c.tif"))
})

test_that("stitch_order_frames stops when a file has no key", {
  # An unmatched path has no defined place in the sequence, and guessing one
  # would put a frame in the wrong part of the mosaic.
  p <- c("scan_1.tif", "scan_2.tif", "notes.txt")
  expect_error(stitch_order_frames(p, order_by = "\\d+(?=\\.)"), "notes.txt")
  expect_error(stitch_order_frames(p, order_by = 7), "regular expression")
})

test_that("stitch_order_frames passes a single frame through", {
  expect_equal(stitch_order_frames("only.tif", order_by = "\\d+"), "only.tif")
})
