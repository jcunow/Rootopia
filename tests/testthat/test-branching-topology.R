# Topology cases the phantom designs do not reach.
#
# The phantoms in test-branching-validation.R are all single, connected,
# acyclic root systems scored against known geometry. The cases here are the
# awkward ones: pruning, two systems in one image, a skeleton with no tips at
# all, and where the classification raster actually puts its values.

# A straight axis with one long lateral and one short spur, plus the filled
# mask that gives both a real thickness.
spur_image <- function() {
  skel <- matrix(0, 61, 61)
  skel[5:55, 31] <- 1                      # main axis
  skel[20, 32:50] <- 1                     # lateral, ~19 px
  skel[40, 32:36] <- 1                     # spur, ~5 px
  mask <- skel
  for (i in 5:55) mask[i, 29:33] <- 1
  for (j in 32:50) mask[19:21, j] <- 1
  for (j in 32:36) mask[39:41, j] <- 1
  list(skel = skel, mask = mask)
}


test_that("strahler_order rises only where two equal orders meet", {
  # The axis carries one lateral and one spur, so under Strahler it stays
  # order 2 above the junction while every tip -- and the basal piece below the
  # last junction -- is 1. tip_order peels from every free end instead, so it
  # scores the same skeleton differently.
  im <- spur_image()
  et <- root_graph_pipeline(im$skel, im$mask, verbose = FALSE)

  expect_true(all(c("strahler_order", "tip_order") %in% names(et)))
  expect_setequal(unique(et$strahler_order), c(1L, 2L))
  expect_equal(sum(et$strahler_order == 1L), 3L)
  expect_equal(sum(et$strahler_order == 2L), 2L)
  expect_equal(sum(et$tip_order == 1L), 4L)
})


test_that("pruning removes short terminal segments and leaves the rest intact", {
  im <- spur_image()
  keep <- root_graph_pipeline(im$skel, im$mask, verbose = FALSE)
  cut  <- root_graph_pipeline(im$skel, im$mask, verbose = FALSE,
                              prune_min_length = 10, prune_iter = 1)

  expect_equal(nrow(cut), nrow(keep) - 1L)             # exactly the spur went
  expect_true(min(keep$length) < 10)
  expect_gte(min(cut$length), 10)
  expect_equal(sum(cut$length), sum(keep$length) - min(keep$length), tolerance = 1e-8)
})


test_that("pruning thresholds the same length the edge table reports", {
  # Pruning and the edge table must measure the same thing: the raw pixel chain
  # plus the rim-to-centroid junction stub (~1 px per junction end). If pruning
  # used the bare chain, a segment could be dropped at a min_length it visibly
  # clears in the reported output.
  skel <- matrix(0, 61, 61); skel[5:55, 31] <- 1; skel[40, 32:44] <- 1
  mask <- skel
  for (i in 5:55) mask[i, 27:35] <- 1
  for (j in 32:44) mask[38:42, j] <- 1

  full <- root_graph_pipeline(skel, mask, verbose = FALSE)
  lat  <- min(full$length)                             # the lateral, stub included
  expect_gt(lat, 12); expect_lt(lat, 13)               # raw chain is 12 px

  # A threshold between the raw length and the reported length must keep it.
  thr <- (12 + lat) / 2
  cut <- root_graph_pipeline(skel, mask, verbose = FALSE,
                             prune_min_length = thr, prune_iter = 1)
  expect_equal(nrow(cut), nrow(full))
  expect_gte(min(cut$length), thr)
})


test_that("prune_skeleton returns the class it was given, in both output modes", {
  im <- spur_image()
  skip_if_not_installed("terra")

  m <- prune_skeleton(im$skel, im$mask, min_length = 10, iter = 1)
  expect_true(is.matrix(m))
  expect_lt(sum(m), sum(im$skel))                      # the spur is gone
  expect_equal(dim(m), dim(im$skel))

  r <- prune_skeleton(terra::rast(im$skel), terra::rast(im$mask),
                      min_length = 10, iter = 1)
  expect_s4_class(r, "SpatRaster")

  # "mask" mode regrows the survivors to their own thickness and clips to the
  # input mask, so it must shrink the mask without ever exceeding it.
  mk <- prune_skeleton(im$skel, im$mask, min_length = 10, iter = 1, output = "mask")
  expect_lt(sum(mk), sum(im$mask))
  expect_true(all(mk[im$mask == 0] == 0))
})


test_that("two disconnected root systems are ordered independently", {
  # branch_order seeds the thickest root of each connected component at 1, so
  # neither system may inherit the other's numbering.
  skel <- matrix(0, 61, 61)
  skel[5:30, 10] <- 1
  for (k in seq(9, 27, by = 6)) skel[k, 10:18] <- 1
  skel[35:58, 40] <- 1
  for (k in seq(39, 55, by = 6)) skel[k, 40:50] <- 1
  mask <- skel
  for (i in 5:30)  mask[i,  9:11] <- 1
  for (i in 35:58) mask[i, 38:42] <- 1                 # second axis deliberately thicker

  et <- root_graph_pipeline(skel, mask, verbose = FALSE, keep_segments = TRUE)
  expect_false(any(is.na(et$branch_order)))

  # Each component must contain a root of order 1: split the segments by which
  # half of the image they sit in and check both halves start at 1.
  segs <- attr(et, "segments")
  half <- vapply(segs, function(s) mean(s$coords[, 1]) < 32, logical(1))
  expect_equal(min(et$branch_order[half]),  1L)
  expect_equal(min(et$branch_order[!half]), 1L)
})


test_that("a skeleton with no tips is reported as unordered, not silently wrong", {
  ring <- matrix(0, 41, 41)
  ring[11:31, 11] <- 1; ring[11:31, 31] <- 1
  ring[11, 11:31] <- 1; ring[31, 11:31] <- 1

  ws <- character()
  et <- withCallingHandlers(
    root_graph_pipeline(ring, ring, verbose = FALSE),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })

  expect_true(any(grepl("cycle", ws)))
  expect_true(all(is.na(et$tip_order)))
  expect_true(all(is.na(et$root_order)))
})


test_that("order_classification_map writes each order at its own pixels", {
  skip_if_not_installed("terra")
  skel <- matrix(0, 41, 41); skel[21, 5:21] <- 1; skel[21:37, 21] <- 1
  res <- branch_order_map(terra::rast(skel), terra::rast(skel),
                          order = "branch_order", unit = "px", verbose = FALSE)
  cm <- terra::as.matrix(res$class_map, wide = TRUE)

  expect_equal(dim(cm), dim(skel))
  expect_true(all(is.na(cm[skel == 0])))               # nothing painted off-root
  expect_true(all(cm[skel == 1] %in% c(1L, NA_integer_)))
  # A row/col transpose would land the values on the mirrored arm, so check a
  # pixel that is on the skeleton but not on its mirror image.
  expect_false(is.na(cm[21, 8]))
  expect_true(is.na(cm[8, 21]))
})


test_that("the classification map skips exactly one pixel per contracted junction", {
  # Junction contraction dissolves the cluster interior, so those pixels belong
  # to no segment and stay NA. This is worth pinning: it is why the map has a
  # pinhole at every branch point, and why pixel counts off the map run slightly
  # below the skeleton.
  skip_if_not_installed("terra")
  ph <- root_phantom("comb", size = 200)
  res <- branch_order_map(terra::rast(ph$skeleton), terra::rast(ph$mask),
                          order = "branch_order", unit = "px", verbose = FALSE)
  cm <- terra::as.matrix(res$class_map, wide = TRUE)

  missing <- sum(ph$skeleton == 1 & is.na(cm))
  expect_equal(missing, sum(res$edges$n_branch_points))
})
