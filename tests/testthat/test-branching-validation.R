# Ground-truth validation of the branching pipeline.
#
# Everything here is scored against synthetic images whose length, width and
# topology are prescribed (see ?root_phantom), so a regression in tracing,
# junction handling, ordering or length integration fails a test instead of
# quietly shifting a number nobody can check.

test_that("root_phantom reports a self-consistent ground truth", {
  ph <- root_phantom("comb", size = 200)
  expect_true(all(c("mask", "skeleton", "truth", "strokes") %in% names(ph)))
  expect_equal(dim(ph$mask), c(200L, 200L))
  expect_gt(sum(ph$mask), sum(ph$skeleton))          # filled roots are thicker
  expect_equal(ph$truth$n_tips, 7L)                  # 2 axis ends + 5 laterals
  expect_equal(ph$truth$n_branch_points, 5L)
  expect_equal(ph$truth$n_roots, 6L)
  expect_equal(sum(ph$truth$by_order$total_length), ph$truth$total_length)
})

test_that("the graph route reproduces length, diameter and topology", {
  skip_if_not_installed("terra")
  for (d in c("comb", "herringbone", "hierarchical", "cross", "fork")) {
    v <- validate_branching(d, size = 200, from = "skeleton", verbose = FALSE)
    expect_true(attr(v, "passed"),
                info = paste0("design '", d, "' failed: ",
                              paste(v$metric[!v$pass], collapse = ", ")))
  }
})

test_that("the end-to-end route survives skeletonisation", {
  skip_if_not_installed("terra")
  # 'fork' is excluded on purpose: its continuation pair is a genuine tie, so
  # the per-root metrics have no ground truth (see ?root_phantom).
  for (d in c("comb", "herringbone", "hierarchical", "cross")) {
    v <- validate_branching(d, size = 200, from = "mask", verbose = FALSE)
    expect_true(attr(v, "passed"),
                info = paste0("design '", d, "' failed: ",
                              paste(v$metric[!v$pass], collapse = ", ")))
  }
})

test_that("a right-angle bend is one root, not a fork", {
  # An 8-connected skeleton makes the pixels either side of a 90-degree corner
  # look like junctions. Contracting them used to leave a 3-px self-loop: a
  # phantom branch point, an NA tip_order and a fragmented root.
  m <- matrix(0, 41, 41); m[21, 5:21] <- 1; m[21:37, 21] <- 1
  et <- root_graph_pipeline(m, m, verbose = FALSE)
  expect_equal(nrow(et), 1L)
  expect_equal(sum(et$n_branch_points), 0)
  expect_false(any(is.na(et$tip_order)))
  expect_equal(sum(et$length), 32, tolerance = 0.05)   # 16 + 16 px of centre line
})

test_that("junction contraction does not eat root length", {
  # Five laterals on a straight axis: contraction dissolves a couple of pixels
  # per branch point, which must be added back as the rim-to-centroid stub.
  m <- matrix(0, 41, 41); m[3:39, 21] <- 1
  for (k in seq(8, 36, by = 7)) m[k, 22:32] <- 1
  et <- root_graph_pipeline(m, m, verbose = FALSE)
  expect_equal(sum(et$length), 36 + 5 * 11, tolerance = 0.03)
  expect_equal(sum(et$n_branch_points), 5)
  expect_equal(sum(et$n_tips), 7)
})

test_that("a crossing resolves into two independent roots", {
  skel  <- matrix(0, 61, 61); skel[5:55, 31] <- 1; skel[31, 11:51] <- 1
  cross <- matrix(0, 61, 61)
  for (i in 5:55) cross[i, 28:34] <- 1
  for (j in 11:51) cross[28:34, j] <- 1
  et <- root_graph_pipeline(skel, cross, verbose = FALSE)
  expect_equal(length(unique(et$root_id)), 2L)
  expect_equal(sum(et$n_branch_points), 0)
})

test_that("crossing_diam_ratio trades a crossing failure for a fork failure", {
  # A thick axis crossed by a thin root and a thick axis with two thin laterals
  # are the same shape in outline. The threshold picks which one is read wrong;
  # it is off by default, so the default must reproduce the crossing reading.
  skel  <- matrix(0, 61, 61); skel[5:55, 31] <- 1; skel[31, 11:51] <- 1
  bilat <- matrix(0, 61, 61)
  for (i in 5:55) bilat[i, 27:35] <- 1
  for (j in 11:51) bilat[30:32, j] <- 1

  off <- root_graph_pipeline(skel, bilat, verbose = FALSE)
  expect_equal(sum(off$n_branch_points), 0)          # read as a crossing
  expect_equal(length(unique(off$root_id)), 2L)

  on <- root_graph_pipeline(skel, bilat, verbose = FALSE, crossing_diam_ratio = 0.5)
  expect_equal(sum(on$n_branch_points), 1)           # read as a bilateral branch
  expect_equal(length(unique(on$root_id)), 3L)
  expect_equal(max(on$branch_order, na.rm = TRUE), 2L)
})

test_that("the continuation pair at a 4-way node is the thicker one", {
  # Both the axis pair and the lateral pair are perfectly straight and perfectly
  # self-similar, so the score ties; the thick pair must win, not whichever the
  # enumeration reaches first.
  skel  <- matrix(0, 61, 61); skel[5:55, 31] <- 1; skel[31, 11:51] <- 1
  bilat <- matrix(0, 61, 61)
  for (i in 5:55) bilat[i, 27:35] <- 1
  for (j in 11:51) bilat[30:32, j] <- 1
  et <- root_graph_pipeline(skel, bilat, verbose = FALSE, crossing_diam_ratio = 0.5)
  thickest <- et$root_id[which.max(et$mean_diameter)]
  expect_equal(unique(et$branch_order[et$root_id == thickest]), 1L)
  expect_equal(max(et$branch_order, na.rm = TRUE), 2L)   # not 3: no phantom generation
})

test_that("convert_root_units is idempotent and reversible", {
  ph <- root_phantom("comb", size = 150)
  et <- root_graph_pipeline(ph$skeleton, ph$mask, verbose = FALSE)
  cm  <- convert_root_units(et, unit = "cm", dpi = 300)
  cm2 <- convert_root_units(cm, unit = "cm", dpi = 300)
  expect_equal(cm$length, cm2$length)
  back <- convert_root_units(cm, unit = "px")
  expect_equal(back$length, et$length, tolerance = 1e-9)
})

test_that("an image with no roots returns an empty result, not an error", {
  skip_if_not_installed("terra")
  z <- terra::rast(matrix(0, 12, 12))
  res <- suppressWarnings(branch_order_map(skel = z, mask = z, unit = "px"))
  expect_s3_class(res, "branchOrderMap")
  expect_null(res$class_map)
  expect_s3_class(res$summary, "data.frame")
})

test_that("a filled mask passed as a skeleton is flagged", {
  blob <- matrix(0, 40, 40); blob[10:30, 10:30] <- 1
  ws <- character()
  withCallingHandlers(
    root_graph_pipeline(blob, blob, verbose = FALSE),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("1-pixel-wide skeleton", ws)))
})
