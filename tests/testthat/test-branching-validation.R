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
  # look like junctions. Contracting them must not leave a 3-px self-loop,
  # which would show up as a phantom branch point, an NA tip_order and a
  # fragmented root.
  m <- matrix(0, 41, 41); m[21, 5:21] <- 1; m[21:37, 21] <- 1
  et <- root_graph_pipeline(m, m, verbose = FALSE)
  expect_equal(nrow(et), 1L)
  expect_equal(sum(et$n_branch_points), 0)
  expect_false(any(is.na(et$tip_order)))
  expect_equal(sum(et$length), 32, tolerance = 0.05)   # 16 + 16 px of center line
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

test_that("the reported length column is Kimura's estimator by default", {
  # validate_branching() deliberately scores the chain-code column, because
  # that is the convention the phantom's truth is written in. This pins the
  # package default separately, so pinning the fixture cannot quietly become a
  # revert of the default.
  expect_equal(eval(formals(convert_root_units)$length_method)[1], "kimura")
  expect_equal(formals(branch_order_map)$length_method, "kimura")

  m <- matrix(0, 11, 41); m[6, 5:37] <- 1
  et <- root_graph_pipeline(m, m, verbose = FALSE)
  conv <- convert_root_units(et, unit = "px", dpi = 300)
  expect_equal(conv$length, conv$length_kimura)
})

test_that("the two length measures agree exactly on pure chains", {
  # Kimura's estimator is sqrt(nd^2 + (nd + no/2)^2) + no/2. With no = 0 that is
  # nd * sqrt(2) and with nd = 0 it is no, so on a chain of only diagonal or
  # only orthogonal steps it reduces to the chain code. That is why the
  # phantom's axis-aligned and 45-degree designs score the same either way on
  # the exact centre line, and it is what makes the fixture's convention choice
  # harmless for the skeleton route.
  horiz <- matrix(0, 11, 41); horiz[6, 5:37] <- 1
  diag_ <- matrix(0, 41, 41); for (k in 0:32) diag_[5 + k, 5 + k] <- 1

  for (m in list(horiz, diag_)) {
    et <- root_graph_pipeline(m, m, verbose = FALSE)
    expect_true(et$n_orth == 0L || et$n_diag == 0L)   # the chain really is pure
    expect_equal(et$length, et$length_kimura)
  }
})

test_that("Kimura is the closer measure once thinning jitters a diagonal", {
  skip_if_not_installed("terra")
  # Thinning a drawn 45-degree lateral does not return a pure diagonal chain:
  # a few orthogonal steps appear, and that is exactly where the two measures
  # part company. The lateral is straight, so the distance between its first
  # and last traced pixel is its true length, and both measures can be scored
  # against it. The chain code reads high because it walks the staircase;
  # Kimura corrects for it and lands closer.
  ph <- root_phantom("herringbone", size = 200)
  skel <- .to_binary_matrix(skeletonize_image(ph$mask, verbose = FALSE))
  et <- root_graph_pipeline(skel, ph$mask, verbose = FALSE, keep_segments = TRUE)
  segs <- attr(et, "segments")

  lat <- which(et$branch_order == 2)
  expect_gt(length(lat), 0)
  i <- lat[1]
  expect_gt(et$n_orth[i], 0)        # jitter is present
  expect_gt(et$n_diag[i], 0)

  p <- segs[[i]]$coords
  span <- sqrt(sum((p[nrow(p), ] - p[1L, ])^2))
  err_poly   <- (et$length[i] - span) / span
  err_kimura <- (et$length_kimura[i] - span) / span

  expect_gt(err_poly, 0)                                 # staircase reads high
  expect_lt(abs(err_kimura), abs(err_poly))              # Kimura reads closer
  expect_lt(abs(err_kimura), 0.05)
})

test_that("Kimura beats the chain code on strokes of known length", {
  skip_if_not_installed("terra")
  # The phantom designs run along the axes or at 45 degrees, the orientations
  # where the two length measures coincide exactly, so validate_branching()
  # cannot tell them apart. This scores them where they do differ: isolated
  # strokes of analytically known length, drawn, filled, thinned and traced, so
  # the figures carry the thinning error as well as the measure. No junctions,
  # so junction contraction cannot confound the comparison.
  N <- 320L
  one <- function(v) {
    true_len <- sum(sqrt(rowSums(diff(v)^2)))
    m  <- .stamp_stroke(matrix(0, N, N), .polyline_px(round(v)), 2L)
    sk <- .to_binary_matrix(skeletonize_image(m, verbose = FALSE))
    et <- root_graph_pipeline(sk, m, verbose = FALSE)
    if (nrow(et) != 1L) return(NULL)        # only score a cleanly traced stroke
    c(true = true_len, poly = et$length[1], kimura = et$length_kimura[1])
  }

  rows <- list()
  for (deg in seq(6, 84, by = 13)) {        # 0, 45 and 90 are the tie cases
    th <- deg * pi / 180; L <- 150
    rows[[length(rows) + 1L]] <-
      one(rbind(c(240, 70), c(240 - L * sin(th), 70 + L * cos(th))))
  }
  for (base in c(10, 35, 70)) {             # gently curved, as a root is
    t <- seq(0, 150, length.out = 150); th <- base * pi / 180
    y <- 8 * sin(2 * pi * t / 60)
    rows[[length(rows) + 1L]] <-
      one(cbind(240 - (t * sin(th) + y * cos(th)),
                 70 + (t * cos(th) - y * sin(th))))
  }
  rows <- do.call(rbind, Filter(Negate(is.null), rows))
  expect_gte(nrow(rows), 8L)

  err  <- function(col) (rows[, col] - rows[, "true"]) / rows[, "true"]
  rms  <- function(x) sqrt(mean(x^2))
  ep <- err("poly"); ek <- err("kimura")

  # the staircase bias is real and upward
  expect_gt(mean(ep), 0.02)
  expect_gt(rms(ep), 0.03)
  # and Kimura removes most of it
  expect_lt(rms(ek), 0.03)
  expect_lt(rms(ek), rms(ep) / 2)
  expect_gte(sum(abs(ek) < abs(ep)), ceiling(0.7 * nrow(rows)))
})
