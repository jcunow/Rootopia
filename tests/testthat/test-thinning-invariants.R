# Thinning is scored against what a medial axis is, not against recorded
# output: a 2-px line must become a 1-px arc with two ends, and a closed ring
# must stay a closed ring. Those are the properties the two-sub-iteration
# structure in lut_thin_fast() exists to preserve.

skip_if_not_installed("terra")

.px     <- function(r) sum(terra::values(r) == 1, na.rm = TRUE)
.mat    <- function(r) matrix(terra::values(r), nrow = terra::nrow(r), byrow = TRUE)
.ncomp  <- function(r) {
  p <- terra::patches(terra::subst(r, 0, NA), directions = 8, zeroAsNA = FALSE)
  length(unique(stats::na.omit(terra::values(p))))
}
.ends   <- function(r) sum(terra::values(detect_skeleton_points(r)$endpoints), na.rm = TRUE)
.branch <- function(r) sum(terra::values(detect_skeleton_points(r)$branching_points), na.rm = TRUE)

test_that("a two-pixel-wide line thins to a one-pixel arc with two ends", {
  m <- matrix(0, 20, 50); m[10:11, 6:45] <- 1
  s <- lut_thin_fast(terra::rast(m))

  # One pixel per occupied column is what "medial axis" means for a
  # horizontal bar, and it is the claim a half-finished thinning would break.
  sm <- .mat(s); cs <- colSums(sm)
  expect_equal(sum(cs > 0), 40L)
  expect_true(all(cs[cs > 0] == 1))
  expect_equal(sum(rowSums(sm) > 0), 1L)

  # A simple open arc: one component, exactly two ends, no junction.
  expect_equal(.ncomp(s), 1L)
  expect_equal(.ends(s), 2)
  expect_equal(.branch(s), 0)
})

test_that("a closed ring survives thinning as a closed ring", {
  # The comment in lut_thin_fast() states that evaluating both deletion sets
  # from one snapshot erodes a loop from both sides and destroys it. A ring has
  # no endpoint to peel from, so it is the case that fails if the two
  # sub-iterations are ever fused into one pass.
  m <- matrix(0, 30, 30); m[6:25, 6:25] <- 1; m[9:22, 9:22] <- 0
  s <- lut_thin_fast(terra::rast(m))

  expect_gt(.px(s), 0)                 # not erased
  expect_lt(.px(s), .px(terra::rast(m)))  # but thinned
  expect_equal(.ncomp(s), 1L)          # still one connected structure
  expect_equal(.ends(s), 0)            # every pixel has two neighbours: a cycle
})

test_that("max_iter stops thinning early instead of running to convergence", {
  m <- matrix(0, 40, 40); m[6:35, 6:35] <- 1
  converged <- lut_thin_fast(terra::rast(m))
  capped    <- lut_thin_fast(terra::rast(m), max_iter = 2L)

  # Two passes cannot thin a 30x30 square, so the capped result must still
  # carry far more foreground than the converged one.
  expect_gt(.px(capped), .px(converged))
  # It is still a usable raster on the input grid, which is what makes an
  # unnoticed cap dangerous: nothing downstream can tell it is unfinished.
  expect_equal(dim(capped), dim(terra::rast(m)))
  expect_setequal(unique(stats::na.omit(terra::values(capped))), c(0, 1))
})

test_that("detect_skeleton_points counts pixels by their eight-neighbourhood", {
  # A one-pixel plus: four arm tips have one neighbour each, so four endpoints.
  m <- matrix(0, 21, 21); m[11, 4:18] <- 1; m[4:18, 11] <- 1
  pts <- detect_skeleton_points(terra::rast(m))

  expect_equal(sum(terra::values(pts$endpoints), na.rm = TRUE), 4)

  # Five pixels, not one, exceed two neighbours: the centre, and each of the
  # four pixels orthogonally adjacent to it -- those see two arm pixels
  # diagonally as well as their own two along the arm. A branch *point* here is
  # a pixel, not a junction, which is why trace_segments() merges adjacent
  # junction pixels into one graph node before building the segment graph.
  bm <- matrix(terra::values(pts$branching_points), nrow = 21, byrow = TRUE)
  expect_equal(sum(bm), 5)
  expect_equal(which(bm == 1, arr.ind = TRUE)[, c("row", "col")],
               cbind(row = c(11L, 10L, 11L, 12L, 11L),
                     col = c(10L, 11L, 11L, 11L, 12L)))
})

test_that("skeletonize_image refuses an image with no foreground", {
  m <- matrix(0, 10, 10)
  expect_error(skeletonize_image(terra::rast(m), verbose = FALSE),
               "No foreground")
})
