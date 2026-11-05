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

test_that("codes at root pixels only give the same skeleton as codes at every cell", {
  # lut_thin_fast() computes the neighborhood code at the current root pixels
  # only. Every other cell has code 0, which the table keeps, so on any mask
  # without NA the result must equal the same rule evaluated at every cell.
  thin_every_cell <- function(m, max_iter = 200L) {
    nr <- nrow(m); nc <- ncol(m)
    deletion_flags <- function(m) {
      mp <- matrix(0L, nr + 2, nc + 2)
      mp[2:(nr + 1), 2:(nc + 1)] <- m
      at <- function(dr, dc) mp[(2:(nr + 1)) + dr, (2:(nc + 1)) + dc]
      code <-   1L*at(-1, -1) +  2L*at(-1, 0) +  4L*at(-1, 1) +
              128L*at( 0, -1) +                  8L*at( 0, 1) +
               64L*at( 1, -1) + 32L*at( 1, 0) + 16L*at( 1, 1)
      code[m != 1L] <- 0L
      .thinning_lut[code + 1L]          # column-major, like m
    }
    for (i in seq_len(max_iter)) {
      d1 <- deletion_flags(m) %in% c(1L, 3L); m[d1] <- 0L
      d2 <- deletion_flags(m) %in% c(2L, 3L); m[d2] <- 0L
      if (!any(d1) && !any(d2)) break
    }
    m
  }

  set.seed(42)
  for (i in 1:60) {
    d <- sample(3:60, 2)
    m <- matrix(rbinom(prod(d), 1, runif(1, 0.1, 0.9)), d[1], d[2])
    if (i %% 4 == 0) m[, 1] <- 1L     # roots on the image border
    if (i %% 9 == 0) m[] <- 1L        # nothing but root
    s <- .mat(lut_thin_fast(terra::rast(m)))
    expect_identical(s == 1, thin_every_cell(m) == 1, info = sprintf("mask %d", i))
  }
})

test_that("masked cells stay NA and act like the image border", {
  # A 5-px bar cut by a masked band. Masked cells read as background, so each
  # side must thin exactly as it would as an image of its own whose edge lies
  # where the mask starts.
  m <- matrix(0, 15, 60); m[6:10, 5:56] <- 1
  masked <- m; masked[, 28:32] <- NA
  s <- .mat(lut_thin_fast(terra::rast(masked)))

  expect_true(all(is.na(s[, 28:32])))
  expect_false(anyNA(s[, -(28:32)]))
  expect_equal(s[, 1:27],  .mat(lut_thin_fast(terra::rast(m[, 1:27]))))
  expect_equal(s[, 33:60], .mat(lut_thin_fast(terra::rast(m[, 33:60]))))
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
