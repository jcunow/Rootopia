# Depth map construction + zoning / slicing utilities.

test_that("create_depthmap returns a raster matching the input grid", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  dm  <- create_depthmap(img, sinoid = TRUE, tube_thicc = 7,
                         angle_from_vertical = 45, dpi = 300)
  expect_s4_class(dm, "SpatRaster")
  # the depth map carries the input grid: same dims, same extent/res/CRS
  expect_equal(dim(dm)[1:2], dim(img)[1:2])
  expect_true(terra::compareGeom(img, dm, stopOnError = FALSE))
})

test_that("create_depthmap validates parameters", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  expect_error(create_depthmap(img, angle_from_vertical = -1))
  expect_error(create_depthmap(img, angle_from_vertical = 91))
  expect_error(create_depthmap(img, tube_thicc = -1))
  expect_error(create_depthmap(img, center_offset = 2))   # must be in [0,1]
  expect_error(create_depthmap(img, geometry = "cylinder"))
})

test_that("both ends of the angle range are real cases, not errors", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  # Upright: one pixel of width is one pixel of depth, so the axis spans
  # ncol * 2.54/dpi. Horizontal: the axis stays at a single depth.
  span <- function(x) unname(diff(terra::minmax(x)[, 1]))

  upright <- create_depthmap(img, geometry = "window",
                             angle_from_vertical = 0, dpi = 300)
  expect_equal(span(upright), (terra::ncol(img) - 1) * 2.54 / 300,
               tolerance = 1e-9)

  horizontal <- create_depthmap(img, geometry = "window",
                                angle_from_vertical = 90, dpi = 300)
  expect_equal(span(horizontal), 0)
})

test_that("the length and circumference terms are complements of each other", {
  skip_if_not_installed("terra")
  # One full turn of the wave, so the spread down a column is the whole
  # curvature term and the spread across a row is the whole axis term.
  tube <- 7; dpi <- 30
  nr   <- round(tube * pi * dpi / 2.54, 0) - 1
  img  <- make_binary_spatraster(nr = nr, nc = 120)

  terms <- function(a) {
    m <- terra::as.matrix(
      create_depthmap(img, geometry = "tube", tube_thicc = tube,
                      angle_from_vertical = a, dpi = dpi, center_offset = 0),
      wide = TRUE)
    c(circumference = diff(range(m[, 1])), length = diff(range(m[1, ])))
  }

  # Upright: all of the depth is along the length, none around the tube.
  up <- terms(0)
  expect_equal(unname(up["circumference"]), 0)
  expect_gt(unname(up["length"]), 0)

  # Horizontal: the reverse, and the circumference spans a full diameter.
  flat <- terms(90)
  expect_equal(unname(flat["length"]), 0)
  expect_equal(unname(flat["circumference"]), tube, tolerance = 1e-3)

  # The ratio has to move with the angle. If both terms carried the same factor
  # it would be constant, and the angle would only rescale the map.
  r <- vapply(c(20, 45, 70), function(a) {
    t <- terms(a); unname(t["circumference"] / t["length"])
  }, numeric(1))
  expect_false(isTRUE(all.equal(r[1], r[3])))
  expect_true(!is.unsorted(r))          # leaning over shifts depth to the tube
})

test_that("tilt is converted from the horizontal convention, with a warning", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  expect_warning(dm <- create_depthmap(img, tilt = 30), "angle_from_vertical")
  expect_equal(terra::as.matrix(dm, wide = TRUE),
               terra::as.matrix(create_depthmap(img, angle_from_vertical = 60),
                                wide = TRUE))
  # 45 is the same angle under either reference.
  expect_equal(
    terra::as.matrix(suppressWarnings(create_depthmap(img, tilt = 45)), wide = TRUE),
    terra::as.matrix(create_depthmap(img, angle_from_vertical = 45), wide = TRUE))
  # Given both and disagreeing, it stops rather than silently picking one.
  expect_error(create_depthmap(img, tilt = 30, angle_from_vertical = 10),
               "different angles")
})

test_that("geometry = 'flat' puts every pixel at one depth", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  # No curvature and no depth axis, so the map is constant at -start_soil
  # whatever the tube arguments say -- and a profile over it has one bin.
  dm <- create_depthmap(img, geometry = "flat", angle_from_vertical = 45,
                        tube_thicc = 7, dpi = 300, start_soil = 2.5)
  expect_true(terra::compareGeom(img, dm, stopOnError = FALSE))
  expect_equal(unname(terra::minmax(dm)[, 1]), c(-2.5, -2.5))
  expect_length(unique(terra::values(binning(dm, nn = 5), mat = FALSE)), 1L)
})

test_that("geometry = 'window' drops the curvature but keeps the depth axis", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  win <- create_depthmap(img, geometry = "window", angle_from_vertical = 45,
                         dpi = 300)
  wmat <- terra::as.matrix(win, wide = TRUE)
  # Curvature varies down the rows; without it every row is identical.
  expect_equal(wmat[1, ], wmat[nrow(wmat), ])
  # The depth axis survives: columns still increase left to right.
  expect_false(is.unsorted(wmat[1, ]))
  expect_gt(diff(range(wmat[1, ])), 0)
})

test_that("geometry supersedes sinoid, and disagreement warns", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  # The old argument keeps working on its own.
  expect_equal(terra::as.matrix(create_depthmap(img, sinoid = FALSE), wide = TRUE),
               terra::as.matrix(create_depthmap(img, geometry = "window"), wide = TRUE))
  expect_equal(terra::as.matrix(create_depthmap(img, sinoid = TRUE), wide = TRUE),
               terra::as.matrix(create_depthmap(img, geometry = "tube"), wide = TRUE))
  # Given both, geometry wins and says so rather than silently picking one.
  expect_warning(create_depthmap(img, sinoid = TRUE, geometry = "window"),
                 "geometry")
})

test_that("center_offset sweeps the tube top over the full image height", {
  skip_if_not_installed("terra")
  # One period of the curvature wave spans the tube circumference,
  # round(tube_thicc * pi * dpi / 2.54) - 1 = 129 px here, so an image 129 rows
  # tall holds exactly one period and the trough row reads the phase directly.
  # The trough is the shallowest point, i.e. the top of the tube, and it sits on
  # row 1 + (nr - 1) * center_offset.
  img <- make_binary_spatraster(nr = 129)
  tube_top_row <- function(offset) {
    dm <- create_depthmap(img, sinoid = TRUE, tube_thicc = 3.5,
                          angle_from_vertical = 45,
                          dpi = 30, center_offset = offset)
    which.min(terra::as.matrix(dm, wide = TRUE)[, 1])
  }
  expect_equal(tube_top_row(0),    1L)
  expect_equal(tube_top_row(0.25), 33L)
  expect_equal(tube_top_row(0.50), 65L)
  expect_equal(tube_top_row(0.75), 97L)
  # Monotone over the whole range, not just the first half.
  offsets <- seq(0, 0.99, by = 0.01)
  expect_true(!is.unsorted(vapply(offsets, tube_top_row, integer(1))))
  expect_equal(tube_top_row(0.99), 128L)
})

test_that("center_offset 0 and 1 are the same point on the tube", {
  skip_if_not_installed("terra")
  # A full turn returns to the image seam, so the two ends agree and the trough
  # sits on the first and last row at once.
  img <- make_binary_spatraster(nr = 129)
  dm <- function(offset) {
    terra::as.matrix(create_depthmap(img, sinoid = TRUE, tube_thicc = 3.5,
                                     angle_from_vertical = 45, dpi = 30,
                                     center_offset = offset), wide = TRUE)
  }
  expect_equal(dm(0), dm(1))
  col1 <- dm(0)[, 1]
  expect_equal(col1[1], col1[129])
  expect_equal(which.min(col1), 1L)
})

test_that("slice_rotation splits into n contiguous slices", {
  skip_if_not_installed("terra")
  img    <- make_binary_spatraster()
  slices <- slice_rotation(img, n = 3)
  expect_length(slices, 3)
  expect_true(all(vapply(slices, inherits, logical(1), "SpatRaster")))
  total_rows <- sum(vapply(slices, terra::nrow, numeric(1)))  # terra::nrow is double
  expect_equal(total_rows, terra::nrow(img))              # exact, no overlap/gap
})

test_that("slice_rotation validates n", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  expect_error(slice_rotation(img, n = 0))
  expect_error(slice_rotation(img, n = terra::nrow(img) + 1))
})

test_that("depth_zoning masks pixels outside the selected depth bin", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster(24, 28)
  dmat <- matrix(10, 24, 28); dmat[1:12, ] <- 5      # top half = 5, bottom = 10
  dm <- terra::rast(dmat)
  z <- depth_zoning(img, depth_map = dm, depth = 5)
  expect_s4_class(z, "SpatRaster")
  expect_equal(terra::nrow(z), terra::nrow(img))      # extent unchanged
  zmat <- matrix(terra::values(z), nrow = 24, byrow = TRUE)
  expect_true(all(is.na(zmat[13:24, ])))              # depth-10 half is masked out
  expect_false(all(is.na(zmat[1:12, ])))              # depth-5 half retained
})

test_that("depth_zoning requires depth_map and depth", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  expect_error(depth_zoning(img, depth = 5))
  expect_error(depth_zoning(img, depth_map = make_binary_spatraster()))
})


test_that("slice_rotation tiles the rotation axis exactly", {
  skip_if_not_installed("terra")
  r  <- terra::rast(matrix(seq_len(20 * 7), 20, 7))
  nr <- terra::nrow(r)
  # The claim is contiguous, non-overlapping slices covering every row. Row
  # counts summing to nr proves no gap and no overlap together, and n values
  # that do not divide nr are the cases a rounding error would show up in.
  for (n in c(1L, 3L, 4L, 7L)) {
    sl   <- slice_rotation(r, n)
    rows <- vapply(sl, function(x) as.numeric(terra::nrow(x)), numeric(1))
    expect_length(sl, n)
    expect_equal(sum(rows), nr, info = paste("n =", n))
    expect_true(all(rows >= 1))
    # Slice i must start where slice i-1 ended: the first cell of each slice
    # is the first cell of the next unused row, counted from the top.
    expect_equal(vapply(sl, function(x) terra::values(x)[1], numeric(1)),
                 cumsum(c(1, head(rows, -1))), info = paste("n =", n))
  }
})

test_that("depth_zoning reads three or more consecutive bins as a range", {
  skip_if_not_installed("terra")
  img  <- terra::rast(matrix(1, 4, 4))
  bins <- terra::rast(matrix(rep(c(0, 5, 10, 15), each = 4), 4, 4, byrow = TRUE))
  kept <- function(d) sort(unique(terra::values(bins)[!is.na(terra::values(
    depth_zoning(img, depth_map = bins, depth = d)))]))

  expect_equal(kept(5), 5)
  # Evenly spaced and longer than two: an inclusive range, so 10 comes along.
  expect_equal(kept(c(5, 10, 15)), c(5, 10, 15))
  # Exactly two values are matched one by one instead, so the bin between them
  # is not kept. The length-2 case is a pair, never a range.
  expect_equal(kept(c(5, 15)), c(5, 15))
})
