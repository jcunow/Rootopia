# clean_image: hole filling / artifact removal. Output-format sweep + a
# behavioral check that an isolated speck is removed and an enclosed hole filled.

test_that("clean_image returns each requested output format", {
  skip_if_not_installed("terra")
  skip_if_not_installed("imager")
  img <- make_binary_spatraster()
  expect_s4_class(clean_image(img, output_format = "spatrast"), "SpatRaster")
  expect_s3_class(clean_image(img, output_format = "cimg"), "cimg")
  expect_true(is.matrix(clean_image(img, output_format = "matrix")))
})

test_that("clean_image removes a small isolated artifact", {
  skip_if_not_installed("terra")
  skip_if_not_installed("imager")
  m <- matrix(0, 30, 30)
  m[10:20, 10:20] <- 1     # solid block (touches nothing else)
  m[3, 3] <- 1             # isolated 1-px artifact
  cleaned <- clean_image(terra::rast(m), drop_objects_upto = 5,
                         output_format = "matrix")
  expect_equal(cleaned[3, 3], 0)
})

test_that("deprecated threshold names still work and warn", {
  skip_if_not_installed("terra")
  skip_if_not_installed("imager")
  m <- matrix(0, 30, 30)
  m[10:20, 10:20] <- 1
  m[3, 3] <- 1
  expect_warning(
    old <- clean_image(terra::rast(m), max_artifact_size = 5,
                       output_format = "matrix"),
    "deprecated")
  new <- clean_image(terra::rast(m), drop_objects_upto = 5,
                     output_format = "matrix")
  expect_equal(old, new)
})

test_that("mm2 thresholds convert to the equivalent pixel thresholds", {
  skip_if_not_installed("terra")
  skip_if_not_installed("imager")
  m <- matrix(0, 30, 30)
  m[10:20, 10:20] <- 1
  m[3, 3] <- 1
  dpi <- 300
  # 5 px at 300 dpi, expressed in mm^2
  mm2 <- 5 / (dpi / 25.4)^2
  a <- clean_image(terra::rast(m), drop_objects_upto = mm2,
                   size_unit = "mm2", dpi = dpi, output_format = "matrix")
  b <- clean_image(terra::rast(m), drop_objects_upto = 5,
                   output_format = "matrix")
  expect_equal(a, b)
  expect_error(clean_image(terra::rast(m), drop_objects_upto = 1,
                           size_unit = "mm2"), "dpi")
})


# ---------------------------------------------------------------------------
# Size thresholds. The holes and specks below are solid, well separated and
# well inside the frame, so their pixel counts are the same under 4- and
# 8-connected labelling and none of them touches the border. The expected
# foreground counts are therefore arithmetic.
# ---------------------------------------------------------------------------

.holey_mask <- function() {
  m <- matrix(0, 24, 24)
  m[4:21, 4:21] <- 1          # 18 x 18 = 324 foreground
  m[7:8,   7:8]  <- 0         # a 2 x 2 hole  (4 px)
  m[13:16, 13:16] <- 0        # a 4 x 4 hole (16 px)
  m
}

test_that("fill_holes_upto fills the holes at or below the threshold only", {
  skip_if_not_installed("imager")
  skip_if_not_installed("terra")
  m    <- .holey_mask()
  base <- sum(m)                                   # 324 - 4 - 16 = 304
  fg   <- function(...) sum(terra::values(
    clean_image(terra::rast(m), output_format = "spatrast", ...)) == 1, na.rm = TRUE)

  expect_equal(base, 304)
  # 0 is "fill none", not "no limit": sizes <= 0 matches no hole.
  expect_equal(fg(fill_holes_upto = 0,  drop_objects_upto = 0), base)
  # Only the 4 px hole is at or below 4.
  expect_equal(fg(fill_holes_upto = 4,  drop_objects_upto = 0), base + 4)
  # Below the smaller hole, nothing qualifies.
  expect_equal(fg(fill_holes_upto = 3,  drop_objects_upto = 0), base)
  # Above both, the whole rectangle comes back.
  expect_equal(fg(fill_holes_upto = 20, drop_objects_upto = 0), base + 20)
  # NULL is the no-limit case, so it matches the threshold that covers both.
  expect_equal(fg(fill_holes_upto = NULL, drop_objects_upto = 0), base + 20)
})

test_that("mm2 thresholds convert to the same pixel counts at a given dpi", {
  skip_if_not_installed("imager")
  skip_if_not_installed("terra")
  m <- .holey_mask()
  fg <- function(...) sum(terra::values(
    clean_image(terra::rast(m), output_format = "spatrast", ...)) == 1, na.rm = TRUE)
  # One pixel is (2.54/dpi)^2 cm^2; at dpi = 25.4 a pixel is exactly 1 mm^2,
  # so a 4 mm^2 threshold must behave exactly like a 4 px one.
  expect_equal(fg(fill_holes_upto = 4, drop_objects_upto = 0,
                  size_unit = "mm2", dpi = 25.4),
               fg(fill_holes_upto = 4, drop_objects_upto = 0))
  expect_error(clean_image(terra::rast(m), fill_holes_upto = 4, size_unit = "mm2"),
               "dpi")
})

test_that("protect_border keeps an edge speck that the size test would remove", {
  skip_if_not_installed("imager")
  skip_if_not_installed("terra")
  m <- matrix(0, 20, 20)
  m[8:13, 8:13] <- 1      # 6 x 6 = 36 px, kept at any threshold below 36
  m[3, 16]      <- 1      # a single interior speck
  m[1, 10]      <- 1      # a single speck on the top border
  fg <- function(...) sum(terra::values(
    clean_image(terra::rast(m), output_format = "spatrast", ...)) == 1, na.rm = TRUE)

  expect_equal(sum(m), 38)
  # Being at the edge confers no protection by default, so both specks go.
  expect_equal(fg(drop_objects_upto = 1, protect_border = FALSE), 36)
  # With protection the border speck is never a candidate.
  expect_equal(fg(drop_objects_upto = 1, protect_border = TRUE), 37)
})

test_that("create_kernel builds the documented shapes", {
  skip_if_not_installed("imager")
  # Asserted as sums and symmetry, which hold whichever way as.cimg() orders
  # the two axes. A disk of radius r keeps cells with sqrt(di^2 + dj^2) <= r,
  # a diamond keeps |di| + |dj| <= r, a square keeps everything.
  ksum <- function(...) sum(as.array(create_kernel(...)))

  expect_equal(ksum("square", 3), 9)     # 3 x 3, all of it
  expect_equal(ksum("diamond", 3), 5)    # centre + 4 orthogonal
  expect_equal(ksum("disk", 3), 5)       # diagonals are sqrt(2) > 1
  expect_equal(ksum("square", 5), 25)
  expect_equal(ksum("diamond", 5), 13)   # 1 + 4 + 8
  expect_equal(ksum("disk", 5), 13)      # 1 + 4 at 1 + 4 at sqrt(2) + 4 at 2

  # An even size is incremented, so it matches the next odd one.
  expect_equal(ksum("disk", 4), ksum("disk", 5))
  expect_error(create_kernel("hexagon", 3), "Invalid kernel shape")
})
