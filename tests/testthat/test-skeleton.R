# Skeletonisation + skeleton point detection.

test_that("skeletonize_image returns a binary single-layer raster", {
  skip_if_not_installed("terra")
  sk <- skeletonize_image(make_binary_spatraster(), verbose = FALSE)
  expect_s4_class(sk, "SpatRaster")
  vals <- terra::values(sk)
  expect_true(all(vals %in% c(0, 1) | is.na(vals)))
  expect_gt(sum(vals == 1, na.rm = TRUE), 0)     # something survived
})

test_that("skeletonize_image errors on an empty image", {
  skip_if_not_installed("terra")
  expect_error(skeletonize_image(terra::rast(matrix(0, 24, 28)), verbose = FALSE))
})

test_that("detect_skeleton_points returns a list", {
  skip_if_not_installed("terra")
  res <- detect_skeleton_points(make_binary_spatraster(), skeletonize = TRUE)
  expect_type(res, "list")
})


test_that("the thinning table never breaks connectivity or eats a root tip", {
  # The 256-entry table inside lut_thin_fast() is otherwise an opaque constant.
  # These are the two properties the whole package leans on: thinning must not
  # split a root, and must not shorten one by deleting its tip. Checking them
  # directly makes the table auditable without needing to trust its provenance.
  lut <- .thinning_lut
  expect_length(lut, 256L)

  # Neighbour offsets in the table's own weight layout.
  wt  <- c(TL = 1, T = 2, TR = 4, R = 8, BR = 16, B = 32, BL = 64, L = 128)
  off <- list(TL = c(-1,-1), T = c(-1,0), TR = c(-1,1), R = c(0,1),
              BR = c(1,1),   B = c(1,0),  BL = c(1,-1), L = c(0,-1))

  # Are the foreground neighbours still 8-connected to each other once the
  # centre pixel is removed?
  stays_connected <- function(code) {
    on <- names(wt)[bitwAnd(code, wt) > 0L]
    if (length(on) <= 1L) return(TRUE)
    p <- do.call(rbind, off[on])
    seen <- rep(FALSE, nrow(p)); seen[1] <- TRUE; q <- 1L; qi <- 1L
    while (qi <= length(q)) {
      v <- q[qi]; qi <- qi + 1L
      for (w in seq_len(nrow(p)))
        if (!seen[w] && max(abs(p[w, ] - p[v, ])) == 1L) { seen[w] <- TRUE; q <- c(q, w) }
    }
    all(seen)
  }

  deletes <- which(lut != 0L) - 1L
  expect_gt(length(deletes), 0L)
  expect_true(all(vapply(deletes, stays_connected, logical(1))))

  n_nb <- vapply(0:255, function(c) sum(bitwAnd(c, wt) > 0L), integer(1))
  expect_true(all(lut[n_nb <= 1L] == 0L))          # tips and isolated pixels stay
  expect_true(all(lut[n_nb >= 7L] == 0L))          # interior pixels are not eroded
})
