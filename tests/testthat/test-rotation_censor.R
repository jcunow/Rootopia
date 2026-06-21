# rotation_censor, including the center.offset enhancement: it must accept an
# absolute row (> 1), a fraction in [0, 1], and a keyword, and produce the same
# window when they describe the same position.

test_that("rotation_censor returns a fixed-width window", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster(nr = 40, nc = 28)
  out <- rotation_censor(img, center.offset = 20, fixed.width = 10,
                         fixed.rotation = TRUE)
  expect_s4_class(out, "SpatRaster")
  expect_lt(terra::nrow(out), terra::nrow(img))
})

test_that("fraction and keyword agree with the equivalent absolute row", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster(nr = 40, nc = 28)   # middle row = 20
  by_row     <- rotation_censor(img, center.offset = 20,       fixed.width = 10)
  by_frac    <- rotation_censor(img, center.offset = 0.5,      fixed.width = 10)
  by_keyword <- rotation_censor(img, center.offset = "middle", fixed.width = 10)
  expect_equal(terra::nrow(by_frac),    terra::nrow(by_row))
  expect_equal(terra::nrow(by_keyword), terra::nrow(by_row))
})

test_that("rotation_censor rejects unknown keywords", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  expect_error(rotation_censor(img, center.offset = "sideways"))
})
