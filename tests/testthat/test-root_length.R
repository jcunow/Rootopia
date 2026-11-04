# root_length: one representative input, parameter sweep over method x unit,
# plus unit-scaling consistency.

test_that("root_length runs across all method x unit combinations", {
  skip_if_not_installed("terra")
  img  <- make_binary_spatraster()
  grid <- expand.grid(
    method = c("kimura2", "kimura1", "freeman_basic", "freeman_corrected"),
    unit   = c("px", "cm", "inch"),
    stringsAsFactors = FALSE
  )
  expect_runs_on_grid(
    function(method, unit)
      root_length(img, method = method, unit = unit, dpi = 300,
                  show_messages = FALSE),
    grid,
    check = function(res, i) expect_true(is.numeric(res) && length(res) == 1)
  )
})

test_that("root_length cm equals px * 2.54/dpi", {
  skip_if_not_installed("terra")
  img <- make_binary_spatraster()
  px  <- root_length(img, unit = "px", method = "kimura2", show_messages = FALSE)
  cm  <- root_length(img, unit = "cm", dpi = 300, method = "kimura2",
                     show_messages = FALSE)
  expect_equal(cm, px * 2.54 / 300, tolerance = 1e-6)
})

test_that("root_length asks for dpi when the unit needs one", {
  skip_if_not_installed("terra")
  m <- matrix(0, 10, 10); m[3, 2:9] <- 1        # 8 px in a row: 7 orthogonal steps
  img <- terra::rast(m)
  expect_error(root_length(img, show_messages = FALSE), "dpi. is required")
  expect_error(root_length(img, dpi = 0, show_messages = FALSE), "dpi")
  expect_equal(root_length(img, dpi = 300, show_messages = FALSE), 7 * 2.54 / 300)
  expect_equal(root_length(img, unit = "px", show_messages = FALSE), 7)
})

test_that("a pixel with two connections of one kind counts both", {
  skip_if_not_installed("terra")
  # A V: five pixels, four diagonal steps, the tip holding two of them.
  v <- matrix(0, 7, 7); v[cbind(1:5, c(1, 2, 3, 2, 1))] <- 1
  # The same five pixels on a straight diagonal.
  s <- matrix(0, 7, 7); s[cbind(1:5, 1:5)] <- 1
  len <- function(m) root_length(terra::rast(m), unit = "px",
                                 method = "freeman_basic", show_messages = FALSE)
  expect_equal(len(v), 4 * sqrt(2))
  expect_equal(len(v), len(s))
  # A plus sign: the center owns its left and upper step.
  p <- matrix(0, 9, 9); p[5, 3:7] <- 1; p[3:7, 5] <- 1
  links <- .skeleton_links(terra::rast(p))
  expect_equal(terra::global(links[["orth"]], "sum")[[1]], 8)
})

test_that("multi-layer input without select_layer warns and names the layer", {
  skip_if_not_installed("terra")
  m <- matrix(0, 10, 10); m[3, 2:9] <- 1
  img <- c(terra::rast(m), terra::rast(m))
  expect_warning(root_length(img, unit = "px", show_messages = FALSE), "layer 1")
  expect_silent(root_length(img, unit = "px", select_layer = 2, show_messages = FALSE))
})
