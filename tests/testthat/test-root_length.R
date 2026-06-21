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
