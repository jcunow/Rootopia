# ---------------------------------------------------------------------------
# Shared test fixtures + generic runners.
#
# Files named helper-*.R are sourced automatically by testthat before any
# test-*.R, so everything here is available in every test file.
#
# The philosophy (see also the package's testing notes):
#   * load_flexible_image() is the input "seam" -- test the FULL input-type
#     matrix there, and at the handful of functions that do their OWN type
#     branching. Everything else gets ONE representative input.
#   * expect_runs_on_all_inputs() / expect_runs_on_grid() make the
#     table-driven parts one-liners.
# ---------------------------------------------------------------------------

# A small synthetic binary "root" image (rows x cols), background 0, root = 1:
# one diagonal main axis plus a short lateral, enough to skeletonise/trace.
make_binary_matrix <- function(nr = 24, nc = 28) {
  m <- matrix(0, nr, nc)
  for (i in 3:18)  m[i, i] <- 1                  # main diagonal root
  for (k in 0:5)   m[10 + k, 12 - k] <- 1        # a lateral off the main axis
  m
}

# Same content as a 3-band (RGB) array, values 0/1.
make_rgb_array <- function(nr = 24, nc = 28) {
  m <- make_binary_matrix(nr, nc)
  array(c(m, m, m), dim = c(nr, nc, 3))
}

# A single-layer binary SpatRaster (the canonical "representative" input).
make_binary_spatraster <- function(nr = 24, nc = 28) {
  skip_if_not_installed("terra")
  terra::rast(make_binary_matrix(nr, nc))
}

# Build the SAME logical image in every input type load_flexible_image()
# supports. Entries that need an optional package (or that don't apply for the
# requested band count) are omitted, so callers can simply iterate the list.
#
# bands = 1 -> matrix + 1-band array/raster; bands = 3 -> RGB array/raster/cimg.
make_all_inputs <- function(nr = 24, nc = 28, bands = 3, include_paths = TRUE) {
  arr <- if (bands == 1) make_binary_matrix(nr, nc) else make_rgb_array(nr, nc)
  out <- list(array = arr)
  if (bands == 1) out$matrix <- make_binary_matrix(nr, nc)

  if (requireNamespace("terra", quietly = TRUE)) {
    out$spatraster <- terra::rast(arr)
  }
  # Plain TIFF (tiff::writeTIFF) avoids the geo-tag read warnings a terra-written
  # GeoTIFF would emit; values must be in [0, 1], which our 0/1 arrays already are.
  if (include_paths && requireNamespace("tiff", quietly = TRUE)) {
    tif <- tempfile(fileext = ".tif")
    suppressWarnings(tiff::writeTIFF(arr, tif))
    out$tif_path <- tif
  }
  if (requireNamespace("imager", quietly = TRUE)) {
    out$cimg <- suppressWarnings(imager::as.cimg(arr))
    if (include_paths) {
      png <- tempfile(fileext = ".png")
      suppressWarnings(imager::save.image(imager::as.cimg(arr), png))
      out$png_path <- png
    }
  }
  out
}

# Assert `fn(input, ...)` runs without error for EVERY input type, with the
# offending type named on failure. `check(res, type_name)` runs extra
# assertions on each successful result (optional).
expect_runs_on_all_inputs <- function(fn, ..., check = NULL, inputs = make_all_inputs()) {
  for (nm in names(inputs)) {
    inp <- inputs[[nm]]
    if (is.null(inp)) next
    res <- tryCatch(fn(inp, ...), error = function(e) e)
    testthat::expect_false(
      inherits(res, "error"),
      info = sprintf("input type '%s' errored: %s", nm,
                     if (inherits(res, "error")) conditionMessage(res) else "")
    )
    if (!is.null(check) && !inherits(res, "error")) check(res, nm)
  }
}

# Run `fn` across every row of a parameter `grid` (a data.frame from
# expand.grid). `fixed` supplies args constant across rows. NA cells are
# dropped so the function's own default applies (handy for select.layer = NA).
# Factors are coerced to character. The failing row is reported on error.
expect_runs_on_grid <- function(fn, grid, fixed = list(), check = NULL) {
  for (i in seq_len(nrow(grid))) {
    row  <- as.list(grid[i, , drop = FALSE])
    row  <- lapply(row, function(x) if (is.factor(x)) as.character(x) else x)
    row  <- row[!vapply(row, function(x) length(x) == 1L && is.na(x), logical(1))]
    args <- c(row, fixed)
    res  <- tryCatch(do.call(fn, args), error = function(e) e)
    testthat::expect_false(
      inherits(res, "error"),
      info = paste0("grid row ", i, " (",
                    paste(names(grid), unlist(grid[i, ]), sep = "=", collapse = ", "),
                    ") errored: ",
                    if (inherits(res, "error")) conditionMessage(res) else "")
    )
    if (!is.null(check) && !inherits(res, "error")) check(res, i)
  }
}
