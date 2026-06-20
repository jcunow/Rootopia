# Stitch two scans into one mosaic with linear feather blending

Core compositing step ported from the Python `ImageStitcher`
`stitch_two_images`. `img2` is placed to the right of `img1` at the
offset estimated by an edge-based FFT phase correlation, so that
`overlap = edge_width - dx`; the overlap band is blended with a
horizontal alpha ramp (1 -\> 0) and the vertical offset is applied by
growing the canvas.

## Usage

``` r
stitch_image_pair(
  img1,
  img2,
  method = "phase",
  edge_width = 250,
  vertical_region = 1000,
  vertical_offset = 300,
  direction = "horizontal",
  preprocess = "none",
  blend = "linear",
  blend_width = NULL
)
```

## Arguments

- img1, img2:

  Image inputs (file path, `SpatRaster`, `Raster*`, `cimg`,
  `magick-image`, matrix, or array). See
  [`load_flexible_image`](https://jcunow.github.io/Rootopia/reference/load_flexible_image.md)
  for supported formats.

- method:

  Alignment method. Currently only `"phase"` (FFT phase correlation) is
  supported. `"feature"` (SIFT/ORB keypoint matching in the Python
  original) requires OpenCV and has no CRAN-compatible R backend; it
  raises an informative error.

- edge_width:

  Width in pixels of the edge band used for alignment. Clamped to the
  image width. Set it close to the true overlap - roughly 1-2x (so the
  overlap is between about half and all of `edge_width`) - for the most
  reliable peak.

- vertical_region:

  Height in pixels of the vertical band used for alignment.

- vertical_offset:

  Starting row (from the top) of the vertical band.

- direction:

  `"horizontal"` (default, frames side by side) or `"vertical"` (frames
  stacked top to bottom, e.g. minirhizotron strips acquired down the
  tube). Vertical transposes the inputs, runs the same horizontal core,
  and transposes the result back.

- preprocess:

  Preprocessing applied to the edge bands before correlation: one of
  `"none"` (default), `"center"` (subtract mean), `"norm"` (divide by
  SD), `"center_norm"`, `"hann"` (demean + Hann window), `"grad"`
  (gradient magnitude) or `"grad_norm"`. `"hann"`, `"grad"` and
  `"center"` help on scans with uneven lighting or sparse texture;
  `"norm"` alone barely changes phase correlation (it is already
  scale-invariant). On well-textured scans `"none"` is usually fine.

- blend:

  How the overlap band is combined: `"linear"` (default, alpha ramp 1
  -\> 0, good for colour scans), `"overlay"` (img2 hides img1),
  `"overlay_first"` (img1 hides img2), `"max"` (lighten / union -
  recommended for segmented/binary masks, where averaging would make
  fractional values and ghost thin roots) or `"min"` (darken).

- blend_width:

  Optional width (px) of the linear ramp, centred in the overlap (hard
  img1 to its left, hard img2 to its right); `NULL` (default) ramps
  across the whole overlap. A smaller value reduces root ghosting on
  colour scans. Ignored unless `blend = "linear"`.

## Value

A numeric `(height, width, channel)` array in 0-255. The estimated
`c(dx, dy, peak)` is attached as `attr(., "offset")` (`dx` the
horizontal shift, with `overlap = edge_width - dx`). Convert for
plotting with e.g.
[`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html)
/
[`terra::plotRGB()`](https://rspatial.github.io/terra/reference/plotRGB.html).

## See also

[`stitch_image_sequence`](https://jcunow.github.io/Rootopia/reference/stitch_image_sequence.md),
[`stitch_root_scans`](https://jcunow.github.io/Rootopia/reference/stitch_root_scans.md)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
img   <- array(runif(80 * 160 * 3) * 255, dim = c(80, 160, 3))
left  <- img[, 1:100, , drop = FALSE]
right <- img[, 81:160, , drop = FALSE]          # 20 px overlap
mosaic <- stitch_image_pair(left, right, edge_width = 30,
                            vertical_region = 80, vertical_offset = 0)
dim(mosaic)
attr(mosaic, "offset")
} # }
```
