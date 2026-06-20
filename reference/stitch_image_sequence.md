# Sequentially stitch a sequence of scans into one mosaic

Frames are stitched in the given order: the first pair is composited,
then each subsequent frame is stitched onto the growing mosaic, in
memory (no temporary PNG round-trips).

## Usage

``` r
stitch_image_sequence(
  images,
  method = "phase",
  edge_width = 250,
  vertical_region = 1000,
  vertical_offset = 300,
  direction = "horizontal",
  preprocess = "none",
  blend = "linear",
  blend_width = NULL,
  return_offsets = FALSE
)
```

## Arguments

- images:

  A character vector of image file paths, or a list of image objects
  accepted by
  [`load_flexible_image`](https://jcunow.github.io/Rootopia/reference/load_flexible_image.md).
  Stitched in the order given, so sort beforehand if needed.

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

- return_offsets:

  Logical. If `TRUE`, return a list with the mosaic and a per-step
  performance data frame (`step`, `dx`, `dy`, `peak`, `overlap`) instead
  of just the mosaic.

## Value

A numeric `(height, width, channel)` array (0-255); or, when
`return_offsets = TRUE`, a list `list(mosaic, offsets)`.

## See also

[`stitch_image_pair`](https://jcunow.github.io/Rootopia/reference/stitch_image_pair.md),
[`stitch_root_scans`](https://jcunow.github.io/Rootopia/reference/stitch_root_scans.md)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
img <- array(runif(70 * 200 * 3) * 255, dim = c(70, 200, 3))
frames <- list(img[, 1:90, , drop = FALSE],
               img[, 71:160, , drop = FALSE],
               img[, 141:200, , drop = FALSE])
res <- stitch_image_sequence(frames, edge_width = 25, vertical_region = 70,
                             vertical_offset = 0, return_offsets = TRUE)
dim(res$mosaic)
res$offsets
} # }
```
