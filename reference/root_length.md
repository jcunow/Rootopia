# Root length estimation from skeleton images

Root length is estimated from a skeletonized binary image using either
Freeman chain-code based estimators or Kimura estimators.

## Usage

``` r
root_length(
  img,
  unit = "cm",
  dpi = NULL,
  select_layer = NULL,
  method = c("kimura2", "kimura1", "freeman_basic", "freeman_corrected"),
  show_messages = TRUE,
  skeletonize = FALSE
)
```

## Arguments

- img:

  Skeletonized binary raster image. If `skeletonize = TRUE`, a segmented
  (non-skeleton) mask can be supplied instead.

- unit:

  Character. Output unit: `"cm"`, `"inch"`, or `"px"`.

- dpi:

  Numeric. Scan resolution (dots per inch). Required when `unit` is
  `"cm"` or `"inch"`; not used for `"px"`.

- select_layer:

  Numeric. Which layer to measure if `img` has multiple layers. Without
  it a multi-layer image is measured on layer 1, with a warning: in the
  bundled RootDetector images layer 1 also carries the tape, and the
  skeleton is layer 2.

- method:

  Character. One of: "freeman_basic", "freeman_corrected", "kimura1",
  "kimura2"

- show_messages:

  Logical. If `TRUE`, prints informational messages during processing.

- skeletonize:

  Logical. If `TRUE`, `img` is treated as a segmented mask and reduced
  to a skeleton internally via
  [`skeletonize_image()`](https://jcunow.github.io/Rootopia/reference/skeletonize_image.md)
  before computing length. Default `FALSE` (assumes `img` is already a
  skeleton).

## Value

Root length in pixels or converted units

## Details

Let Nd be the number of diagonal pixel connections and No the number of
orthogonal pixel connections in the skeleton. Every connection between
two 8-connected skeleton pixels is counted once, so a pixel at a bend or
a junction contributes each of its connections.

Freeman methods treat the skeleton as a discrete chain-code path:

- freeman_basic: L = sqrt(2) \* Nd + No

- freeman_corrected: L = 0.948 \* (sqrt(2) \* Nd + No)

Kimura methods treat the skeleton as a discretized representation of an
underlying continuous curve and reduce orientation bias:

- kimura1: L = sqrt(Nd^2 + (Nd + No)^2)

- kimura2 (default): L = sqrt(Nd^2 + (Nd + No/2)^2) + No/2

The Kimura2 estimator is generally preferred due to improved stability
across object orientations and curvature distributions.

## Examples

``` r
data(skl_Oulanka2023_Session01_T067)
skl <- terra::rast(skl_Oulanka2023_Session01_T067)
root_length(skl, unit = "cm", dpi = 150, select_layer = 2)
#> Diagonal: 14811 | Orthogonal: 41399
#> [1] 1002.029
```
