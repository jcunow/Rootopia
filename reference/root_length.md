# Root length estimation from skeleton images

Root length is estimated from a skeletonized binary image using either
Freeman chain-code based estimators or Kimura estimators.

## Usage

``` r
root_length(
  img,
  unit = "cm",
  dpi = 300,
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

  Numeric. Scan resolution (dots per inch); required for cm/inch
  conversion.

- select_layer:

  Numeric. Which layer to select if `img` has multiple layers.

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
orthogonal pixel connections in the skeleton.

Freeman methods treat the skeleton as a discrete chain-code path: -
freeman_basic: L = sqrt(2) \* Nd + No

\- freeman_corrected: L = 0.948 \* (sqrt(2) \* Nd + No)

Kimura methods treat the skeleton as a discretized representation of an
underlying continuous curve and reduce orientation bias:

\- kimura1: L = sqrt(Nd^2 + (Nd + No)^2)

\- kimura2 (default): L = sqrt(Nd^2 + (Nd + No/2)^2) + No/2

The Kimura2 estimator is generally preferred due to improved stability
across object orientations and curvature distributions.
