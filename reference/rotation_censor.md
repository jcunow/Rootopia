# Censor image edges based on rotation

Crops the rotation axis (rows) of a root scan. In fixed mode it returns
a fixed-width window centred on a given row (e.g. the rotation centre
from
[`estimate_rotation_center()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_center.md));
in variable mode it trims a band sized by a measured offset. Optionally
previews what is kept versus cut.

## Usage

``` r
rotation_censor(
  img,
  center.offset = 0,
  cut.buffer = 0.02,
  fixed.rotation = TRUE,
  fixed.width = 500,
  select.layer = NULL,
  overlay = FALSE,
  ...
)
```

## Arguments

- img:

  Input image as raster, file name, or array.

- center.offset:

  Numeric. Meaning depends on `fixed.rotation`: when `TRUE`, the row to
  centre the kept window on (an absolute row, e.g. from
  [`estimate_rotation_center()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_center.md));
  when `FALSE`, the rotation shift in rows to trim (e.g. from
  [`estimate_rotation_shift()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_shift.md)).

- cut.buffer:

  Extra proportion of the rotation axis to trim (variable mode).

- fixed.rotation:

  Logical. If `TRUE`, return a fixed-width window.

- fixed.width:

  Output width in rows when `fixed.rotation = TRUE`.

- select.layer:

  Integer or `NULL`. Layer to use for multi-band inputs.

- overlay:

  Logical. If `TRUE`, plot the full image with the kept window (green
  outline) and discarded margins (red shading) before cropping. Default
  `FALSE`.

- ...:

  Passed to the underlying `terra` plotting call when `overlay = TRUE`.

## Value

A cropped `SpatRaster` (returned invisibly when `overlay = TRUE`), or
`NULL` if the result is empty.

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
img <- terra::rast(seg_Oulanka2023_Session01_T067)
r0  <- estimate_rotation_center(img)
rotation_censor(img, center.offset = r0, fixed.width = 800,
                fixed.rotation = TRUE, overlay = TRUE)
```
