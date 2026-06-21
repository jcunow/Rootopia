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

  Numeric or character. Where to centre the kept window (in fixed mode),
  given in one of three forms:

  - **Absolute row** - a number `> 1`: the exact row to centre on (e.g.
    from
    [`estimate_rotation_center()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_center.md)).

  - **Fraction** - a number in `[0, 1]`: a fraction of the image height.
    `0` = top, `0.25` = a quarter down, `0.5` = middle, `1` = bottom (so
    the centre row is `center.offset * nrow`).

  - **Keyword** - `"top"` (= 0), `"middle"` / `"center"` / `"centre"` (=
    0.5), or `"bottom"` (= 1).

  The default `0` centres on the top row. When `fixed.rotation = FALSE`
  the resolved value is used as the rotation shift in rows to trim (e.g.
  from
  [`estimate_rotation_shift()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_shift.md));
  pass an absolute number there.

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


# Same window centred on the middle of the image, via a fraction or keyword
rotation_censor(img, center.offset = 0.5,      fixed.width = 800)
#> class       : SpatRaster
#> size        : 800, 4900, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 4900, 0, 800  (xmin, xmax, ymin, ymax)
#> coord. ref. : 
#> source(s)   : memory
#> names       : lyr.1, lyr.2, lyr.3
#> min values  :     0,     0,     0
#> max values  :   255,   255,   255
rotation_censor(img, center.offset = "middle", fixed.width = 800)
#> class       : SpatRaster
#> size        : 800, 4900, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 4900, 0, 800  (xmin, xmax, ymin, ymax)
#> coord. ref. : 
#> source(s)   : memory
#> names       : lyr.1, lyr.2, lyr.3
#> min values  :     0,     0,     0
#> max values  :   255,   255,   255
# Top of the tube a quarter of the way down
rotation_censor(img, center.offset = 0.25, fixed.width = 800)
#> fixed.width = 800 cannot be centred symmetrically on row 286 (image is 1144 rows). Max symmetric width here is 570 px; clamping to image bounds.
#> class       : SpatRaster
#> size        : 685, 4900, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 4900, 0, 685  (xmin, xmax, ymin, ymax)
#> coord. ref. : 
#> source(s)   : memory
#> names       : lyr.1, lyr.2, lyr.3
#> min values  :     0,     0,     0
#> max values  :   255,   255,   255
```
