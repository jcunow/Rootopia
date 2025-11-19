# Censor image edges based on rotation

Crops image edges to handle non-overlapping regions between sequential
scans.

## Usage

``` r
rotation_censor(
  img,
  center.offset = 0,
  cut.buffer = 0.02,
  fixed.rotation = TRUE,
  fixed.width = 500,
  select.layer = NULL
)
```

## Arguments

- img:

  Input image to censor

- center.offset:

  Rotation shift in rows (from estimate_rotation_shift())

- cut.buffer:

  Proportion of image to cut when fixed_rotation=FALSE

- fixed.rotation:

  Use fixed output dimensions

- fixed.width:

  Output width when fixed_rotation=TRUE

- select.layer:

  Integer. Specifies which layer to use if the input is a multi-band
  image. Default is \`NULL\`.

## Value

Cropped raster image

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
img = terra::rast(seg_Oulanka2023_Session01_T067)
censored.raster = rotation_censor(img,
                         center.offset = 120,
                         cut.buffer = 0.02,
                         fixed.rotation = FALSE)
                         
censored.raster = rotation_censor(img,
                         center.offset = 220,
                         cut.buffer = 0.02,
                         fixed.width = 1000,
                         fixed.rotation = TRUE)
#> New image dimension: 720 is smaller than specified fixed.width: 1000. Too strong offset for this fixed.width. Consider adjusting the fixed.width.
#> Max. possible fixed.width with this offset is: 1848(+-1 rounding error)
```
