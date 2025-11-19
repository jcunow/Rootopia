# Estimates rotation from tape coverage

This function analyzes image data to determine rotation based on tape
coverage, assuming more tape is present on the upper side of the tube.

## Usage

``` r
estimate_rotation_center(
  img,
  tape.brightness = 0.66,
  extra.rows = 100,
  search.area = 0.45,
  tape.quantile = 0.98,
  nclasses = 3,
  select.layer = NULL
)
```

## Arguments

- img:

  Input image as raster, file name, or array

- tape.brightness:

  Brightness threshold for tape detection (0-1)

- extra.rows:

  Additional rows to add for analysis

- search.area:

  Proportion of image to analyze (0-1)

- tape.quantile:

  Quantile used to align brightness with tape (0-1)

- nclasses:

  Number of classes for pixel clustering

- select.layer:

  Integer. Specifies which layer to use if the input is a multi-band
  image. Default is \`NULL\`.

## Value

numeric Position of the center of extruding tape

## Examples

``` r
img = seg_Oulanka2023_Session01_T067
r0 = estimate_rotation_center(img)
```
