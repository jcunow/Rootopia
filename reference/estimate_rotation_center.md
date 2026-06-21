# Estimates rotation from tape coverage

This function analyzes image data to determine rotation based on tape
coverage, assuming more tape is present on the upper side of the tube.

## Usage

``` r
estimate_rotation_center(
  img,
  tape_brightness = 0.66,
  extra_rows = 100,
  search_area = 0.45,
  tape_quantile = 0.98,
  nclasses = 3,
  select_layer = NULL
)
```

## Arguments

- img:

  Input image as raster, file name, or array

- tape_brightness:

  Brightness threshold for tape detection (0-1)

- extra_rows:

  Additional rows to add for analysis

- search_area:

  Proportion of image to analyze (0-1)

- tape_quantile:

  Quantile used to align brightness with tape (0-1)

- nclasses:

  Number of classes for pixel clustering

- select_layer:

  Integer. Specifies which layer to use if the input is a multi-band
  image. Default is \`NULL\`.

## Value

numeric Position of the center of extruding tape

## Examples

``` r
img = seg_Oulanka2023_Session01_T067
r0 = estimate_rotation_center(img)
```
