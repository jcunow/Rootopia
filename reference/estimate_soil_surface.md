# Estimate soil surface position using tape markers

Detects the soil surface by analyzing tape coverage patterns in the
image.

## Usage

``` r
estimate_soil_surface(
  img,
  search_area = 0.45,
  tape_thresh = 0.33,
  dpi = 150,
  nclasses = 3,
  inverse = FALSE,
  tape_overlap = 0.5,
  tape_brightness = 0.6,
  extra_rows = 100,
  tape_quantile = 0.98,
  select_layer = NULL
)
```

## Arguments

- img:

  Input image (raster, filename, or array)

- search_area:

  Proportion of image to analyze

- tape_thresh:

  Minimum tape coverage ratio

- dpi:

  Image resolution

- nclasses:

  Number of clustering classes

- inverse:

  Invert detection for dark markers

- tape_overlap:

  Safety margin for tape (cm)

- tape_brightness:

  Brightness threshold for tape

- extra_rows:

  Additional analysis rows

- tape_quantile:

  Brightness alignment quantile

- select_layer:

  Integer. Specifies which layer to use if the input is a multi-band
  image. Default is \`NULL\`.

## Value

data.frame with soil surface and tape end positions

## Examples

``` r
img = rgb_Oulanka2023_Session03_T067
Soil0Estimates = estimate_soil_surface(img)
```
