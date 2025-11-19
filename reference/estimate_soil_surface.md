# Estimate soil surface position using tape markers

Detects the soil surface by analyzing tape coverage patterns in the
image.

## Usage

``` r
estimate_soil_surface(
  img,
  search.area = 0.45,
  tape.tresh = 0.33,
  dpi = 150,
  nclasses = 3,
  inverse = FALSE,
  tape.overlap = 0.5,
  tape.brightness = 0.6,
  extra.rows = 100,
  tape.quantile = 0.98,
  select.layer = NULL
)
```

## Arguments

- img:

  Input image (raster, filename, or array)

- search.area:

  Proportion of image to analyze

- tape.tresh:

  Minimum tape coverage ratio

- dpi:

  Image resolution

- nclasses:

  Number of clustering classes

- inverse:

  Invert detection for dark markers

- tape.overlap:

  Safety margin for tape (cm)

- tape.brightness:

  Brightness threshold for tape

- extra.rows:

  Additional analysis rows

- tape.quantile:

  Brightness alignment quantile

- select.layer:

  Integer. Specifies which layer to use if the input is a multi-band
  image. Default is \`NULL\`.

## Value

data.frame with soil surface and tape end positions

## Examples

``` r
img = rgb_Oulanka2023_Session03_T067
Soil0Estimates = estimate_soil_surface(img)
```
