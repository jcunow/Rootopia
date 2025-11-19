# Calculate Image Coloration Metrics

Calculate Image Coloration Metrics

## Usage

``` r
tube_coloration(img, r = 0.2126, g = 0.7152, b = 0.0722)
```

## Arguments

- img:

  Three-band raster (RGB) or path to image

- r:

  Red channel weight

- g:

  Green channel weight

- b:

  Blue channel weight

## Value

Data frame of color metrics

## Examples

``` r
data(rgb_Oulanka2023_Session03_T067)
img = terra::rast(rgb_Oulanka2023_Session03_T067)
colorvector = tube_coloration(img)
#> Some pixels have zero intensity, which may affect color calculations
```
