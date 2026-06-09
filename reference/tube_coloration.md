# Calculate Image Coloration Metrics

Calculate Image Coloration Metrics

## Usage

``` r
tube_coloration(img, r = 0.2126, g = 0.7152, b = 0.0722)
```

## Arguments

- img:

  Three-band raster (RGB) or path to image.

- r:

  Red channel luminosity weight. Default follows ITU-R BT.709.

- g:

  Green channel luminosity weight.

- b:

  Blue channel luminosity weight.

## Value

A data frame with columns: rcc, gcc, bcc, hue, saturation, luminosity,
red, green, blue.

## Examples

``` r
data(rgb_Oulanka2023_Session03_T067)
img <- terra::rast(rgb_Oulanka2023_Session03_T067)
colorvector <- tube_coloration(img)
#> Some pixels have zero intensity, which may affect colour calculations
```
