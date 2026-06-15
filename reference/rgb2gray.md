# Convert RGB image to grayscale with optimized memory management and parallel processing

Convert RGB image to grayscale with optimized memory management and
parallel processing

## Usage

``` r
rgb2gray(img, r = 0.21, g = 0.72, b = 0.07)
```

## Arguments

- img:

  SpatRaster RGB image

- r:

  Weight for red channel

- g:

  Weight for green channel

- b:

  Weight for blue channel

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
img = seg_Oulanka2023_Session01_T067
gray.raster = rgb2gray(img)
#> Error in rgb2gray(img): could not find function "rgb2gray"
```
