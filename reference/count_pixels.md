# Count all pixels in a segmented image

Count all pixels in a segmented image

## Usage

``` r
count_pixels(img)
```

## Arguments

- img:

  A single-layer raster image (SpatRaster or compatible format).

## Value

A numeric value — the sum of all non-NA pixel values.

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
img <- terra::rast(seg_Oulanka2023_Session01_T067)[[2]]
rootpixel <- count_pixels(img)
```
