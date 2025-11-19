# counts all pixels in a segmented image

counts all pixels in a segmented image

## Usage

``` r
count_pixels(img)
```

## Arguments

- img:

  one layer image

## Value

a numeric value

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
img = terra::rast(seg_Oulanka2023_Session01_T067)[[2]]
rootpixel  = count_pixels(img)
```
