# Create a buffer halo) around non-zero pixels

Create a buffer halo) around non-zero pixels

## Usage

``` r
create_root_buffer(img, width = 2, halo.only = TRUE, kernel = "circle")
```

## Arguments

- img:

  SpatRaster/matrix/array - segmented image

- width:

  numeric - buffer width in pixels (default: 2)

- halo.only:

  logical - if TRUE, returns only the buffer zone (default: TRUE)

- kernel:

  character - shape of the thickening kernel: "circle" or "diamond"

## Value

SpatRast - buffer zone around non-zero pixels

## Examples

``` r
data(seg_Oulanka2023_Session03_T067)
img <- terra::rast(seg_Oulanka2023_Session03_T067)
create_root_buffer(img, width = 2)
#> class       : SpatRaster
#> size        : 1161, 4900, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 4900, 0, 1161  (xmin, xmax, ymin, ymax)
#> coord. ref. : 
#> source(s)   : memory
#> names       : lyr.1, lyr.2, lyr.3
#> min values  :  -254,  -254,  -254
#> max values  :     1,     1,     1
```
