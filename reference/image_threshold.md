# Threshold or deblur an image to binarize features

Applies global or adaptive thresholding to a raster image, with optional
masking and simple deblurring based on structural layer separation.

## Usage

``` r
image_threshold(
  img,
  threshold = 0.4,
  method = "global",
  window_size = 15,
  select.layer = 2,
  mask.layer = 1,
  binary_01 = FALSE,
  deblur = FALSE
)
```

## Arguments

- img:

  SpatRaster object

- threshold:

  Numeric value (0–1). For global thresholding, it's a fraction of the
  global max. For adaptive thresholding, it's a fraction of local mean.
  For deblurring, it's the fraction of max used to recover structure.

- method:

  "global" or "adaptive". Ignored if \`deblur = TRUE\`.

- window_size:

  Integer (odd), only used for adaptive thresholding.

- select.layer:

  Integer or NULL. Which layer to use for thresholding or deblurring. If
  NULL and multilayer, the mean of all layers is used.

- mask.layer:

  Integer or NULL. If set, used to preserve masked regions or enhance
  structure.

- binary_01:

  Logical. If TRUE, binarized output uses 0/1. If FALSE, it retains max
  value of input.

- deblur:

  Logical. If TRUE, applies deblurring logic instead of standard
  thresholding.

## Value

SpatRaster object

## Examples

``` r
img = terra::rast(seg_Oulanka2023_Session03_T067)
image_threshold(img, threshold = 0.3, method = "global")
#> class       : SpatRaster 
#> size        : 1161, 4900, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 4900, 0, 1161  (xmin, xmax, ymin, ymax)
#> coord. ref. :  
#> source(s)   : memory
#> names       : lyr.1, lyr.2, lyr.2 
#> min values  :     0,     0,     0 
#> max values  :   255,   255,   255 
image_threshold(img, threshold = 0.9, method = "adaptive", window_size = 15, binary_01 = TRUE)
#> class       : SpatRaster 
#> size        : 1161, 4900, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 4900, 0, 1161  (xmin, xmax, ymin, ymax)
#> coord. ref. :  
#> source(s)   : memory
#> names       : lyr.1, lyr.2, lyr.2 
#> min values  :     0,     0,     0 
#> max values  :     1,     1,     1 
image_threshold(img, threshold = 0.4, select.layer = 2, mask.layer = 1, deblur = TRUE)
#> class       : SpatRaster 
#> size        : 1161, 4900, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 4900, 0, 1161  (xmin, xmax, ymin, ymax)
#> coord. ref. :  
#> source(s)   : memory
#> names       : lyr.1, lyr.2, lyr.2 
#> min values  :     0,     0,     0 
#> max values  :   255,   255,   255 
```
