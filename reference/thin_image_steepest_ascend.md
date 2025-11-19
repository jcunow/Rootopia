# Thin binary image using steepest ascent skeletonization

This function performs skeletonization of a binary image using a
steepest ascent algorithm on the distance transform. Each foreground
pixel traces its path to the ridge (skeleton) by following the steepest
ascent gradient in the distance transform field.

## Usage

``` r
thin_image_steepest_ascend(img, verbose = TRUE, select.layer = 2)
```

## Arguments

- img:

  Input image. Can be a file path, image object, or any format supported
  by `load_flexible_image`.

- verbose:

  Logical. If `TRUE`, prints progress information including the final
  skeleton pixel count. Default is `TRUE`.

- select.layer:

  Integer. Which layer/channel to select from the input image. Default
  is 2. This parameter is passed to `load_flexible_image`.

## Value

A binary matrix of the same dimensions as the input image, where 1
indicates skeleton pixels and 0 indicates background.

## Details

The algorithm works by:

1.  Computing the Euclidean distance transform of the binary image

2.  For each foreground pixel, tracing the steepest ascent path in the
    distance transform until reaching a local maximum (ridge point)

3.  Marking all ridge points as skeleton pixels

The function uses padding to handle edge cases during the steepest
ascent tracing process.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load and skeletonize an image
raster <- terra::rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0,0), nrow = 3))
skeleton <- thin_image_steepest_ascend(raster)

# Visualize the result
plot(as.cimg(skeleton))
} # }
```
