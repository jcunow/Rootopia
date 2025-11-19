# Detect Skeleton Points: Branching Points and Endpoints

Identifies the branching points and endpoints of a skeletonized binary
image.

## Usage

``` r
detect_skeleton_points(img, select.layer = 2)
```

## Arguments

- img:

  A matrix, data frame, or \`SpatRaster\` object representing the
  skeletonized binary image.

- select.layer:

  Integer. Specifies which layer to use if the input is a multi-band
  image. Default is \`2\`.

## Value

A named list containing two \`SpatRaster\` objects:

- `endpoints`: A binary raster where endpoints are marked as `1`.

- `branching_points`: A binary raster where branching points are marked
  as `1`.

## Details

This function detects key points in a skeletonized binary image:

- **Endpoints**: Pixels with exactly one neighbor in the skeleton.

- **Branching Points**: Pixels with more than two neighbors in the
  skeleton.

The function uses a 3x3 neighborhood kernel to count the number of
neighbors for each foreground pixel (`1`) in the image. Based on the
neighbor count, points are classified as endpoints or branching points.

The input image should be skeletonized (thin and connected) before using
this function. If not already binary, the input image will be binarized
internally.

## See also

[`skeletonize_image`](https://jcunow.github.io/RootScanR/reference/skeletonize_image.md),
[`thin_image_zhangsuen`](https://jcunow.github.io/RootScanR/reference/thin_image_zhangsuen.md),
[`thin_image_guohall`](https://jcunow.github.io/RootScanR/reference/thin_image_guohall.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(terra)

# Example skeletonized image
skeleton <- rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0), nrow = 4))

# Detect endpoints and branching points
points <- detect_skeleton_points(skeleton)

# Access results
endpoints <- points$endpoints
branching_points <- points$branching_points
} # }
```
