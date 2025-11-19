# Medial Axis Transform (Internal)

This internal function computes the medial axis transform of a binary
image, identifying the set of skeleton points equidistant to the
object's boundaries.

## Usage

``` r
medial_axis_transform(img, verbose = FALSE, select.layer = NULL)
```

## Arguments

- img:

  A matrix, data frame, or \`SpatRaster\` object representing the binary
  image for transformation.

- verbose:

  Logical. If \`TRUE\`, outputs diagnostic information such as image
  dimensions, progress of computation, and final skeleton size. Default
  is \`FALSE\`.

- select.layer:

  Integer indicating the layer to use if \`img\` is a multi-layer
  \`SpatRaster\`. Default is 2.

## Value

A binary matrix where \`1\` represents skeleton pixels and \`0\`
represents the background.

## Details

\- The input image is first processed using
[`load_flexible_image`](https://jcunow.github.io/RootScanR/reference/load_flexible_image.md)
to ensure it is a binary matrix. - The algorithm proceeds through the
following steps: 1. \*\*Distance Transform\*\*: Computes the distance of
each foreground pixel to the nearest background pixel using a two-pass
algorithm. 2. \*\*Local Maxima Detection\*\*: Identifies local maxima in
the distance transform to mark potential skeleton points. 3.
\*\*Skeleton Refinement\*\*: Ensures connectivity by connecting skeleton
points within an 8-neighborhood. - The result is a binary image
representing the medial axis of the input object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage
raster <- terra::rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0,1), nrow = 3))
skeleton <- medial_axis_transform(raster, verbose = TRUE)
} # }
```
