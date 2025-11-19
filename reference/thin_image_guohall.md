# Enhanced Guo-Hall Thinning Algorithm Thin Binary Image using Guo-Hall Algorithm (Internal)

This internal function applies the Guo-Hall thinning algorithm to reduce
binary images to their skeletons while preserving connectivity and
structure.

## Usage

``` r
thin_image_guohall(img, verbose = FALSE, select.layer = 2)
```

## Arguments

- img:

  A matrix, data frame, or \`SpatRaster\` object representing the binary
  image to be thinned.

- verbose:

  Logical. If \`TRUE\`, outputs diagnostic information such as image
  dimensions, pixel removal counts, and iteration progress. Default is
  \`FALSE\`.

- select.layer:

  Integer indicating the layer to use if \`img\` is a multi-layer
  \`SpatRaster\`. Default is 2.

## Value

A binary matrix representing the thinned image (skeleton).

## Details

\- The input image is first processed using
[`load_flexible_image`](https://jcunow.github.io/RootScanR/reference/load_flexible_image.md)
to ensure it is a binary matrix. - Thinning is performed in an iterative
process consisting of two subiterations per cycle: 1. In the first
subiteration, pixels are marked for removal based on specific Guo-Hall
conditions. 2. In the second subiteration, a different set of conditions
is applied to mark additional pixels for removal. - The process
continues until no pixels are removed in an iteration or the maximum
number of iterations (default: 1000) is reached. - The Guo-Hall
algorithm ensures that the skeleton of the image is preserved.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage
raster <- terra::rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0,0), nrow = 3))
thinned_image <- thin_image_guohall(raster, verbose = TRUE)
} # }
```
