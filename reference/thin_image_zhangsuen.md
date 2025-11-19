# Thin Binary Image using Zhang-Suen Algorithm (Internal)

This internal function performs image thinning using the Zhang-Suen
thinning algorithm. It reduces binary images to their skeleton while
preserving the structure and connectivity of the foreground pixels.

## Usage

``` r
thin_image_zhangsuen(img, verbose = TRUE, select.layer = 2)
```

## Arguments

- img:

  A matrix, data frame, or \`SpatRaster\` object representing the binary
  image to be thinned.

- verbose:

  Logical. If \`TRUE\`, prints diagnostic information such as iteration
  progress and pixel removal counts. Default is \`TRUE\`.

- select.layer:

  Integer indicating the layer to use if \`img\` is a multi-layer
  \`SpatRaster\`. Default is 2.

## Value

A binary matrix representing the thinned image (skeleton).

## Details

\- The function first prepares the image using the
[`load_flexible_image`](https://jcunow.github.io/RootScanR/reference/load_flexible_image.md)
function, ensuring binary matrix format. - Thinning is performed
iteratively in two subiterations per cycle: 1. Identifying pixels to be
removed based on Zhang-Suen conditions (first subiteration). 2. Refining
removal decisions in the second subiteration. - The algorithm continues
until no pixels are removed in an iteration or a maximum number of
iterations is reached (default: 1000).

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage
raster <- terra::rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0,0), nrow = 3))
thinned_image <- thin_image_zhangsuen(raster, verbose = TRUE, select.layer = NULL)
} # }
```
