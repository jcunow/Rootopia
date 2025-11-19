# Skeletonization Wrapper Function

This function serves as a wrapper for applying different skeletonization
methods to a binary image, including the Zhang-Suen, Guo-Hall, Medial
Axis Transform (MAT), and Steepest Ascent algorithms.

## Usage

``` r
skeletonize_image(
  img,
  methods = c("ZhangSuen", "GuoHall", "MAT", "SteepestAscend"),
  verbose = TRUE,
  select.layer = NULL
)
```

## Arguments

- img:

  A matrix, data frame, or \`SpatRaster\` object representing the binary
  image to be skeletonized.

- methods:

  A character vector specifying the skeletonization methods to apply.
  Valid options are `"ZhangSuen"`, `"GuoHall"`, `"MAT"`, and
  `"SteepestAscend"`. Defaults to all four methods.

- verbose:

  Logical. If `TRUE`, displays progress and diagnostic messages during
  processing. Defaults to `TRUE`.

- select.layer:

  Integer specifying the layer to use if `img` is a multi-layer
  \`SpatRaster\`. Defaults to `NULL`, which may use package-specific
  defaults for each method.

## Value

If a single method is selected, the function returns a \`SpatRaster\`
object representing the skeletonized image. If multiple methods are
selected, a named list of \`SpatRaster\` objects is returned, where each
element is named according to the method used.

## Details

This function allows for flexible and streamlined skeletonization of
binary images using one or more supported algorithms:

- `"ZhangSuen"`: Implements the Zhang-Suen thinning algorithm, a
  parallel iterative method that preserves connectivity while reducing
  binary objects to their skeletal representation.

- `"GuoHall"`: Implements the Guo-Hall thinning algorithm, an improved
  parallel thinning method that often produces cleaner skeletons with
  better preservation of shape characteristics.

- `"MAT"`: Computes the Medial Axis Transform to extract the skeleton
  based on the distance transform of the binary image.

- `"SteepestAscend"`: Uses a steepest ascent algorithm on the distance
  transform, where each foreground pixel traces its path to ridge points
  (skeleton) by following the steepest gradient.

The function processes the input image with the specified methods and
returns the results. If multiple methods are chosen, the results are
returned as a named list, with each element corresponding to a method.
Each method may have different computational complexity and produce
slightly different skeletal representations.

## Note

Different skeletonization methods may produce varying results depending
on the input image characteristics. The Zhang-Suen and Guo-Hall methods
are iterative thinning approaches, while MAT and SteepestAscend are
based on distance transforms. Consider the specific requirements of your
analysis when choosing methods.

## See also

[`thin_image_zhangsuen`](https://jcunow.github.io/RootScanR/reference/thin_image_zhangsuen.md),
[`thin_image_guohall`](https://jcunow.github.io/RootScanR/reference/thin_image_guohall.md),
[`medial_axis_transform`](https://jcunow.github.io/RootScanR/reference/medial_axis_transform.md),
[`thin_image_steepest_ascend`](https://jcunow.github.io/RootScanR/reference/thin_image_steepest_ascend.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load a binary image as a SpatRaster
binary_image <- terra::rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0, 0), nrow = 3))

# Apply all skeletonization methods
skeletons <- skeletonize_image(binary_image, verbose = TRUE)

# Apply only Zhang-Suen method
zhang_skeleton <- skeletonize_image(binary_image, methods = "ZhangSuen")

# Apply only Zhang-Suen method
SA_skeleton <- skeletonize_image(binary_image, methods = "SteepestAscend")

# Apply multiple specific methods
selected_skeletons <- skeletonize_image(binary_image, 
                                       methods = c("ZhangSuen", "MAT"),
                                       verbose = FALSE)

# Access results from multiple methods
zhang_result <- selected_skeletons$ZhangSuen
mat_result <- selected_skeletons$MAT
} # }
```
