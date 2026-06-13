# Skeletonize binary image

Applies iterative morphological thinning using a LUT-based
implementation of the Zhang–Suen skeletonization algorithm.

## Usage

``` r
skeletonize_image(
  img,
  verbose = TRUE,
  select.layer = NULL,
  overlay_png_path = NULL
)
```

## Arguments

- img:

  Input image (SpatRaster, matrix, or compatible format)

- verbose:

  Logical. Print summary statistics

- select.layer:

  Layer index for multi-layer rasters

- overlay_png_path:

  Optional file path for saving overlay visualization

## Value

Binary SpatRaster representing skeletonized image

## Details

Processing steps:

1\. Input is converted to a single-layer binary SpatRaster using
\`load_flexible_image()\`. 2. Foreground pixel count is computed. 3.
Skeletonization is performed using \`lut_thin_fast()\`: - iterative
removal of pixels based on 3x3 neighbourhood codes - lookup table
determines pixel deletions in two sub-steps per iteration 4. Output is
the final thinned binary raster. 5. Optionally, an overlay image is
generated: - original image marked as base layer - skeleton pixels
overlaid in a separate class - saved using base R PNG plotting

## Examples

``` r
library(terra)
#> terra 1.9.27

# Load example binary segmentation
data(seg_Oulanka2023_Session01_T067)

# Ensure single-layer raster if needed
img <- seg_Oulanka2023_Session01_T067

# Skeletonize
skel <- skeletonize_image(img, select.layer = 2, verbose = FALSE)

# Visual check
if (FALSE) { # \dontrun{
skeletonize_image(img, select.layer = 2, overlay_png_path = "overlay.png")
} # }
```
