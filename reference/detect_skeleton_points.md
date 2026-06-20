# Detect endpoints and branching points in a skeleton image

Computes local connectivity of each foreground pixel using an
8-neighbourhood (Moore neighbourhood).

## Usage

``` r
detect_skeleton_points(img, select.layer = NULL, skeletonize = FALSE)
```

## Arguments

- img:

  Binary skeleton image. If `skeletonize = TRUE`, a segmented
  (non-skeleton) mask can be supplied instead.

- select.layer:

  Layer index for multi-layer rasters

- skeletonize:

  Logical. If `TRUE`, `img` is treated as a segmented mask and reduced
  to a skeleton internally via
  [`skeletonize_image()`](https://jcunow.github.io/Rootopia/reference/skeletonize_image.md)
  before detecting points. Default `FALSE` (assumes `img` is already a
  skeleton).

## Value

List with:

- endpoints:

  SpatRaster marking pixels with exactly one neighbor

- branching_points:

  SpatRaster marking pixels with more than two neighbors

## Details

The computation proceeds as follows:

1\. Raster values are converted into a numeric matrix. 2. A 1-pixel
zero-padding border is added around the matrix. 3. For each pixel in the
original image: - A 3x3 window is extracted from the padded matrix - The
center pixel is excluded - Neighbor count is computed as the sum of
remaining 8 values

Classification rules: - Endpoint: pixel == 1 AND neighbor count == 1 -
Branch point: pixel == 1 AND neighbor count \> 2

Outputs are converted back into SpatRaster objects.

## Examples

``` r
# Load example binary segmentation
data(seg_Oulanka2023_Session01_T067)

# Ensure single-layer raster if needed
img <- seg_Oulanka2023_Session01_T067

# Skeletonize
skel <- skeletonize_image(img, select.layer = 2, verbose = FALSE)
skel.points <- detect_skeleton_points(skel)
```
