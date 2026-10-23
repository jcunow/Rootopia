# Detect endpoints and branching points in a skeleton image

Computes local connectivity of each foreground pixel using an
8-neighborhood (Moore neighborhood).

## Usage

``` r
detect_skeleton_points(img, select_layer = NULL, skeletonize = FALSE)
```

## Arguments

- img:

  Binary skeleton image. If `skeletonize = TRUE`, a segmented
  (non-skeleton) mask can be supplied instead.

- select_layer:

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

Both are rebuilt without the input's extent, so carry it across with
[`terra::ext()`](https://rspatial.github.io/terra/reference/ext.html)
before overlaying them on the skeleton.

## Details

The computation proceeds as follows:

1.  Raster values are converted into a numeric matrix.

2.  A 1-pixel zero-padding border is added around the matrix.

3.  For each pixel in the original image:

    - A 3x3 window is extracted from the padded matrix

    - The center pixel is excluded

    - Neighbor count is computed as the sum of remaining 8 values

Classification rules:

- Endpoint: pixel == 1 AND neighbor count == 1

- Branch point: pixel == 1 AND neighbor count \> 2

Outputs are converted back into SpatRaster objects.

## What this counts, and what it over-counts

The test is applied to each pixel on its own, which makes the two
rasters the right input for an overlay – they say *where* the points are
– and the wrong input for a count. Branch points are over-counted for
two reasons, and tips for a third:

- **A junction is several pixels wide after thinning.** Each of those
  pixels can have more than two neighbors, so one real fork marks
  several branch-point pixels.

- **Crossings are not resolved.** Two roots that cross in projection
  touch in the image, and the skeleton joins them into a four-armed
  node. Nothing branches there, but the node has more than two neighbors
  and is counted as a branch point.

- **Thinning spurs end in a tip.** A short terminal stub left by
  thinning has a degree-1 end, so it is counted as a root tip unless it
  is pruned first with
  [`prune_skeleton`](https://jcunow.github.io/Rootopia/reference/prune_skeleton.md).

For counts, use the segment graph instead:
[`branch_order_map()`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md)
merges adjacent junction pixels into one node, separates four-armed
crossings where both candidate through-roots are straight enough, and
optionally prunes weak terminal segments. Its `$edges` table carries
`n_tips` and `n_branch_points` per segment, summed per order class in
`$summary`. The two will not agree, and the graph's numbers are the
defensible ones.

## See also

[`branch_order_map`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md)
for counts from the segment graph,
[`prune_skeleton`](https://jcunow.github.io/Rootopia/reference/prune_skeleton.md)
for removing thinning spurs first.

## Examples

``` r
# Load example binary segmentation
data(seg_Oulanka2023_Session01_T067)

# Ensure single-layer raster if needed
img <- seg_Oulanka2023_Session01_T067

# Skeletonize
skel <- skeletonize_image(img, select_layer = 2, verbose = FALSE)
skel.points <- detect_skeleton_points(skel)
```
