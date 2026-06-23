# Prune short or thin spurs from a skeleton or segmentation image

Standalone clean-up step that removes short/thin terminal branches
("spurs") that survive segmentation. Spurs are defined on the skeleton
graph (that is where connectivity lives), but the result can be returned
either as the cleaned skeleton (`output = "skeleton"`) or propagated
back into segmentation space (`output = "mask"`).

## Usage

``` r
prune_skeleton(
  skel,
  mask = NULL,
  min_length = 0,
  min_diameter = 0,
  iter = 1L,
  output = c("skeleton", "mask"),
  verbose = FALSE
)
```

## Arguments

- skel:

  Binary skeleton: single-layer `SpatRaster` or 0/1 matrix.

- mask:

  Filled root mask on the same grid, for the distance transform and
  (when `output = "mask"`) reconstruction. Optional but recommended.

- min_length:

  Minimum terminal-segment length (px) to keep.

- min_diameter:

  Minimum terminal-segment diameter (px) to keep.

- iter:

  Number of pruning passes.

- output:

  `"skeleton"` (default) returns the pruned skeleton; `"mask"` returns
  the pruned binary segmentation.

- verbose:

  Print progress.

## Value

A cleaned image matching the class of `skel` (`SpatRaster` in,
`SpatRaster` out; matrix in, matrix out).

## Details

The skeleton is traced into segments and
[`prune_terminal_segments`](https://jcunow.github.io/Rootopia/reference/prune_terminal_segments.md)
removes terminal segments below `min_length` / `min_diameter`, iterated
`iter` times. For `output = "mask"`, each surviving skeleton pixel is
regrown by its local distance-transform radius (the medial axis
reconstruction) and intersected with the original `mask`, so the spur's
body is removed from the segmentation while the real roots are restored
to their original thickness. A `mask` is therefore strongly recommended
(and required for trustworthy `output = "mask"` and the `min_diameter`
test); without it diameters collapse to ~1 px.

## See also

[`prune_terminal_segments`](https://jcunow.github.io/Rootopia/reference/prune_terminal_segments.md),
[`clean_image`](https://jcunow.github.io/Rootopia/reference/clean_image.md)

## Examples

``` r
if (FALSE) { # \dontrun{
skl <- skeletonize_image(mask)
# remove spurs shorter than 15 px, return a cleaned skeleton
skl2 <- prune_skeleton(skl, mask, min_length = 15, iter = 2)
# ... or clean the segmentation itself
seg2 <- prune_skeleton(skl, mask, min_length = 15, output = "mask")
} # }
```
