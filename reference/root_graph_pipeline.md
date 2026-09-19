# Skeleton-to-graph root ordering pipeline (pixel units)

Core engine: crops to the foreground, computes the distance transform,
traces segments with junction contraction, resolves crossings,
optionally prunes weak tips, and assigns all three order schemes.
Returns a per-segment edge table in *pixels*; use
[`branch_order_map`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md)
for a unit-aware wrapper.

## Usage

``` r
root_graph_pipeline(
  skel = NULL,
  mask = NULL,
  verbose = TRUE,
  dt_backend = "auto",
  crop = TRUE,
  overlay_png = NULL,
  max_side = 2000,
  keep_segments = FALSE,
  resolve_overlaps = TRUE,
  splice_passthrough = TRUE,
  crossing_straight = -0.5,
  crossing_diam_ratio = 0,
  color_by = c("strahler_order", "branch_order", "root_order", "tip_order"),
  diam_weight = 0.5,
  prune_min_length = 0,
  prune_min_diameter = 0,
  prune_iter = 0L
)
```

## Arguments

- skel:

  Binary skeleton: single-layer `SpatRaster` (preferred) or 0/1 matrix.
  Do not pre-convert a raster with
  [`as.matrix()`](https://rspatial.github.io/terra/reference/coerce.html).
  If `NULL`, it is computed from `mask` via
  [`skeletonize_image()`](https://jcunow.github.io/Rootopia/reference/skeletonize_image.md).

- mask:

  Filled root mask on the same grid for the distance transform. Required
  if `skel` is `NULL`.

- verbose:

  Print progress.

- dt_backend:

  Distance-transform backend: `"auto"`, `"imager"`, or `"baseR"`.

- crop:

  Crop to the foreground bounding box (with a 1-px pad) first.

- overlay_png:

  Optional path for the validation PNG.

- max_side:

  Long-side cap (px) for the overlay.

- keep_segments:

  Attach the traced segments as `attr(., "segments")` (needed for
  re-plotting and classification maps).

- resolve_overlaps:

  Resolve degree-4 crossings by continuity.

- splice_passthrough:

  Dissolve contracted junctions that carry only two segment ends (bends
  and thinning artefacts, not branch points) by splicing the two arms
  into one segment. Leave `TRUE` unless you specifically want one edge
  per traced skeleton chain.

- crossing_straight:

  Straightness threshold for crossing resolution.

- crossing_diam_ratio:

  Opt-in thickness test at a degree-4 node, off by default (`0` =
  resolve on geometry alone, the historical behaviour). When set (0.5 is
  a reasonable value), a node is treated as a crossing only if the two
  candidate through-roots are within this thickness ratio of each other;
  otherwise it is read as a bilateral branch and left intact. This is a
  genuine trade-off, not a free win: a thick axis crossed by a thin root
  and a thick axis with two thin laterals are the same shape in outline,
  so raising the threshold fixes the second case and breaks the first.
  Set it only if your material has more bilateral branching than
  fine-over-coarse crossing.

- color_by:

  Which order column the overlay colors by.

- diam_weight:

  Diameter-vs-angle weight for the root-continuation choice.

- prune_min_length, prune_min_diameter, prune_iter:

  Optional terminal-segment pruning (off when `prune_iter = 0`); see
  [`prune_terminal_segments`](https://jcunow.github.io/Rootopia/reference/prune_terminal_segments.md).

## Value

A per-segment edge table (data.frame) with order and diameter columns in
pixels, carrying `attr`s `crop_offset`, `dims`, and (if requested)
`segments`.

## Details

Three order columns are produced (see
[`branch_order_map`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md)
for the full rules): `tip_order` (per-segment leaf-peeling),
`root_order` (per-root max tip_order), `branch_order` (per-root
centrifugal generation from the thickest root). `color_by` selects which
one the overlay uses.

## See also

[`branch_order_map`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md)
