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
  crossing_straight = -0.5,
  color_by = c("branch_order", "root_order", "tip_order"),
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

- crossing_straight:

  Straightness threshold for crossing resolution.

- color_by:

  Which order column the overlay colours by.

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
