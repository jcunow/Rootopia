# Branch-order classification of a root skeleton

Main entry point. Converts a binary root skeleton into a per-segment
graph, assigns each segment a branching order, and returns the
per-segment table, a per-order summary, and a classification raster
aligned to the input. Lengths and diameters are reported in real units.

## Usage

``` r
branch_order_map(
  skel = NULL,
  mask = NULL,
  order = c("branch_order", "root_order", "tip_order"),
  unit = "cm",
  dpi = 300,
  length_method = "polyline",
  template = NULL,
  overlay_png = NULL,
  return_map = TRUE,
  ...
)
```

## Arguments

- skel:

  Binary skeleton: single-layer `SpatRaster` (preferred) or 0/1 matrix.
  Pass the raster directly; do not pre-convert with
  [`as.matrix()`](https://rspatial.github.io/terra/reference/coerce.html).
  If `NULL`, it is computed from `mask` via
  [`skeletonize_image()`](https://jcunow.github.io/RootScanR/reference/skeletonize_image.md).

- mask:

  Filled (un-thinned) root mask on the same grid, used for the distance
  transform / diameters. If `NULL`, the skeleton is used and diameters
  collapse to ~1 px. Required if `skel` is `NULL`.

- order:

  Which scheme labels `$class_map`/`$summary`: `"branch_order"`,
  `"root_order"`, or `"tip_order"`.

- unit:

  Reporting unit: `"cm"`, `"inch"`, or `"px"`.

- dpi:

  Scan resolution (dots per inch); required for cm/inch.

- length_method:

  `"polyline"` (sqrt(2) chain code, follows curves) or `"kimura"`
  (per-segment Kimura correction, better for straight segments).

- template:

  `SpatRaster` defining the `$class_map` grid; defaults to `skel` when
  it is a raster.

- overlay_png:

  Optional path; writes the order-coloured validation PNG.

- return_map:

  Build `$class_map`.

- ...:

  Passed to
  [`root_graph_pipeline`](https://jcunow.github.io/RootScanR/reference/root_graph_pipeline.md)
  (e.g. `dt_backend`, `crossing_straight`, `prune_iter`, and
  `diam_weight` — the diameter-vs-angle weight for the continuation
  rule, \>= 0).

## Value

An object of class `"branchOrderMap"`: a list with `$edges`, `$summary`,
`$class_map` (`SpatRaster`; chosen order per pixel, `NA` off-root), and
`$order`, `$unit`, `$dpi`, `$length_method`.

## Details

All three order schemes are always computed and stored on `$edges`;
`order` only selects which labels `$class_map` and `$summary`.

- `tip_order` (per segment):

  Topological leaf-peeling (Strahler-like). Every terminal segment is
  order 1; peeling terminals away round by round, a segment's order is
  `1 + max(child orders)`. Order rises toward the interior, so the
  distal end of even a thick root is 1.

- `root_order` (per root):

  Segments are grouped into continuous roots (continuation rule below);
  each root takes the *maximum* `tip_order` along it, so a thick main
  axis keeps its high order out to its tip. Sensitive to how deep the
  deepest subtree runs.

- `branch_order` (per root):

  Centrifugal generation. The thickest root (length-weighted mean
  diameter) in each connected component is order 1; counting branching
  hops outward, its laterals are 2, theirs 3, and so on. Independent of
  subtree depth, so one heavily branched root does not inflate the rest.

**Continuation rule** (groups segments into roots for `root_order` and
`branch_order`): at each junction the two arms forming the same root are
chosen by `straightness + diam_weight * diameter_similarity`.
`diam_weight = 0` uses angle only; larger values let thickness dominate.
The odd arm out is a lateral. `tip_order` does not use this rule.

## See also

[`root_graph_pipeline`](https://jcunow.github.io/RootScanR/reference/root_graph_pipeline.md),
[`order_classification_map`](https://jcunow.github.io/RootScanR/reference/order_classification_map.md),
[`summarize_orders`](https://jcunow.github.io/RootScanR/reference/summarize_orders.md),
[`convert_root_units`](https://jcunow.github.io/RootScanR/reference/convert_root_units.md)

## Examples

``` r
if (FALSE) { # \dontrun{
skel <- skeletonize_image(mask)
res  <- branch_order_map(skel, mask, order = "branch_order",
                         unit = "cm", dpi = 300, dt_backend = "imager")
res$summary
terra::plot(res$class_map)
} # }
```
