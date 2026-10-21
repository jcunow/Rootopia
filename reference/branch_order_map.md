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
  order = c("strahler_order", "order_from_base", "root_max_order_from_tip",
    "order_from_tip"),
  unit = "cm",
  dpi = 300,
  length_method = "kimura",
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
  [`skeletonize_image()`](https://jcunow.github.io/Rootopia/reference/skeletonize_image.md).

- mask:

  Filled (un-thinned) root mask on the same grid, used for the distance
  transform / diameters. If `NULL`, the skeleton is used and diameters
  collapse to ~1 px. Required if `skel` is `NULL`.

- order:

  Which scheme labels `$class_map`/`$summary`: `"strahler_order"`,
  `"order_from_tip"`, `"root_max_order_from_tip"` or
  `"order_from_base"`. The 1.1.0 spellings `"tip_order"`, `"root_order"`
  and `"branch_order"` are accepted and resolved to these.

- unit:

  Reporting unit: `"cm"`, `"inch"`, or `"px"`.

- dpi:

  Scan resolution (dots per inch); required for cm/inch.

- length_method:

  `"kimura"` (default) or `"polyline"`. The sqrt(2) chain code walks the
  skeleton pixel to pixel, so it traces a staircase and reports more
  length than the line it represents; Kimura's estimator corrects for
  that and is the better estimate. The two agree exactly on a perfectly
  horizontal, vertical or 45-degree segment and differ only where a
  chain mixes orthogonal and diagonal steps. See
  [`convert_root_units`](https://jcunow.github.io/Rootopia/reference/convert_root_units.md)
  for the definition. Use `"polyline"` to reproduce output from Rootopia
  1.1.0 and earlier.

- template:

  `SpatRaster` defining the `$class_map` grid; defaults to `skel` when
  it is a raster.

- overlay_png:

  Optional path; writes the order-colored validation PNG.

- return_map:

  Build `$class_map`.

- ...:

  Passed to
  [`root_graph_pipeline`](https://jcunow.github.io/Rootopia/reference/root_graph_pipeline.md)
  (e.g. `dt_backend`, `crossing_straight`, `prune_iter`, and
  `diam_weight` — the diameter-vs-angle weight for the continuation
  rule, \>= 0).

## Value

An object of class `"branchOrderMap"`: a list with `$edges`, `$summary`,
`$class_map` (`SpatRaster`; chosen order per pixel, `NA` off-root), and
`$order`, `$unit`, `$dpi`, `$length_method`.

## Details

All four order schemes are always computed and stored on `$edges`;
`order` only selects which labels `$class_map` and `$summary`. Two label
a single graph segment and two label a whole continuous root, which is
what the column names say: a per-root column holds the same value on
every segment of that root, so averaging it over segments weights long
roots more.

- `order_from_tip` (per segment):

  Leaf-peeling count. Every terminal segment is order 1; peeling
  terminals away round by round, a segment's order is
  `1 + max(child orders)`. It rises toward the interior, so the distal
  end of even a thick root is 1. The count runs from the tips inward
  (centripetal), the direction root ecology's morphometric
  classification uses, but it is not strict Strahler ordering: it
  increments at *every* junction, so an axis carrying ten single
  laterals reaches 11.

- `strahler_order` (per segment):

  The same peeling with the strict rule: a segment's order rises only
  where two of its children share the highest order, otherwise it
  inherits that order unchanged. A lone thin lateral joining a thick
  axis leaves the axis where it was.

- `root_max_order_from_tip` (per root):

  No count of its own. Segments are grouped into continuous roots
  (continuation rule below) and every segment of a root takes the
  *largest* `order_from_tip` found along it – which, since that column
  never falls as you move basally, is its value at the root's basal end.
  A thick main axis therefore keeps its high value out to its tip
  instead of dropping to 1. Sensitive to how deep the deepest subtree
  runs.

- `order_from_base` (per root):

  Generation counting, the direction the developmental classification
  uses. The thickest root (length-weighted mean diameter) in each
  connected component is order 1; counting branching hops outward, its
  laterals are 2, theirs 3, and so on. Independent of subtree depth, so
  one heavily branched root does not inflate the rest. Order 1 is chosen
  by thickness, not by knowing which root grew first.

Up to Rootopia 1.1.0 these columns were called `tip_order`, `root_order`
and `branch_order`. Those names are still returned, as copies of the
columns above, and are still accepted by `order`, `order_col` and
`value`.

**Continuation rule** (groups segments into roots, for the two per-root
columns): at each junction the two arms forming the same root are chosen
by `straightness + diam_weight * diameter_similarity`. `diam_weight = 0`
uses angle only; larger values let thickness dominate. The odd arm out
is a lateral. The per-segment columns do not use this rule.

## See also

[`root_graph_pipeline`](https://jcunow.github.io/Rootopia/reference/root_graph_pipeline.md),
[`order_classification_map`](https://jcunow.github.io/Rootopia/reference/order_classification_map.md),
[`summarize_orders`](https://jcunow.github.io/Rootopia/reference/summarize_orders.md),
[`convert_root_units`](https://jcunow.github.io/Rootopia/reference/convert_root_units.md),
[`root_transport`](https://jcunow.github.io/Rootopia/reference/root_transport.md)
and
[`network_metrics`](https://jcunow.github.io/Rootopia/reference/network_metrics.md),
which both take the `$edges` table returned here.

## Examples

``` r
if (FALSE) { # \dontrun{
skel <- skeletonize_image(mask)
res  <- branch_order_map(skel, mask, order = "order_from_base",
                         unit = "cm", dpi = 300, dt_backend = "imager")
res$summary
terra::plot(res$class_map, maxcell = Inf)
} # }
```
