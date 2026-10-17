# Changelog

## Rootopia (development version)

- New
  [`network_metrics()`](https://jcunow.github.io/Rootopia/reference/network_metrics.md).
  Summarizes any segment network given as an edge table: nodes, edges,
  components, tips, branch points, and the alpha, beta and gamma
  connectivity coefficients. The `topology` argument decides whether
  loop metrics are reported as traits (`"reticulate"`, e.g. mycelium or
  leaf venation) or under `$qc` (`"tree"`, where they are identities,
  and `"projected_tree"`, where a loop is two roots crossing in
  projection and the count becomes a segmentation diagnostic).

- New
  [`root_transport()`](https://jcunow.github.io/Rootopia/reference/root_transport.md).
  For every segment in an edge table it reports the path length and
  Hagen-Poiseuille path resistance from the base of its root system, and
  the number of tips, the root length and the lateral surface area that
  lie beyond it. Each component is rooted at its thickest free end, the
  same convention the Strahler ordering uses. Segments on a cycle have
  no subtree and get `NA` with a warning.

- [`render_order_overlay()`](https://jcunow.github.io/Rootopia/reference/render_order_overlay.md)
  paints continuous per-segment values, not only branch orders. It picks
  the legend from the data – one labelled swatch per order for small
  integers, a viridis ramp with a scale bar otherwise – and
  `continuous`, `value_label` and `n_bins` control it. The first
  argument is renamed from `tip_order` to `values`; the old name still
  works and warns.

- [`convert_root_units()`](https://jcunow.github.io/Rootopia/reference/convert_root_units.md)
  and
  [`branch_order_map()`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md)
  default to `length_method = "kimura"`. The previous default, the
  sqrt(2) chain code, walks the skeleton pixel to pixel and so measures
  the staircase rather than the line it represents, which reads long.
  Kimura’s estimator corrects for that and is the better estimate of
  root length. The two agree exactly on a perfectly horizontal, vertical
  or 45-degree segment and differ only where a chain mixes orthogonal
  and diagonal steps. This changes reported root length, and everything
  derived from it, by a few percent. Pass `length_method = "polyline"`
  to reproduce earlier output or to compare against software that
  reports the chain code.

- `image_threshold(method = "adaptive")` gains `adaptive_k`, the
  multiple of the local mean a pixel must reach, defaulting to 1.15. The
  adaptive rule compares a pixel against the mean of its own
  neighborhood, and a window holding only background averages to the
  background level, so about half the background sits above it: a
  multiplier at or below 1 marks half or more of a root-free area as
  root, which the previous default of 0.4 did. Values at or below 1 now
  warn, and `threshold` – a fraction of the global maximum, used by the
  global and deblur modes – warns rather than being silently ignored
  when passed alongside `method = "adaptive"`.

- [`image_threshold()`](https://jcunow.github.io/Rootopia/reference/image_threshold.md)
  raises the `window_size` default from 15 to 51. The local mean
  includes the pixel under test together with its root, so a root
  filling much of its own window raises the level it is judged against
  and is rejected. The window has to be large compared with the widest
  root; the right value scales with scan resolution. A window wider than
  the image is clamped to the largest odd window that fits, with a
  warning, instead of reaching
  [`terra::focal()`](https://rspatial.github.io/terra/reference/focal.html)
  and failing there.

- `image_threshold(method = "adaptive")` computes the local mean as two
  1-D
  [`terra::focal()`](https://rspatial.github.io/terra/reference/focal.html)
  passes rather than one square kernel. A box mean is separable, so the
  work per pixel grows with the window side instead of its square, which
  makes a usefully large window affordable. Results match the square
  kernel to floating-point precision when no cell is `NA`.

- [`?image_threshold`](https://jcunow.github.io/Rootopia/reference/image_threshold.md)
  states what neither method can do. Both mark pixels brighter than a
  reference, so where root and soil are about equally bright no setting
  separates them: the margin the background noise demands exceeds what
  the contrast allows. The bundled minirhizotron scan is that case,
  which is why Rootopia expects a segmentation from a trained model
  rather than offering to produce one. The adaptive example runs on a
  small synthetic gradient, where both conditions hold.

- [`clean_image()`](https://jcunow.github.io/Rootopia/reference/clean_image.md)
  renames `max_hole_size` to `fill_holes_upto` and `max_artifact_size`
  to `drop_objects_upto`. The old names said how big a thing was but not
  what happened to it; the new ones state that anything larger is left
  alone. The old names still work and warn.

- [`clean_image()`](https://jcunow.github.io/Rootopia/reference/clean_image.md)
  gains `size_unit` and `dpi`. With `size_unit = "mm2"` the thresholds
  are physical areas rather than pixel counts, so the same call cleans
  the same physical size whether a scan is 300 or 1200 dpi.

- [`root_depth_metrics()`](https://jcunow.github.io/Rootopia/reference/root_depth_metrics.md)
  /
  [`batch_root_traits()`](https://jcunow.github.io/Rootopia/reference/root_depth_metrics.md)
  rename `clean_max_hole_size` to `clean_fill_holes_upto` and
  `clean_max_artifact_size` to `clean_drop_objects_upto`, and gain
  `clean_size_unit`. The old names still work and warn.

- [`root_depth_metrics()`](https://jcunow.github.io/Rootopia/reference/root_depth_metrics.md)
  validates the cleaning thresholds before the image loop. A non-numeric
  or `NA` threshold now stops with a message naming the argument, where
  it previously failed inside the per-image `tryCatch` and left the
  image silently uncleaned. `NULL` means no limit and switches cleaning
  on, matching
  [`clean_image()`](https://jcunow.github.io/Rootopia/reference/clean_image.md).

## Rootopia 1.1.0

Beta release, available from the GitHub repository only. Not submitted
to CRAN.

- Requires R \>= 4.1 (the package uses the native `|>` pipe) and `terra`
  \>= 1.8.
