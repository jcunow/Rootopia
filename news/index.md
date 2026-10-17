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
  sqrt(2) chain code, walks the skeleton pixel to pixel and so traces a
  staircase rather than the line it represents. Measured against digital
  straight lines of known length it runs about 5.7% long on average, up
  to 8.2% near 22.5 degrees, and zero only for exactly horizontal or
  vertical segments; Kimura’s estimator reduces that to about 1.8%. This
  changes reported root length, and everything derived from it, by a few
  percent. Pass `length_method = "polyline"` to reproduce earlier
  output.

- `image_threshold(method = "adaptive")` gains `adaptive_k`, the
  multiple of the local mean a pixel must reach, defaulting to 1.15. The
  adaptive rule compares a pixel against the mean of its own
  neighborhood, and a window holding only background averages to the
  background level, so a multiplier at or below 1 marks half or more of
  a root-free area as root. Measured on synthetic scans with a 5%
  background coefficient of variation, a multiplier of 1.15 leaves 0.1%
  of a root-free patch marked where 1.05 leaves 16% and 0.4 leaves all
  of it. Across root/background contrasts of 1.3 to 2.2 and background
  CVs of 2% to 10% the best value fell between 1.05 and 1.30. Values at
  or below 1 warn, and `threshold` – a fraction of the global maximum,
  used by the global and deblur modes – warns rather than being silently
  ignored when passed alongside `method = "adaptive"`.

- [`image_threshold()`](https://jcunow.github.io/Rootopia/reference/image_threshold.md)
  raises the `window_size` default from 15 to 51. The local mean
  includes the pixel under test together with its root, so a root
  filling much of its own window raises the level it is judged against
  and is rejected. On synthetic scans with roots up to 9 px wide, the
  best F1 against a known mask was 0.54 at a window of 15, 0.90 at 51
  and 0.96 at 101. A window wider than the image is clamped to the
  largest odd window that fits, with a warning, instead of reaching
  [`terra::focal()`](https://rspatial.github.io/terra/reference/focal.html)
  and failing there.

- [`?image_threshold`](https://jcunow.github.io/Rootopia/reference/image_threshold.md)
  states what neither method can do. Both mark pixels brighter than a
  reference, so where root and soil are equally bright the useful range
  for `adaptive_k` is empty and no setting works. Measured against its
  own segmentation, the bundled minirhizotron scan has roots averaging
  1.07 times the soil brightness while the soil’s coefficient of
  variation runs from 0.23 to 0.93 across the frame; adaptive
  thresholding reaches an F1 of about 0.1 there. The adaptive example
  runs on a small synthetic gradient instead, where both conditions hold
  and the three roots are recovered exactly.

- `image_threshold(method = "adaptive")` computes the local mean as two
  1-D
  [`terra::focal()`](https://rspatial.github.io/terra/reference/focal.html)
  passes rather than one square kernel. A box mean is separable, so the
  work per pixel is proportional to the window side instead of its
  square: on a 1000 x 1000 raster the local mean takes 0.73 s at a
  window of 101 where a square kernel takes 29 s. Results are identical
  to machine precision when no cell is `NA`.

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
