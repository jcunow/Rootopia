# Changelog

## RootScanR 1.1.0

- [`load_flexible_image()`](https://jcunow.github.io/RootScanR/reference/load_flexible_image.md):
  `SpatRaster` outputs with 3 or 4 layers now get
  [`terra::RGB()`](https://rspatial.github.io/terra/reference/RGB.html)
  metadata set automatically, so
  [`terra::plotRGB()`](https://rspatial.github.io/terra/reference/plotRGB.html)
  works without an extra manual call. The orphaned, unintegrated
  `cimg_to_spatrast()` patch (and its man page) have been removed since
  the fix is now applied directly.

- [`root_diameter()`](https://jcunow.github.io/RootScanR/reference/root_diameter.md):
  fixed `root_volume`, which previously summed `diameter^2 * pi` (off by
  a factor of 4, and missing the per-pixel length term entirely). Volume
  is now `sum(pi * (diameter/2)^2)` per skeleton pixel (cylinder of
  radius `diameter/2`, length 1 px), converted to `unit^3`. Also adds
  `root_surface_area` (`sum(pi * diameter)` per skeleton pixel,
  converted to `unit^2`) — the lateral surface area of the same
  cylinders.

- [`root_thickness()`](https://jcunow.github.io/RootScanR/reference/root_thickness.md)
  is now deprecated
  ([`.Deprecated()`](https://rdrr.io/r/base/Deprecated.html)); it is a
  naive area/length estimator that ignores branching and local width
  variation. Use
  [`root_diameter()`](https://jcunow.github.io/RootScanR/reference/root_diameter.md)’s
  `mean_diameter`, `median_diameter`, `root_volume`, and new
  `root_surface_area` instead.

- [`clean_image()`](https://jcunow.github.io/RootScanR/reference/clean_image.md)
  gains optional `pre_threshold`, `pre_threshold_method`, and
  `pre_threshold_window_size` arguments. When `pre_threshold` is set,
  [`image_threshold()`](https://jcunow.github.io/RootScanR/reference/image_threshold.md)
  binarizes the input *before* hole-filling and artifact removal, so
  non-binary (e.g. probability/grayscale) inputs can be cleaned
  directly.

- `zoning(mode = "rotation")`: `rotation_slices = c(i, i)` now selects a
  single circumferential slice `i` (previously required
  `rotation_slices[1] < rotation_slices[2]` and errored on equal
  values). This enables looping over individual rotation-axis slices,
  e.g. for per-slice trait extraction.

- Minirhizotron vignette: added a “Cumulative root accumulation with
  depth” example using
  [`root_accumulation()`](https://jcunow.github.io/RootScanR/reference/root_accumulation.md).

- Rotation Bias vignette: replaced the temporally-framed “Rhythmicity
  Analysis” section with a “Circumferential zoning and rhythmicity”
  section that uses `zoning(mode = "rotation")` to split the tube
  surface into circumferential slices and applies
  [`rhythmicity()`](https://jcunow.github.io/RootScanR/reference/rhythmicity.md)
  /
  [`fit_sine_curve()`](https://jcunow.github.io/RootScanR/reference/fit_sine_curve.md)
  to test for a systematic (top-down / side-to-side) pattern in root
  distribution around the tube. Notes planned future work: a combined
  rotation-shift + censor helper, and estimating the rotation center
  from the amplitude peak of this circumferential fit.

- New root branching-order pipeline:
  [`root_graph_pipeline()`](https://jcunow.github.io/RootScanR/reference/root_graph_pipeline.md),
  [`branch_order_map()`](https://jcunow.github.io/RootScanR/reference/branch_order_map.md),
  [`order_metrics()`](https://jcunow.github.io/RootScanR/reference/order_metrics.md),
  [`summarize_orders()`](https://jcunow.github.io/RootScanR/reference/summarize_orders.md),
  [`convert_root_units()`](https://jcunow.github.io/RootScanR/reference/convert_root_units.md),
  [`order_classification_map()`](https://jcunow.github.io/RootScanR/reference/order_classification_map.md),
  [`prune_terminal_segments()`](https://jcunow.github.io/RootScanR/reference/prune_terminal_segments.md),
  [`render_order_overlay()`](https://jcunow.github.io/RootScanR/reference/render_order_overlay.md),
  and
  [`plot_order_window()`](https://jcunow.github.io/RootScanR/reference/plot_order_window.md)
  classify a root skeleton into a segment graph and assign
  tip/root/branch order (main axis vs. lateral roots).

- [`root_depth_metrics()`](https://jcunow.github.io/RootScanR/reference/root_depth_metrics.md)
  gains `calc_root_order_metrics`, which runs the new branching-order
  pipeline per image and adds per-depth-bin (`mean.branch_order`,
  `max.branch_order`, `mean.root_order`, `lateral_root_fraction`) and
  per-tube (`main_root.*`, `lateral_roots.*`, `n_root_orders`) columns.

- Vignettes updated with a “Root branching order” section in the
  Minirhizotron, Flatbed, and Batch Processing tutorials.

- [`root_depth_metrics()`](https://jcunow.github.io/RootScanR/reference/root_depth_metrics.md):
  `path.skl` is now optional. If it is `NULL` or a skeleton file is
  missing for an image, the skeleton is computed internally via
  [`skeletonize_image()`](https://jcunow.github.io/RootScanR/reference/skeletonize_image.md).
  The “Minimum requirements” table and documentation in the Batch
  Processing vignette have been corrected to reflect that both
  `path.skl` and `soil_starts` are optional.

- Fixed
  [`skeletonize_image()`](https://jcunow.github.io/RootScanR/reference/skeletonize_image.md),
  which previously called a non-existent `lut_thin()` function and would
  error on every call; it now correctly uses the LUT-based Zhang-Suen
  thinning implementation
  ([`lut_thin_fast()`](https://jcunow.github.io/RootScanR/reference/lut_thin_fast.md)).

- Fixed
  [`root_diameter()`](https://jcunow.github.io/RootScanR/reference/root_diameter.md)’s
  internal skeletonization fallback (used when `skeleton.img` is not
  supplied), which called
  [`skeletonize_image()`](https://jcunow.github.io/RootScanR/reference/skeletonize_image.md)
  with a non-existent `methods = skeleton_method` argument and errored
  on every call. `skeleton_method` is now documented as currently
  unused.

- All remaining functions that consume a skeleton can now compute one
  internally via
  [`skeletonize_image()`](https://jcunow.github.io/RootScanR/reference/skeletonize_image.md)
  if none is supplied:

  - `root_length(img, skeletonize = TRUE)` and
    `detect_skeleton_points(img, skeletonize = TRUE)` treat `img` as a
    segmented mask and skeletonize it first.
  - [`root_graph_pipeline()`](https://jcunow.github.io/RootScanR/reference/root_graph_pipeline.md)
    and
    [`branch_order_map()`](https://jcunow.github.io/RootScanR/reference/branch_order_map.md)
    now accept `skel = NULL`; if omitted, the skeleton is computed from
    `mask` via
    [`skeletonize_image()`](https://jcunow.github.io/RootScanR/reference/skeletonize_image.md).

## RootScanR 1.0.0

- Initial CRAN submission.
