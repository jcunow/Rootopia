# RootScanR 1.1.0

* `load_flexible_image()`: `SpatRaster` outputs with 3 or 4 layers now get
  `terra::RGB()` metadata set automatically, so `terra::plotRGB()` works
  without an extra manual call. The orphaned, unintegrated
  `cimg_to_spatrast()` patch (and its man page) have been removed since the
  fix is now applied directly.
* `root_diameter()`: fixed `root_volume`, which previously summed
  `diameter^2 * pi` (off by a factor of 4, and missing the per-pixel length
  term entirely). Volume is now `sum(pi * (diameter/2)^2)` per skeleton pixel
  (cylinder of radius `diameter/2`, length 1 px), converted to `unit^3`. Also
  adds `root_surface_area` (`sum(pi * diameter)` per skeleton pixel, converted
  to `unit^2`) — the lateral surface area of the same cylinders.
* `root_thickness()` is now deprecated (`.Deprecated()`); it is a naive
  area/length estimator that ignores branching and local width variation. Use
  `root_diameter()`'s `mean_diameter`, `median_diameter`, `root_volume`, and
  new `root_surface_area` instead.
* `clean_image()` gains optional `pre_threshold`, `pre_threshold_method`, and
  `pre_threshold_window_size` arguments. When `pre_threshold` is set,
  `image_threshold()` binarizes the input *before* hole-filling and artifact
  removal, so non-binary (e.g. probability/grayscale) inputs can be cleaned
  directly.
* `zoning(mode = "rotation")`: `rotation_slices = c(i, i)` now selects a
  single circumferential slice `i` (previously required
  `rotation_slices[1] < rotation_slices[2]` and errored on equal values).
  This enables looping over individual rotation-axis slices, e.g. for
  per-slice trait extraction.
* Minirhizotron vignette: added a "Cumulative root accumulation with depth"
  example using `root_accumulation()`.
* Rotation Bias vignette: replaced the temporally-framed "Rhythmicity
  Analysis" section with a "Circumferential zoning and rhythmicity" section
  that uses `zoning(mode = "rotation")` to split the tube surface into
  circumferential slices and applies `rhythmicity()` / `fit_sine_curve()` to
  test for a systematic (top-down / side-to-side) pattern in root
  distribution around the tube. Notes planned future work: a combined
  rotation-shift + censor helper, and estimating the rotation center from the
  amplitude peak of this circumferential fit.

* New root branching-order pipeline: `root_graph_pipeline()`,
  `branch_order_map()`, `order_metrics()`, `summarize_orders()`,
  `convert_root_units()`, `order_classification_map()`,
  `prune_terminal_segments()`, `render_order_overlay()`, and
  `plot_order_window()` classify a root skeleton into a segment graph and
  assign tip/root/branch order (main axis vs. lateral roots).
* `root_depth_metrics()` gains `calc_root_order_metrics`, which runs the new
  branching-order pipeline per image and adds per-depth-bin
  (`mean.branch_order`, `max.branch_order`, `mean.root_order`,
  `lateral_root_fraction`) and per-tube (`main_root.*`, `lateral_roots.*`,
  `n_root_orders`) columns.
* Vignettes updated with a "Root branching order" section in the
  Minirhizotron, Flatbed, and Batch Processing tutorials.
* `root_depth_metrics()`: `path.skl` is now optional. If it is `NULL` or a
  skeleton file is missing for an image, the skeleton is computed internally
  via `skeletonize_image()`. The "Minimum requirements" table and
  documentation in the Batch Processing vignette have been corrected to
  reflect that both `path.skl` and `soil_starts` are optional.
* Fixed `skeletonize_image()`, which previously called a non-existent
  `lut_thin()` function and would error on every call; it now correctly uses
  the LUT-based Zhang-Suen thinning implementation (`lut_thin_fast()`).
* Fixed `root_diameter()`'s internal skeletonization fallback (used when
  `skeleton.img` is not supplied), which called `skeletonize_image()` with a
  non-existent `methods = skeleton_method` argument and errored on every
  call. `skeleton_method` is now documented as currently unused.
* All remaining functions that consume a skeleton can now compute one
  internally via `skeletonize_image()` if none is supplied:
  - `root_length(img, skeletonize = TRUE)` and
    `detect_skeleton_points(img, skeletonize = TRUE)` treat `img` as a
    segmented mask and skeletonize it first.
  - `root_graph_pipeline()` and `branch_order_map()` now accept
    `skel = NULL`; if omitted, the skeleton is computed from `mask` via
    `skeletonize_image()`.

# RootScanR 1.0.0

* Initial CRAN submission.
