# RootScanR 1.1.0

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
