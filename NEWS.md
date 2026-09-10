# Rootopia 1.1.0

* Beta Version, Github repo only.

## Bug fixes

* `root_depth_metrics()` joined directory paths with `paste0()`, so a path
  without a trailing separator produced `<dir><file>` and every image failed to
  load. The wrapper's fault tolerance reported this as "skipped" and returned an
  empty result. Directories without a trailing slash now work.

* `root_depth_metrics()`: the distribution-indices and advanced-metrics blocks
  referenced columns as `dplyr::.data$col`, which raises an error rather than
  resolving. Both groups are on by default, so every run silently returned `NA`
  for `mrd`, `total.length.density`, `rootlength.fraction`, `ent_per_rootpx`,
  `patch_density_norm` and `mean.var.diameter`. They now compute.

* `order_metrics()` divided by zero when computing `branching_frequency` for a
  group of zero total length; it now returns `NA`, matching the guard already
  applied to `mean_diameter`.

## Behaviour changes

* `prune_terminal_segments()` (and `prune_iter` in `root_graph_pipeline()` /
  `branch_order_map()`) previously thresholded the raw pixel chain, while the
  edge table reported that length plus the rim-to-centroid junction stub
  (roughly 1 px per junction end). A segment could therefore be pruned at a
  `min_length` it visibly cleared in the returned table. Both now use the same
  definition, so slightly fewer segments are pruned at a given threshold.
  Pruning is off by default, so only callers who set `prune_iter > 0` are
  affected.

## Performance

* The branching pipeline is roughly twice as fast on large root systems
  (10.3 s to 4.8 s on a 4700-segment skeleton), with an identical edge table.
  `build_edge_table()` no longer builds one `data.frame` per segment, and the
  node-incidence table is built once by `split()` instead of by `rbind()` in
  three separate functions.

## Documentation and API

* `list_tubes()` and `list_scan_files()` are now exported. The stitching
  vignette instructed readers to call them, and `R/stitch_scans.R` lists both as
  part of its public API, but they were tagged `@keywords internal`.

* `tip_order` is documented as leaf-peeling (centrifugal) depth rather than
  "Strahler-like"; it increments at every junction, which strict Strahler
  ordering does not.

* `order_classification_map()` documents that junction contraction leaves one
  unpainted pixel per branch point, so pixel counts taken off the map run below
  the skeleton. Lengths should be read from `et$length`.

* The `select_layer` default is documented as `NULL` in `root_diameter()` and
  `deep_drive()`, matching the signatures (both said `2`).
