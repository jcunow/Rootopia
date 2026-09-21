# Changelog

## Rootopia 1.1.0

Beta release, available from the GitHub repository only. Not submitted
to CRAN.

- Requires R \>= 4.1 (the package uses the native `|>` pipe) and `terra`
  \>= 1.8.
- [`root_depth_metrics()`](https://jcunow.github.io/Rootopia/reference/root_depth_metrics.md)
  /
  [`batch_root_traits()`](https://jcunow.github.io/Rootopia/reference/root_depth_metrics.md)
  run a full depth-profile analysis over a directory of images; flatbed
  scans and flat rhizotron windows are the default geometry, cylindrical
  minirhizotron tubes are opt-in via the tube parameters.
- [`stitch_root_scans()`](https://jcunow.github.io/Rootopia/reference/stitch_root_scans.md)
  combines overlapping scans of one tube into a mosaic.
- [`branch_order_map()`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md)
  and
  [`root_graph_pipeline()`](https://jcunow.github.io/Rootopia/reference/root_graph_pipeline.md)
  assign tip, branch and root order from a skeleton.
