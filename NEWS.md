# Rootopia 1.1.0

Beta release, available from the GitHub repository only. Not submitted to CRAN.

* Requires R >= 4.1 (the package uses the native `|>` pipe) and `terra` >= 1.8.
* `root_depth_metrics()` / `batch_root_traits()` run a full depth-profile
  analysis over a directory of images; flatbed scans and flat rhizotron windows
  are the default geometry, cylindrical minirhizotron tubes are opt-in via the
  tube parameters.
* `stitch_root_scans()` combines overlapping scans of one tube into a mosaic.
* `branch_order_map()` and `root_graph_pipeline()` assign tip, branch and root
  order from a skeleton.
