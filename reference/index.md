# Package index

## Batch processing

High-level wrapper that runs a full depth-profile analysis over a
directory of images.

- [`root_depth_metrics()`](https://jcunow.github.io/Rootopia/reference/root_depth_metrics.md)
  : Compute root traits over a depth profile from segmented
  (mini)rhizotron images

## Stitching scans from the same Tube together

High-level wrapper that enables to combine multiple depths from a tube
into a single large mosaic iamge

- [`stitch_root_scans()`](https://jcunow.github.io/Rootopia/reference/stitch_root_scans.md)
  : Batch-stitch grouped scan sequences (tubes) into mosaics

## Image loading and preprocessing

Functions for loading, cleaning, thresholding, and transforming root
images.

- [`load_flexible_image()`](https://jcunow.github.io/Rootopia/reference/load_flexible_image.md)
  : Load an image flexibly from file or convert from memory
- [`clean_image()`](https://jcunow.github.io/Rootopia/reference/clean_image.md)
  : Clean a binary root image
- [`report_image_components()`](https://jcunow.github.io/Rootopia/reference/report_image_components.md)
  : Summarise the component-size distribution of a binary image
- [`image_threshold()`](https://jcunow.github.io/Rootopia/reference/image_threshold.md)
  : Threshold or deblur an image to binarize features
- [`rgb2gray()`](https://jcunow.github.io/Rootopia/reference/rgb2gray.md)
  : Convert RGB image to grayscale with optimized memory management and
  parallel processing
- [`prune_skeleton()`](https://jcunow.github.io/Rootopia/reference/prune_skeleton.md)
  : Prune short or thin spurs from a skeleton or segmentation image
- [`rotation_censor()`](https://jcunow.github.io/Rootopia/reference/rotation_censor.md)
  : Censor image edges based on rotation
- [`skeletonize_image()`](https://jcunow.github.io/Rootopia/reference/skeletonize_image.md)
  : Skeletonize binary image

## Root traits

Per-image trait extraction: length, diameter, architecture, and spatial
metrics.

- [`root_length()`](https://jcunow.github.io/Rootopia/reference/root_length.md)
  : Root length estimation from skeleton images
- [`root_diameter()`](https://jcunow.github.io/Rootopia/reference/root_diameter.md)
  : Estimate Root Diameters
- [`detect_skeleton_points()`](https://jcunow.github.io/Rootopia/reference/detect_skeleton_points.md)
  : Detect endpoints and branching points in a skeleton image
- [`count_pixels()`](https://jcunow.github.io/Rootopia/reference/count_pixels.md)
  : Count all pixels in a segmented image
- [`deep_drive()`](https://jcunow.github.io/Rootopia/reference/deep_drive.md)
  : Assess Root Growth Direction Relative to Depth Gradient
- [`root_scape_metrics()`](https://jcunow.github.io/Rootopia/reference/root_scape_metrics.md)
  : RootScapeMetric relies on Landscapemetrics to extract 'Root Scape'
  Features akin to landscape analysis.

## Root branching order

Convert a root skeleton into a segment graph, assign tip/root/branch
order, and summarise root architecture by order class.

- [`branch_order_map()`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md)
  : Branch-order classification of a root skeleton
- [`root_graph_pipeline()`](https://jcunow.github.io/Rootopia/reference/root_graph_pipeline.md)
  : Skeleton-to-graph root ordering pipeline (pixel units)
- [`order_metrics()`](https://jcunow.github.io/Rootopia/reference/order_metrics.md)
  : Aggregate root architecture by order, or split focal-vs-rest
- [`summarize_orders()`](https://jcunow.github.io/Rootopia/reference/summarize_orders.md)
  : Per-order summary of length and diameter
- [`convert_root_units()`](https://jcunow.github.io/Rootopia/reference/convert_root_units.md)
  : Convert edge-table lengths and diameters to real units
- [`order_classification_map()`](https://jcunow.github.io/Rootopia/reference/order_classification_map.md)
  : Rasterise a per-segment value onto the image grid
- [`prune_terminal_segments()`](https://jcunow.github.io/Rootopia/reference/prune_terminal_segments.md)
  : Prune short or thin terminal segments
- [`render_order_overlay()`](https://jcunow.github.io/Rootopia/reference/render_order_overlay.md)
  : Write an order-colored validation overlay (PNG)
- [`plot_order_window()`](https://jcunow.github.io/Rootopia/reference/plot_order_window.md)
  : Native-resolution validation of a sub-window

## Root turnover

Estimate root production, decay, and turnover from temporal image
comparison or RootDetector DPC images.

- [`root_turnover()`](https://jcunow.github.io/Rootopia/reference/root_turnover.md)
  : Unified Root Turnover Analysis

## Depth mapping

Functions for constructing and binning depth maps from minirhizotron
geometry.

- [`create_depthmap()`](https://jcunow.github.io/Rootopia/reference/create_depthmap.md)
  : Create A Phase-Shifted, Tilt-Amplitude Sine Depth Map
- [`binning()`](https://jcunow.github.io/Rootopia/reference/binning.md)
  : Bin continuous depth values into discrete intervals
- [`depth_zoning()`](https://jcunow.github.io/Rootopia/reference/depth_zoning.md)
  : Mask a scan to a depth zone (depth zone masking)
- [`estimate_soil_surface()`](https://jcunow.github.io/Rootopia/reference/estimate_soil_surface.md)
  : Estimate soil surface position using tape markers
- [`estimate_rotation_center()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_center.md)
  : Estimates rotation from tape coverage
- [`estimate_rotation_shift()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_shift.md)
  : Estimate rotational/depth shift between two root scans

## Soil and colour characterisation

Colour metrics, soil classification, and texture analysis.

- [`create_root_buffer()`](https://jcunow.github.io/Rootopia/reference/create_root_buffer.md)
  : Create a buffer halo) around non-zero pixels
- [`tube_coloration()`](https://jcunow.github.io/Rootopia/reference/tube_coloration.md)
  : Calculate Image Coloration Metrics
- [`analyze_soil_texture()`](https://jcunow.github.io/Rootopia/reference/analyze_soil_texture.md)
  : Texture calculation using Gray-Level Co-occurrence Matrix (GLCM)
- [`classify_soil_rgb()`](https://jcunow.github.io/Rootopia/reference/classify_soil_rgb.md)
  : Classify soil material classes from a minirhizotron RGB raster
- [`build_soil_centroids()`](https://jcunow.github.io/Rootopia/reference/build_soil_centroids.md)
  : Build soil class centroids from manual RGB color picks
- [`plot_soil_classification()`](https://jcunow.github.io/Rootopia/reference/plot_soil_classification.md)
  : Plot the output of classify_soil_rgb

## Distribution analysis

Statistical indices summarising root distribution with depth.

- [`MRD()`](https://jcunow.github.io/Rootopia/reference/MRD.md) :
  Calculate Mean Rooting Depth
- [`RPI()`](https://jcunow.github.io/Rootopia/reference/RPI.md) :
  Calculate Root Penetration Index
- [`rhythmicity()`](https://jcunow.github.io/Rootopia/reference/rhythmicity.md)
  : Assess rhythmicity via sine curve fitting and model comparison
- [`modal_peaks()`](https://jcunow.github.io/Rootopia/reference/modal_peaks.md)
  : Detect and Classify Modes in a Distribution Using Prominence or
  Mclust
- [`circular_mean()`](https://jcunow.github.io/Rootopia/reference/circular_mean.md)
  : Calculate a circular mean to determine average Directionality
- [`root_accumulation()`](https://jcunow.github.io/Rootopia/reference/root_accumulation.md)
  : Calculate root accumulation
- [`compare_depth_distribution()`](https://jcunow.github.io/Rootopia/reference/compare_depth_distribution.md)
  : Compare depth distributions using multiple metrics

## Rotation bias

Detect and correct for rotational artefacts in minirhizotron tube scans.

- [`slice_rotation()`](https://jcunow.github.io/Rootopia/reference/slice_rotation.md)
  : Slice a scan along the rotation (circumferential) axis
- [`fit_sine_curve()`](https://jcunow.github.io/Rootopia/reference/fit_sine_curve.md)
  : Fit a sine curve to data with optional fixed period
- [`rhythmicity()`](https://jcunow.github.io/Rootopia/reference/rhythmicity.md)
  : Assess rhythmicity via sine curve fitting and model comparison

## Visualization

Helper to zoom into an image.

- [`zoom_plot()`](https://jcunow.github.io/Rootopia/reference/zoom_plot.md)
  : Plot a SpatRaster with a magnified native-resolution inset

## Example data

Built-in example images from the Oulanka 2023 field campaign.

- [`seg_Oulanka2023_Session01_T067`](https://jcunow.github.io/Rootopia/reference/seg_Oulanka2023_Session01_T067.md)
  : Segmented Minirhizotron Root Scan - Session 1, Tube 67
- [`seg_Oulanka2023_Session03_T067`](https://jcunow.github.io/Rootopia/reference/seg_Oulanka2023_Session03_T067.md)
  : Segmented Minirhizotron Root Scan - Session 3, Tube 67
- [`skl_Oulanka2023_Session01_T067`](https://jcunow.github.io/Rootopia/reference/skl_Oulanka2023_Session01_T067.md)
  : Skeletonized Root Scan - Session 1, Tube 67
- [`skl_Oulanka2023_Session03_T067`](https://jcunow.github.io/Rootopia/reference/skl_Oulanka2023_Session03_T067.md)
  : Skeletonized Root Scan - Session 3, Tube 67
- [`rgb_Oulanka2023_Session03_T067`](https://jcunow.github.io/Rootopia/reference/rgb_Oulanka2023_Session03_T067.md)
  : Original Minirhizotron Root Scan - Session 3, Tube 67
- [`TurnoverDPC_data`](https://jcunow.github.io/Rootopia/reference/TurnoverDPC_data.md)
  : Root Turnover Analysis Data
- [`flatbed_scan_example`](https://jcunow.github.io/Rootopia/reference/flatbed_scan_example.md)
  : Example Flatbed Root Scan (segmented)
