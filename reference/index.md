# Package index

## Batch processing

High-level wrapper that runs a full depth-profile analysis over a
directory of images.

- [`root_depth_metrics()`](https://jcunow.github.io/RootScanR/reference/root_depth_metrics.md)
  : Compute root traits over a depth profile from segmented
  (mini)rhizotron images

## Image loading and preprocessing

Functions for loading, cleaning, thresholding, and transforming root
images.

- [`load_flexible_image()`](https://jcunow.github.io/RootScanR/reference/load_flexible_image.md)
  : Load an image flexibly from file or convert from memory
- [`clean_image()`](https://jcunow.github.io/RootScanR/reference/clean_image.md)
  : Clean a binary root image
- [`image_threshold()`](https://jcunow.github.io/RootScanR/reference/image_threshold.md)
  : Threshold or deblur an image to binarize features
- [`rgb2gray()`](https://jcunow.github.io/RootScanR/reference/rgb2gray.md)
  : Convert RGB image to grayscale with optimized memory management and
  parallel processing
- [`rotation_censor()`](https://jcunow.github.io/RootScanR/reference/rotation_censor.md)
  : Censor image edges based on rotation
- [`skeletonize_image()`](https://jcunow.github.io/RootScanR/reference/skeletonize_image.md)
  : Skeletonize binary image

## Root traits

Per-image trait extraction: length, diameter, architecture, and spatial
metrics.

- [`root_length()`](https://jcunow.github.io/RootScanR/reference/root_length.md)
  : Root length estimation from skeleton images
- [`root_diameter()`](https://jcunow.github.io/RootScanR/reference/root_diameter.md)
  : Estimate Root Diameters
- [`root_thickness()`](https://jcunow.github.io/RootScanR/reference/root_thickness.md)
  : Approximate average root thickness (deprecated)
- [`detect_skeleton_points()`](https://jcunow.github.io/RootScanR/reference/detect_skeleton_points.md)
  : Detect endpoints and branching points in a skeleton image
- [`count_pixels()`](https://jcunow.github.io/RootScanR/reference/count_pixels.md)
  : Count all pixels in a segmented image
- [`px.sum()`](https://jcunow.github.io/RootScanR/reference/px.sum.md) :
  Count pixels (deprecated alias for count_pixels)
- [`deep_drive()`](https://jcunow.github.io/RootScanR/reference/deep_drive.md)
  : Assess Root Growth Direction Relative to Depth Gradient
- [`root_scape_metrics()`](https://jcunow.github.io/RootScanR/reference/root_scape_metrics.md)
  : RootScapeMetric relies on Landscapemetrics to extract 'Root Scape'
  Features akin to landscape analysis.
- [`root_accumulation()`](https://jcunow.github.io/RootScanR/reference/root_accumulation.md)
  : Calculate root accumulation

## Root branching order

Convert a root skeleton into a segment graph, assign tip/root/branch
order, and summarise root architecture by order class.

- [`branch_order_map()`](https://jcunow.github.io/RootScanR/reference/branch_order_map.md)
  : Branch-order classification of a root skeleton
- [`root_graph_pipeline()`](https://jcunow.github.io/RootScanR/reference/root_graph_pipeline.md)
  : Skeleton-to-graph root ordering pipeline (pixel units)
- [`order_metrics()`](https://jcunow.github.io/RootScanR/reference/order_metrics.md)
  : Aggregate root architecture by order, or split focal-vs-rest
- [`summarize_orders()`](https://jcunow.github.io/RootScanR/reference/summarize_orders.md)
  : Per-order summary of length and diameter
- [`convert_root_units()`](https://jcunow.github.io/RootScanR/reference/convert_root_units.md)
  : Convert edge-table lengths and diameters to real units
- [`order_classification_map()`](https://jcunow.github.io/RootScanR/reference/order_classification_map.md)
  : Rasterise a per-segment value onto the image grid
- [`prune_terminal_segments()`](https://jcunow.github.io/RootScanR/reference/prune_terminal_segments.md)
  : Prune short or thin terminal segments
- [`render_order_overlay()`](https://jcunow.github.io/RootScanR/reference/render_order_overlay.md)
  : Write an order-coloured validation overlay (PNG)
- [`plot_order_window()`](https://jcunow.github.io/RootScanR/reference/plot_order_window.md)
  : Native-resolution validation of a sub-window

## Depth mapping

Functions for constructing and binning depth maps from minirhizotron
geometry.

- [`create_depthmap()`](https://jcunow.github.io/RootScanR/reference/create_depthmap.md)
  : Create A Phase-Shifted, Tilt-Amplitude Sine Depth Map
- [`binning()`](https://jcunow.github.io/RootScanR/reference/binning.md)
  : Bin continuous depth values into discrete intervals
- [`slice_rotation()`](https://jcunow.github.io/RootScanR/reference/slice_rotation.md)
  : Slice a scan along the rotation (circumferential) axis
- [`estimate_soil_surface()`](https://jcunow.github.io/RootScanR/reference/estimate_soil_surface.md)
  : Estimate soil surface position using tape markers
- [`estimate_rotation_center()`](https://jcunow.github.io/RootScanR/reference/estimate_rotation_center.md)
  : Estimates rotation from tape coverage
- [`estimate_rotation_shift()`](https://jcunow.github.io/RootScanR/reference/estimate_rotation_shift.md)
  : Estimate rotational/depth shift between two root scans

## Soil and colour characterisation

Colour metrics, peat classification, and texture analysis.

- [`tube_coloration()`](https://jcunow.github.io/RootScanR/reference/tube_coloration.md)
  : Calculate Image Coloration Metrics
- [`analyze_soil_texture()`](https://jcunow.github.io/RootScanR/reference/analyze_soil_texture.md)
  : Texture calculation using Gray-Level Co-occurrence Matrix (GLCM)
- [`classify_peat_rgb()`](https://jcunow.github.io/RootScanR/reference/classify_peat_rgb.md)
  : Classify peat material classes from a minirhizotron RGB raster
- [`build_peat_centroids()`](https://jcunow.github.io/RootScanR/reference/build_peat_centroids.md)
  : Build peat class centroids from manual RGB colour picks
- [`plot_peat_classification()`](https://jcunow.github.io/RootScanR/reference/plot_peat_classification.md)
  : Plot the output of classify_peat_rgb

## Distribution analysis

Statistical indices summarising root distribution with depth.

- [`MRD()`](https://jcunow.github.io/RootScanR/reference/MRD.md) :
  Calculate Mean Rooting Depth
- [`RPI()`](https://jcunow.github.io/RootScanR/reference/RPI.md) :
  Calculate Root Penetration Index
- [`root_turnover()`](https://jcunow.github.io/RootScanR/reference/root_turnover.md)
  : Unified Root Turnover Analysis
- [`rhythmicity()`](https://jcunow.github.io/RootScanR/reference/rhythmicity.md)
  : Assess rhythmicity via sine curve fitting and model comparison
- [`modal_peaks()`](https://jcunow.github.io/RootScanR/reference/modal_peaks.md)
  : Detect and Classify Modes in a Distribution Using Prominence or
  Mclust
- [`circular_mean()`](https://jcunow.github.io/RootScanR/reference/circular_mean.md)
  : Calculate a circular mean to determine average Directionality
- [`tail_weighted_js_divergence()`](https://jcunow.github.io/RootScanR/reference/tail_weighted_js_divergence.md)
  : Calculate tail-weighted Jensen-Shannon divergence
- [`tail_weighted_kl_divergence()`](https://jcunow.github.io/RootScanR/reference/tail_weighted_kl_divergence.md)
  : Calculate tail-weighted KL divergence for discrete distributions
- [`tail_weighted_wasserstein_distance()`](https://jcunow.github.io/RootScanR/reference/tail_weighted_wasserstein_distance.md)
  : A tailweighted Version of 1 dimensional Wasserstein distance betwwen
  two probability vectors

## Rotation bias

Detect and correct for rotational artefacts in minirhizotron tube scans.

- [`fit_sine_curve()`](https://jcunow.github.io/RootScanR/reference/fit_sine_curve.md)
  : Fit a sine curve to data with optional fixed period
- [`create_root_buffer()`](https://jcunow.github.io/RootScanR/reference/create_root_buffer.md)
  : Create a buffer halo) around non-zero pixels
- [`create_depthmap()`](https://jcunow.github.io/RootScanR/reference/create_depthmap.md)
  : Create A Phase-Shifted, Tilt-Amplitude Sine Depth Map

## Example data

Built-in example images from the Oulanka 2023 field campaign.

- [`seg_Oulanka2023_Session01_T067`](https://jcunow.github.io/RootScanR/reference/seg_Oulanka2023_Session01_T067.md)
  : Segmented Minirhizotron Root Scan - Session 1, Tube 67
- [`seg_Oulanka2023_Session03_T067`](https://jcunow.github.io/RootScanR/reference/seg_Oulanka2023_Session03_T067.md)
  : Segmented Minirhizotron Root Scan - Session 3, Tube 67
- [`skl_Oulanka2023_Session01_T067`](https://jcunow.github.io/RootScanR/reference/skl_Oulanka2023_Session01_T067.md)
  : Skeletonized Root Scan - Session 1, Tube 67
- [`skl_Oulanka2023_Session03_T067`](https://jcunow.github.io/RootScanR/reference/skl_Oulanka2023_Session03_T067.md)
  : Skeletonized Root Scan - Session 3, Tube 67
- [`rgb_Oulanka2023_Session03_T067`](https://jcunow.github.io/RootScanR/reference/rgb_Oulanka2023_Session03_T067.md)
  : Original Minirhizotron Root Scan - Session 3, Tube 67
- [`TurnoverDPC_data`](https://jcunow.github.io/RootScanR/reference/TurnoverDPC_data.md)
  : Root Turnover Analysis Data
