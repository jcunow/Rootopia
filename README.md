# Rootopia: Extract Root Traits from Flatbed and Minirhizotron Images

[![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://jcunow.github.io/Rootopia/)

Rootopia extracts quantitative belowground traits from **flatbed root scans** and
**minirhizotron image sequences**. It provides flexible image preprocessing, trait
estimation such as root length and diameter, depth mapping, and distribution tools.

Some features are still experimental. This package is an insect's playground and we
love bugs here — please call them out if their behaviour is undesirable.

## What Rootopia does

- Cleans and processes root scans from (mini)rhizotrons and flatbed scanners
- Estimates root length, depth profiles, and turnover between sessions
- Measures root diameters and diameter distributions
- Analyses root architecture: branching points, root tips, branch and root order,
  and root angle distribution
- Computes spatial indices for root distribution with depth and for rotation bias
- Stitches overlapping scan sequences into one mosaic per tube

## Tutorials

New to Rootopia? Start with one of the step-by-step tutorials — they walk through
loading, depth mapping, and trait extraction one function at a time, so you learn
what each stage does. Once you are comfortable, the batch-processing tutorial wraps
the whole minirhizotron pipeline and flatbed scan analysis over a folder of images
into a single function call.

- [Start here: Minirhizotron Workflow (step-by-step)](articles/MinirhizotronScans_vignettes.html)
- [Flatbed Scan Workflow (step-by-step)](articles/FlatBedScans_vignettes.html)
- [Stitching Scan Sequences into Mosaics](articles/Stitching_vignette.html)
- [Rotation Bias Correction](articles/Rotation_Bias_vignettes.html)
- [Special Topics: soil colour, texture, halos, turnover](articles/SpecialTopics_vignette.html)
- [Batch Processing (whole folders at once)](articles/BatchProcessing_vignette.html)

Each tutorial includes code, images, overlays, and tips for interpretation.

## Quick start

```r
library(Rootopia)
library(terra)

# A segmented minirhizotron scan ships with the package.
# select_layer = 2 is the root channel; scale = "binary" maps it to 0/1.
data("seg_Oulanka2023_Session01_T067")
img <- load_flexible_image(seg_Oulanka2023_Session01_T067,
                           output_format = "spatrast",
                           scale         = "binary",
                           select_layer  = 2)

# Preprocess and skeletonise
cleaned <- clean_image(img, max_artifact_size = 5, max_hole_size = 5)
skel    <- skeletonize_image(cleaned, verbose = FALSE)

# Length and diameter
root_length(skel, unit = "cm", dpi = 150)
diam_map <- root_diameter(cleaned, unit = "cm", dpi = 150)
hist(diam_map$diameters)

# Architecture
points        <- detect_skeleton_points(skel)
root_tips     <- count_pixels(points$endpoints)
root_branches <- count_pixels(points$branching_points)

# Branch order (most informative for flatbed scans)
order_res <- branch_order_map(skel, mask = cleaned, order = "root_order",
                              unit = "cm", dpi = 150)
order_res$summary
order_metrics(order_res, focal = "thinnest")

# Visualise
zoom_plot(diam_map$diameter_rast)
zoom_plot(skel)
zoom_plot(order_res$class_map)
```

## Documentation

- [Function reference](reference/index.html)
- [Vignettes](articles/index.html)
- [Source code on GitHub](https://github.com/jcunow/Rootopia)

## Contributions and issues

Found a bug or have a feature request? Please
[open an issue](https://github.com/jcunow/Rootopia/issues) or submit a pull request.
