

---

## Rootopia: R Tool to Extract Root Traits from Flatbed and Minirhizotron Images


[![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://jcunow.github.io/Rootopia/)



Welcome to **Rootopia**, an R package designed to extract quantitative belowground traits from **flatbed root scans** and **minirhizotron image sequences**. It provides flexible image preprocessing, trait estimation like root length & diameter, depth mapping, and distribution tools. 

 - some features are still experimental. This package is an insect's playground and we love bugs here. Just call them out if their behaviour is undesirable.


## 💡 What can Rootopia do?

- ️ Clean and process root scans from (mini)rhizotrons and flatbed scanners
-  Estimate root length, depth profiles, and turnover between sessions
-  Measure root diameters and diameter distributions
- ️ Analyze root architecture: branching points, root tips, branch/root order, and root angle distribution
-  Spatial indices (root distribution with depth, rotation bias)
-  Stitch overlapping scan sequences into one mosaic per tube

---

## 📘 Tutorials

New to Rootopia? Start with one of the step-by-step tutorials — they walk
through loading, depth mapping, and trait extraction one function at a time,
so you learn what each stage does. Once you are comfortable, the
batch-processing tutorial wraps the whole minirhizotron pipeline and flatbed scan analysis over a folder of images into a single function call.

* 📄 [Start here: Minirhizotron Workflow (step-by-step)](articles/MinirhizotronScans_vignettes.html)
* 📄 [Flatbed Scan Workflow (step-by-step)](articles/FlatBedScans_vignettes.html)
* 📄 [Stitching Scan Sequences into Mosaics](articles/Stitching_vignette.html)
* 📄 [Rotation Bias Correction](articles/Rotation_Bias_vignettes.html)
* 📄 [Batch Processing (whole folders at once)](articles/BatchProcessing_vignette.html)

Each tutorial includes code, images, overlays, and tips for interpretation.

---

<br />

## 🚀 Quick Start

```r
# Load example (segmented) minirhizotron scan
library(Rootopia)
library(terra)
data("seg_Oulanka2023_Session01_T067")
img <- load_flexible_image(img, output_format = "spatrast", scale = "binary")

# Preprocess and extract traits
tape_mask <- img[[2]] - img[[1]] 
terra::values(tape_mask)[terra::values(tape_mask) != -1] <- NA
binary <- load_flexible_image(img, select_layer = 2)
cleaned <- clean_image(binary, max_artifact_size = 5, max_hole_size = 5)
skel <- skeletonize_image(cleaned, verbose = F)

# Length and diameter
root_length(skel, unit = "cm", dpi = 150)
diam_map <- root_diameter(cleaned, unit = "cm", dpi = 150 )
hist(diam_map$diameters)


# Architecture
metrics <- detect_skeleton_points(skel)
root_tips = count_pixels(metrics$endpoints)
root_branches = count_pixels(metrics$branching_points)


# Branch order (more useful for flatbed scans):
order_res <- branch_order_map(skel, mask = cleaned, order = "root_order", unit = "cm", dpi = 150)
order_metrics(order_res, focal = "thinnest")

# Visualize
zoom_plot(diam_map$diameter_rast)
zoom_plot(skel)
zoom_plot(order_res$class_map)
```

---





## 🔍 Documentation

* [Function reference](reference/index.html)
* [Vignettes](articles/index.html)
* [Source code on GitHub](https://github.com/jcunow/Rootopia)

---

### 🙌 Contributions & Issues

Found a bug or have a feature request? Please [open an issue](https://github.com/jcunow/Rootopia/issues) or submit a pull request.

