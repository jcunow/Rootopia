

---

## RootScanR: R Tool to Extract Root Traits from Flatbed and Minirhizotron Images


[![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://jcunow.github.io/RootScanR/)



Welcome to **RootScanR**, an R package designed to extract quantitative belowground traits from **flatbed root scans** and **minirhizotron image sequences**. It provides flexible image preprocessing, trait estimation like root length & diameter, depth mapping, and distribution tools. 

 - some features are still experimental. This package is an insect's playground and we love bugs here. Just call them out if their behaviour is undesirable.


## 💡 What can RootScanR do?

- 🖼️ Clean and process root scans from (mini)rhizotrons and flatbed scanners
- 🌱 Estimate root length, depth profiles, and turnover between sessions
- 📏 Measure root diameters and diameter distributions
- 🕸️ Analyze root architecture: branching points, root tips, branch/root order
  (main axis vs. laterals), and root angle distribution
- 🥕 Spatial indices (root distribution with depth, rotation bias)
---

## 📘 Tutorials

New to RootScanR? Start with the batch-processing tutorial — it runs the
whole pipeline (loading, depth mapping, trait extraction) over a folder of
images with a single function call. The step-by-step tutorials below are for
when you need full control over an individual processing stage.

* 🚀 [Start here: Batch Processing](articles/BatchProcessing_vignette.html)
* 📄 [Minirhizotron Workflow (step-by-step)](articles/MinirhizotronScans_vignettes.html)
* 📄 [Flatbed Scan Workflow (step-by-step)](articles/FlatBedScans_vignettes.html)
* 📄 [Rotation Bias Correction](articles/Rotation_Bias_vignettes.html)

Each tutorial includes code, images, overlays, and tips for interpretation.

---

<br />

## 🚀 Quick Start

```r
# Load example flatbed scan
library(RootScanR)
library(terra)
data("seg_Oulanka2023_Session01_T067")
img <- seg_Oulanka2023_Session01_T067

# Preprocess and extract traits
binary <- load_flexible_image(img)
cleaned <- clean_image(binary, max_artifact_size = 5, max_hole_size = 5, select.layer = 2)
skel <- skeletonize_image(cleaned,  methods = "GuoHall", verbose = F)

# Length and diameter
root_length(skel, unit = "cm", dpi = 150)
diam_map <- root_diameter(cleaned, unit = "cm", dpi = 150, select.layer = NULL )
modal_peaks(diam_map$diameters, prominence_threshold = 10, mclust = F, adjust = 5, display_type = "density")


# Architecture
metrics <- detect_skeleton_points(skel)
print(metrics)

# Branch order: main axis (order 1) vs. lateral roots (order >= 2)
order_res <- branch_order_map(skel, mask = cleaned, order = "branch_order", unit = "cm", dpi = 150)
order_metrics(order_res, focal = "thickest")

# Visualize
terra::plot(diam_map$diameter_rast)
terra::plot(skel)
```

---





## 🔍 Documentation

* [Function reference](reference/index.html)
* [Vignettes](articles/index.html)
* [Source code on GitHub](https://github.com/jcunow/RootScanR)

---

### 🙌 Contributions & Issues

Found a bug or have a feature request? Please [open an issue](https://github.com/jcunow/RootScanR/issues) or submit a pull request.

