

---

## RootScanR: Extract Root Traits from Images
output: distill::distill\_article
description: "R tools for root trait extraction from flatbed and minirhizotron images."
---

[![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://jcunow.github.io/RootScanR/)

---

Welcome to **RootScanR**, an R package designed to extract quantitative belowground traits from **flatbed root scans** and **minirhizotron image sequences**. It provides flexible image preprocessing, trait estimation like root length & diameter, depth mapping, and distribution tools. 

 - some features are still experimental. At the moment, you will still need to provide a skeleton created elewhere (for example RootDetector). 


## 💡 What can RootScanR do?

- 🖼️ Clean and process root scans from (mini)rhizotrons and flatbed scanners
- 🌱 Estimate root length
- 📏 Measure root diameters and distributions
- 🕸️ Analyze root architecture (branching point, root tips, root angle distribution)
- 🥕 Spatial-spatial indices (root angle distribution, rotation bias)
---

## 📘 Tutorials

Jump straight into hands-on tutorials:

* 📄 [Flatbed Scan Workflow](articles/FlatBedScans_vignettes.html)
* 📄 [Minirhizotron Workflow](articles/MinirhizotronScans_vignettes.html)
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

