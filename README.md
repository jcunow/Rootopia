

---

## RootScanR: Extract Root Traits from Images
output: distill::distill\_article
description: "R tools for root trait extraction from flatbed and minirhizotron images."
---

[![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://jcunow.github.io/RootScanR/)

---

Welcome to **RootScanR**, an R package designed to extract quantitative belowground traits from **flatbed root scans** and **minirhizotron image sequences**. It provides flexible image preprocessing, trait estimation like root length & diameter, depth mapping, and distribution tools.




## 💡 What can RootScanR do?

- 🖼️ Clean and process binary root images
- 🌱 Estimate root length (skeletonize and Kimura Root Length)
- 📏 Measure root diameters and modal distributions
- 🕸️ Analyze root architecture (branching points, MST, loopiness)
- 🥕 Spatial-spatial indecies (root angle distribution, rotation bias)
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

## 📚 Reference & Citation

To cite `RootScanR`, please use the relevant methods cited in the vignettes:

* Kimura et al. (1999) – Pixel-based root length estimation
* Zhang & Suen (1984) – Skeletonization
* Otsu (1979) – Thresholding
* Scrucca et al. (2016) – Model-based clustering (Mclust)

You can also generate the full citation in R:

```r
citation("RootScanR")
```

---

## 📢 News & Updates

See the [News page](news/index.html) for the latest changes and version notes.

---

## 🔍 Documentation

* [Function reference](reference/index.html)
* [Vignettes](articles/index.html)
* [Source code on GitHub](https://github.com/jcunow/RootScanR)

---

### 🙌 Contributions & Issues

Found a bug or have a feature request? Please [open an issue](https://github.com/jcunow/RootScanR/issues) or submit a pull request.

