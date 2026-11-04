# Skeletonized Root Scan - Session 3, Tube 67

Skeletonized root scan image from a sedge fen in northern Finland
(October 2023). The image was processed using RootDetector for
segmentation and skeletonization.

## Usage

``` r
data(skl_Oulanka2023_Session03_T067)
```

## Format

A numeric array (rows x columns x layers), values 0-255. Rebuild a
SpatRaster with
[`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html).
Dimensions:

- 4900 columns (width)

- 1161 rows (height)

- 3 layers (channels)

## Source

Original minirhizotron images acquired by Johannes Cunow and Gesche
Blume-Werry at Oulanka Research Station, Finland (October 2023). Images
were processed using RootDetector.

## Details

Binary mask, stored as 0 (background) and a single non-zero value
(foreground):

- Layer 1: foreign objects such as tape, filled

- Layers 2 and 3 (identical): root skeleton

At the stored resolution the skeleton's one-pixel lines are broken into
many short pieces. For traits that depend on connectivity (tips, branch
points, branching order), skeletonize the matching segmentation with
[`skeletonize_image()`](https://jcunow.github.io/Rootopia/reference/skeletonize_image.md).

## References

Peters B, Blume-Werry G, Gillert A, et al. (2023) RootDetector: a
convolutional neural network for root detection. Scientific Reports
13:1399. https://doi.org/10.1038/s41598-023-28400-x

## Author

Johannes Cunow

## Examples

``` r
if (FALSE) { # \dontrun{
  data(skl_Oulanka2023_Session03_T067)
  skl_Oulanka2023_Session03_T067 = terra::rast(skl_Oulanka2023_Session03_T067)
  terra::plot(skl_Oulanka2023_Session03_T067, maxcell = Inf)
} # }
```
