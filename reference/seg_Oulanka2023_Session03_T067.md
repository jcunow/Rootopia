# Segmented Minirhizotron Root Scan - Session 3, Tube 67

Segmented root scan image from a sedge fen in northern Finland (October
2023). The image was processed using RootDetector for root segmentation.

## Usage

``` r
data(seg_Oulanka2023_Session03_T067)
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

Binary mask, stored as 0 (background) and 255 (foreground):

- Layer 1: foreign objects such as tape

- Layers 2 and 3 (identical): roots

Spatial dimensions correspond to physical tube measurements: columns =
tube length, rows = tube rotation.

## References

Peters B, Blume-Werry G, Gillert A, et al. (2023) RootDetector: a
convolutional neural network for root detection. Scientific Reports
13:1399. https://doi.org/10.1038/s41598-023-28400-x

## Author

Johannes Cunow

## Examples

``` r
if (FALSE) { # \dontrun{
  data(seg_Oulanka2023_Session03_T067)
  seg_Oulanka2023_Session03_T067 = terra::rast(seg_Oulanka2023_Session03_T067)
  terra::plot(seg_Oulanka2023_Session03_T067, maxcell = Inf)
} # }
```
