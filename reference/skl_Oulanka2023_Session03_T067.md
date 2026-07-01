# Skeletonized Root Scan - Session 3, Tube 67

Skeletonized root scan image from a sedge fen in northern Finland
(October 2023). The image was processed using RootDetector for
segmentation and skeletonization.

## Usage

``` r
data(skl_Oulanka2023_Session03_T067)
```

## Format

A RasterBrick object with dimensions:

- 4900 columns (width)

- 1161 rows (height)

- 3 layers (channels)

## Source

Original minirhizotron images acquired by Johannes Cunow and Gesche
Blume-Werry at Oulanka Research Station, Finland (October 2023). Images
were processed using RootDetector.

## Details

Binary mask representation where:

- Root skeletons = 1

- Background = 0

- Layer 1 includes foreign objects (e.g., tape) marked as 1

Skeletonization reduces root width to single-pixel lines while
preserving the root system topology.

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
