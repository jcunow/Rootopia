# Skeletonized Root Scan - Session 1, Tube 67

Skeletonized root scan image from a sedge fen in northern Finland (June
2023). The image was processed using RootDetector for segmentation and
skeletonization.

## Usage

``` r
data(skl_Oulanka2023_Session01_T067)
```

## Format

A RasterBrick object with dimensions:

- 4900 columns (width)

- 1144 rows (height)

- 3 layers (channels)

## Source

Images by J.Cunow

## Details

Binary mask representation where:

- Root skeletons = 1

- Background = 0

- Layer 1 includes foreign objects (e.g., tape) marked as 1

Skeletonization reduces root width to single-pixel lines while
preserving the root system topology.

## Author

Johannes Cunow <johannes.cunow@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
  data(skl_Oulanka2023_Session01_T067)
  skl_Oulanka2023_Session01_T067 = terra::rast(skl_Oulanka2023_Session01_T067)
  terra::plot(skl_Oulanka2023_Session01_T067)
} # }
```
