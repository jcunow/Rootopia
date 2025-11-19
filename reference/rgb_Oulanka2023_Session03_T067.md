# Original Minirhizotron Root Scan - Session 3, Tube 67

Original RGB root scan image from a sedge fen in northern Finland
(October 2023). This image represents a composite of multiple stitched
scans.

## Usage

``` r
data(rgb_Oulanka2023_Session03_T067)
```

## Format

A RasterBrick object with dimensions:

- 4900 columns (width)

- 1161 rows (height)

- 3 layers (RGB channels)

## Source

Images by J.Cunow

## Details

The image represents a complete tube scan where:

- Columns correspond to tube length

- Rows correspond to tube rotation

- RGB channels represent true color information

## Author

Johannes Cunow <johannes.cunow@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
  data(rgb_Oulanka2023_Session03_T067)
  rgb_Oulanka2023_Session03_T067 = terra::rast(rgb_Oulanka2023_Session03_T067)
  terra::plotRGB(rgb_Oulanka2023_Session03_T067)
} # }
```
