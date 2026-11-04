# Root Turnover Analysis Data

Root turnover analysis from a sedge fen in northern Finland, comparing
root presence between two time points (June 2023 and October 2023) using
RootDetector root tracking feature.

## Usage

``` r
data(TurnoverDPC_data)
```

## Format

A numeric array (rows x columns x layers), values 0-255. Rebuild a
SpatRaster with
[`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html).
Dimensions:

- 2273 rows (height)

- 2550 columns (width)

- 3 layers (red, green, blue)

## Source

Original minirhizotron images acquired by Johannes Cunow and Gesche
Blume-Werry at Oulanka Research Station, Finland (June and October
2023). Images were processed using RootDetector.

## Details

An RGB image in which each turnover class has its own color:

|                                     |            |               |
|-------------------------------------|------------|---------------|
| **Class**                           | **Color**  | **R, G, B**   |
| background                          | dark slate | 39, 54, 59    |
| foreign objects (e.g. tape)         | red        | 255, 0, 0     |
| constant (present at both dates)    | white      | 255, 255, 255 |
| production (new at the second date) | green      | 96, 209, 130  |
| decay (gone at the second date)     | rose       | 226, 106, 116 |

`root_turnover(method = "dpc")` reads these classes with its default
`decay_layer = 1` and `product_layer = 2`.

## References

Peters B, Blume-Werry G, Gillert A, et al. (2023) RootDetector: a
convolutional neural network for root detection. Scientific Reports
13:1399. https://doi.org/10.1038/s41598-023-28400-x

## Author

Johannes Cunow

## Examples

``` r
if (FALSE) { # \dontrun{
  data(TurnoverDPC_data)
  img <- terra::rast(TurnoverDPC_data)
  terra::plotRGB(img, maxcell = Inf)
  root_turnover(img, method = "dpc")
} # }
```
