# Root Turnover Analysis Data

Root turnover analysis from a sedge fen in northern Finland, comparing
root presence between two time points (June 2023 and October 2023) using
RootDetector root tracking feature.

## Usage

``` r
data(TurnoverDPC_data)
```

## Format

A RasterBrick object with dimensions:

- 2550 columns (width)

- 2273 rows (height)

- 3 layers (decay, growth, persistent)

## Source

Original minirhizotron images acquired by Johannes Cunow and Gesche
Blume-Werry at Oulanka Research Station, Finland (June 2023). Images
were processed using RootDetector.

## Details

Multi-layer representation of root dynamics where:

- Layer 1: Root decay (value = 1)

  - Roots that disappeared between time points

  - Foreign objects (e.g., tape)

  - Persistent roots

- Layer 2: Root growth (value = 1)

  - New roots that appeared

  - Persistent roots

- Layer 3: Persistent roots (value = 1)

  - Roots present in both time points

Background is represented as 0 in all layers.

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
  img = terra::rast(TurnoverDPC_data)
  # Plot individual layers
  terra::plot(img, maxcell = Inf)
  # Calculate turnover statistics
  decay <- sum(terra::values(img[[1]]) == 1)
  growth <- sum(terra::values(img[[2]]) == 1)
  persistent <- sum(terra::values(img[[3]]) == 1)
} # }
```
