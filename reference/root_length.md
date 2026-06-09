# Calculate Root Length using Kimura's Method with optimizations

Calculate Root Length using Kimura's Method with optimizations

## Usage

``` r
root_length(img, unit = "cm", dpi = 300, select.layer = 2)
```

## Arguments

- img:

  A skeletonized root image raster

- unit:

  Output unit ("px", "cm", or "inch")

- dpi:

  Image resolution (required when unit = "cm" or "inch")

- select.layer:

  Integer. Specifies which layer to use if the input is a multi-band
  image. Default is \`2\`, matching the RootDetector output format where
  layer 2 contains the root channel.

## Value

Numeric value representing root length in specified unit

## Examples

``` r
data(skl_Oulanka2023_Session01_T067)
img <- terra::rast(skl_Oulanka2023_Session01_T067)
RL <- root_length(img = img, unit = "cm", dpi = 300, select.layer = 2)
```
