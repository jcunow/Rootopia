# Texture calculation

Texture calculation

## Usage

``` r
analyze_soil_texture(
  img.color,
  grays = 7,
  window = c(9, 9),
  metrics = c("variance", "second_moment")
)
```

## Arguments

- img.color:

  Three-band raster or path to image

- grays:

  Number of gray levels

- window:

  Window size for GLCM

- metrics:

  Texture metrics to calculate

## Value

Raster with texture metrics

## Examples

``` r
data(rgb_Oulanka2023_Session03_T067)
img = raster::brick(rgb_Oulanka2023_Session03_T067)
analyze_soil_texture(img, 7, c(9,9), metrics = "second_moment")
#> class      : RasterLayer 
#> dimensions : 1160, 4899, 5682840  (nrow, ncol, ncell)
#> resolution : 0.0002041233, 0.000862069  (x, y)
#> extent     : 0, 1, 0, 1  (xmin, xmax, ymin, ymax)
#> crs        : NA 
#> source     : memory
#> names      : glcm_second_moment 
#> values     : 0.1111111, 1  (min, max)
#> 
```
