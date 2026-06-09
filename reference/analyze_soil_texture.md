# Texture calculation using Gray-Level Co-occurrence Matrix (GLCM)

Texture calculation using Gray-Level Co-occurrence Matrix (GLCM)

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

  Three-band raster (RGB) or path to image. Internally converted to a
  \`raster::RasterBrick\` as required by the \`glcm\` package.

- grays:

  Number of gray levels for GLCM quantization. Must be between 2
  and 255. Default is 7.

- window:

  Window size for GLCM calculation as a length-2 vector of odd positive
  integers, e.g. \`c(9, 9)\`. Default is \`c(9, 9)\`.

- metrics:

  Character vector of GLCM texture statistics to calculate. Valid
  options: "mean", "variance", "homogeneity", "contrast",
  "dissimilarity", "entropy", "second_moment", "correlation".

## Value

A RasterLayer (or RasterBrick for multiple metrics) with texture values.

## Examples

``` r
data(rgb_Oulanka2023_Session03_T067)
img <- raster::brick(rgb_Oulanka2023_Session03_T067)
analyze_soil_texture(img, grays = 7, window = c(9, 9),
                     metrics = "second_moment")
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
