# Estimate Root Diameters

This function estimates root diameters and root volume from an input
image using skeletonization and distance transform methods. The input
can be a file path, raster, image object, or array, which is converted
to a binary image before processing.

## Usage

``` r
root_diameter(
  img,
  skeleton_method = "GuoHall",
  select.layer = NULL,
  diagnostics = FALSE,
  unit = "cm",
  dpi = 300
)
```

## Arguments

- img:

  A character string (file path), \`SpatRaster\`, \`RasterBrick\`,
  \`RasterLayer\`, \`cimg\`, \`magick-image\`, or array. The input image
  to process.

- skeleton_method:

  Character. The method to use for skeletonization. Default is
  \`"Guo-Hall"\`.

- select.layer:

  Integer. Specifies which layer to use if the input is a multi-band
  image. Default is \`2\`.

- diagnostics:

  Logical. If \`TRUE\`, enables diagnostic plots and logging. Default is
  \`FALSE\`.

- unit:

  output in pixel 'px', 'inch' or in 'cm'

- dpi:

  scan resolution. Only used if unit = 'cm' or 'inch'

## Value

A list containing:

- quantiles:

  Numeric vector of diameter quantiles (10th to 100th percentile).

- mean_diameter:

  Numeric. The mean root diameter.

- median_diameter:

  Numeric. The median root diameter.

- diameters:

  Numeric vector of all diameter values in the skeletonized regions.

- skeleton_rast:

  \`SpatRaster\`. Binary raster mask of skeletonized regions.

- diameter_rast:

  \`SpatRaster\`. Raster showing diameters in the skeletonized regions.

- distance_map_rast:

  \`SpatRaster\`. Raster showing the distance transform values.

- root_volume:

  Numeric. The sum of root volume - assuming cylindrical roots

## Details

The function works as follows: - Converts the input image to a binary
format (\`cimg\`). - Applies a distance transform to compute the
Euclidean distance for the foreground (root) pixels. - Skeletonizes the
binary image to identify root centerlines. - Filters distance values to
retain only those corresponding to the skeletonized regions. - Computes
diameter statistics, including quantiles, mean, and median diameters.

The function supports various input formats and normalizes image values
to the range \[0, 1\] if needed. It uses the \`terra\` package for
raster operations and the \`imager\` package for image processing.

## Examples

``` r
# Example usage:
data(seg_Oulanka2023_Session01_T067)
result <- root_diameter(img = seg_Oulanka2023_Session01_T067,
  skeleton_method = "GuoHall", select.layer = 2, unit = "px",
  diagnostics = TRUE)
#> 
#> Applying method: GuoHall 
#> Image dimensions: 1144 x 4900 
#> Initial foreground pixels: 200781 
#> Iteration 1 : Removed 95916 pixels
#> Iteration 2 : Removed 2872 pixels
#> Iteration 3 : Removed 822 pixels
#> Iteration 4 : Removed 326 pixels
#> Iteration 5 : Removed 149 pixels
#> Iteration 6 : Removed 66 pixels
#> Iteration 7 : Removed 29 pixels
#> Iteration 8 : Removed 19 pixels
#> Iteration 9 : Removed 16 pixels
#> Iteration 10 : Removed 8 pixels
#> Iteration 11 : Removed 6 pixels
#> Iteration 12 : Removed 3 pixels
#> Final foreground pixels: 100549 
#> Total iterations: 13 
#> Processing complete. Summary statistics:
#> Mean diameter: 2.21
#> Median diameter: 2.00
#> Number of valid measurements: 100549

# Access results:
print(result$mean_diameter)
#> [1] 2.20623
terra::plot(result$skeleton_rast)

```
