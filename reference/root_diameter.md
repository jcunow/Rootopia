# Estimate Root Diameters

This function estimates root diameters and root volume from an input
image using skeletonization and distance transform methods. The input
can be a file path, raster, image object, or array, which is converted
to a binary image before processing.

## Usage

``` r
root_diameter(
  img,
  skeleton_method = "MAT",
  skeleton.img = NULL,
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
  \`"MAT"\`. Will be skipped if skeleton \`SpatRaster\`is provided.

- skeleton.img:

  A character string (file path), \`SpatRaster\`, \`RasterBrick\`,
  \`RasterLayer\`, \`cimg\`, \`magick-image\`, or array. Uses this
  object instead of computing it from scratch.

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
  skeleton_method = "MAT", select.layer = 2, unit = "px",
  diagnostics = TRUE)
#> 
#> Applying method: MAT 
#> 
#> Distance transform computed
#> Processing complete. Summary statistics:
#> Mean diameter: 2.14
#> Median diameter: 2.00
#> Number of valid measurements: 161904

# Access results:
print(result$mean_diameter)
#> [1] 2.14191
terra::plot(result$skeleton_rast)

```
