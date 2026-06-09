# Count pixels (deprecated alias for count_pixels)

\`px.sum()\` is a deprecated alias for \[count_pixels()\]. Please update
your code to use \`count_pixels()\` instead.

## Usage

``` r
px.sum(img, layer = NULL)
```

## Arguments

- img:

  A single-layer raster image.

- layer:

  Ignored. Provided for backward compatibility only; use
  \[terra::subset()\] to select a layer before calling this function.

## Value

A numeric value — the sum of all non-NA pixel values.
