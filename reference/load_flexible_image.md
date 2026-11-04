# Load an image flexibly from file or convert from memory

Load an image flexibly from file or convert from memory

## Usage

``` r
load_flexible_image(
  input,
  output_format = "cimg",
  scale = c("to_01", "to_255", "binary", "none"),
  select_layer = NULL,
  normalize = NULL,
  binarize = NULL,
  denormalize = NULL
)
```

## Arguments

- input:

  File path or image object

- output_format:

  Character, one of "cimg", "spatrast", "matrix", "array", "brick",
  "raster", "magick-image". Capitalization does not matter, and
  "SpatRaster", "RasterBrick", "RasterLayer" and "magick" are accepted
  as well.

- scale:

  Character, the value rescaling to apply. One of `"to_01"` (0-255 -\>
  0-1, the default), `"to_255"` (0-1 -\> 0-255), `"binary"` (every
  positive value becomes 1), or `"none"` (leave values untouched).
  `"to_01"` and `"to_255"` leave data that is already in the target
  range untouched.

- select_layer:

  Numeric, which layer to select if input has multiple layers

- normalize, binarize, denormalize:

  Deprecated logical flags kept for backward compatibility; use `scale`
  instead. `normalize = TRUE` maps to `scale = "to_01"`,
  `denormalize = TRUE` to `"to_255"`, and `binarize = TRUE` to
  `"binary"`. Supplying these together with `scale`, or setting both
  `normalize` and `denormalize`, is an error.

## Value

The image in the requested format. Arrays and matrices are (rows,
columns, layers); a SpatRaster built from an array spans `0..ncol` by
`0..nrow`, while a SpatRaster passed in with `scale = "none"` and no
`select_layer` is returned unchanged. Reading `.tif`/`.tiff` files needs
the tiff package.
