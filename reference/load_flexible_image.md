# Load an image flexibly from file or convert from memory

Load an image flexibly from file or convert from memory

## Usage

``` r
load_flexible_image(
  input,
  output_format = "cimg",
  normalize = TRUE,
  select.layer = NULL,
  binarize = FALSE
)
```

## Arguments

- input:

  File path or image object

- output_format:

  Character, one out of "cimg", "spatrast", "matrix", "array", "brick",
  "raster", "spatrast", "magick-image". Other spellings are accepted.

- normalize:

  Logical, whether to normalize values to 0-1 range if they're in 0-255

- select.layer:

  Numeric, which layer to select if input has multiple layers

- binarize:

  Logical, whether the output is strictly 0 and 1. Overwrites normalize
