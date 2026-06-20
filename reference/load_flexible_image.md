# Load an image flexibly from file or convert from memory

Load an image flexibly from file or convert from memory

## Usage

``` r
load_flexible_image(
  input,
  output_format = "cimg",
  scale = c("to_01", "to_255", "binary", "none"),
  select.layer = NULL,
  normalize = NULL,
  binarize = NULL,
  denormalize = NULL
)
```

## Arguments

- input:

  File path or image object

- output_format:

  Character, one out of "cimg", "spatrast", "matrix", "array", "brick",
  "raster", "spatrast", "magick-image". Other spellings are accepted.

- scale:

  Character, the value rescaling to apply. One of \`"to_01"\` (0-255 -\>
  0-1, the default), \`"to_255"\` (0-1 -\> 0-255), \`"binary"\`
  (strictly 0/1), or \`"none"\` (leave values untouched). Each
  conversion is a no-op if the data is already in the target range.

- select.layer:

  Numeric, which layer to select if input has multiple layers

- normalize, binarize, denormalize:

  Deprecated logical flags kept for backward compatibility; use
  \`scale\` instead. \`normalize = TRUE\` maps to \`scale = "to_01"\`,
  \`denormalize = TRUE\` to \`"to_255"\`, and \`binarize = TRUE\` to
  \`"binary"\`. Supplying these together with \`scale\`, or setting both
  \`normalize\` and \`denormalize\`, is an error.
