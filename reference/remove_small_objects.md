# Remove small isolated white artifacts from a binary image

Sets to 0 all white regions (value = 1) that are not connected to the
image border.

## Usage

``` r
remove_small_objects(img, max_size = NULL)
```

## Arguments

- img:

  A \`cimg\` binary image (values 0 and 1).

- max_size:

  Maximum artifact size in pixels to remove. If \`NULL\` (default), all
  isolated white regions are removed.

## Value

A \`cimg\` with artifacts removed.
