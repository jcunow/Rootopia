# Remove small isolated white artifacts from a binary image

Sets to 0 all white regions (value = 1) that are not connected to the
image border.

## Usage

``` r
remove_small_objects(img, max_size = NULL, protect_border = FALSE)
```

## Arguments

- img:

  A \`cimg\` binary image (values 0 and 1).

- max_size:

  Maximum artifact size in pixels to remove. If \`NULL\`, all candidate
  white regions are removed.

- protect_border:

  Logical. If \`TRUE\`, white regions touching the image border are
  never removed (they are assumed to be roots leaving the frame). If
  \`FALSE\` (default), border-touching regions are subject to the same
  \`max_size\` test as any other region, so edge specks are removed too.
  Border location confers no protection on its own.

## Value

A \`cimg\` with artifacts removed.
