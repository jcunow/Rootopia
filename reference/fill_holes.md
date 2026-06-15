# Fill internal black holes in a binary image

Sets to 1 all black regions (value = 0) that are completely surrounded
by white (value = 1) – i.e. not connected to the image border.

## Usage

``` r
fill_holes(img, max_size = NULL)
```

## Arguments

- img:

  A \`cimg\` binary image (values 0 and 1).

- max_size:

  Maximum hole size in pixels to fill. If \`NULL\` (default), all holes
  are filled regardless of size.

## Value

A \`cimg\` with holes filled.
