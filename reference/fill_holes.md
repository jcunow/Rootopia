# Fill internal black holes in a binary image

Sets to 1 every black region (value 0) of at most `max_size` pixels that
is enclosed by white (value 1), i.e. does not touch the image border.
Holes are 4-connected: two background pixels that touch only at a corner
belong to different holes, because the 8-connected foreground passes
between them.

## Usage

``` r
fill_holes(img, max_size = NULL)
```

## Arguments

- img:

  A 0/1 matrix.

- max_size:

  Maximum hole size in pixels to fill. `NULL` (default) fills nothing;
  `Inf` fills every hole.

## Value

The matrix with holes filled.
