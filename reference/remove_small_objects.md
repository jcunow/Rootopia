# Remove small white objects from a binary image

Sets to 0 every white region (value 1) of at most `max_size` pixels.
Regions are 8-connected: pixels that touch at a corner belong to the
same object, so a diagonal run of root pixels is one object. With
`protect_border = TRUE`, regions touching the image border are kept
whatever their size.

## Usage

``` r
remove_small_objects(img, max_size = NULL, protect_border = FALSE)
```

## Arguments

- img:

  A 0/1 matrix.

- max_size:

  Maximum object size in pixels to remove. `NULL` (default) removes
  nothing; `Inf` removes every candidate region.

- protect_border:

  Logical. If `TRUE`, white regions touching the image border are never
  removed (they are assumed to be roots leaving the frame). If `FALSE`
  (default), border-touching regions are subject to the same `max_size`
  test as any other region, so edge specks are removed too.

## Value

The matrix with objects removed.
