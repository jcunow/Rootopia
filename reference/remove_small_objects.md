# Remove small white artifacts from binary images

Identifies small internal white regions (objects) in a binary image and
removes them by setting their pixel values to 0. Artifacts are defined
as white areas (value = 1) not connected to the image border.

## Usage

``` r
remove_small_objects(img, max_size = NULL)
```

## Arguments

- img:

  A \`cimg\` object representing a binary image (values 0 and 1).

- max_size:

  Optional maximum size (in pixels) of white objects to remove. If
  \`NULL\`, all isolated objects are removed.

## Value

A \`cimg\` object with small white artifacts removed (set to 0).
