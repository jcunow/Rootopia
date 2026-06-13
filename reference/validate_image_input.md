# Validate image input

Checks structural validity and (optionally) enforces binary constraints
for image processing routines.

## Usage

``` r
validate_image_input(
  img,
  allow_empty = FALSE,
  min_dim = c(3, 3),
  require_binary = TRUE,
  select.layer = NULL
)
```

## Arguments

- img:

  Input image (SpatRaster or compatible object)

- allow_empty:

  Logical. Allow all-zero images.

- min_dim:

  Minimum allowed image dimensions.

- require_binary:

  Logical. Enforce binary values (0/1).

- select.layer:

  Layer index for multi-layer inputs.

## Value

List containing validated image and metadata.
