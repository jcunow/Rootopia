# Validate Image Input Parameters

Internal function to validate input parameters for image processing
functions Used for image skeletonizing

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

  The input image to validate

- allow_empty:

  Logical, whether to allow empty images

- min_dim:

  Minimum required dimensions

- require_binary:

  Logical, whether to require binary values

- select.layer:

  Layer to validate for multi-layer images

## Value

List containing validated and processed image data
