# Extract Root Decay, New Root Production, and No-Change Roots (only 'RootDetector' images)

Extract Root Decay, New Root Production, and No-Change Roots (only
'RootDetector' images)

## Usage

``` r
turnover_dpc(
  img,
  product_layer = 2,
  decay_layer = 1,
  blur_capture = 0.95,
  im_return = FALSE,
  include_virtualroots = FALSE
)
```

## Arguments

- img:

  SpatRaster with three layers for production, decay, and stagnation

- product_layer:

  Integer indicating the production layer index (1-3)

- decay_layer:

  Integer indicating the decay & tape layer index (1-3)

- blur_capture:

  Threshold for pixel inclusion (0-1). Default: 0.95

- im_return:

  Logical: return images instead of values? Default: FALSE

- include_virtualroots:

  Logical: consider all roots present at any timepoint? Default: FALSE

## Value

If im_return = FALSE: tibble with pixel sums and ratios If im_return =
TRUE: list of SpatRaster layers for tape, constant, production, and
decay

## Examples

``` r
if (FALSE) { # \dontrun{
data(TurnoverDPC_data)
img = terra::rast(TurnoverDPC_data)
DPCs = turnover_dpc(img = img, im_return = FALSE)
} # }
```
