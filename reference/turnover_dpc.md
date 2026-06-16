# Extract Root Decay, New Root Production, and No-Change Roots (only 'RootDetector' images)

Extract Root Decay, New Root Production, and No-Change Roots (only
'RootDetector' images)

## Usage

``` r
turnover_dpc(
  img,
  product.layer = 2,
  decay.layer = 1,
  blur.capture = 0.95,
  im.return = FALSE,
  include.virtualroots = FALSE
)
```

## Arguments

- img:

  SpatRaster with three layers for production, decay, and stagnation

- product.layer:

  Integer indicating the production layer index (1-3)

- decay.layer:

  Integer indicating the decay & tape layer index (1-3)

- blur.capture:

  Threshold for pixel inclusion (0-1). Default: 0.95

- im.return:

  Logical: return images instead of values? Default: FALSE

- include.virtualroots:

  Logical: consider all roots present at any timepoint? Default: FALSE

## Value

If im.return = FALSE: tibble with pixel sums and ratios If im.return =
TRUE: list of SpatRaster layers for tape, constant, production, and
decay

## Examples

``` r
if (FALSE) { # \dontrun{
data(TurnoverDPC_data)
img = terra::rast(TurnoverDPC_data)
DPCs = turnover_dpc(img = img, im.return = FALSE)
} # }
```
