# Unified Root Turnover Analysis

Performs root turnover analysis for either a single multi-layer image or
two separate images

## Usage

``` r
root_turnover(
  img1,
  img2 = NULL,
  method = "kimura",
  unit = "cm",
  dpi = 300,
  select.layer = NULL,
  product.layer = 2,
  decay.layer = 1,
  blur.capture = 0.95,
  im.return = FALSE,
  include.virtualroots = FALSE
)
```

## Arguments

- img1:

  Primary SpatRaster input (either multi-layer or first timepoint image)

- img2:

  Optional second timepoint image (if img1 is single timepoint)

- method:

  Analysis method: "kimura", "rootpx", or "dpc" (root decomposition)

- unit:

  Unit of root length measurement (only for method = "kimura"). Default:
  "cm"

- dpi:

  Image resolution (only for method = "kimura"). Default: 300

- select.layer:

  Integer or NULL. When two images are provided with multiple layers,
  specifies which layer to use. When img1 is multi-layer, ignored for
  DPC method.

- product.layer:

  Integer indicating the production layer index for DPC method (1-3)

- decay.layer:

  Integer indicating the decay & tape layer index for DPC method (1-3)

- blur.capture:

  Threshold for pixel inclusion in DPC method (0-1). Default: 0.95

- im.return:

  Logical: return images instead of values for DPC method? Default:
  FALSE

- include.virtualroots:

  Logical: consider all roots present at any timepoint in DPC method?
  Default: FALSE

## Value

Depends on method and parameters: - For temporal comparison: data.frame
with root production and turnover - For DPC method: tibble with pixel
sums and ratios or list of SpatRaster layers

## See also

[`turnover_tc`](https://jcunow.github.io/RootScanR/reference/turnover_tc.md),
[`turnover_dpc`](https://jcunow.github.io/RootScanR/reference/turnover_dpc.md)
