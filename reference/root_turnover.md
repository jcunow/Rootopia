# Unified Root Turnover Analysis

Wrapper around the two root-turnover methods. `method` selects which one
runs:

- `"tc"` (Temporal Comparison):

  Compares two timepoint images (`img1`, `img2`) and reports standing
  roots, their net change, and that change as a fraction of each date's
  standing roots. Dispatches to
  [`turnover_tc`](https://jcunow.github.io/Rootopia/reference/turnover_tc.md).
  The `tc_method` argument chooses how root amount is measured:
  `"kimura"` (root length) or `"rootpx"` (root pixel count).

- `"dpc"` (Decay, Production, Constant):

  Decomposes a single multi-layer 'RootDetector' image into decayed,
  newly produced, and unchanged (constant) root fractions. Dispatches to
  [`turnover_dpc`](https://jcunow.github.io/Rootopia/reference/turnover_dpc.md);
  `img2` is not used.

## Usage

``` r
root_turnover(
  img1,
  img2 = NULL,
  method = c("tc", "dpc"),
  tc_method = c("kimura", "rootpx"),
  unit = "cm",
  dpi = 300,
  select_layer = NULL,
  product_layer = 2,
  decay_layer = 1,
  blur_capture = 0.95,
  im_return = FALSE,
  include_virtualroots = FALSE
)
```

## Arguments

- img1:

  Primary SpatRaster input. For `method = "tc"` this is the first
  timepoint image; for `method = "dpc"` this is the multi-layer DPC
  image.

- img2:

  Second timepoint SpatRaster. Required for `method = "tc"`; ignored
  (with a warning) for `method = "dpc"`.

- method:

  Which turnover method to run: `"tc"` (temporal comparison of two
  images) or `"dpc"` (Decay/Production/Constant decomposition of one
  multi-layer image). Default `"tc"`.

- tc_method:

  Measurement sub-method for `method = "tc"`: `"kimura"` (root length)
  or `"rootpx"` (root pixel count). Ignored for `method = "dpc"`.
  Default `"kimura"`.

- unit:

  Unit of root length measurement (only for `tc_method = "kimura"`).
  Default: "cm"

- dpi:

  Image resolution (only for `tc_method = "kimura"`). Default: 300

- select_layer:

  Integer or NULL. For `method = "tc"` with multi-layer images, selects
  which layer to compare. Ignored for `method = "dpc"`.

- product_layer:

  Integer indicating the production layer index for the DPC method (1-3)

- decay_layer:

  Integer indicating the decay & tape layer index for the DPC method
  (1-3)

- blur_capture:

  Threshold for pixel inclusion in the DPC method (0-1). Default: 0.95

- im_return:

  Logical: return images instead of values for the DPC method? Default:
  FALSE

- include_virtualroots:

  Logical: consider all roots present at any timepoint in the DPC
  method? Default: FALSE

## Value

Depends on the method:

- `"tc"`: data.frame with standing roots, their net change
  (`production`), and that change as a fraction of each date's standing
  roots.

- `"dpc"`: data.frame of pixel sums and ratios, or (if
  `im_return = TRUE`) a list of SpatRaster layers.

## See also

[`turnover_tc`](https://jcunow.github.io/Rootopia/reference/turnover_tc.md),
[`turnover_dpc`](https://jcunow.github.io/Rootopia/reference/turnover_dpc.md)

## Examples

``` r
# DPC: single multi-layer 'RootDetector' image
data(TurnoverDPC_data)
img <- terra::rast(TurnoverDPC_data)
root_turnover(img, method = "dpc")
#>      tape constant production  decay newgrowth.ratio decay.ratio constant.ratio
#> 1 1763826   177422     246436 180452          0.5814      0.5042         0.2936

# TC: two timepoint images compared by root length (kimura)
data(skl_Oulanka2023_Session01_T067)
data(skl_Oulanka2023_Session03_T067)
t1 <- terra::rast(skl_Oulanka2023_Session01_T067)
t2 <- terra::rast(skl_Oulanka2023_Session03_T067)
root_turnover(t1, t2, method = "tc", tc_method = "kimura")
#> Warning: img has 3 layers and select_layer is not set; measuring layer 1. Layer 1 is often not the skeleton: set select_layer to the layer that holds it.
#> Diagonal: 914793 | Orthogonal: 916350
#> Warning: img has 3 layers and select_layer is not set; measuring layer 1. Layer 1 is often not the skeleton: set select_layer to the layer that holds it.
#> Diagonal: 861596 | Orthogonal: 863048
#>   standingroot_t1 standingroot_t2 production newroot.per_t1 newroot.per_t2
#> 1        17847.64        16809.66  -1037.987        -0.0582        -0.0617
```
