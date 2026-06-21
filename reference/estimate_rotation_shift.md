# Estimate rotational/depth shift between two root scans

Estimate rotational/depth shift between two root scans

## Usage

``` r
estimate_rotation_shift(
  img1,
  img2,
  cor_type = "phase",
  fixed_depth_pixel = NULL,
  fixed_width = NULL,
  select_layer = NULL,
  window = TRUE,
  overlay = FALSE,
  overlay_layer = 2
)
```

## Arguments

- img1, img2:

  Image inputs (paths, arrays, or rasters).

- cor_type:

  "phase" (phase correlation) or "ccf" (normalized cross-corr).

- fixed_depth_pixel:

  Depth band along COLUMNS. Length-2 = range start:end; longer =
  explicit column indices; NULL = use full width.

- fixed_width:

  Optional: restrict the ROTATION axis (rows), centered.

- select_layer:

  Layer to use for multi-band inputs.

- window:

  Demean + Hann-window before FFT to suppress edge artifacts.

- overlay:

  If TRUE, also draw a before/after magenta-green overlay.

- overlay_layer:

  Layer to display in the overlay (root mask, default 2).

## Value

Named numeric vector: depth (column lag), rotation (row lag), peak.

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
data(seg_Oulanka2023_Session03_T067)
img1 <- terra::rast(seg_Oulanka2023_Session01_T067)
img2 <- terra::rast(seg_Oulanka2023_Session03_T067)
estimate_rotation_shift(img1, img2, cor_type = "phase", select_layer = 2)
#> Warning: Image size mismatch detected; cropping to common extent
#>        depth     rotation         peak 
#> -18.00000000  -9.00000000   0.03476661 
```
