# Estimate rotational/depth shift between two root scans

Estimate rotational/depth shift between two root scans

## Usage

``` r
estimate_rotation_shift(
  img1,
  img2,
  cor.type = "phase",
  fixed.depth.pixel = NULL,
  fixed.width = NULL,
  select.layer = NULL,
  window = TRUE,
  overlay = FALSE,
  overlay.layer = 2
)
```

## Arguments

- img1, img2:

  Image inputs (paths, arrays, or rasters).

- cor.type:

  "phase" (phase correlation) or "ccf" (normalized cross-corr).

- fixed.depth.pixel:

  Depth band along COLUMNS. Length-2 = range start:end; longer =
  explicit column indices; NULL = use full width.

- fixed.width:

  Optional: restrict the ROTATION axis (rows), centered.

- select.layer:

  Layer to use for multi-band inputs.

- window:

  Demean + Hann-window before FFT to suppress edge artefacts.

- overlay:

  If TRUE, also draw a before/after magenta-green overlay.

- overlay.layer:

  Layer to display in the overlay (root mask, default 2).

## Value

Named numeric vector: depth (column lag), rotation (row lag), peak.

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
data(seg_Oulanka2023_Session03_T067)
img1 <- terra::rast(seg_Oulanka2023_Session01_T067)
img2 <- terra::rast(seg_Oulanka2023_Session03_T067)
estimate_rotation_shift(img1, img2, cor.type = "phase", select.layer = 2)
#> Warning: Image size mismatch detected; cropping to common extent
#>        depth     rotation         peak 
#> -18.00000000  -9.00000000   0.03476661 
```
