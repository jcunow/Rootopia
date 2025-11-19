# Detect rotation shift between two images

Calculates the rotation shift between two sequential images using either
cross-correlation or phase correlation methods.

## Usage

``` r
estimate_rotation_shift(
  img1,
  img2,
  cor.type = "phase",
  fixed.depth.pixel = c(1000, 4000),
  fixed.width = NULL,
  select.layer = NULL
)
```

## Arguments

- img1:

  Reference image (3-channel RGB)

- img2:

  Subsequent image to compare (3-channel RGB)

- cor.type:

  Correlation type: "ccf" (cross) or "phase" (frequency domain)

- fixed.depth.pixel:

  Depth range to analyze c(start, end)

- fixed.width:

  Width of analysis region in pixels

- select.layer:

  Integer. Specifies which layer to use if the input is a multi-band
  image. Default is \`NULL\`.

## Value

Vector of shifts (x,y) in pixels

## Examples

``` r
img1 = seg_Oulanka2023_Session01_T067
img2 = seg_Oulanka2023_Session03_T067
y.lag = estimate_rotation_shift(img1,img2,"phase")
#> Warning: Warning in estimate_rotation_shift: Images differ in size by -17 x 0 pixels
```
