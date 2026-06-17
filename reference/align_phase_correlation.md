# Edge-based phase-correlation alignment of two scans

Estimates the relative pixel shift between the right edge of `img1` and
the left edge of `img2` using FFT phase correlation - the alignment step
behind
[`stitch_image_pair`](https://jcunow.github.io/RootScanR/reference/stitch_image_pair.md).

## Usage

``` r
align_phase_correlation(
  img1,
  img2,
  edge_width = 250,
  vertical_region = 1000,
  vertical_offset = 300,
  preprocess = "none"
)
```

## Arguments

- img1, img2:

  Image inputs accepted by
  [`load_flexible_image`](https://jcunow.github.io/RootScanR/reference/load_flexible_image.md).

- edge_width:

  Width in pixels of the edge band used for alignment (clamped to image
  width).

- vertical_region:

  Height in pixels of the vertical band used for alignment. The actual
  sampled region is clipped to image bounds if needed.

- vertical_offset:

  Starting row (from the top) of the vertical band. Used as a reference
  position for selecting the vertical window.

- preprocess:

  Preprocessing applied to edge bands before FFT. One of: `"none"`,
  `"center"`, `"norm"`, `"center_norm"`, `"hann"`, `"grad"`, or
  `"grad_norm"`.

## Value

Named numeric vector `c(dx, dy, peak)`:

- `dx`: horizontal placement shift (`overlap = edge_width - dx`)

- `dy`: vertical shift (pixels, applied to `img2`)

- `peak`: normalized correlation peak height

## Details

A vertical band of the edge region is used for alignment. The band is
defined by a requested starting position (`vertical_offset`) and height
(`vertical_region`). If the requested region exceeds image bounds, it is
automatically truncated to fit within the image. Only the outermost
`edge_width` columns are used.

The cross-power spectrum is \\F_1 \bar{F_2} / (\|F_1\|\|F_2\| +
\epsilon)\\ and the shift is obtained from the correlation peak (with
wrap-around past half the transform dimensions).

The returned `dx` is the placement shift: `overlap = edge_width - dx`.
`dy` is the vertical shift applied to `img2`.

## See also

[`stitch_image_pair`](https://jcunow.github.io/RootScanR/reference/stitch_image_pair.md),
[`estimate_rotation_shift`](https://jcunow.github.io/RootScanR/reference/estimate_rotation_shift.md)
