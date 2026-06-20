# Edge-based phase-correlation alignment of two scans

Estimates the relative pixel shift between the right edge of `img1` and
the left edge of `img2` using FFT phase correlation - the alignment step
behind
[`stitch_image_pair`](https://jcunow.github.io/Rootopia/reference/stitch_image_pair.md).

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
  [`load_flexible_image`](https://jcunow.github.io/Rootopia/reference/load_flexible_image.md).

- edge_width:

  Width in pixels of the edge band used for alignment (clamped to the
  image width).

- vertical_region:

  Height in pixels of the vertical band used for alignment.

- vertical_offset:

  Starting row (from the top) of the vertical band.

- preprocess:

  Preprocessing of the edge bands: one of `"none"`, `"center"`,
  `"norm"`, `"center_norm"`, `"hann"`, `"grad"` or `"grad_norm"`.

## Value

Named numeric vector `c(dx, dy, peak)`: `dx` horizontal placement shift
(`overlap = edge_width - dx`), `dy` vertical shift (pixels) and `peak`
the normalised correlation peak height.

## Details

Only a vertical band of the edges is used (rows `vertical_offset` ..
`vertical_offset + vertical_region`, clamped to the image) and only the
outermost `edge_width` columns, optionally `preprocess`ed. The
cross-power spectrum is \\F_1 \bar{F_2} / (\|F_1\|\|F_2\| + \epsilon)\\
and the shift is read from the correlation peak (with wrap-around past
the half dimension). The returned `dx` is the placement shift:
`overlap = edge_width - dx`. `dy` is the vertical shift to apply to
`img2`.

## See also

[`stitch_image_pair`](https://jcunow.github.io/Rootopia/reference/stitch_image_pair.md),
[`estimate_rotation_shift`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_shift.md)
