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
  preprocess = "none",
  min_overlap = 0,
  max_overlap = Inf
)
```

## Arguments

- img1, img2:

  Image inputs accepted by
  [`load_flexible_image`](https://jcunow.github.io/Rootopia/reference/load_flexible_image.md).

- edge_width:

  Width in pixels of the band searched for the overlap, not the overlap
  itself (clamped to the narrower frame). It must be at least the true
  overlap; beyond that it is only an upper bound, and a band many times
  the overlap still locates it.

- vertical_region:

  Height in pixels of the vertical band used for alignment.

- vertical_offset:

  Starting row (from the top) of the vertical band.

- preprocess:

  Preprocessing of the edge bands: one of `"none"`, `"center"`,
  `"norm"`, `"center_norm"`, `"hann"`, `"grad"` or `"grad_norm"`.

- min_overlap, max_overlap:

  Bounds in pixels on the overlap to look for. The peak is read only
  from the columns of the correlation surface that correspond to an
  overlap of `min_overlap` to `max_overlap` pixels, so an overlap
  outside the bracket is never returned. The defaults `0` and `Inf`
  leave both ends open and search the whole band. `max_overlap` is
  clamped to the band width, since the band cannot express a wider
  overlap; `min_overlap` wider than the band, or a `max_overlap` below
  one pixel, is an error rather than a silently empty search.

## Value

Named numeric vector `c(dx, dy, peak)`: `dx` horizontal placement shift
(`overlap = edge_width - dx`), `dy` vertical shift (pixels) and `peak`
the normalized correlation peak height.

## Details

Only a vertical band of the edges is used (rows `vertical_offset` ..
`vertical_offset + vertical_region`, clamped to the image) and only the
outermost `edge_width` columns, optionally `preprocess`ed. The
cross-power spectrum is \\F_1 \bar{F_2} / (\|F_1\|\|F_2\| + \epsilon)\\
and the shift is read from the correlation peak. The returned `dx` is
the placement shift: `overlap = edge_width - dx`, and it is the peak
column as found, because the lag between the two bands is
`edge_width - overlap` and so already lies in `[0, edge_width)`. `dy` is
the vertical shift to apply to `img2`, and does wrap past the half
height, since a vertical drift can go either way.

## See also

[`stitch_image_pair`](https://jcunow.github.io/Rootopia/reference/stitch_image_pair.md),
[`estimate_rotation_shift`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_shift.md)
