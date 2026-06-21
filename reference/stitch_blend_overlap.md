# Combine the overlap band of two aligned regions

Combine the overlap band of two aligned regions

## Usage

``` r
stitch_blend_overlap(reg1, reg2, blend = "linear", blend_width = NULL)
```

## Arguments

- reg1, reg2:

  Numeric `(h, w, c)` arrays - the aligned overlap of img1 (left) and
  img2 (right).

- blend:

  One of `"linear"` (alpha ramp 1 -\> 0), `"overlay"` (img2 on top),
  `"overlay_first"` (img1 on top), `"max"` (lighten / union) or `"min"`
  (darken).

- blend_width:

  Optional ramp width (px) for `"linear"`, centered in the overlap;
  `NULL` ramps across the whole overlap.

## Value

The combined `(h, w, c)` array.
