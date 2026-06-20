# Compose two already-loaded (H, W, C) arrays into a mosaic

The compositing core shared by
[`stitch_image_pair`](https://jcunow.github.io/Rootopia/reference/stitch_image_pair.md)
and
[`stitch_image_sequence`](https://jcunow.github.io/Rootopia/reference/stitch_image_sequence.md).

## Usage

``` r
stitch_compose_pair(
  a1,
  a2,
  method,
  edge_width,
  vertical_region,
  vertical_offset,
  direction = "horizontal",
  preprocess = "none",
  blend = "linear",
  blend_width = NULL
)
```

## Arguments

- a1, a2:

  Already-loaded `(H, W, C)` numeric arrays (0-255).

- method, edge_width, vertical_region, vertical_offset, direction,
  preprocess, blend, blend_width:

  As described for
  [`stitch_image_pair`](https://jcunow.github.io/Rootopia/reference/stitch_image_pair.md).

## Value

A list `list(mosaic, dx, dy, peak, overlap)`.
