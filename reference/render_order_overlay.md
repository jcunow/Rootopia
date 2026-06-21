# Write an order-colored validation overlay (PNG)

Paints each segment by its order value over the mask; unordered (cyclic)
segments are drawn in red so failures are visible.

## Usage

``` r
render_order_overlay(segs, tip_order, dims, file, mask = NULL, max_side = 2000)
```

## Arguments

- segs:

  Segment list (kept via `keep_segments = TRUE`).

- tip_order:

  Numeric vector of the order value per segment to color by.

- dims:

  Integer length-2 vector, the (cropped) raster dimensions.

- file:

  Output PNG path.

- mask:

  Optional 0/1 mask drawn as a gray background.

- max_side:

  Long-side cap (px) for the output image.

## Value

Invisibly, `file`.
