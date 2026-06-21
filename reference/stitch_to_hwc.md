# Load any supported image to a (height, width, channel) array in 0-255

Normalizes orientation so that dimension 1 is rows (Y), dimension 2 is
columns (X) and dimension 3 is channels - matching the (H, W, C) layout
used by the Python reference (OpenCV) so the stitching math ports
one-to-one.

## Usage

``` r
stitch_to_hwc(input)
```

## Arguments

- input:

  An image in any format accepted by
  [`load_flexible_image`](https://jcunow.github.io/Rootopia/reference/load_flexible_image.md).

## Value

A numeric `(H, W, C)` array scaled to 0-255.
