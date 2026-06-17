# Rec. 601 luma projection of an (H, W, C) array

OpenCV's COLOR_BGR2GRAY uses the BT.601 weights (0.299, 0.587, 0.114);
the reference is ported with the same weights so the correlation peak
matches.

## Usage

``` r
stitch_luma(a)
```

## Arguments

- a:

  A numeric `(H, W, C)` array.

## Value

A numeric `(H, W)` luma matrix.
