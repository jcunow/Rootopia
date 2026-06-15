# Slice a scan along the rotation (circumferential) axis

Splits a minirhizotron strip into `n` equal, contiguous, non-overlapping
slices along the rotation axis (image rows; depth runs along the
columns). Slice 1 is the top row band. Use to measure a trait per
circumferential position, e.g. as input to
[`rhythmicity()`](https://jcunow.github.io/RootScanR/reference/rhythmicity.md).

## Usage

``` r
slice_rotation(img, n)
```

## Arguments

- img:

  A terra::SpatRaster.

- n:

  Number of slices around the circumference.

## Value

A list of `n` SpatRasters, in rotation order.

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
seg <- terra::rast(seg_Oulanka2023_Session01_T067)
slices <- slice_rotation(seg, 16)
```
