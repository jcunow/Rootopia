# Mask a scan to a depth zone (depth zone masking)

`depth_zoning()` performs *depth zone masking*: it returns the input
image with every pixel *outside* the requested depth bin(s) set to `NA`,
keeping the original grid and extent. This lets you run any per-pixel
trait function (length, diameter, color, landscape metrics) on a single
depth slice without splitting the raster into separate objects.

## Usage

``` r
depth_zoning(img, depth_map, depth, select_layer = NULL, crop_extent = NULL)
```

## Arguments

- img:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html).

- depth_map:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  of binned depth values (see
  [`binning`](https://jcunow.github.io/Rootopia/reference/binning.md)),
  aligned to `img`.

- depth:

  Numeric. Target depth bin(s) to keep.

- select_layer:

  Integer. Optionally return a single layer of the result.

- crop_extent:

  Numeric length-4 `c(xmin, xmax, ymin, ymax)`. Optional crop applied
  before masking.

## Value

A
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
masked (extent unchanged) to the requested depth zone.

## Details

Depth is matched against the binned values in `depth_map`: a single
value selects the closest available bin; a consecutive sequence selects
an inclusive range; discrete values select each closest bin.

To split a tube into circumferential (rotation) slices, use
[`slice_rotation`](https://jcunow.github.io/Rootopia/reference/slice_rotation.md)
instead.

## See also

[`binning`](https://jcunow.github.io/Rootopia/reference/binning.md),
[`slice_rotation`](https://jcunow.github.io/Rootopia/reference/slice_rotation.md)

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
img <- terra::rast(seg_Oulanka2023_Session01_T067)
mask <- img[[1]] - img[[2]]; mask[mask == 255] <- NA
depth_map  <- create_depthmap(img, mask, start_soil = 2.9, dpi = 150 )
depth_bins <- binning(depth_map, nn = 5)
# depth_bins is already aligned with img -- no transpose/flip needed
# Keep only the root layer pixels that fall in the 10 cm depth bin
slice_10cm <- depth_zoning(img[[2]], depth_map = depth_bins, depth = 10)
```
