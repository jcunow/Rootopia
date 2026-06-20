# Mask a scan to a depth and/or rotation zone (zone masking)

`zoning()` performs *zone masking*: it returns the input image with
every pixel *outside* the requested zone set to `NA`, keeping the
original grid and extent. This lets you run any per-pixel trait function
(length, diameter, colour, landscape metrics) on a single depth slice
without splitting the raster into separate objects.

## Usage

``` r
zoning(
  img,
  depth_map = NULL,
  depth = NULL,
  select_layer = NULL,
  crop_extent = NULL,
  mode = c("rotation", "depth", "both"),
  rotation_slices = NULL,
  rotation_total_slices = NULL
)
```

## Arguments

- img:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html).

- depth_map:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  of binned depth values (see
  [`binning`](https://jcunow.github.io/RootScanR/reference/binning.md)),
  aligned to `img`. Required for `mode = "depth"` or `"both"`.

- depth:

  Numeric. Target depth bin(s) to keep. Required for `mode = "depth"` or
  `"both"`.

- select_layer:

  Integer. Optionally return a single layer of the result.

- crop_extent:

  Numeric length-4 `c(xmin, xmax, ymin, ymax)`. Optional crop applied
  before masking.

- mode:

  One of `"rotation"`, `"depth"`, or `"both"`.

- rotation_slices:

  Numeric length-2 `c(from, to)` band of slices to keep (1 = top).
  Required for `mode = "rotation"` or `"both"`.

- rotation_total_slices:

  Numeric. Total number of slices the circumference is divided into.
  Required for `mode = "rotation"` or `"both"`.

## Value

A
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
masked (and/or cropped) to the requested zone.

## Details

The zone is defined by one or both of:

- **depth** (`mode = "depth"`): keep only pixels whose binned depth
  (from `depth_map`) matches `depth`. A single value selects the closest
  available bin; a consecutive sequence selects an inclusive range;
  discrete values select each closest bin.

- **rotation** (`mode = "rotation"`): crop to a contiguous band of
  circumferential slices. For simply splitting a tube into `n` equal
  slices, prefer
  [`slice_rotation`](https://jcunow.github.io/RootScanR/reference/slice_rotation.md);
  use the rotation mode here only when you need to combine a rotation
  band with depth masking in one call (`mode = "both"`).

Note that depth masking sets out-of-zone pixels to `NA` (extent
unchanged), whereas rotation masking *crops* the extent to the selected
band.

## See also

[`binning`](https://jcunow.github.io/RootScanR/reference/binning.md),
[`slice_rotation`](https://jcunow.github.io/RootScanR/reference/slice_rotation.md)

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
img <- terra::rast(seg_Oulanka2023_Session01_T067)
mask <- img[[1]] - img[[2]]; mask[mask == 255] <- NA
depth_map  <- create_depthmap(img, mask, start.soil = 2.9, dpi = 150 )
depth_bins <- binning(depth_map, nn = 5)
depth_bins <- terra::flip(terra::t(depth_bins))
# Keep only the root layer pixels that fall in the 10 cm depth bin
slice_10cm <- zoning(img[[2]], mode = "depth",
                     depth_map = depth_bins, depth = 10)
```
