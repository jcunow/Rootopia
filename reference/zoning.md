# Zone Image Data by Depth and/or Rotation Slices with Optional Spatial Cropping

This function extracts and zones image data based on depth bins and/or
rotation slices. It supports three modes of operation:

- `"depth"`:

  Zones image pixels according to binned depth values using a provided
  depth map and selected depth indices.

- `"rotation"`:

  Zones the image by slicing it according to specified rotation slices
  along the x-axis (rows).

- `"both"`:

  Applies depth zoning first, then applies rotation slicing on the
  depth-zoned image.

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
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object representing the input image. Can be multi-layer.

- depth_map:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with binned depth values corresponding spatially to `img`.
  Required if `mode` is `"depth"` or `"both"`. Must have compatible
  spatial properties with `img`.

- depth:

  Numeric vector or range specifying depth values to select from
  `depth_map`. Supports single values, sequences (e.g. `3:6`), or
  arbitrary numeric vectors. Required if `mode` is `"depth"` or
  `"both"`.

- select_layer:

  Integer scalar. Selects which layer to extract from `img` if
  multi-layer. Use `NULL` to keep all layers.

- crop_extent:

  Numeric vector of length 4 in format `c(xmin, xmax, ymin, ymax)` for
  spatial cropping applied before zoning operations. Use `NULL` for no
  cropping.

- mode:

  Character string specifying the zoning mode. One of `"depth"`,
  `"rotation"`, or `"both"`. Default is `"rotation"`.

- rotation_slices:

  Numeric vector of length 2 specifying start and end slice indices for
  rotation zoning (e.g., `c(2, 4)`). Values must be between 1 and
  `rotation_total_slices`. Required if `mode` is `"rotation"` or
  `"both"`.

- rotation_total_slices:

  Integer scalar specifying the total number of conceptual slices along
  the x-axis for rotation zoning. Required if `mode` is `"rotation"` or
  `"both"`.

## Value

A
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object containing the zoned and optionally cropped image data. Pixels
not matching the zoning criteria are set to `NA`.

## Details

The function supports flexible cropping of the input image using a
specified spatial extent, which is applied before any zoning operations.
When depth mode is used, the function warns if cropping would remove
pixels containing desired depth values.

The function processes operations in the following order:

1.  **Spatial cropping**: If `crop_extent` is provided, both `img` and
    `depth_map` are cropped first.

2.  **Depth zoning**: If mode includes "depth", pixels are masked based
    on depth values in `depth_map`.

3.  **Rotation zoning**: If mode includes "rotation", the image is
    sliced along the x-axis (rows) according to `rotation_slices`.

4.  **Layer selection**: If `select_layer` is specified, only that layer
    is retained.

**Depth Processing**:

- If `depth` is a single value, selects the closest matching depth in
  `depth_map`

- If `depth` is a sequence (consecutive values), selects all depths
  within the range

- If `depth` is an arbitrary vector, selects the closest unique depths
  for each value

**Rotation Processing**:

- Divides the image rows into `rotation_total_slices` conceptual slices

- Extracts rows corresponding to slices `rotation_slices[1]` through
  `rotation_slices[2]`

- Useful for analyzing angular segments in circular or rotational
  patterns

**Performance Optimization**:

- Uses vectorized operations with terra::ifel() for fast depth filtering

- Direct logical operations for creating depth masks

- Eliminates intermediate raster creation for better performance

**Warnings**: The function issues warnings when:

- Cropping would remove pixels containing desired depth values

- The resulting image contains only NA values

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example data
data(seg_Oulanka2023_Session01_T067)
img = terra::rast(seg_Oulanka2023_Session01_T067)
depth_map = terra::t(create_depthmap(img))

# Depth zoning example - select depths 3 through 6
zone_img <- zone_image(
  img = img,
  mode = "depth",
  depth_map = depth_map,
  depth = 3:6
)

# Rotation zoning example - extract slices 2-4 out of 8 total
zone_img <- zone_image(
  img = img,
  mode = "rotation",
  rotation_slices = c(2, 4),
  rotation_total_slices = 8
)

# Combined zoning with cropping
zone_img <- zone_image(
  img = img,
  mode = "both",
  depth_map = depth_map,
  depth = c(4.5, 5.2, 6.0),
  rotation_slices = c(1, 3),
  rotation_total_slices = 5,
  crop_extent = c(100, 500, 200, 400)
)

# Select only the first layer and crop spatially
zone_img <- zone_image(
  img = img,
  mode = "depth",
  depth_map = depth_map,
  depth = 5,
  select_layer = 1,
  crop_extent = c(0, 1000, 0, 1000)
)

# Mixed depth selection: range plus individual values
zone_img <- zone_image(
  img = img,
  mode = "depth",
  depth_map = depth_map,
  depth = c(3:8, 13, 15.5)  # Range 3-8 plus closest to 13 and 15.5
)
} # }

```
