# Assign a soil depth in centimeters to every pixel

Builds a per-pixel depth map from the imaging geometry. `geometry`
selects the shape being imaged, and the shape decides which of the two
corrections apply: foreshortening of the depth axis by the insertion
angle, and the depth spread across a curved surface.

## Usage

``` r
create_depthmap(
  img,
  mask = NULL,
  sinoid = TRUE,
  tube_thicc = 7,
  tilt = NULL,
  dpi = 300,
  start_soil = 0,
  center_offset = 0.5,
  geometry = NULL,
  angle_from_vertical = 45
)
```

## Arguments

- img:

  Input image (accepts terra SpatRaster, matrix, array, or file path).
  Only its size is used, so any layer of a multi-layer image will do. A
  SpatRaster passes its extent, resolution and CRS on to the depth map.

- mask:

  Raster mask on the same grid as `img`, marking foreign objects such as
  tape with 1 (or TRUE); those cells are `NA` in the depth map. 0 and NA
  cells are kept.

- sinoid:

  Logical. Superseded by `geometry`, and kept so existing code runs
  unchanged: `TRUE` means `geometry = "tube"`, `FALSE` means
  `geometry = "window"`. Ignored when `geometry` is given.

- tube_thicc:

  Numeric; diameter of minirhizotron tube in cm. Used by
  `geometry = "tube"` only

- tilt:

  Numeric. Superseded by `angle_from_vertical`, which states the
  convention in its name. Degrees from *horizontal*, the opposite
  reference, so `tilt = t` is `angle_from_vertical = 90 - t`. Given it,
  the conversion is applied and a warning says so. `45` is the same
  angle either way

- dpi:

  Numeric; image resolution in dots per inch

- start_soil:

  Numeric; soil surface boundary in cm (0 = surface)

- center_offset:

  Numeric; the image row of the top of the tube – the shallowest line of
  the curvature wave – given in one of two forms:

  - **Fraction** - a number in `[0, 1]`: the tube top is on row
    `1 + center_offset * (nrow - 1)`, so `0` puts it on the first row,
    `0.5` halfway down and `1` on the last row.

  - **Absolute row** - a number `> 1`: the row itself, at most `nrow`,
    e.g. from
    [`estimate_rotation_center()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_center.md).

  The wave itself always spans one full tube circumference
  (`tube_thicc * pi` cm at `dpi`); an image that covers only part of the
  circumference shows the matching part of the wave. Used by
  `geometry = "tube"` only

- geometry:

  Character; the shape being imaged. One of `"tube"` (default),
  `"window"` or `"flat"` – see **Geometries**.

- angle_from_vertical:

  Numeric; insertion angle in degrees measured **from vertical**, in
  `[0, 90]`. `0` is an upright tube and `90` a horizontal one, so the
  number is how far the tube leans over – the convention field notes
  usually record ("installed at 45 degrees"). Used by
  `geometry = "tube"` and `"window"`.

  The angle splits a tube's depth between its two axes, and the two
  terms are complements. The depth gained along the tube's length scales
  by `cos(angle_from_vertical)`: an upright tube gains a full centimeter
  of depth per centimeter of tube, a horizontal one gains none. The
  depth spread around its circumference scales by
  `sin(angle_from_vertical)`: an upright tube's cross-section is
  horizontal, so its whole circumference sits at one depth, while a
  horizontal tube's upper side is a full diameter above its lower side.

## Value

terra raster object containing the depth map

## Geometries

- `"tube"`:

  A cylindrical minirhizotron tube, unrolled by the scanner. Depth
  increases along the image width and is scaled by
  `cos(angle_from_vertical)`; a curvature wave spanning one
  circumference is added across the image height, scaled by
  `sin(angle_from_vertical)`. Needs `tube_thicc`, `angle_from_vertical`
  and `center_offset`.

- `"window"`:

  A planar rhizotron window. Depth increases along the image width and
  is scaled by `cos(angle_from_vertical)`, with no curvature term.
  `tube_thicc` and `center_offset` are unused.

- `"flat"`:

  A surface with no depth axis at all – a flatbed scan of washed roots
  on a tray. Every pixel is at the same depth, `-start_soil`, so a
  profile over this map collapses to one depth bin. `tube_thicc`,
  `angle_from_vertical` and `center_offset` are unused.

## Author

Johannes Cunow, Robert Weigel

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
seg_Oulanka2023_Session01_T067 = terra::rast(seg_Oulanka2023_Session01_T067)
img = seg_Oulanka2023_Session01_T067
mask = seg_Oulanka2023_Session01_T067[[1]] > 0   # layer 1 marks the tape
map = create_depthmap(img, mask, start_soil = 0.1,
  geometry = "tube",
  tube_thicc = 7,
  angle_from_vertical = 45,
  dpi = 300,
  center_offset = 0.1 )

# A tray of washed roots: no depth axis, so every pixel gets one depth
flat = create_depthmap(img, geometry = "flat", dpi = 300)
terra::minmax(flat)
#>     lyr.1
#> min     0
#> max     0
```
