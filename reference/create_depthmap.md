# Create A Phase-Shifted, Tilt-Amplitude Sine Depth Map

This function generates a depth map for minirhizotron images, accounting
for tube geometry and insertion angle.

## Usage

``` r
create_depthmap(
  img,
  mask = NULL,
  sinoid = TRUE,
  tube_thicc = 7,
  tilt = 45,
  dpi = 300,
  start_soil = 0,
  center_offset = 0.5
)
```

## Arguments

- img:

  Input image (accepts terra SpatRaster, matrix, array, or file path).
  For multi-band images, specify band_index parameter

- mask:

  Raster mask indicating foreign objects (1 = mask, 0 or NA = keep)

- sinoid:

  Logical; if TRUE, accounts for tube curvature in depth calculation

- tube_thicc:

  Numeric; diameter of minirhizotron tube in cm

- tilt:

  Numeric; insertion angle in degrees measured from horizontal, in (0,
  90\]. The depth axis is scaled by `sin(tilt)`, so 90 is the head-on
  case (vertical tube, or a flatbed scan) with no foreshortening.
  Minirhizotron tubes are commonly installed 30-45 degrees from
  vertical, which is 45-60 in this convention

- dpi:

  Numeric; image resolution in dots per inch

- start_soil:

  Numeric; soil surface boundary in cm (0 = surface)

- center_offset:

  Numeric in `[0, 1]`; where the top of the tube – the shallowest line
  of the curvature wave – sits across the image height, as a fraction of
  one full tube circumference. `0` puts it on the first row, `0.5`
  halfway down, `1` on the last row. `0` and `1` are the same point on
  the tube, the image seam

## Value

terra raster object containing the depth map

## Author

Johannes Cunow, Robert Weigel

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
seg_Oulanka2023_Session01_T067 = terra::rast(seg_Oulanka2023_Session01_T067)
img = seg_Oulanka2023_Session01_T067
mask = seg_Oulanka2023_Session01_T067[[1]] - seg_Oulanka2023_Session01_T067[[2]]
mask[mask == 255] <- NA
map = create_depthmap(img,mask,start_soil = 0.1,
  sinoid = TRUE,
  tube_thicc = 7,
  tilt = 45,
  dpi = 300,
  center_offset = 0.1 )
```
