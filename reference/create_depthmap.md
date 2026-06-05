# Create A Phase-Shifted, Tilt-Amplitude Sine Depth Map

This function generates a depth map for minirhizotron images, accounting
for tube geometry and insertion angle.

## Usage

``` r
create_depthmap(
  img,
  mask = NULL,
  sinoid = TRUE,
  tube.thicc = 7,
  tilt = 45,
  dpi = 300,
  start.soil = 0,
  center.offset = 0.5,
  progress = FALSE
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

- tube.thicc:

  Numeric; diameter of minirhizotron tube in cm

- tilt:

  Numeric; minirhizotron tube insertion angle in degrees (typically
  30-45)

- dpi:

  Numeric; image resolution in dots per inch

- start.soil:

  Numeric; soil surface boundary in cm (0 = surface)

- center.offset:

  Numeric; rotational center offset (0 = centered, 1 = edge)

- progress:

  Message; indicates how mny rows have been processed

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
map = create_depthmap(img,mask,start.soil = 0.1,
  sinoid = TRUE,
  tube.thicc = 7,
  tilt = 45,
  dpi = 300,
  center.offset = 0.1 )
```
