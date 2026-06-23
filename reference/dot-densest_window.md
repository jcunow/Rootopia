# Locate the densest frac-sized window in a raster

Returns the map-coordinate centre \`c(cx, cy)\` of the \`frac\`-sized
window containing the most non-zero / non-\`NA\` cells. Uses a coarse
block scan so it stays cheap on large rasters.

## Usage

``` r
.densest_window(x, frac, layer = NULL)
```
