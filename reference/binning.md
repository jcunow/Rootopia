# Bin continuous depth values into discrete intervals

Bin continuous depth values into discrete intervals

## Usage

``` r
binning(depthmap, nn, round.option = "rounding")
```

## Arguments

- depthmap:

  SpatRaster/matrix/array - continuous depth values

- nn:

  numeric - bin width

- round.option:

  character - binning method: "rounding", "ceiling", or "floor"

## Value

SpatRaster - binned depth values

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
img = terra::rast(seg_Oulanka2023_Session01_T067)
mask = img[[1]] - img[[2]]
mask[mask == 255] <- NA
img = img
depthmap = create_depthmap(img,mask,start.soil = 2.9 )
binned.map = binning(depthmap,nn = 5)
```
