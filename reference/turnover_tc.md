# Calculate global root production and root turnover from temporal comparison

Calculate global root production and root turnover from temporal
comparison

## Usage

``` r
turnover_tc(
  im.t1,
  im.t2,
  method = "kimura",
  unit = "cm",
  dpi = 300,
  select.layer = 2
)
```

## Arguments

- im.t1:

  SpatRaster object for timepoint 1

- im.t2:

  SpatRaster object for timepoint 2

- method:

  Analysis method: "kimura" or "rootpx"

- unit:

  Unit of root length measurement (only for method = "kimura"). Default:
  "cm"

- dpi:

  Image resolution (only for method = "kimura"). Default: 300

- select.layer:

  Integer specifying the layer to use in both timesteps if `img` is a
  multi-layer \`SpatRaster\`. Defaults to 2.

## Value

data.frame containing: - standingroot_t1: Standing roots at first
timepoint - standingroot_t2: Standing roots at second timepoint -
production: Root production between timepoints - newroot - newroot

## Examples

``` r
  data(skl_Oulanka2023_Session01_T067)
  data(skl_Oulanka2023_Session03_T067)
  time1 <- terra::rast(skl_Oulanka2023_Session01_T067)
  time2 <- terra::rast(skl_Oulanka2023_Session03_T067)
  turnover.values <- turnover_tc(
    im.t1 = time1,
    im.t2 = time2,
    method = "kimura")
```
