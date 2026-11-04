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
  select_layer = 2
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

- select_layer:

  Integer specifying the layer to use in both timesteps if `img` is a
  multi-layer `SpatRaster`. Defaults to 2.

## Value

data.frame containing:

- standingroot_t1: Standing roots at first timepoint

- standingroot_t2: Standing roots at second timepoint

- production: Net change in standing roots,
  `standingroot_t2 - standingroot_t1`. Roots lost between the dates
  offset roots gained, so the value is negative when more roots
  disappeared than appeared.

- newroot.per_t1: `production` as a fraction of `standingroot_t1` (0.25
  = 25 percent)

- newroot.per_t2: `production` as a fraction of `standingroot_t2`

The last two names carry a dot where the percent sign would read more
naturally: [`data.frame()`](https://rdrr.io/r/base/data.frame.html)
makes the names syntactic, so index them as `newroot.per_t1`.

## Examples

``` r
if (FALSE) { # \dontrun{
  data(skl_Oulanka2023_Session01_T067)
  data(skl_Oulanka2023_Session03_T067)
  time1 <- terra::rast(skl_Oulanka2023_Session01_T067)
  time2 <- terra::rast(skl_Oulanka2023_Session03_T067)
  turnover.values <- turnover_tc(
    im.t1 = time1,
    im.t2 = time2,
    method = "kimura")
    } # }
```
