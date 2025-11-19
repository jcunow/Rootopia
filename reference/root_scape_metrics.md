# RootScapeMetric relies on Landscapemetrics to extract 'Root Scape' Features akin to landscape analysis.

RootScapeMetric relies on Landscapemetrics to extract 'Root Scape'
Features akin to landscape analysis.

## Usage

``` r
root_scape_metrics(
  img,
  indexD = NA,
  select.layer = NULL,
  metrics = c("lsm_c_ca", "lsm_l_ent", "lsm_c_pd", "lsm_c_np", "lsm_c_pland",
    "lsm_c_area_mn", "lsm_c_area_cv", "lsm_c_enn_mn", "lsm_c_enn_cv")
)
```

## Arguments

- img:

  segmented raster (values = 0,1). Consider whether skeletonized raster
  is appropriate.

- indexD:

  please specify depth. Will only affect the output column = "depth".
  Useful when used in a loop.

- select.layer:

  Integer. Specifies which layer to use if the input is a multi-band
  image. Default is \`2\`.

- metrics:

  which ,metrics should be calculated from the available ones in
  'landscapemetrics::calculate_lsm()'.

## Value

a bunch of metric values

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
img = terra::rast(seg_Oulanka2023_Session01_T067)
RootScapeObject  = root_scape_metrics(img,indexD = 80, select.layer = 2,  metrics = c("lsm_c_ca"))
#> Warning: Please use 'check_landscape()' to ensure the input data is valid.
```
