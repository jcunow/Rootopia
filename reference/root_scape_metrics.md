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

  Segmented raster (values = 0, 1). Consider whether a skeletonized
  raster is more appropriate for your use case.

- indexD:

  Depth index for the output column "depth". Useful when called inside a
  loop over depth bins.

- select.layer:

  Integer. Specifies which layer to use if the input is a multi-band
  image. Default is \`NULL\` (single-layer expected).

- metrics:

  Character vector of metrics to calculate; must be valid names from
  \`landscapemetrics::list_lsm()\`.

## Value

A data frame of metric values with columns: metric, value, object,
depth.

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
img <- terra::rast(seg_Oulanka2023_Session01_T067)
RootScapeObject <- root_scape_metrics(img, indexD = 80, select.layer = 2,
                                       metrics = c("lsm_c_ca"))
#> Warning: Please use 'check_landscape()' to ensure the input data is valid.
```
