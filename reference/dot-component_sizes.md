# Connected-component sizes (px) of a logical matrix

Labels the connected `TRUE` regions of `x` with the same connectivity
clean_image() uses – 8 for root objects, 4 for holes – and returns one
size per region. With `exclude_border = TRUE`, regions touching the
image edge are dropped (used for enclosed holes).

## Usage

``` r
.component_sizes(x, directions, exclude_border)
```
