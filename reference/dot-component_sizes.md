# Connected-component sizes (px) of a binary cimg

Labels the connected \`1\`-regions of \`bin_cimg\` and returns one size
per component. With \`exclude_border = TRUE\`, components touching the
image edge are dropped (used for enclosed holes).

## Usage

``` r
.component_sizes(bin_cimg, exclude_border)
```
