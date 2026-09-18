# Score the branching pipeline against a phantom with known properties

Runs
[`branch_order_map`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md)
on a synthetic root image whose geometry is known exactly and reports,
metric by metric, what was expected, what was measured and whether the
difference is within tolerance. This is the check to run after touching
the tracing, ordering or length code, and the one to cite when reporting
what the package's numbers mean.

## Usage

``` r
validate_branching(
  design = c("comb", "herringbone", "hierarchical", "cross", "fork"),
  phantom = NULL,
  from = c("skeleton", "mask"),
  size = 400,
  dpi = 300,
  tolerance = NULL,
  overlay_png = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- design:

  Phantom design, passed to
  [`root_phantom`](https://jcunow.github.io/Rootopia/reference/root_phantom.md);
  ignored when `phantom` is supplied. For `"fork"` the per-root metrics
  (`n_roots`, `max_branch_order`, per-order length and diameter) have no
  ground truth and are omitted from the report.

- phantom:

  A phantom from
  [`root_phantom`](https://jcunow.github.io/Rootopia/reference/root_phantom.md)
  (or any list with `$mask`, `$skeleton` and `$truth`).

- from:

  `"skeleton"` to score the graph alone, `"mask"` to score the whole
  pipeline including skeletonisation.

- size, dpi:

  Passed to
  [`root_phantom`](https://jcunow.github.io/Rootopia/reference/root_phantom.md)
  when building a phantom.

- tolerance:

  Named list with `count`, `length` and `diameter` relative tolerances;
  defaults depend on `from`.

- overlay_png:

  Optional path for the order-coloured validation image.

- verbose:

  Print the report.

- ...:

  Passed to
  [`branch_order_map`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md).

## Value

A data.frame with one row per metric (`metric`, `expected`, `observed`,
`abs_error`, `rel_error`, `tolerance`, `pass`). `attr(., "passed")` is
`TRUE` when every row passes; `attr(., "result")` holds the
`branchOrderMap` object.

## Details

`from = "skeleton"` feeds the exact one-pixel centre line, so the score
isolates the graph: tracing, junction contraction, crossing resolution,
ordering and length integration. `from = "mask"` skeletonises the filled
image first and therefore also carries the thinning error – chiefly the
erosion of about one root radius at every tip, which shortens the total
by a few percent. Both are legitimate; they answer different questions,
and the default tolerances differ accordingly.

Counts (`n_tips`, `n_branch_points`, `n_roots`, `max_branch_order`,
`n_unordered`) are required to be exact. Lengths and diameters are
scored as relative error.

## See also

[`root_phantom`](https://jcunow.github.io/Rootopia/reference/root_phantom.md),
[`branch_order_map`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md)

## Examples

``` r
# \donttest{
v <- validate_branching("comb", size = 200, verbose = FALSE)
v[, c("metric", "expected", "observed", "pass")]
#>                  metric   expected   observed pass
#> 1          total_length 530.000000 529.057764 TRUE
#> 2                n_tips   7.000000   7.000000 TRUE
#> 3       n_branch_points   5.000000   5.000000 TRUE
#> 4           n_unordered   0.000000   0.000000 TRUE
#> 5               n_roots   6.000000   6.000000 TRUE
#> 6      max_branch_order   2.000000   2.000000 TRUE
#> 7   order1_total_length 180.000000 180.307764 TRUE
#> 8  order1_mean_diameter   9.968956   9.968804 TRUE
#> 9   order2_total_length 350.000000 348.750000 TRUE
#> 10 order2_mean_diameter   6.142957   6.087856 TRUE
attr(v, "passed")
#> [1] TRUE
# }
```
