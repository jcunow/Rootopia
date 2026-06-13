# Aggregate root architecture by order, or split focal-vs-rest

Flexible query over a `branch_order_map` result (or its `$edges` table):
length, diameter, and branching architecture, either per order class or
split into a focal group versus the rest. Diameter is length-weighted so
it reflects how much root length sits at each width. Lengths/diameters
are in the unit of `x`; counts are integers; `branching_frequency` is
per unit length.

## Usage

``` r
order_metrics(x, order_col = NULL, focal = NULL)
```

## Arguments

- x:

  A `"branchOrderMap"` object or an `edges` data.frame.

- order_col:

  Order column to aggregate by; defaults to the object's `$order`, else
  `"branch_order"`.

- focal:

  Controls the split:

  - `NULL` (default): one row per order class.

  - `"thinnest"` / `"thickest"` (aliases `"finest"` / `"coarsest"`): the
    order class with the smallest / largest length-weighted diameter,
    versus all others. Selected by *diameter*, so it is independent of
    the order numbering; the `orders` column reports which order
    number(s) fell in each group.

  - a numeric order value or vector (e.g. `1` or `c(4, 5)`): those
    orders as the focal group versus the rest.

## Value

A data.frame with one row per group:

- `group`:

  Order value (no split) or `"focal"`/`"rest"`.

- `orders`:

  Comma-separated order numbers contained in the group (so a
  thinnest/thickest split shows which orders it picked).

- `n_segments`:

  Graph segments in the group.

- `n_tips`:

  Root tips/apices (degree-1 endpoints).

- `n_branch_points`:

  Branching junctions on roots of the group.

- `total_length`:

  Summed root length.

- `length_fraction`:

  Share of total ordered length.

- `mean_segment_length`:

  Mean length per segment.

- `branching_frequency`:

  Branch points per unit length.

- `mean_diameter`:

  Length-weighted mean diameter.

- `median_diameter`:

  Median of segment median diameters.

`attr(., "n_unordered")` counts NA-order segments; for a split,
`attr(., "focal_orders")` lists the focal order number(s).

## Examples

``` r
if (FALSE) { # \dontrun{
res <- branch_order_map(skel, mask, order = "branch_order", unit = "cm", dpi = 300)
order_metrics(res)                       # per-order table
order_metrics(res, focal = "thinnest")   # thinnest order class vs all others
order_metrics(res, focal = c(1))         # order 1 vs the rest
} # }
```
