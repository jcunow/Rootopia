# Per-order summary of length and diameter

Per-order summary of length and diameter

## Usage

``` r
summarize_orders(et, order_col = "branch_order")
```

## Arguments

- et:

  An `edges` table (ideally already in real units via
  [`convert_root_units`](https://jcunow.github.io/RootScanR/reference/convert_root_units.md)
  or
  [`branch_order_map`](https://jcunow.github.io/RootScanR/reference/branch_order_map.md)).

- order_col:

  Which order column to group by.

## Value

A data.frame with one row per order:

- `order`:

  The order class value.

- `n_segments`:

  Number of graph segments (edges between nodes).

- `n_tips`:

  Number of root tips/apices (degree-1 endpoints).

- `n_branch_points`:

  Number of branching junctions on roots of this order (= number of
  laterals departing, for ordinary 3-way forks).

- `total_length`:

  Summed root length.

- `mean_segment_length`:

  Mean length per segment.

- `branching_frequency`:

  Branch points per unit length (`n_branch_points / total_length`).

- `mean_diameter`:

  Length-weighted mean diameter.

- `median_diameter`:

  Median of segment median diameters.

Unordered (NA) segments are excluded and counted in
`attr(., "n_unordered")`.
