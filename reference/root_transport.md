# Per-segment transport quantities

Adds, for every segment, how far and through how much resistance water
would travel from the base of its root system to that segment, and how
much root the segment serves. The result is a per-segment vector, so it
can be painted onto the image with
[`render_order_overlay`](https://jcunow.github.io/Rootopia/reference/render_order_overlay.md).

## Usage

``` r
root_transport(et, eta = 1, min_radius = NULL)
```

## Arguments

- et:

  Edge table from
  [`root_graph_pipeline`](https://jcunow.github.io/Rootopia/reference/root_graph_pipeline.md)
  or
  [`branch_order_map`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md),
  needing `from`, `to`, `length` and `mean_diameter`.

- eta:

  Dynamic viscosity. Cancels in any ratio; change it only if you want
  resistance in physical units.

- min_radius:

  Radius floor, applied before the fourth power so a zero-radius segment
  cannot give infinite resistance. `NULL` (default) uses the smallest
  positive radius in the table, so the floor carries whatever unit the
  table is in; a segment with no usable diameter is then treated as the
  thinnest root present. Give a number to set the floor yourself, in the
  table's length unit.

## Value

`et` with the six columns above added. Segments on a cycle get `NA` for
the downstream columns, since they have no subtree, and a warning names
how many there were. A component with no free end at all – a pure ring,
which a root system does not produce – has no proximal end to root at,
so an arbitrary node is used and its `path_*` values are measured from a
meaningless start.

## What the columns mean

- `resistance`:

  Hagen-Poiseuille resistance of the segment alone, \\8 \eta L / (\pi
  r^4)\\. A geometric proxy: it treats the root as an open pipe, where
  real resistance is set by xylem anatomy. The viscosity cancels in
  every comparison, so the default `eta = 1` leaves the values relative.

- `path_length`, `path_resistance`:

  Cumulative from the base of the component to the far end of this
  segment.

- `downstream_tips`, `downstream_length`, `downstream_surface`:

  How many tips, how much root length, and how much lateral surface area
  lie beyond this segment – what it serves. These peak at the base.

Lengths and diameters are used in whatever unit the edge table carries,
so run
[`convert_root_units`](https://jcunow.github.io/Rootopia/reference/convert_root_units.md)
first if you want centimetres.

## See also

[`render_order_overlay`](https://jcunow.github.io/Rootopia/reference/render_order_overlay.md)
to paint any of these columns.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- branch_order_map(skel, mask, unit = "cm", dpi = 300)
et  <- root_transport(res$edges)

# a transport map: colour every segment by its path resistance
render_order_overlay(attr(et, "segments"), et$path_resistance,
                     attr(et, "dims"), "transport.png",
                     value_label = "path resistance")
} # }
```
