# Prune short or thin terminal segments

Iteratively removes terminal (degree-1) segments below a length or
diameter threshold. Segments are only deleted, never re-routed, so the
ordering stays valid. Skeleton-level pruning before the pipeline is an
alternative, fully modular approach.

## Usage

``` r
prune_terminal_segments(
  segs,
  DT,
  min_length = 0,
  min_diameter = 0,
  iter = 1L,
  splice = TRUE
)
```

## Arguments

- segs:

  Segment list from `trace_segments`.

- DT:

  Distance-transform matrix (for the diameter test).

- min_length:

  Minimum segment length (px) to keep a terminal segment.

- min_diameter:

  Minimum segment diameter (px) to keep a terminal segment.

- iter:

  Number of pruning passes.

- splice:

  Re-join the parent across a junction that a removed lateral has left
  with only two arms, before the next pass.

## Value

The pruned segment list.

## Details

Between passes the junctions left behind are re-checked: once a lateral
is gone its branch point carries only two arms, so the parent's two
pieces are spliced back into one segment before the next pass measures
them. Without that step a multi-pass prune eats healthy axes, because
each parent piece between two former branch points is short enough on
its own to look like a spur. Set `splice = FALSE` to delete laterals
without splicing.
