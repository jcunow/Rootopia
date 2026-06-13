# Prune short or thin terminal segments

Iteratively removes terminal (degree-1) segments below a length or
diameter threshold. Only deletes segments (never rewires), so ordering
remains valid. Skeleton-level pruning before the pipeline is an
alternative, fully modular approach.

## Usage

``` r
prune_terminal_segments(segs, DT, min_length = 0, min_diameter = 0, iter = 1L)
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

## Value

The pruned segment list.
