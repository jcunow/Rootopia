# Topology metrics for a spatial network

Summarizes an edge table as a graph: node and edge counts, components,
tips, branch points, and the connectivity coefficients used in network
ecology. Works on any network whose segments are listed as node pairs,
not only root systems.

## Usage

``` r
network_metrics(
  edges,
  topology = c("auto", "projected_tree", "tree", "reticulate")
)
```

## Arguments

- edges:

  An edge table with character columns `from` and `to`, as returned by
  [`root_graph_pipeline`](https://jcunow.github.io/Rootopia/reference/root_graph_pipeline.md),
  or the object returned by
  [`branch_order_map`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md),
  whose `$edges` element is used. `length` is used for `total_length`
  and for `artefact_rate` when present.

- topology:

  One of `"auto"`, `"projected_tree"`, `"tree"`, `"reticulate"`.
  `"auto"` calls a network with no cycles a tree, and anything else a
  projected tree, because an edge table alone cannot distinguish a real
  loop from a crossing. Declare the topology whenever you know it.

## Value

A list with `topology`, `traits`, `qc` and `note`. An empty or zero-row
input returns the same structure with empty `traits` rather than an
error.

## Why the topology has to be declared

The same quantity means three different things depending on whether
loops in the data are real.

- `"reticulate"`:

  Loops are biology: fungal mycelium with anastomoses, leaf venation,
  the endoplasmic reticulum, vasculature. The cyclomatic number and the
  alpha coefficient measure redundancy, so they are returned as traits.

- `"tree"`:

  Loops are impossible: a root system, a bronchial tree, a river
  network. Then \\E = V - G\\ holds by definition, so the cyclomatic
  number is identically zero and the alpha coefficient carries no
  information. They are returned under `$qc`, where a value above zero
  means the input is not actually a forest.

- `"projected_tree"`:

  A tree photographed in two dimensions, which is what a root scan is.
  Every loop is two roots crossing in projection, so the cyclomatic
  number counts unresolved overlaps and belongs under `$qc`.
  `artefact_rate` normalizes that count by network length so it is
  comparable between images.

Reporting a meshedness of 0.002 as a biological result, when it is a
count of root crossings the segmentation failed to separate, is the
mistake this separation is meant to prevent.

## See also

[`root_graph_pipeline`](https://jcunow.github.io/Rootopia/reference/root_graph_pipeline.md),
[`branch_order_map`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md)

## Examples

``` r
# A three-segment tree: one basal edge splitting into two laterals
et <- data.frame(from = c("A", "B", "B"), to = c("B", "C", "D"),
                 length = c(10, 6, 7), stringsAsFactors = FALSE)
network_metrics(et, topology = "tree")$traits
#> $n_nodes
#> [1] 4
#> 
#> $n_edges
#> [1] 3
#> 
#> $n_components
#> [1] 1
#> 
#> $n_tips
#> [1] 3
#> 
#> $n_branch_points
#> [1] 1
#> 
#> $tip_ratio
#> [1] 0.75
#> 
#> $beta
#> [1] 0.75
#> 
#> $total_length
#> [1] 23
#> 
#> $mean_edge_length
#> [1] 7.666667
#> 

# The same edges plus one closing a loop: a crossing on a root scan
et2 <- rbind(et, data.frame(from = "C", to = "D", length = 4))
network_metrics(et2, topology = "projected_tree")$qc
#> $cyclomatic
#> [1] 1
#> 
#> $alpha_meshedness
#> [1] 0.3333333
#> 
#> $gamma
#> [1] 0.6666667
#> 
#> $n_self_loops
#> [1] 0
#> 
#> $n_duplicate_edges
#> [1] 0
#> 
#> $is_forest
#> [1] FALSE
#> 
#> $artefact_cycles
#> [1] 1
#> 
#> $artefact_rate
#> [1] 0.03703704
#> 
```
