# Synthetic root image with known geometry

Draws a root system whose length, width and topology are prescribed, so
the branching pipeline can be scored against ground truth. Returns both
the filled mask (as a scanner would see it) and the exact one-pixel
centre line, which lets graph errors be told apart from thinning errors.

## Usage

``` r
root_phantom(
  design = c("comb", "herringbone", "hierarchical", "cross", "fork"),
  size = 400,
  dpi = 300,
  as = c("matrix", "spatraster")
)
```

## Arguments

- design:

  One of `"comb"` (main axis + 5 perpendicular laterals),
  `"herringbone"` (45-degree laterals), `"hierarchical"` (three
  generations of T-attachments), `"cross"` (two roots overlapping in
  an X) or `"fork"` (symmetric dichotomous Y).

- size:

  Image side length in pixels (square image).

- dpi:

  Nominal scan resolution recorded in `$truth$dpi`; used to express the
  ground truth in cm as well as pixels.

- as:

  `"matrix"` (default) or `"spatraster"` for the returned images.

## Value

A list with `$mask` (filled binary root image), `$skeleton` (exact 1-px
centre line), `$truth` (ground-truth list, including `$by_order`) and
`$strokes` (the stroke definitions).

## Details

Every design is built from strokes – polylines stamped with a disk of
radius `r`, giving a root `2r + 1` px wide – and each lateral starts
exactly on its parent's centre line. The ground truth is therefore
analytic:

- `total_length`:

  Sum of the stroke centre-line lengths (px).

- `n_tips`, `n_branch_points`:

  Counted from the stroke attachment graph, so they hold for
  T-attachments and forks alike.

- `n_roots`, `max_branch_order`, `by_order`:

  Per-root truth, defined only where the continuation rule is
  unambiguous (see `continuation_defined`).

Two conventions are worth stating, because the truth is written to match
the package and not the other way round:

- **Diameter.** Rootopia reports `2 * EDT`, the inscribed diameter
  measured to the nearest *background* pixel, exactly as
  [`root_diameter`](https://jcunow.github.io/Rootopia/reference/root_diameter.md)
  does. For a stroke `W` px wide that is `W + 1`, so
  `by_order$mean_diameter` carries the `+1`.

- **Length.** The designs run along the axes or at 45 degrees, so the
  drawn skeleton's chain-code length equals its Euclidean length. At
  other angles the pipeline would *correctly* report the well-known
  chain-code overestimate, which is a property of the measure rather
  than a defect.

`"fork"` is a symmetric dichotomous Y. Which arm continues the parent is
genuinely undefined there, so `truth$continuation_defined` is `FALSE`
and
[`validate_branching`](https://jcunow.github.io/Rootopia/reference/validate_branching.md)
scores only length and topology for it.

## See also

[`validate_branching`](https://jcunow.github.io/Rootopia/reference/validate_branching.md),
[`branch_order_map`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md)

## Examples

``` r
ph <- root_phantom("comb", size = 200)
ph$truth$n_tips
#> [1] 7
ph$truth$total_length
#> [1] 530
```
