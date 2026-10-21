# Convert edge-table lengths and diameters to real units

Pixels convert as `inch = px / dpi` and `cm = px * 2.54 / dpi`. `length`
is set to the chosen method; `length_poly` and `length_kimura` are both
retained in the same unit.

## Usage

``` r
convert_root_units(
  et,
  unit = c("cm", "inch", "px"),
  dpi = 300,
  length_method = c("kimura", "polyline")
)
```

## Arguments

- et:

  Edge table in pixels (from
  [`root_graph_pipeline`](https://jcunow.github.io/Rootopia/reference/root_graph_pipeline.md)).

- unit:

  One of `"cm"`, `"inch"`, `"px"`.

- dpi:

  Scan resolution (dots per inch).

- length_method:

  `"kimura"` (default) or `"polyline"`. `"polyline"` is the sqrt(2)
  chain code: it walks the skeleton pixel to pixel, counting 1 for an
  orthogonal step and sqrt(2) for a diagonal one. A digital line is a
  staircase, so that walk is longer than the line it represents and the
  chain code reports the staircase. Kimura's estimator, \\\sqrt{n_d^2 +
  (n_d + n_o/2)^2} + n_o/2\\ over the diagonal and orthogonal step
  counts, corrects for that and is the better estimate of root length.

  Both reduce to the same value when a chain is all orthogonal or all
  diagonal, so a perfectly horizontal, vertical or 45-degree segment
  reports the same length either way. They differ only where a chain
  mixes the two, which is where the staircase is, and there the chain
  code reads high. Kimura's correction is derived for straight digital
  lines, so it is most exact on segments that are nearly straight; both
  measures follow whatever path the skeleton takes, so neither can
  recover length the skeleton does not have.

  Use `"polyline"` to compare against software that reports the chain
  code.

## Value

`et` with length/diameter columns in `unit`; records `unit`, `dpi`,
`length_method` as attributes.

## Examples

``` r
if (FALSE) { # \dontrun{
et <- root_graph_pipeline(skel, mask)
et_cm <- convert_root_units(et, unit = "cm", dpi = 300)
} # }
```
