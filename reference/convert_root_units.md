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
  staircase, so that path is longer than the line it represents –
  measured against straight lines of known length, the chain code runs
  about 5.7% long on average and up to 8.2%, worst near 22.5 degrees and
  zero only for perfectly horizontal or vertical segments. Kimura's
  estimator corrects for the staircase and reduces that to about 1.8% on
  average. Use `"polyline"` to reproduce output from Rootopia 1.1.0 and
  earlier, which defaulted to it.

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
