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
  length_method = c("polyline", "kimura")
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

  `"polyline"` (sqrt(2) chain code) or `"kimura"`.

## Value

`et` with length/diameter columns in `unit`; records `unit`, `dpi`,
`length_method` as attributes.
