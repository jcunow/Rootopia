# Approximate average root thickness (deprecated)

\`root_thickness()\` is a naive estimator that back-calculates an
average diameter from a total root length and a total root pixel count
(\`area / length\`), assuming every pixel belongs to a single root of
uniform width. It ignores branching, overlapping roots, and the actual
local width distribution.

\[root_diameter()\] computes per-pixel diameters directly from the
distance transform of the segmented mask, and returns \`mean_diameter\`,
\`median_diameter\`, \`quantiles\`, \`root_volume\`, and
\`root_surface_area\` (lateral surface area, assuming cylindrical root
segments). Use those instead of \`root_thickness()\`.

## Usage

``` r
root_thickness(kimuralength, rootpx, dpi = 300)
```

## Arguments

- kimuralength:

  Total root length in cm (e.g. from \[root_length()\]).

- rootpx:

  Total number of root pixels in the image section.

- dpi:

  Image resolution in dots per inch. Default is 300.

## Value

A numeric value in cm representing approximate average root diameter.

## Examples

``` r
root.thicc <- root_thickness(kimuralength = 300, rootpx = 9500, dpi = 300)
#> Warning: root_thickness() is a naive area/length estimator and is deprecated. Use root_diameter() for per-pixel diameter, root_volume, and root_surface_area computed from the distance transform.
```
