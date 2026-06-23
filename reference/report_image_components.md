# Summarise the component-size distribution of a binary image

Reports the size distribution of the two component types that
\[clean_image()\] acts on, to help you choose \`max_hole_size\` and
\`max_artifact_size\`:

- \*\*holes\*\* — \`background (0)\` regions fully enclosed by \`root
  (1)\` (segmentation gaps inside roots). Regions touching the image
  border are background, not holes, and are excluded.

- \*\*root components\*\* — connected \`root (1)\` regions. The whole
  root system is usually one large component; genuine artifacts are the
  small components at the low end of this distribution.

## Usage

``` r
report_image_components(img, plot = TRUE, breaks = 30, select_layer = NULL)
```

## Arguments

- img:

  A binary image (\`SpatRaster\`, \`cimg\`, matrix, or file path); any
  format accepted by \[load_flexible_image()\].

- plot:

  Logical. Draw size-distribution histograms (log10 x-axis, with mean
  and median marked). Default \`TRUE\`.

- breaks:

  Number of histogram breaks. Default \`30\`.

- select_layer:

  Integer or \`NULL\`. Layer to use for multi-layer inputs.

## Value

Invisibly, a list with numeric vectors \`holes\` and \`objects\` giving
the size (in pixels) of every enclosed hole and every root component.

## No size threshold here

This function applies \*\*no\*\* size cutoff — it characterises
\*every\* component. There is therefore nothing special separating a
"big root" from an "artifact": both are root components, distinguished
only by size. Use the printed summary and the histograms to pick the
\`max\_\*\_size\` values you then pass to \[clean_image()\]. (Values are
in pixels; scale with your scanner DPI.)

## Examples

``` r
if (FALSE) { # \dontrun{
data(seg_Oulanka2023_Session01_T067)
sizes <- report_image_components(seg_Oulanka2023_Session01_T067)
quantile(sizes$objects)   # inspect the small end to set max_artifact_size
} # }
```
