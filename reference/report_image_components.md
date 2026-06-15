# Report sizes of holes and isolated artifacts in a binary image

Prints a human-readable summary of all internal black holes and isolated
white artifacts in \`img\`, together with their pixel counts. Use this
\*\*before\*\* calling \[clean_image()\] to decide on appropriate values
for \`max_hole_size\` and \`max_artifact_size\`.

## Usage

``` r
report_image_components(img)
```

## Arguments

- img:

  A \`cimg\` binary image (values 0 and 1), or any format accepted by
  \[load_flexible_image()\].

## Value

Invisibly \`NULL\`. Prints to the console.

## Choosing thresholds

At 300 DPI a single root cross-section in a minirhizotron scan is
roughly 5–150 px^2 depending on root diameter. As a starting point:

- \`max_artifact_size = 10\` removes single-pixel noise and tiny
  segmentation specks while preserving fine roots.

- \`max_hole_size = 50\` fills small gaps inside roots without merging
  genuinely separate objects.

Scale these linearly if your scanner DPI differs (e.g. at 150 DPI, halve
both values).

## Examples

``` r
if (FALSE) { # \dontrun{
img <- imager::as.cimg(matrix(0, 50, 50))
img[10:20, 10:20] <- 1   # white square
img[13:15, 13:15] <- 0   # hole inside it
img[40, 40]        <- 1  # isolated artifact
report_image_components(img)
} # }
```
