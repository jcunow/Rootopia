# Clean a binary root image

Performs three sequential cleaning operations on a binary segmented
image: 1. \*\*Hole filling\*\* — fills black regions enclosed by white
(segmentation gaps inside roots). 2. \*\*Artifact removal\*\* — removes
isolated white specks not connected to the image border (false-positive
root detections). 3. \*\*Edge smoothing\*\* \*(optional, off by
default)\* — applies morphological closing to smooth jagged root edges.

## Usage

``` r
clean_image(
  img,
  max_hole_size = NULL,
  max_artifact_size = NULL,
  edge_smooth = FALSE,
  kernel_shape = "disk",
  kernel_size = 3,
  iterations = 1,
  select.layer = NULL,
  output_format = "spatrast",
  report = FALSE
)
```

## Arguments

- img:

  A \`cimg\` object, \`SpatRaster\`, matrix, or file path.

- max_hole_size:

  Maximum hole size in pixels to fill. If \`NULL\`, all enclosed holes
  are filled. See \*\*Choosing thresholds\*\* above.

- max_artifact_size:

  Maximum artifact size in pixels to remove. If \`NULL\`, all isolated
  white regions are removed.

- edge_smooth:

  Logical. Apply morphological closing after hole/artifact cleaning.
  Default \`FALSE\`.

- kernel_shape:

  Structuring element shape for edge smoothing: \`"disk"\` (default),
  \`"square"\`, or \`"diamond"\`.

- kernel_size:

  Structuring element size (odd integer). Default \`3\`.

- iterations:

  Number of closing iterations for edge smoothing. Default \`1\`.

- select.layer:

  Integer or \`NULL\`. Which layer to use for multi-layer inputs.

- output_format:

  Character. Format of the returned object. One of \`"spatrast"\`
  (default), \`"cimg"\`, or \`"matrix"\`. Using \`"spatrast"\` means the
  result can be passed directly to \`terra::plot()\`,
  \`skeletonize_image()\`, \`root_length()\`, etc. without any further
  conversion.

- report:

  Logical. If \`TRUE\`, also calls \[report_image_components()\] on the
  \*original\* (uncleaned) image before cleaning. When \`output_format =
  "spatrast"\` (default), the cleaned raster is returned directly even
  when \`report = TRUE\`; the report is printed as a side effect.
  Default \`FALSE\`.

## Value

A cleaned image in the format specified by \`output_format\`
(\`SpatRaster\`, \`cimg\`, or matrix).

## Why clean before skeletonisation

\`skeletonize_image()\` uses the Medial Axis Transform, which is driven
by the distance transform. Small holes inside a root inflate the local
distance values and force the medial axis to bifurcate around the hole,
producing spurious branching points. Isolated artifact pixels produce
phantom skeleton segments. Cleaning first yields a much cleaner
skeleton.

## Choosing thresholds

Call \[report_image_components()\] on your image first to see the actual
pixel counts of all holes and artifacts. At 300 DPI, sensible starting
values are \`max_hole_size = 50\` and \`max_artifact_size = 10\`.

## Edge smoothing caution

\`edge_smooth = TRUE\` applies a morphological closing that slightly
dilates then erodes root edges. This can merge closely adjacent roots
and alter root diameter measurements. Only use it when the segmentation
output has very jagged edges; leave it off (\`FALSE\`, the default)
otherwise.

## See also

\[report_image_components()\], \[skeletonize_image()\]

## Examples

``` r
data(seg_Oulanka2023_Session01_T067)
img <- terra::rast(seg_Oulanka2023_Session01_T067)


# Clean: fill small holes, remove tiny artifacts — returns SpatRaster
cleaned <- clean_image(img,
                       max_hole_size     = 50,
                       max_artifact_size = 10,
                       select.layer      = 2)

# If you need a cimg for further imager operations:
cleaned_cimg <- clean_image(img, max_hole_size = 50,
                            output_format = "cimg", select.layer = 2)
```
