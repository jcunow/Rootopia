# Smooth object edges with morphological closing

Applies a morphological closing (dilation then erosion) to smooth the
edges of binary objects. Intended as a post-cleaning step before
skeletonisation; do \*\*not\*\* apply after skeletonisation.

## Usage

``` r
smooth_root_edges(img, kernel_shape = "disk", kernel_size = 3, iterations = 1)
```

## Arguments

- img:

  A \`cimg\` binary image or any format accepted by
  \[load_flexible_image()\].

- kernel_shape:

  One of \`"disk"\` (default), \`"square"\`, \`"diamond"\`.

- kernel_size:

  Structuring element size (odd integer). At 300 DPI, \`kernel_size =
  3\` is a good starting point.

- iterations:

  Number of closing iterations.

## Value

A \`cimg\` binary image.
