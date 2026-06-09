# Create a morphological structuring element

Returns a \`cimg\` kernel of a given shape for use with
\`imager::dilate()\` / \`imager::erode()\`.

## Usage

``` r
create_kernel(shape = "disk", size = 3)
```

## Arguments

- shape:

  One of \`"disk"\` (default), \`"square"\`, or \`"diamond"\`.

- size:

  Kernel size in pixels (odd integer). Even values are silently
  incremented by 1.

## Value

A \`cimg\` kernel.
