# Convert a cimg object to a SpatRaster, preserving RGB metadata

Internal helper that replaces the bare \`terra::rast(as.array(img))\`
pattern used throughout the package. The key addition is \`terra::RGB(r)
\<- ...\`, which registers the channel indices so that
\`terra::plotRGB()\` works without the user having to know about this
terra requirement.

## Usage

``` r
cimg_to_spatrast(img, normalize = FALSE)
```

## Arguments

- img:

  A \`cimg\` object (any number of channels).

- normalize:

  Logical. If \`TRUE\`, rescale values to \[0, 1\].

## Value

A \`SpatRaster\` with RGB metadata set when the image has 3 or 4
channels.
