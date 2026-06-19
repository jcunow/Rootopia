# Write an (H, W, C) 0-255 array to PNG or TIFF (chosen by file extension)

Write an (H, W, C) 0-255 array to PNG or TIFF (chosen by file extension)

## Usage

``` r
stitch_write_image(a, path)
```

## Arguments

- a:

  A numeric `(H, W, C)` array (0-255).

- path:

  Output path; the extension (`.png`, `.tif`/`.tiff`) selects the
  format.

## Value

The path, invisibly.
