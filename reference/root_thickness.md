# Approximate average Root Thickness

Estimates mean root diameter from total root length and root pixel
count.

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
```
