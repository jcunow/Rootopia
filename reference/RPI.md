# Calculate Root Penetration Index

Calculate Root Penetration Index

## Usage

``` r
RPI(roots, w)
```

## Arguments

- roots:

  Numeric vector of root coverage values

- w:

  Numeric vector of weights (typically depths)

## Value

Numeric RPI value between -1 and 1

## Examples

``` r
w <- seq(5, 25, 5)
roots <- c(0, 10, 7, 3, 1)
rpi <- RPI(roots, w)
```
