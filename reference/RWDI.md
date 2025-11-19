# Calculate Root Weight Depth Index

Calculate Root Weight Depth Index

## Usage

``` r
RWDI(w, roots)
```

## Arguments

- w:

  Numeric vector of weights (typically depths)

- roots:

  Numeric vector of root coverage values

## Value

Numeric RWDI value

## Examples

``` r
w <- seq(5, 25, 5)
roots <- c(0, 10, 7, 3, 1)
rwdi <- RWDI(w, roots)
```
