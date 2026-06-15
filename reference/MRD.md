# Calculate Mean Rooting Depth

Calculate Mean Rooting Depth

## Usage

``` r
MRD(w, roots)
```

## Arguments

- w:

  Numeric vector of weights (typically depths)

- roots:

  Numeric vector of root coverage values

## Value

Numeric MRD value

## Examples

``` r
w <- seq(5, 25, 5)
roots <- c(0, 10, 7, 3, 1)
mean_rooting_depth <- MRD(w, roots)
```
