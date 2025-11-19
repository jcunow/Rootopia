# Calculate a circular mean to determine average Directionality

Calculate a circular mean to determine average Directionality

## Usage

``` r
circular_mean(angles, input_units = "degrees", output_units = "degrees")
```

## Arguments

- angles:

  Numeric vector of input angles

- input_units:

  Character string specifying input units ("radians" or "degrees")

- output_units:

  Character string specifying output units ("radians" or "degrees")

## Value

Numeric value representing the average angle

## Examples

``` r
circular_mean(angles = c(360,90,0), input_units = "degrees", output_units = "degrees")
#> [1] 26.56505
```
