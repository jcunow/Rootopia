# Assess rhythmicity via sine curve fitting and model comparison

Fits a sine curve to time-series data and computes rhythmicity
statistics. The amplitude significance is tested via various statistical
methods.

## Usage

``` r
rhythmicity(
  x,
  y,
  fix_period = 24,
  method = "F",
  parStart = list(amp = 3, phase = 0, offset = 0, period = 12)
)
```

## Arguments

- x:

  Numeric vector of time or phase values.

- y:

  Numeric vector of observed values (same length as `x`).

- fix_period:

  Either numeric (e.g., 24) to fix the period, or `NULL` to estimate it.

- method:

  Type of test to use: `"F"` for F-test (default), `"FLR"` for finite
  sample likelihood ratio test, or `"LR"` for large sample likelihood
  ratio test.

- parStart:

  Named list with starting values for parameters: `amp`, `phase`,
  `offset`, and `period`. Only used when `fix_period` is `NULL`.

## Value

A named list with the following elements:

- `amplitude`: Fitted amplitude

- `phase`: Fitted phase

- `offset`: Fitted offset

- `period`: Fitted or fixed period

- `peakTime`: Peak time (phase + period/4)

- `R2`: Coefficient of determination

- `residual_se`: Residual standard error (F-test only)

- `df_residual`: Residual degrees of freedom (F-test only)

- `sigma02`: Variance under null (LR methods only)

- `sigmaA2`: Variance under alternative (LR methods only)

- `l0`: Log-likelihood under null (LR methods only)

- `l1`: Log-likelihood under alternative (LR methods only)

- `df`: Degrees of freedom for the test

- `stat`: Test statistic

- `p_value`: P-value for the test

## Examples

``` r
set.seed(123)
x <- seq(0, 48, length.out = 100)
y <- 1.5 * sin(2 * pi / 24 * (x + 4)) + 4 + rnorm(100, 0, 1.5)

# Fixed period with F-test
rhythmicity(x, y, fix_period = 24, method = "F")
#> $amplitude
#> [1] 1.694882
#> 
#> $phase
#> [1] 4.786127
#> 
#> $offset
#> [1] 4.132499
#> 
#> $period
#> [1] 24
#> 
#> $peakTime
#> [1] 1.213873
#> 
#> $R2
#> [1] 0.4480709
#> 
#> $residual_se
#> [1] 1.35584
#> 
#> $df_residual
#> [1] 97
#> 
#> $df
#>   df1 df2
#> 1   2  97
#> 
#> $stat
#>     F_stat
#> 1 39.37361
#> 
#> $p_value
#> [1] 3.02928e-13
#> 

# Fixed period with finite sample LR test
rhythmicity(x, y, fix_period = 24, method = "FLR")
#> $amplitude
#> [1] 1.694882
#> 
#> $phase
#> [1] 4.786127
#> 
#> $offset
#> [1] 4.132499
#> 
#> $period
#> [1] 24
#> 
#> $peakTime
#> [1] 1.213873
#> 
#> $R2
#> [1] 0.4480709
#> 
#> $sigma02
#> [1] 3.230765
#> 
#> $sigmaA2
#> [1] 1.783153
#> 
#> $l0
#> [1] -200.5298
#> 
#> $l1
#> [1] -170.813
#> 
#> $df
#>   df1 df2
#> 1   2  97
#> 
#> $stat
#>   FLR_stat
#> 1 39.37361
#> 
#> $p_value
#> [1] 3.02928e-13
#> 

# Estimate period with F-test
rhythmicity(x, y, fix_period = NULL, method = "F")
#> $amplitude
#>      amp 
#> 0.731638 
#> 
#> $phase
#>     phase 
#> 0.1453872 
#> 
#> $offset
#>   offset 
#> 4.082082 
#> 
#> $period
#>   period 
#> 13.82691 
#> 
#> $peakTime
#>  period 
#> 3.31134 
#> 
#> $R2
#> [1] 0.08131785
#> 
#> $residual_se
#> [1] 1.758327
#> 
#> $df_residual
#> [1] 96
#> 
#> $df
#>   df1 df2
#> 1   3  96
#> 
#> $stat
#>     F_stat
#> 1 2.832505
#> 
#> $p_value
#> [1] 0.04235487
#> 
```
