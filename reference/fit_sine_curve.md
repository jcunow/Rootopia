# Fit a sine curve to data with optional fixed period

Fits a sine curve model: \\ y = amp \* sin(2 \* pi / period \* (x +
phase)) + offset \\ to data points (x, y) using nonlinear least squares.

## Usage

``` r
fit_sine_curve(
  x,
  y,
  parStart = list(amp = 3, phase = 0, offset = 0, period = 24),
  fix_period = NULL
)
```

## Arguments

- x:

  Numeric vector of independent variable values (e.g., time).

- y:

  Numeric vector of dependent variable values.

- parStart:

  Named list of starting values for parameters: - amp: Amplitude
  (default 3) - phase: Phase shift (default 0) - offset: Vertical offset
  (default 0) - period: Period (default 24)

- fix_period:

  Numeric or NULL. If numeric, period is fixed to this value and not
  estimated. If NULL (default), period is estimated.

## Value

A list containing: - amp: Estimated amplitude - phase: Estimated phase
(modulo period) - offset: Estimated offset - period: Estimated or fixed
period - tss: Total sum of squares - rss: Residual sum of squares - R2:
Coefficient of determination - residual_se: Residual standard error -
df_residual: Residual degrees of freedom - predicted: Fitted values -
residuals: Residuals (y - predicted)

## Examples

``` r
set.seed(1)
x <- seq(0, 48, length.out = 100)
period = 30.7
y <- 2 * sin(2 * pi / period * (x + 3)) + 5 + rnorm(100, 0, 0.8)
fit_3parameters <- fit_sine_curve(x, y, fix_period = 24)
fit_3parameters$R2
#> [1] 0.4241515
fit_3parameters$period
#> [1] 24
# estimate period
fit_4parameters <- fit_sine_curve(x, y, fix_period = NULL)
fit_4parameters$R2 
#> [1] 0.7948954
fit_4parameters$period
#>   period 
#> 30.79209 

 plot(x,y)
 lines(x,fit_3parameters$predicted, col = "steelblue")
 lines(x,fit_4parameters$predicted, col = "coral")
```
