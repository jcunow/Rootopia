# Calculate weights for probability distribution comparison

Calculate weights for probability distribution comparison

## Usage

``` r
tail_weight_function(
  index = NULL,
  parameter = list(lambda = 0.2, x0 = 5),
  index.spacing = "equal",
  method = "sigmoid",
  baseline.weight = 0,
  inverse = FALSE
)
```

## Arguments

- index:

  Numeric vector specifying the positions for weight calculation

- parameter:

  List containing: - lambda: Shape parameter (0 = constant weighting) -
  x0: Curve offset/inflection point

- index.spacing:

  Character, either "equal" or "custom" for index spacing type

- method:

  Character, weighting function type: "constant", "asymptotic",
  "linear", "exponential", "sigmoid", "gompertz", "step"

- baseline.weight:

  Numeric between 0-1, minimum weight value

- inverse:

  Logical, if TRUE reverses weight distribution (left vs right tail)

## Value

Normalized weight vector
