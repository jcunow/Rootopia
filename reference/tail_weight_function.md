# Calculate weights for probability distribution comparison

Constructs a flexible weighting vector used to emphasize or de-emphasize
portions of a probability distribution along an index (e.g., depth).

## Usage

``` r
tail_weight_function(
  index = NULL,
  parameter = list(lambda = 0.2, x0 = 5),
  index.spacing = "equal",
  weighting = "sigmoid",
  baseline.weight = 0,
  inverse = FALSE
)
```

## Arguments

- index:

  Numeric vector specifying positions.

- parameter:

  List with:

  - `lambda`: shape parameter controlling steepness

  - `x0`: inflection or threshold parameter

- index.spacing:

  Character. Either `"equal"` or `"custom"`.

- weighting:

  Character. Weighting function: `"constant"`, `"asymptotic"`,
  `"linear"`, `"exponential"`, `"sigmoid"`, `"gompertz"`, `"step"`.

- baseline.weight:

  Numeric in \[0,1\]. Minimum weight applied.

- inverse:

  Logical. If TRUE, reverses weighting direction.

## Value

Numeric vector of normalized weights.
