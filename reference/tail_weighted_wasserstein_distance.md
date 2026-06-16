# Tail-weighted Wasserstein distance (1D)

Computes the 1D Wasserstein distance between two probability
distributions with optional depth-dependent weighting.

## Usage

``` r
tail_weighted_wasserstein_distance(
  P,
  Q,
  index = 1:min(length(P), length(Q)),
  index.spacing = "equal",
  parameter = list(lambda = 0.2, x0 = 10),
  weighting = "constant",
  baseline.weight = 0,
  inverse = FALSE
)
```

## Arguments

- P:

  Numeric probability vector.

- Q:

  Numeric probability vector.

- index:

  Numeric vector defining position (e.g., depth layers).

- index.spacing:

  Character. `"equal"` or `"custom"`.

- parameter:

  List with `lambda` and `x0`.

- weighting:

  Character. Weighting function applied before comparison.

- baseline.weight:

  Numeric in \[0,1\]. Minimum weight.

- inverse:

  Logical. Reverse weighting direction.

## Value

Numeric Wasserstein distance.

## Details

Wasserstein distance measures the minimal "transport cost" required to
transform one distribution into another along a one-dimensional index.

This implementation uses
[`transport::wasserstein1d()`](https://rdrr.io/pkg/transport/man/wasserstein1d.html)
after applying optional weighting.
