# Tail-weighted Kullback-Leibler divergence

Computes KL divergence between two probability distributions with
optional depth-dependent weighting along an index (e.g., soil depth).

## Usage

``` r
tail_weighted_kl_divergence(
  P,
  Q,
  index = 1:min(length(P), length(Q)),
  index.spacing = "equal",
  parameter = list(lambda = 0.2, x0 = 30),
  weighting = "constant",
  baseline.weight = 0,
  inverse = FALSE,
  alignPQ = TRUE,
  cut = FALSE
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

- alignPQ:

  Logical. If TRUE, aligns vectors of unequal length.

- cut:

  Logical. If TRUE, truncates longer vector instead of padding.

## Value

Numeric KL divergence value.

## Details

KL divergence is asymmetric: \$\$D\_{KL}(P \parallel Q) = \sum P \log(P
/ Q)\$\$

Weighting modifies contribution of each index position prior to
normalization.
