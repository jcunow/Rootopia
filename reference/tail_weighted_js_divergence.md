# Tail-weighted Jensen-Shannon divergence

Computes a symmetric divergence between two distributions using a
weighted KL decomposition.

## Usage

``` r
tail_weighted_js_divergence(
  P,
  Q,
  index = 1:min(length(P), length(Q)),
  index.spacing = "equal",
  parameter = list(lambda = 0.2, x0 = 30),
  weighting = "constant",
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

- inverse:

  Logical. Reverse weighting direction.

- alignPQ:

  Logical. If TRUE, aligns vectors of unequal length.

- cut:

  Logical. If TRUE, truncates longer vector instead of padding.

## Value

Numeric JS divergence value.

## Details

Jensen-Shannon divergence is defined as: \$\$JS(P,Q) = 1/2 KL(P \|\|
M) + 1/2 KL(Q \|\| M)\$\$ where \\M = (P + Q)/2\\.

Weighting is applied prior to normalization of P, Q, and M.
