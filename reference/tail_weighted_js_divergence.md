# Tail-weighted Jensen-Shannon divergence

Computes a symmetric divergence between two distributions using a
weighted KL decomposition.

## Usage

``` r
tail_weighted_js_divergence(
  P,
  Q,
  index = NULL,
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

  Numeric positions of the bins (e.g. depths), one per bin after
  alignment. `NULL` (default) uses `1, 2, ..., n`.

- index.spacing:

  Character. `"equal"` or `"custom"`.

- parameter:

  List with `lambda` and `x0`.

- weighting:

  Character. Weighting function applied before comparison.

- baseline.weight:

  Numeric in `[0,1]`. Minimum weight.

- inverse:

  Logical. Reverse weighting direction.

- alignPQ:

  Logical. If TRUE, aligns vectors of unequal length; if FALSE, unequal
  lengths are an error.

- cut:

  Logical. If TRUE, truncates the longer vector instead of padding the
  shorter one with zeros.

## Value

Numeric JS divergence value.

## Details

Jensen-Shannon divergence is defined as: \$\$JS(P,Q) = 1/2 KL(P \|\|
M) + 1/2 KL(Q \|\| M)\$\$ where \\M = (P + Q)/2\\.

Weighting is applied prior to normalization of P, Q, and M.
