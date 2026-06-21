# Compare depth distributions using multiple metrics

Unified interface for comparing two probability distributions using
Wasserstein distance, Jensen-Shannon divergence, or KL divergence, with
optional depth-dependent weighting.

## Usage

``` r
compare_depth_distribution(
  P,
  Q,
  metric = "wasserstein",
  tail_weight = FALSE,
  weighting = "sigmoid",
  parameter = list(lambda = 0.2, x0 = 30),
  inverse = FALSE,
  index = 1:min(length(P), length(Q)),
  index.spacing = "equal",
  alignPQ = TRUE,
  cut = FALSE,
  baseline.weight = 0
)
```

## Arguments

- P:

  Numeric probability vector.

- Q:

  Numeric probability vector.

- metric:

  Character. One of `"wasserstein"` (default), `"js"`, or `"kl"`.

- tail_weight:

  Logical. If TRUE, enables depth-dependent weighting.

- weighting:

  Character. Weighting function used when `tail_weight = TRUE`.

- parameter:

  List with weighting parameters `lambda` and `x0`.

- inverse:

  Logical. Reverse weighting direction.

- index:

  Numeric vector defining depth or ordering.

- index.spacing:

  Character. `"equal"` or `"custom"`.

- alignPQ:

  Logical. Aligns unequal-length vectors.

- cut:

  Logical. If TRUE, truncates longer vector when aligning.

- baseline.weight:

  Numeric in \[0,1\]. Minimum weight applied.

## Value

Numeric distance or divergence value.

## Details

This function wraps:

- [`tail_weighted_wasserstein_distance()`](https://jcunow.github.io/Rootopia/reference/tail_weighted_wasserstein_distance.md)

- [`tail_weighted_js_divergence()`](https://jcunow.github.io/Rootopia/reference/tail_weighted_js_divergence.md)

- [`tail_weighted_kl_divergence()`](https://jcunow.github.io/Rootopia/reference/tail_weighted_kl_divergence.md)

If `tail_weight = FALSE`, uniform weighting is used.

## Data requirements

`P` and `Q` are non-negative vectors representing discrete probability
mass functions over a shared ordered index (e.g. depth bins).

They are internally normalized to sum to 1 if necessary.

Typical construction:


    counts_P <- c(5, 10, 20, 15)
    counts_Q <- c(2, 8, 25, 20)

    P <- counts_P / sum(counts_P)
    Q <- counts_Q / sum(counts_Q)

Continuous data must be discretised (e.g. binning) before use.

## Alignment

If `P` and `Q` differ in length, they are made comparable as follows:

**cut = TRUE:** both vectors are truncated to the shared minimum length.

**cut = FALSE:** the shorter vector is padded with zeros.

No interpolation or re-binning is performed. Alignment only enforces
equal vector length; consistent bin definitions are assumed.

## Interpretation

Elements of `P` and `Q` correspond to probability mass at positions
given by `index`. The same index must apply to both distributions for
valid comparison.

## Examples

``` r
counts_P <- c(5, 10, 20, 15)
counts_Q <- c(2, 8, 25, 20)

P <- counts_P / sum(counts_P)
Q <- counts_Q / sum(counts_Q) 
compare_depth_distribution(P, Q)
#> [1] 0.05909091
compare_depth_distribution(P, Q, metric = "js")
#> [1] 0.01228139
compare_depth_distribution(P, Q, tail_weight = TRUE, weighting = "sigmoid")
#> [1] 0.04693082
```
