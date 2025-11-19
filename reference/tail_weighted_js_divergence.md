# Calculate tail-weighted Jensen-Shannon divergence

Calculate tail-weighted Jensen-Shannon divergence

## Usage

``` r
tail_weighted_js_divergence(
  P,
  Q,
  parameter = list(lambda = 0.2, x0 = 30),
  method = "constant",
  inverse = FALSE,
  alignPQ = TRUE,
  cut = FALSE,
  index = 1:min(c(length(Q), length(P))),
  index.spacing = "equal"
)
```

## Arguments

- P:

  probability vector 2

- Q:

  probability vector 1

- parameter:

  list with lambda -\> shape parameter (0 = constant weighting) & x0 -\>
  curve offset (= inflexion point )

- method:

  weighting function along index. Available options are: c("constant",
  "asymptotic", "linear, "exponential", "sigmoid", "gompertz","step")

- inverse:

  changes from right tail to left tail if TRUE

- alignPQ:

  if TRUE, index end values will be cut off in case of unequal length of
  P & Q so that length of P & Q is equal

- cut:

  if FALSE, 0 will be added to the shorter vector. If TRUE, the longer
  vector will be shortened at the end.

- index:

  a positive numeric vector containing probability spacing e.g., depth

- index.spacing:

  whether index intervals are equally distant i.e., c(1,2,3,4....n), if
  "equal" than index is c(1,n)

## Value

Numeric JS divergence value

## Examples

``` r
P <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)  # Distribution P
Q <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)**6  # Distribution Q

# Ensure the distributions are valid (non-negative and sum to 1)
P <- P / sum(P)
Q <- Q / sum(Q)

tail_weighted_js_divergence(P,Q,parameter = list(lambda = 0.2,x0=30))
#> [1] 0.197
```
