# Calculate tail-weighted KL divergence for discrete distributions

Calculate tail-weighted KL divergence for discrete distributions

## Usage

``` r
tail_weighted_kl_divergence(
  P,
  Q,
  index = 1:min(c(length(Q), length(P))),
  index.spacing = "equal",
  parameter = list(lambda = 0.2, x0 = 30),
  cut = FALSE,
  inverse = FALSE,
  method = "constant",
  alignPQ = TRUE
)
```

## Arguments

- P:

  probability vector 2

- Q:

  probability vector 1

- index:

  a positive numeric vector containing probability spacing e.g., depth

- index.spacing:

  whether index intervals are equally distant i.e., c(1,2,3,4....n), if
  "equal" than index is c(1,n)

- parameter:

  list with lambda -\> shape parameter (0 = constant weighting) & x0 -\>
  curve offset (= inflexion point )

- cut:

  if FALSE, 0 will be added to the shorter vector. If TRUE, the longer
  vector will be shortened at the end.

- inverse:

  changes from right tail to left tail if TRUE

- method:

  weighting function along index. Available options are: c("constant",
  "asymptotic", "linear, "exponential", "sigmoid", "gompertz","step")

- alignPQ:

  if TRUE, index end values will be cut off in case of unequal length of
  P & Q so that length of P & Q is equal

## Value

KL divergence, not symmetrical - changing the input order will change
the result

## Details

Kullback-Leibler Divergence
