# A tailweighted Version of 1 dimensional Wasserstein distance betwwen two probability vectors

A tailweighted Version of 1 dimensional Wasserstein distance betwwen two
probability vectors

## Usage

``` r
tail_weighted_wasserstein_distance(
  Q,
  P,
  inverse = FALSE,
  parameter = list(lambda = 0.2, x0 = 10),
  method = "step",
  baseline.weight = 0,
  index = 1:min(c(length(Q), length(P))),
  index.spacing = "equal"
)
```

## Arguments

- Q:

  probability vector 1

- P:

  probability vector 2

- inverse:

  changes from right tail to left tail if TRUE

- parameter:

  list with lambda -\> shape parameter (0 = constant weighting) & x0 -\>
  curve offset (= inflexion point )

- method:

  weighting function along index. Available options are: c("constant",
  "asymptotic", "linear, "exponential", "sigmoid", "gompertz","step")

- baseline.weight:

  Numeric between 0-1

- index:

  a positive numeric vector containing probability spacing e.g., depth

- index.spacing:

  whether index intervals are equally distant i.e., c(1,2,3,4....n), if
  "equal" than index is c(1,n)

## Value

Numeric Wasserstein distance

## Examples

``` r
P <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)  # Distribution P
Q <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)**6  # Distribution Q

# Ensure the distributions are valid (non-negative and sum to 1)
P <- P / sum(P)
Q <- Q / sum(Q)

tail_weighted_wasserstein_distance(P,Q,
  inverse=FALSE,method="constant",parameter = list(lambda = 0.2,x0=3))
#> [1] 0.1008476
```
