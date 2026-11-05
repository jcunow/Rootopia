# Tail-weighted Wasserstein distance (1D)

Computes the 1D Wasserstein distance between two probability
distributions with optional depth-dependent weighting.

## Usage

``` r
tail_weighted_wasserstein_distance(
  P,
  Q,
  index = NULL,
  index.spacing = "equal",
  parameter = list(lambda = 0.2, x0 = 10),
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

Numeric Wasserstein distance, in the units of `index`.

## Details

Wasserstein distance measures the minimal transport cost required to
transform one distribution into another along a one-dimensional index:
the mass-weighted mean distance the probability has to travel, in the
units of `index`. For depth profiles indexed in centimeters it is
therefore a depth, and a profile shifted `d` cm returns exactly `d`.

`P` and `Q` are the mass at the positions in `index`. On a shared
one-dimensional support the distance is the area between the two
cumulative distributions: after sorting by `index`, the absolute
difference of the cumulative masses times the gap to the next position,
summed over all gaps. A profile with all its mass in one bin is a valid
input.
