# Preprocess an edge band before correlation

Preprocess an edge band before correlation

## Usage

``` r
stitch_preprocess(m, method = "none")
```

## Arguments

- m:

  A numeric matrix (edge band).

- method:

  One of `"none"`, `"center"` (subtract mean), `"norm"` (divide by SD),
  `"center_norm"`, `"hann"` (demean + separable Hann window), `"grad"`
  (central-difference gradient magnitude) or `"grad_norm"`.

## Value

The preprocessed matrix.
