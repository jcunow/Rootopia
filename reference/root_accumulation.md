# Calculate root accumulation

Calculate root accumulation

## Usage

``` r
root_accumulation(x, group, depth, variable, stdrz = "counts")
```

## Arguments

- x:

  Data frame containing group, depth, and variable columns

- group:

  Character vector naming one or more grouping columns. With several,
  each combination of their values is accumulated separately, e.g.
  `group = c("Plot", "Site")`.

- depth:

  Character string specifying depth column name

- variable:

  Character string specifying accumulating values column

- stdrz:

  Character string naming the denominator applied to the cumulative sum:
  `"counts"` (default) returns it raw, `"additive"` divides by its
  maximum so each group's curve ends at 1, and `"relative"` divides by
  the sum of the cumulative values. Use `"additive"` to compare the
  shape of profiles between groups

## Value

Numeric vector of accumulated values, one per row of `x` and in the same
row order.

## Examples

``` r
df = data.frame(depth = c(seq(0,80,20),seq(0,80,20)),
               Plot = c(rep("a",5),rep("b",5)), rootpx = c(5,50,20,15,5,10,40,30,10,5) )
accum_root = root_accumulation(df,group = "Plot", depth = "depth", variable = "rootpx")

# Several grouping columns: one curve per Plot x Session combination
df2 <- rbind(transform(df, Session = "June"), transform(df, Session = "Oct"))
root_accumulation(df2, group = c("Plot", "Session"), depth = "depth",
                  variable = "rootpx")
#>  [1]  5 55 75 90 95 10 50 80 90 95  5 55 75 90 95 10 50 80 90 95
```
