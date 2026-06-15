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

  Character vector specifying grouping variable(s)

- depth:

  Character string specifying depth column name

- variable:

  Character string specifying accumulating values column

- stdrz:

  Character string specifying standardization method

## Value

Numeric vector of accumulated values

## Examples

``` r
df = data.frame(depth = c(seq(0,80,20),seq(0,80,20)),
               Plot = c(rep("a",5),rep("b",5)), rootpx = c(5,50,20,15,5,10,40,30,10,5) )
accum_root = root_accumulation(df,group = "Plot", depth = "depth", variable = "rootpx")
```
