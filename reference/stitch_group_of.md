# Extract a group id (e.g. tube label) from each path

Extract a group id (e.g. tube label) from each path

## Usage

``` r
stitch_group_of(x, group_regex)
```

## Arguments

- x:

  Character vector of paths.

- group_regex:

  Regular expression for the group id, or `NULL` to place everything in
  a single group.

## Value

Character vector of group ids (`NA` where no match).
