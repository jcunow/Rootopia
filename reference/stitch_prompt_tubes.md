# Interactively prompt for a tube selection

Interactively prompt for a tube selection

## Usage

``` r
stitch_prompt_tubes(unique_groups, groups)
```

## Arguments

- unique_groups:

  Sorted character vector of tube ids.

- groups:

  Character vector of every kept file's tube id (for counts).

## Value

Integer indices into `unique_groups`, or `NULL` for all.
