# Order the frames of one tube

Frames are stitched in the order given, so this decides the sequence.
The default key is the file name, which sorts lexicographically and
therefore puts `_10` before `_2` when frame numbers are not zero-padded.
`order_by` extracts a key from each path instead: a number if the match
parses as one, otherwise the matched text.

## Usage

``` r
stitch_order_frames(paths, order_by = NULL, decreasing = FALSE)
```

## Arguments

- paths:

  Character vector of file paths for one tube.

- order_by:

  Optional PCRE pattern whose match in each path is the sort key, e.g.
  `"\d+(?=\.)"` for the number before the extension. `NULL` sorts by the
  whole path.

- decreasing:

  Logical; reverse the order.

## Value

`paths`, reordered.
