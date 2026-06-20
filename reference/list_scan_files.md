# List scan files with index and group id

Discovers the image files that
[`stitch_root_scans`](https://jcunow.github.io/Rootopia/reference/stitch_root_scans.md)
would process and returns them as an indexed table, so you can see what
is in a folder and pick a range of files to stitch (e.g.
`select = 1:36`). For tube-level selection see
[`list_tubes`](https://jcunow.github.io/Rootopia/reference/list_tubes.md).
Uses the same discovery and sort order as `stitch_root_scans`, so the
indices line up.

## Usage

``` r
list_scan_files(input, pattern = NULL, group_regex = "T0\\d{2}")
```

## Arguments

- input:

  Either a directory (searched recursively) or a character vector of
  image file paths.

- pattern:

  Optional substring used to keep only matching file names (e.g.
  `".tiff"`). `NULL` keeps all files.

- group_regex:

  Regular expression identifying the group id within each path. Default
  `"T0\d{2}"` matches tube labels such as `T067`.

## Value

A data frame with columns `index` (1-based position), `file` (full path)
and `group` (matched id, or `NA`).

## See also

[`list_tubes`](https://jcunow.github.io/Rootopia/reference/list_tubes.md),
[`stitch_root_scans`](https://jcunow.github.io/Rootopia/reference/stitch_root_scans.md)
