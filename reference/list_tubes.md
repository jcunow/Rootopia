# List the tubes (groups) found in a scan folder

Summarises the unique groups (tubes) that
[`stitch_root_scans`](https://jcunow.github.io/RootScanR/reference/stitch_root_scans.md)
would build, so you can see the tube names and pick a range to stitch
(e.g. `tubes = 1:36`).

## Usage

``` r
list_tubes(input, pattern = NULL, group_regex = "T0\\d{2}")
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

A data frame with columns `index` (1-based, the value to pass to
`stitch_root_scans(tubes = ...)`), `tube` (group id) and `n_frames`.

## See also

[`stitch_root_scans`](https://jcunow.github.io/RootScanR/reference/stitch_root_scans.md),
[`list_scan_files`](https://jcunow.github.io/RootScanR/reference/list_scan_files.md)

## Examples

``` r
d <- file.path(tempdir(), "tubes_demo")
dir.create(d, showWarnings = FALSE)
file.create(file.path(d, c("S_T037_01.tiff", "S_T037_02.tiff",
                           "S_T038_01.tiff", "S_T039_01.tiff")))
#> [1] TRUE TRUE TRUE TRUE
list_tubes(d, pattern = ".tiff")
#> Error in list_tubes(d, pattern = ".tiff"): could not find function "list_tubes"
unlink(d, recursive = TRUE)
```
