# Batch-stitch grouped scan sequences (tubes) into mosaics

Stitches sets of related scan images into continuous mosaics. Input
files are optionally filtered, grouped by an identifier extracted from
filenames, and stitched sequentially within each group. By default,
files are stitched in lexicographic filename order, but custom
within-group ordering can be specified.

## Usage

``` r
stitch_root_scans(
  input,
  pattern = NULL,
  group_regex = "T0\\d{2}",
  select = NULL,
  tubes = NULL,
  order_by = NULL,
  decreasing = FALSE,
  out_dir = NULL,
  out_prefix = "",
  method = "phase",
  edge_width = 250,
  vertical_region = 1000,
  vertical_offset = 300,
  direction = "horizontal",
  preprocess = "none",
  blend = "overlay_first",
  blend_width = NULL,
  report = FALSE,
  verbose = TRUE
)
```

## Arguments

- input:

  Either a directory (searched recursively) or a character vector of
  image file paths.

- pattern:

  Optional filter applied to filenames (e.g. `".tiff"`). If `NULL`, all
  files are used.

- group_regex:

  Regular expression used to define grouping keys (e.g. tube identifiers
  such as `T067`). If `NULL`, all files are treated as a single group.

- select:

  Optional integer index vector selecting a subset of input files prior
  to grouping.

- tubes:

  Optional selection of groups to process. Can be integer indices into
  the sorted group list, character group names, or `"ask"` for
  interactive selection. If `NULL`, all groups are processed.

- order_by:

  Optional regular expression used to extract an ordering key from
  filenames within each group. If `NULL` (default), files are stitched
  in lexicographic filename order. When a match is found, numeric
  components are ordered numerically; otherwise ordering is
  lexicographic.

- decreasing:

  Logical. If `TRUE`, reverses the within-group order. Useful when
  acquisition order is reversed (e.g. bottom-to-top scans or
  deepest-to-shallowest sequences).

- out_dir:

  Optional output directory for writing mosaics as PNG files. If `NULL`,
  results are returned only.

- out_prefix:

  Filename prefix for exported mosaics.

- method:

  Alignment method. Currently only `"phase"` (FFT phase correlation) is
  supported. `"feature"` (SIFT/ORB keypoint matching in the Python
  original) requires OpenCV and has no CRAN-compatible R backend; it
  raises an informative error.

- edge_width:

  Width in pixels of the edge band used for alignment. Clamped to the
  image width. Set it close to the true overlap - roughly 1-2x (so the
  overlap is between about half and all of `edge_width`) - for the most
  reliable peak.

- vertical_region:

  Height in pixels of the vertical band used for alignment.

- vertical_offset:

  Starting row (from the top) of the vertical band.

- direction:

  `"horizontal"` (default, frames side by side) or `"vertical"` (frames
  stacked top to bottom, e.g. minirhizotron strips acquired down the
  tube). Vertical transposes the inputs, runs the same horizontal core,
  and transposes the result back.

- preprocess:

  Preprocessing applied to the edge bands before correlation: one of
  `"none"` (default), `"center"` (subtract mean), `"norm"` (divide by
  SD), `"center_norm"`, `"hann"` (demean + Hann window), `"grad"`
  (gradient magnitude) or `"grad_norm"`. `"hann"`, `"grad"` and
  `"center"` help on scans with uneven lighting or sparse texture;
  `"norm"` alone barely changes phase correlation (it is already
  scale-invariant). On well-textured scans `"none"` is usually fine.

- blend:

  How the overlap band is combined: `"linear"` (default, alpha ramp 1
  -\> 0, good for colour scans), `"overlay_second"` (img2 hides img1),
  `"overlay_first"` (img1 hides img2), `"max"` (lighten / union -
  recommended for segmented/binary masks, where averaging would make
  fractional values and ghost thin roots) or `"min"` (darken).

- blend_width:

  Optional width (px) of the linear ramp, centred in the overlap (hard
  img1 to its left, hard img2 to its right); `NULL` (default) ramps
  across the whole overlap. A smaller value reduces root ghosting on
  colour scans. Ignored unless `blend = "linear"`.

- report:

  Logical. If `TRUE`, returns additional alignment diagnostics per
  stitching step.

- verbose:

  Logical; enables progress output during processing.

## Value

If `report = FALSE`, a named list of mosaics (one array per group). If
`report = TRUE`, a list:

- mosaics:

  Named list of stitched images.

- report:

  Data frame with per-step diagnostics: `tube`, `step`, `dx`, `dy`,
  `peak`, `overlap`.

## Details

Each group (tube) produces one stitched image. Groups containing a
single frame are returned unchanged.

## See also

[`list_tubes`](https://jcunow.github.io/RootScanR/reference/list_tubes.md),
[`list_scan_files`](https://jcunow.github.io/RootScanR/reference/list_scan_files.md),
[`stitch_image_sequence`](https://jcunow.github.io/RootScanR/reference/stitch_image_sequence.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Default: lexicographic filename order within each tube
res <- stitch_root_scans(
  "path/to/scans",
  pattern = ".tiff"
)

# Order by scan depth code (e.g. L001, L002, ...)
res <- stitch_root_scans(
  "path/to/scans",
  pattern = ".tiff",
  order_by = "L\\d{3}"
)

# Reverse acquisition order (e.g. bottom-to-top scans)
res <- stitch_root_scans(
  "path/to/scans",
  pattern = ".tiff",
  order_by = "L\\d{3}",
  decreasing = TRUE
)

# Order by date-stamped filenames (YYYYMMDD)
res <- stitch_root_scans(
  "path/to/scans",
  pattern = ".tiff",
  order_by = "\\d{8}"
)
} # }
```
