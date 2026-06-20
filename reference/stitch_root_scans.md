# Batch-stitch grouped scan sequences (tubes) into mosaics

High-level driver ported from the Python `ImageStitcher` `main` /
`process_subset`. Files are discovered, optionally subset, grouped by an
id pattern (one group per tube), sorted within each group, and stitched
into one mosaic per group. Single-frame groups are passed through
unchanged. Mosaics are returned and, optionally, written to disk as
PNGs.

## Usage

``` r
stitch_root_scans(
  input,
  pattern = NULL,
  group_regex = "T0\\d{2}",
  select = NULL,
  tubes = NULL,
  out_dir = NULL,
  out_prefix = "",
  out_format = "png",
  method = "phase",
  edge_width = 250,
  vertical_region = 1000,
  vertical_offset = 300,
  direction = "horizontal",
  preprocess = "none",
  blend = "linear",
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

  Optional substring used to keep only matching file names (e.g.
  `".tiff"`). `NULL` keeps all files.

- group_regex:

  Regular expression identifying the group id within each path. Default
  `"T0\d{2}"` matches tube labels such as `T067`. Use `NULL` to stitch
  every file into a single mosaic.

- select:

  Optional integer vector of indices into the (sorted) *file* list, e.g.
  `1:36`. See
  [`list_scan_files`](https://jcunow.github.io/Rootopia/reference/list_scan_files.md).
  `NULL` uses all files. Applied before grouping.

- tubes:

  Optional *tube* selection: integer indices into the sorted tube list
  (e.g. `1:36`, see
  [`list_tubes`](https://jcunow.github.io/Rootopia/reference/list_tubes.md)),
  a character vector of tube names (e.g. `c("T037", "T040")`), or the
  string `"ask"` to print the tubes and choose a range interactively in
  one call (interactive sessions only). `NULL` keeps all tubes.

- out_dir:

  Optional directory to write one mosaic per tube to (named
  `<out_prefix><tube>.<ext>`). Created if needed. `NULL` (default)
  returns mosaics only.

- out_prefix:

  Filename prefix for written mosaics.

- out_format:

  Output image format when `out_dir` is set: `"png"` (default) or
  `"tiff"`. Requires the corresponding package (png or tiff).

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
  -\> 0, good for colour scans), `"overlay"` (img2 hides img1),
  `"overlay_first"` (img1 hides img2), `"max"` (lighten / union -
  recommended for segmented/binary masks, where averaging would make
  fractional values and ghost thin roots) or `"min"` (darken).

- blend_width:

  Optional width (px) of the linear ramp, centred in the overlap (hard
  img1 to its left, hard img2 to its right); `NULL` (default) ramps
  across the whole overlap. A smaller value reduces root ghosting on
  colour scans. Ignored unless `blend = "linear"`.

- report:

  Logical. If `TRUE`, return a list with both the mosaics and a per-step
  performance table instead of just the mosaics (see Value).

- verbose:

  Logical; print per-tube progress and a mean/min alignment peak.

## Value

If `report = FALSE` (default), invisibly a named list of mosaics (one
numeric `(H, W, C)` array per tube). If `report = TRUE`, a list
`list(mosaics, report)` where `report` is a data frame with columns
`tube`, `step`, `dx`, `dy`, `peak` (confidence; higher is better) and
`overlap` (`= edge_width - dx`).

## Details

Call
[`list_tubes`](https://jcunow.github.io/Rootopia/reference/list_tubes.md)
first to see the tube names, then pass a range to `tubes` (e.g.
`tubes = 1:36`) to stitch just those tubes.

## See also

[`list_tubes`](https://jcunow.github.io/Rootopia/reference/list_tubes.md),
[`list_scan_files`](https://jcunow.github.io/Rootopia/reference/list_scan_files.md),
[`stitch_image_sequence`](https://jcunow.github.io/Rootopia/reference/stitch_image_sequence.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# 1) See the tubes (names + frame counts)
list_tubes("path/to/scans", pattern = ".tiff")

# 2) Stitch the first 36 tubes, with a performance report and a preprocess
res <- stitch_root_scans("path/to/scans", pattern = ".tiff",
                         tubes = 1:36, preprocess = "grad", report = TRUE)
res$report
aggregate(peak ~ tube, res$report, mean)

# 3) A named subset, written straight to PNG
stitch_root_scans("path/to/scans", pattern = ".tiff",
                  tubes = c("T037", "T040"), out_dir = "path/to/output")
} # }
```
