# Example Flatbed Root Scan (segmented)

A segmented flatbed root scan used to demonstrate the depth-independent
flatbed workflow (skeletonisation, root length, diameter, branching, and
branch order). Roots are foreground, background is 0.

## Usage

``` r
data(flatbed_scan_example)
```

## Format

A 3-dimensional numeric array (rows x columns x layers), matching the
other bundled datasets. Rebuild a SpatRaster with
[`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html).
Dimensions:

- 2000 rows (height)

- 2000 columns (width)

- 3 layers (segmentation channels)

## Details

Flatbed scans capture the full root system in a 2D plane, so traits are
computed globally rather than per depth bin. The image is a binary
segmentation (root vs. background) stored across three channels; layer 2
is used as the root channel in the vignette.

## References

Cunow J, Pijcke F, Olofsson J, Väisänen M, Blume-Werry G (accepted,
2026). Reindeer grazing induces spatial and functional shifts in root
systems of boreal pine forests. Oikos. https://doi.org/10.1002/oik.12211

## Author

Johannes Cunow <johannes.cunow@gmail.com>

Johannes Cunow Original flatbed scans from boreal pine forests in
northern Finland, acquired by Johannes Cunow and others (field campaign
2022).

## Examples

``` r
if (FALSE) { # \dontrun{
  data(flatbed_scan_example)
  seg <- terra::rast(flatbed_scan_example)
  seg <- terra::ifel(seg[[2]] > 0, 1, 0)
  terra::plot(seg, maxcell = Inf)
} # }
```
