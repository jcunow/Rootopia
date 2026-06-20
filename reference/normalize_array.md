# Rescale the array according to \`scale\`

All conversions use fixed factors (255), never a per-image max, and are
guarded so a conversion is a no-op when the data is already in the
target range:

- \`"none"\` leave values untouched

- \`"to_01"\` 0-255 -\> 0-1 (divide by 255; skipped if already \<= 1)

- \`"to_255"\` 0-1 -\> 0-255 (multiply by 255; skipped if already \> 1)

- \`"binary"\` strictly 0/1 via ceiling(arr / max)

## Usage

``` r
normalize_array(arr, scale = "none")
```
