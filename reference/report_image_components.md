# Report sizes of holes and white artifacts in binary images

Analyzes a binary image and prints the sizes of internal black holes
(0-valued regions enclosed by 1s) and isolated white artifacts (1-valued
regions not connected to the border). This function is useful for
diagnosing what would be affected by \`fill_holes()\` and
\`remove_small_objects()\` operations.

## Usage

``` r
report_image_components(img)
```

## Arguments

- img:

  A \`cimg\` object representing a binary image with values 0 and 1.

## Value

None (invisible \`NULL\`). The function prints human-readable summaries
to the console.

## Details

Holes are detected by inverting the image and applying connected
component labeling. Any region that does not touch the image border is
considered a candidate hole or artifact.

\- Holes are black regions (0-valued pixels) completely enclosed by
white (1-valued) areas. - Artifacts are small white regions (1-valued
pixels) that are not connected to the image border. - Pixel counts are
printed for each detected region.

## Examples

``` r
#' # Create a complex test image with holes and artifacts

  img <- imager::as.cimg(matrix(0, 150, 150))  # Start with black background

  # Create multiple white objects with black holes
  img[20:50, 20:50] <- 1       # White square 1
  img[30:35, 30:35] <- 0       # Small black hole in square 1

  img[70:120, 70:120] <- 1     # White square 2
  img[80:85, 80:85] <- 0       # Small black hole 1 in square 2
  img[100:115, 100:115] <- 0   # Large black hole 2 in square 2
  
# Add small artifacts (1-pixel specks)
img[10, 140] <- 1
img[145, 15] <- 1

# Add a 2×2 speck
img[130:131, 40:41] <- 1

# Add an irregular blob
img[100:102, 10] <- 1
img[101:102, 11] <- 1
img[101, 12] <- 1

  # Create a white ring (donut shape)
  center_x <- 40
  center_y <- 100
for (i in 1:150) {
 for (j in 1:150) {
   dist <- sqrt((i - center_x)^2 + (j - center_y)^2)
   if (dist <= 20 && dist >= 10) {
     img[i, j,,] <- 1  
 }}}
 
report_image_components(img)
#> === HOLES (black regions inside objects) ===
#> Hole 1: 0 pixels
#> Hole 2: 0 pixels
#> Hole 3: 0 pixels
#> Hole 4: 36 pixels
#> Hole 5: 0 pixels
#> Hole 6: 0 pixels
#> Hole 7: 0 pixels
#> Hole 8: 36 pixels
#> Hole 9: 305 pixels
#> Hole 10: 256 pixels
#> Hole 11: 0 pixels
#> 
#> === ARTIFACTS (small isolated white objects) ===
#> Artifact 1: 6 pixels
#> Artifact 2: 1 pixels
#> Artifact 3: 925 pixels
#> Artifact 4: 0 pixels
#> Artifact 5: 4 pixels
#> Artifact 6: 2309 pixels
#> Artifact 7: 952 pixels
#> Artifact 8: 0 pixels
#> Artifact 9: 0 pixels
#> Artifact 10: 0 pixels
#> Artifact 11: 1 pixels
```
