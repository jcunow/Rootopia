# Threshold or deblur an image to binarize features

Turns a grayscale image into a binary one, either against a single
image-wide level (`method = "global"`) or against the local brightness
around each pixel (`method = "adaptive"`). A separate `deblur = TRUE`
mode instead separates two layers of a segmented image.

## Usage

``` r
image_threshold(
  img,
  threshold = 0.4,
  method = "global",
  adaptive_k = 1.15,
  window_size = 51,
  select_layer = 2,
  mask_layer = 1,
  binary_01 = FALSE,
  deblur = FALSE,
  invert = FALSE
)
```

## Arguments

- img:

  SpatRaster object.

- threshold:

  Fraction of the global maximum, between 0 and 1, used by
  `method = "global"` and by `deblur = TRUE`. Not used by
  `method = "adaptive"`, which takes `adaptive_k` instead.

- method:

  `"global"` or `"adaptive"`. Ignored when `deblur = TRUE`.

- adaptive_k:

  Multiple of the local mean a pixel must reach to count as foreground,
  for `method = "adaptive"` only. Must be above 1 for the method to
  separate anything; see the section above.

- window_size:

  Odd integer, the side of the local-mean box in pixels, for
  `method = "adaptive"` only.

- select_layer:

  Integer or NULL. The layer the threshold is computed on. The output
  has as many layers as `img`, and every layer except `mask_layer`
  receives this layer's result. If NULL and the image is multilayer, the
  threshold is computed on the mean of all layers and that one result
  goes into every layer.

- mask_layer:

  Integer or NULL. A layer that keeps its own result instead of
  receiving `select_layer`'s: it is thresholded on its own pixel values
  and against its own level (its own maximum, or its own local mean for
  `method = "adaptive"`). It excludes no pixels. The default 1 suits
  RootDetector output, where layer 1 marks foreign objects such as tape
  and layer 2 the roots, so the tape stays in its own layer after
  thresholding. NULL treats every layer alike.

- binary_01:

  Logical. Coding of the output: TRUE gives 0/1; FALSE gives 0 and the
  maximum of the thresholded layer, e.g. 0/255 for 8-bit input.

- deblur:

  Logical. If TRUE, runs the deblurring mode instead of thresholding.

- invert:

  Logical. FALSE, the default, keeps the pixels at or above the cut, so
  the brighter side is the foreground: roots on a dark tray. TRUE keeps
  the pixels at or below it, for the commoner flatbed case of dark roots
  on a bright tray. `threshold` remains a fraction of the image maximum
  and `adaptive_k` a factor on the local mean; with `invert = TRUE` the
  local cut is the mean divided by `adaptive_k` rather than multiplied
  by it. Ignored when `deblur = TRUE`.

## Value

SpatRaster object.

## Global thresholding

A pixel is foreground when it is at least `threshold` times the
brightest value in the image. One level for the whole image, so it only
works when the illumination is even: a flatbed scan with a bright center
and dark corners has no single level that is right everywhere.

## Adaptive thresholding

A pixel is foreground when it is at least `adaptive_k` times the mean of
the `window_size` x `window_size` box around it. The reference level
follows the illumination, so a brightness gradient leaves the decision
unchanged.

Two constraints decide whether it works, and both are properties of the
image rather than preferences:

**`adaptive_k` must exceed 1.** In a window holding only background the
local mean *is* the background level, so about half the background sits
above it. At `adaptive_k = 1` roughly half of a root-free patch
therefore comes out as root, and below 1 nearly all of it does. The
method needs a dead band above 1, and the band has to clear the
background's own noise: the noisier the background, the larger
`adaptive_k` has to be. Its upper bound is the root-to-background
brightness ratio – past that the roots fail their own test and nothing
is found. The default suits a reasonably clean scan; a noisy one wants
more, and a scan with little root-to-background contrast leaves no
useful range at all.

**`window_size` must be large compared with the widest root.** The
window average includes the root itself, so a root that fills much of
its own window raises its own reference level and is rejected. Recall is
what fails first: the roots are suppressed rather than the background. A
window several times the widest root is a reasonable starting point;
measure the widest root in pixels on one scan and scale from there,
since a root's pixel width depends on the scan resolution. A larger
window also leaves the full root-to-background contrast available, so it
allows a larger `adaptive_k`; set the window first, then the multiplier.
Going large costs little, because the local mean is computed as two 1-D
passes and so scales with the window side rather than its square.

## What neither method can do

Both methods mark pixels that are *brighter* than a reference, so they
separate bright features on a darker background and nothing else. Where
root and soil are about equally bright, the useful range for
`adaptive_k` closes: the margin the background noise demands is larger
than the root-to-background contrast allows, and no setting works. The
bundled minirhizotron scan is such an image – scored against the
segmentation shipped with it, its roots are not reliably brighter than
the surrounding soil, and the soil's own texture varies more than root
differs from soil. That is a property of the image rather than of the
implementation, and tuning does not rescue it. Scans like it need a
trained segmenter such as RootDetector, whose binary output the rest of
Rootopia reads directly.

## Deblurring

`deblur = TRUE` ignores `method`. It subtracts `select_layer` from
`mask_layer` and thresholds both at `threshold` times the maximum of
`select_layer`, which separates in-focus structure from an out-of-focus
halo in a two-layer segmented image.

## Examples

``` r
img = terra::rast(seg_Oulanka2023_Session03_T067)
image_threshold(img, threshold = 0.3, method = "global")
#> class       : SpatRaster
#> size        : 1161, 4900, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 4900, 0, 1161  (xmin, xmax, ymin, ymax)
#> coord. ref. : 
#> source(s)   : memory
#> names       : lyr.1, lyr.2, lyr.2
#> min values  :     0,     0,     0
#> max values  :   255,   255,   255
image_threshold(img, threshold = 0.4, select_layer = 2, mask_layer = 1, deblur = TRUE)
#> class       : SpatRaster
#> size        : 1161, 4900, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 4900, 0, 1161  (xmin, xmax, ymin, ymax)
#> coord. ref. : 
#> source(s)   : memory
#> names       : lyr.1, lyr.2, lyr.2
#> min values  :     0,     0,     0
#> max values  :   255,   255,   255

# Adaptive thresholding, on a small image that meets both conditions above:
# the background falls from 0.40 to 0.11 across the frame, a 3.6-fold range
# no single global level can follow, and three 1 px roots sit at twice their
# own local background. The window is far wider than the roots, and the
# multiplier is above 1.
m <- outer(seq(1, 0.4, length.out = 120), seq(1, 0.7, length.out = 120)) * 0.4
m[, c(20, 50, 90)] <- m[, c(20, 50, 90)] * 2
bw <- image_threshold(terra::rast(m), method = "adaptive",
                      adaptive_k = 1.15, window_size = 31, binary_01 = TRUE)
terra::plot(bw)
```
