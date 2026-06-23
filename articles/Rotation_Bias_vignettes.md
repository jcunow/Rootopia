# Rotation Bias and Rhythmicity Analysis with Rootopia

## Rotation Bias Analysis

### Introduction

This vignette demonstrates how to detect and correct for rotational bias
in minirhizotron setups using the **Rootopia** package. In minirhizotron
studies, the scanner tube is inserted at an angle into the soil and can
rotate slightly between sessions. Because the CI-600 and similar
scanners do not cover a full 360° arc, sequential images from the same
tube may not perfectly overlap. Left uncorrected, this rotation
introduces a systematic spatial bias — roots near the scan edges are
observed in some sessions but not others.

Rootopia addresses this through four functions:

- [`estimate_rotation_center()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_center.md)
  — locates the rotational zero point from tape coverage in a single
  image
- [`estimate_rotation_shift()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_shift.md)
  — quantifies the pixel offset between two sessions using image
  correlation
- [`rotation_censor()`](https://jcunow.github.io/Rootopia/reference/rotation_censor.md)
  — crops images to the shared, overlap region
- `zoning(mode = "rotation")` — splits the tube surface into slices
  along the rotation axis (circumference), so that root traits can be
  summarized separately for each slice

The last point matters beyond rotation correction itself: once the tube
circumference is split into slices, the
[`rhythmicity()`](https://jcunow.github.io/Rootopia/reference/rhythmicity.md)
/
[`fit_sine_curve()`](https://jcunow.github.io/Rootopia/reference/fit_sine_curve.md)
family — normally introduced for time-series data — can be applied to
the sequence of per-slice root traits to test whether roots are
distributed **evenly around the tube circumference**, or whether there
is a systematic top-down / side-to-side pattern (e.g. more roots on the
underside of the tube). This is a *spatial*, not temporal, application
of those functions — see the [Circumferential zoning and
rhythmicity](#circumferential-zoning-and-rhythmicity) section below.

### Installation

``` r

# install.packages("remotes")
# remotes::install_github("jcunow/Rootopia")

library(Rootopia)
library(terra)
```

    ## terra 1.9.34

------------------------------------------------------------------------

### Workflow

#### 1. Estimate the rotational center from a single image

[`estimate_rotation_center()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_center.md)
detects the white adhesive tape attached to the upper side of the tube.
Because more tape is visible on one side, the function infers where the
top of the tube (rotational zero) lies. It returns a pixel row index.

Since a tube has a fixed geometry, this only needs to be done once per
tube — the rotational center does not change between sessions.

``` r

data(seg_Oulanka2023_Session01_T067)
img <- terra::rast(seg_Oulanka2023_Session01_T067)

r0 <- estimate_rotation_center(img)
print(paste("Estimated rotation center (row):", r0))
```

Key parameters:

| Parameter | Description |
|----|----|
| `tape_brightness` | Brightness threshold (0–1) for classifying tape pixels |
| `search_area` | Fraction of image width to analyze (tape is near one edge) |
| `nclasses` | Number of unsupervised clustering classes |
| `tape_quantile` | Quantile used to align the brightness scale |

------------------------------------------------------------------------

#### 2. Estimate the rotation shift between two sessions

When the same tube is scanned in two different sessions, the scanner may
have rotated.
[`estimate_rotation_shift()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_shift.md)
uses either cross-correlation (`"ccf"`) or phase correlation (`"phase"`)
on a shared depth window to find the pixel offset. Phase correlation is
generally more robust to brightness differences between sessions.

``` r

data(seg_Oulanka2023_Session01_T067)
data(seg_Oulanka2023_Session03_T067)

shift <- estimate_rotation_shift(
  seg_Oulanka2023_Session01_T067,
  seg_Oulanka2023_Session03_T067,
  cor_type          = "phase",
  fixed_depth_pixel = c(1000, 4000)
)
```

    ## Warning in doTryCatch(return(expr), name, parentenv, handler): Image size
    ## mismatch detected; cropping to common extent

``` r

cat("Rotation shift (depth px, rotation px): ", shift[1:2], "\n")
```

    ## Rotation shift (depth px, rotation px):  -18 -9

``` r

# Visual inspection
estimate_rotation_shift(seg_Oulanka2023_Session01_T067, 
                        seg_Oulanka2023_Session03_T067, 
                        cor_type = "phase", 
                        select_layer = 2, 
                        overlay = T)
```

    ## Warning in doTryCatch(return(expr), name, parentenv, handler): Image size
    ## mismatch detected; cropping to common extent

![](Rotation_Bias_vignettes_files/figure-html/unnamed-chunk-3-1.png)

The returned vector is `c(x.lag, y.lag)` — the horizontal (depth axis)
and vertical (rotation axis) pixel shifts. A large vertical shift means
the tube rotated substantially between sessions.

> **Planned**: a single helper that estimates
> [`estimate_rotation_shift()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_shift.md)
> between two sessions and then applies
> [`rotation_censor()`](https://jcunow.github.io/Rootopia/reference/rotation_censor.md)
> using that shift directly, so that aligning session coverage becomes a
> one-step operation. Until then, use the two functions as shown in
> steps 2 and 3.

------------------------------------------------------------------------

#### 3. Censor image edges to the shared overlap region

[`rotation_censor()`](https://jcunow.github.io/Rootopia/reference/rotation_censor.md)
crops each image so that only the rows present in every session are
retained. This eliminates the non-overlapping margins and makes root
counts directly comparable across sessions.

Two modes are available:

- **`fixed_rotation = FALSE`** — cuts proportionally based on the
  measured offset; output width varies
- **`fixed_rotation = TRUE`** — centers the crop on a specified row and
  forces a fixed output width; recommended when comparing multiple
  sessions

``` r

data(seg_Oulanka2023_Session01_T067)
img <- load_flexible_image(seg_Oulanka2023_Session01_T067, "spatrast")

# rotation center (absolute row) to center the crop on
r0 <- estimate_rotation_center(img)

# preview the crop (green = kept, red = cut) and apply it in one call.
# fixed_width must be <= image height (1144 rows here) to sit symmetrically.
censored <- rotation_censor(
  img,
  center_offset  = r0,
  cut_buffer     = 0.02,
  fixed_width    = 800,
  fixed_rotation = TRUE,
  overlay        = TRUE,
  main = "Rotation censor: kept window (green) vs cut (red)"
)
```

![](Rotation_Bias_vignettes_files/figure-html/unnamed-chunk-4-1.png)

> **Note on tube geometry.** The inner and outer tube diameters differ,
> so the observed root length slightly underestimates the true root
> length in the soil. A resize coefficient may be applied separately;
> [`rotation_censor()`](https://jcunow.github.io/Rootopia/reference/rotation_censor.md)
> does not handle this correction.

------------------------------------------------------------------------

### Circumferential zoning and rhythmicity

A minirhizotron image is a long, narrow strip that represents a slice of
the tube’s circumference, with depth running along one axis and the
rotation (circumferential) position running along the other.
`zoning(mode = "rotation")` divides the image along this rotation axis
into `rotation_total_slices` equal slices, and
`rotation_slices = c(i, i)` extracts a single slice `i`.

Looping over slices gives a sequence of root-trait values indexed by
circumferential position. The
[`rhythmicity()`](https://jcunow.github.io/Rootopia/reference/rhythmicity.md)
/
[`fit_sine_curve()`](https://jcunow.github.io/Rootopia/reference/fit_sine_curve.md)
functions — which fit
$`y = A \sin\!\left(\frac{2\pi}{P}(x + \phi)\right) + c`$ and test
whether $`A \neq 0`$ — are agnostic to what `x` represents. Applied to
this sequence with `x` = slice index and the period `P` fixed to
`rotation_total_slices` (one full turn around the tube), they test
whether roots are distributed **evenly around the tube** or concentrated
on one side (e.g. a systematic top-down bias from gravitropism, light
incidence on one face of the installation, or an installation-angle
artifact).

#### Extract a root trait per circumferential slice

``` r

data(seg_Oulanka2023_Session01_T067)

root_layer <- load_flexible_image(seg_Oulanka2023_Session01_T067, scale = "binary", select_layer = 2,
                                  output_format = "spatrast")


# rotation trim: drop the top 33% of rows
e <- terra::ext(root_layer)
root_layer <- terra::crop(root_layer, terra::ext(e[1]+ e[2] / 4, e[2], e[3], e[4])) 
terra::plot(root_layer, maxcell = Inf)
```

![](Rotation_Bias_vignettes_files/figure-html/unnamed-chunk-5-1.png)

``` r

# --> the trimming avoids noisy deep layer and tape artifacts at the top

# create circumferential zones
n_slices <- 48
slices <- slice_rotation(root_layer, n_slices)

# tallies from the binary mask: roots = 1, soil/background = 0
n_root <- function(s) { v <- terra::values(s, mat = FALSE); sum(v > 0,  na.rm = TRUE) }
n_soil <- function(s) { v <- terra::values(s, mat = FALSE); sum(v == 0, na.rm = TRUE) }

slice_traits <- data.frame(
  slice   = seq_len(n_slices),
  rootlen = sapply(slices, root_length, unit = "cm", dpi = 150),  # cm
  rootpx  = sapply(slices, n_root),
  soilpx  = sapply(slices, n_soil)
)
```

    ## Diagonal: 539 | Orthogonal: 655

    ## Diagonal: 1001 | Orthogonal: 1242

    ## Diagonal: 1266 | Orthogonal: 1552

    ## Diagonal: 1125 | Orthogonal: 1436

    ## Diagonal: 1254 | Orthogonal: 1545

    ## Diagonal: 1366 | Orthogonal: 1606

    ## Diagonal: 1290 | Orthogonal: 1516

    ## Diagonal: 1116 | Orthogonal: 1380

    ## Diagonal: 1543 | Orthogonal: 1873

    ## Diagonal: 1361 | Orthogonal: 1732

    ## Diagonal: 1765 | Orthogonal: 1883

    ## Diagonal: 2158 | Orthogonal: 2282

    ## Diagonal: 1666 | Orthogonal: 1738

    ## Diagonal: 1913 | Orthogonal: 2041

    ## Diagonal: 2333 | Orthogonal: 2597

    ## Diagonal: 2963 | Orthogonal: 3306

    ## Diagonal: 2361 | Orthogonal: 2845

    ## Diagonal: 2362 | Orthogonal: 2801

    ## Diagonal: 2265 | Orthogonal: 2483

    ## Diagonal: 2539 | Orthogonal: 2757

    ## Diagonal: 2547 | Orthogonal: 2875

    ## Diagonal: 3009 | Orthogonal: 3339

    ## Diagonal: 2260 | Orthogonal: 2654

    ## Diagonal: 2932 | Orthogonal: 3418

    ## Diagonal: 1936 | Orthogonal: 2234

    ## Diagonal: 2642 | Orthogonal: 3073

    ## Diagonal: 2362 | Orthogonal: 2847

    ## Diagonal: 2988 | Orthogonal: 3400

    ## Diagonal: 2584 | Orthogonal: 3010

    ## Diagonal: 2292 | Orthogonal: 2686

    ## Diagonal: 1713 | Orthogonal: 2028

    ## Diagonal: 2044 | Orthogonal: 2479

    ## Diagonal: 2298 | Orthogonal: 2772

    ## Diagonal: 1715 | Orthogonal: 1976

    ## Diagonal: 1145 | Orthogonal: 1433

    ## Diagonal: 1107 | Orthogonal: 1345

    ## Diagonal: 1041 | Orthogonal: 1243

    ## Diagonal: 846 | Orthogonal: 998

    ## Diagonal: 769 | Orthogonal: 937

    ## Diagonal: 472 | Orthogonal: 528

    ## Diagonal: 533 | Orthogonal: 621

    ## Diagonal: 639 | Orthogonal: 844

    ## Diagonal: 740 | Orthogonal: 959

    ## Diagonal: 648 | Orthogonal: 814

    ## Diagonal: 647 | Orthogonal: 803

    ## Diagonal: 1098 | Orthogonal: 1378

    ## Diagonal: 1522 | Orthogonal: 1898

    ## Diagonal: 866 | Orthogonal: 1433

``` r

# root length per observed area (root + soil pixels)
slice_traits$rld <- slice_traits$rootlen /
                    ((slice_traits$rootpx + slice_traits$soilpx) * (2.54/150)^2)
#slice_traits$rld.fraction <- slice_traits$rld / sum(slice_traits$rld, na.rm = TRUE)
```

#### Test for a circumferential (top-down) pattern

``` r

res <- rhythmicity(
  x          = slice_traits$slice,
  y          = slice_traits$rld,
  fix_period = n_slices,
  method     = "F",
  par_start   = list(amp = 0.01, phase = 0, offset = mean(slice_traits$rootlen),
                     period = n_slices)
)

cat(sprintf("Mean: %.4f | Amplitude: %.4f | Phase (slice): %.2f | R²: %.3f | p: %.4f\n",
            res$offset, res$amplitude, res$phase, res$R2, res$p_value))
```

    ## Mean: 2.7540 | Amplitude: 1.5502 | Phase (slice): 38.31 | R²: 0.801 | p: 0.0000

``` r

cat("The top-down difference is", round((res$amplitude*2) / res$offset,2), " -> circa twice as large as the mean"    )
```

    ## The top-down difference is 1.13  -> circa twice as large as the mean

A significant, non-trivial amplitude indicates that root abundance
varies systematically with circumferential position — i.e. there is more
root material on one side of the tube than the other. `res$phase` (and
`res$peakTime`) identify *which* slice the pattern peaks at.

#### Visualize the fit

``` r

slice_seq <- seq(1, n_slices, length.out = 200)
fit_seq   <- res$amplitude *
              sin(2 * pi / n_slices * (slice_seq + res$phase)) +
              res$offset

plot(slice_traits$slice, slice_traits$rld, pch = 16, col = "gray50",
     xlab = "Circumferential slice (1 = top, going around the tube)",
     ylab = "Root Length Density (cm / cm2)",
     main = "Root distribution around tube circumference")
lines(slice_seq, fit_seq, col = "coral", lwd = 2)
abline( a = res$offset, b = 0, col  = "gray40", lty = 3)

segments(x0 = res$peakTime,
         x1 = res$peakTime,
         y0 = res$offset,
         y1 = res$offset + res$amplitude,
         lwd = 2, lty = 2,
         col = "gray40")
segments(x0 = res$peakTime + res$period/2,
         x1 = res$peakTime +  res$period/2,
         y0 = res$offset,
         y1 = res$offset -  res$amplitude,
         lwd = 2, lty = 2,
         col = "gray40")
legend("topright",
       legend = sprintf("A = %.4f, R² = %.3f, p = %.4f",
                        res$amplitude, res$R2, res$p_value),
       bty = "n")
```

![](Rotation_Bias_vignettes_files/figure-html/unnamed-chunk-7-1.png)

``` r

# cross section view
library(ggplot2)
inner <- 6
ggplot(slice_traits, aes(slice, rld + inner)) +
    geom_col(fill = "steelblue") +
    coord_polar(start = pi) +
  xlab("") + ylab ("")+ ggtitle("Circumferential Root Length Distribution")+
  theme_light()+
  theme(axis.text = element_blank()) + 
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 0, ymax = inner,
           fill = "white", color = NA)
```

![](Rotation_Bias_vignettes_files/figure-html/unnamed-chunk-7-2.png)

> **Going further**: the same slice loop can be combined with
> `slice_rotation(mode = "both")` to additionally split each
> circumferential slice by depth bin, producing a slice x depth grid.
> Fitting
> [`rhythmicity()`](https://jcunow.github.io/Rootopia/reference/rhythmicity.md)
> separately within each depth bin (x = slice, y = trait per slice)
> tests whether the circumferential pattern is consistent with depth or
> changes from shallow to deep — i.e. a genuine *top-down amplitude*
> gradient.

------------------------------------------------------------------------

### Future directions

Two extensions are planned and not yet implemented:

1.  **Combined rotation-shift + censor helper** (mentioned in step 2):
    estimate
    [`estimate_rotation_shift()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_shift.md)
    between two sessions and automatically apply
    [`rotation_censor()`](https://jcunow.github.io/Rootopia/reference/rotation_censor.md)
    with that shift, so aligning session coverage is a single call.
2.  **Estimating the rotation center from the amplitude peak**: if the
    circumferential analysis above reveals a strong, consistent
    sinusoidal pattern, the fitted phase (`res$phase` / `res$peakTime`)
    marks the circumferential position of peak root density. This could
    be used as a data-driven alternative or complement to
    [`estimate_rotation_center()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_center.md),
    which currently relies only on tape coverage. Though, field
    calibrations are absolutely the preferred option.

------------------------------------------------------------------------

### Further reading

- [Batch
  Processing](https://jcunow.github.io/Rootopia/articles/BatchProcessing_vignette.md)
  — the recommended starting point for processing multiple images
- [Minirhizotron
  Scans](https://jcunow.github.io/Rootopia/articles/MinirhizotronScans_vignettes.md)
  — step-by-step depth analysis workflow
- [Function
  reference](https://jcunow.github.io/Rootopia/reference/index.md)
- Source and issues: <https://github.com/jcunow/Rootopia>
