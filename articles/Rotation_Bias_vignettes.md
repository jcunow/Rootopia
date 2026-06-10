# Rotation Bias and Rhythmicity Analysis with RootScanR

## Rotation Bias Analysis

### Introduction

This vignette demonstrates how to detect and correct for rotational bias
in minirhizotron setups using the **RootScanR** package. In
minirhizotron studies, the scanner tube is inserted at an angle into the
soil and can rotate slightly between sessions. Because the CI-600 and
similar scanners do not cover a full 360° arc, sequential images from
the same tube may not perfectly overlap. Left uncorrected, this rotation
introduces a systematic spatial bias — roots near the scan edges are
observed in some sessions but not others.

RootScanR addresses this through three calibration functions:

- [`estimate_rotation_center()`](https://jcunow.github.io/RootScanR/reference/estimate_rotation_center.md)
  — locates the rotational zero point from tape coverage in a single
  image
- [`estimate_rotation_shift()`](https://jcunow.github.io/RootScanR/reference/estimate_rotation_shift.md)
  — quantifies the pixel offset between two sessions using image
  correlation
- [`rotation_censor()`](https://jcunow.github.io/RootScanR/reference/rotation_censor.md)
  — crops images to the shared, overlap region

### Installation

``` r

# install.packages("remotes")
# remotes::install_github("jcunow/RootScanR")

library(RootScanR)
library(terra)
```

------------------------------------------------------------------------

### Workflow

#### 1. Estimate the rotational center from a single image

[`estimate_rotation_center()`](https://jcunow.github.io/RootScanR/reference/estimate_rotation_center.md)
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
| `tape.brightness` | Brightness threshold (0–1) for classifying tape pixels |
| `search.area` | Fraction of image width to analyse (tape is near one edge) |
| `nclasses` | Number of unsupervised clustering classes |
| `tape.quantile` | Quantile used to align the brightness scale |

------------------------------------------------------------------------

#### 2. Estimate the rotation shift between two sessions

When the same tube is scanned in two different sessions, the scanner may
have rotated.
[`estimate_rotation_shift()`](https://jcunow.github.io/RootScanR/reference/estimate_rotation_shift.md)
uses either cross-correlation (`"ccf"`) or phase correlation (`"phase"`)
on a shared depth window to find the pixel offset. Phase correlation is
generally more robust to brightness differences between sessions.

``` r

data(seg_Oulanka2023_Session01_T067)
data(seg_Oulanka2023_Session03_T067)

shift <- estimate_rotation_shift(
  seg_Oulanka2023_Session01_T067,
  seg_Oulanka2023_Session03_T067,
  cor.type          = "phase",
  fixed.depth.pixel = c(1000, 4000)
)

cat("Rotation shift (depth px, rotation px):", shift, "\n")
```

The returned vector is `c(x.lag, y.lag)` — the horizontal (depth axis)
and vertical (rotation axis) pixel shifts. A large vertical shift means
the tube rotated substantially between sessions.

------------------------------------------------------------------------

#### 3. Censor image edges to the shared overlap region

[`rotation_censor()`](https://jcunow.github.io/RootScanR/reference/rotation_censor.md)
crops each image so that only the rows present in every session are
retained. This eliminates the non-overlapping margins and makes root
counts directly comparable across sessions.

Two modes are available:

- **`fixed.rotation = FALSE`** — cuts proportionally based on the
  measured offset; output width varies
- **`fixed.rotation = TRUE`** — centres the crop on a specified row and
  forces a fixed output width; recommended when comparing multiple
  sessions

``` r

data(seg_Oulanka2023_Session01_T067)
img <- terra::rast(seg_Oulanka2023_Session01_T067)

# Fixed mode — crop to a 1800-pixel-wide window centred on the rotation center
censored <- rotation_censor(
  img,
  center.offset  = r0,
  cut.buffer     = 0.02,
  fixed.width    = 1800,
  fixed.rotation = TRUE
)

terra::plot(censored, main = "Censored image")
```

> **Note on tube geometry.** The inner and outer tube diameters differ,
> so the observed root length slightly underestimates the true root
> length in the soil. A resize coefficient may be applied separately;
> [`rotation_censor()`](https://jcunow.github.io/RootScanR/reference/rotation_censor.md)
> does not handle this correction.

------------------------------------------------------------------------

### Rhythmicity Analysis

In some deployments, minirhizotron observations are made at high
temporal frequency (e.g. multiple scans per day over a short window). In
these datasets it can be useful to test whether root observations show a
rhythmic pattern across the sampling cycle.

[`rhythmicity()`](https://jcunow.github.io/RootScanR/reference/rhythmicity.md)
fits a sine curve of the form

``` math
y = A \sin\!\left(\frac{2\pi}{P}(x + \phi)\right) + c
```

and tests whether the amplitude $`A`$ is significantly different from
zero, using your choice of three statistical tests.

> **On period choice**: the period $`P`$ should reflect a biologically
> or methodologically meaningful cycle in your sampling design — for
> example the length of your monitoring campaign, or the interval at
> which conditions repeat. Use `fix_period = NULL` only if you genuinely
> want to estimate the period from the data and have enough observations
> to support a four-parameter fit.

#### Fitting and testing

``` r

set.seed(42)
x <- seq(0, 48, length.out = 120)
y <- 1.8 * sin(2 * pi / 24 * (x + 5)) + 4 + rnorm(120, 0, 1.2)

# F-test (default) — most interpretable for small samples
res_F <- rhythmicity(x, y, fix_period = 24, method = "F")
cat(sprintf("Amplitude: %.2f | Phase: %.1f | R²: %.3f | p: %.4f\n",
            res_F$amplitude, res_F$phase, res_F$R2, res_F$p_value))

# Finite-sample likelihood ratio test
res_FLR <- rhythmicity(x, y, fix_period = 24, method = "FLR")

# Large-sample likelihood ratio test
res_LR  <- rhythmicity(x, y, fix_period = 24, method = "LR")
```

#### Estimating an unknown period

``` r

res_free <- rhythmicity(x, y,
                        fix_period = NULL,
                        parStart   = list(amp = 2, phase = 0,
                                          offset = 4, period = 20))

cat(sprintf("Estimated period: %.1f | R²: %.3f | p: %.4f\n",
            res_free$period, res_free$R2, res_free$p_value))
```

> **Caution with free-period fits.** With four free parameters the
> likelihood surface can be multimodal. Provide starting values grounded
> in your sampling design and always inspect the fitted curve before
> interpreting results.

#### Plotting the fitted curve

``` r

x_seq <- seq(min(x), max(x), length.out = 300)
y_fit <- res_F$amplitude *
         sin(2 * pi / res_F$period * (x_seq + res_F$phase)) +
         res_F$offset

plot(x, y, pch = 16, col = "grey50",
     xlab = "Observation index", ylab = "Root measurement",
     main = "Rhythmicity fit")
lines(x_seq, y_fit, col = "coral", lwd = 2)
legend("topright",
       legend = sprintf("A = %.2f, R² = %.3f, p = %.4f",
                        res_F$amplitude, res_F$R2, res_F$p_value),
       bty = "n")
```

#### Interpreting the output

| Element     | Description                          |
|-------------|--------------------------------------|
| `amplitude` | Estimated amplitude $`A`$            |
| `phase`     | Phase shift $`\phi`$ (modulo period) |
| `offset`    | Baseline $`c`$                       |
| `period`    | Fixed or estimated period $`P`$      |
| `peakTime`  | Time of peak within one cycle        |
| `R2`        | Coefficient of determination         |
| `stat`      | Test statistic (F, FLR, or LR)       |
| `p_value`   | P-value for $`A \neq 0`$             |

A significant p-value indicates a rhythmic component. Always inspect
$`R^2`$ alongside it — a significant but low-$`R^2`$ result may simply
reflect a noisy dataset rather than a meaningful pattern.

------------------------------------------------------------------------

### Further reading

- [Batch
  Processing](https://jcunow.github.io/RootScanR/articles/BatchProcessing_vignette.md)
  — the recommended starting point for processing multiple images
- [Minirhizotron
  Scans](https://jcunow.github.io/RootScanR/articles/MinirhizotronScans_vignettes.md)
  — step-by-step depth analysis workflow
- [Function
  reference](https://jcunow.github.io/RootScanR/reference/index.md)
- Source and issues: <https://github.com/jcunow/RootScanR>
