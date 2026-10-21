# How the insertion angle splits a tube's depth between its two axes.
#
# Depth along the tube's length varies across columns; depth around its
# circumference varies down the rows. The angle trades one for the other, so
# their ratio is the geometry and everything else is scale.
#
# Nothing is written to disk.

library(Rootopia)
library(terra)

tube_cm <- 7
dpi     <- 150
nc      <- 400                                       # columns = depth axis
nr      <- round(tube_cm * pi * dpi / 2.54, 0) - 1   # 1298 = exactly one turn

img   <- terra::rast(matrix(0, nrow = nr, ncol = nc))
angles <- seq(0, 90, by = 5)   # degrees from vertical

# One map at a time: keep the two scalars, drop the raster. Down one column the
# axis term is constant, so that spread is the curvature alone; across one row
# the curvature is constant, so that spread is the axis term alone.
res <- lapply(angles, function(tl) {
  dm <- tryCatch(
    create_depthmap(img, sinoid = TRUE, tube_thicc = tube_cm,
                    angle_from_vertical = tl,
                    dpi = dpi, start_soil = 0, center_offset = 0),
    error = function(e) NULL)
  if (is.null(dm))
    return(list(curvature = NA_real_, axis = NA_real_, prof = rep(NA_real_, nr)))
  m <- terra::as.matrix(dm, wide = TRUE)
  list(curvature = diff(range(m[, 1])),
       axis      = diff(range(m[1, ])),
       prof      = m[, 1])          # one column = the curvature wave itself
})

d <- data.frame(angle     = angles,
                curvature = vapply(res, `[[`, numeric(1), "curvature"),
                axis      = vapply(res, `[[`, numeric(1), "axis"))
prof <- lapply(res, `[[`, "prof")
skipped <- d$angle[is.na(d$curvature)]
if (length(skipped))
  message("angle rejected by this build: ", paste(skipped, collapse = ", "))

# Measured ratio of the two terms.
d$ratio        <- ifelse(d$axis > 0, d$curvature / d$axis, NA_real_)
axis_px        <- (nc - 1) * 2.54 / dpi
d$ratio_expected <- tube_cm * sin(d$angle * pi/180) /
                    (axis_px * cos(d$angle * pi/180))

print(round(d, 3), row.names = FALSE)

# --- 1. the two terms across the sweep --------------------------------------
op <- par(mfrow = c(2, 1), mar = c(4, 4, 2.5, 1))

plot(d$angle, d$curvature, type = "o", pch = 16, col = "coral3", ylim = c(0, tube_cm),
     xlab = "insertion angle (degrees from vertical)", ylab = "cm",
     main = "The two terms of the depth map")
lines(d$angle, d$axis, type = "o", pch = 16, col = "steelblue4")
lines(d$angle, tube_cm * sin(d$angle * pi/180), lty = 2, col = "coral3")
legend("left", bty = "n", lwd = c(2, 2, 1), lty = c(1, 1, 2),
       col = c("coral3", "steelblue4", "coral3"),
       legend = c("around the circumference (measured)",
                  "along the tube length (measured)",
                  "expected: tube * sin(angle)"))

# --- 2. the ratio: the part that is geometry rather than scale ---------------
ok <- is.finite(d$ratio) & is.finite(d$ratio_expected) & d$ratio_expected > 0
plot(d$angle[ok], d$ratio[ok], type = "o", pch = 16, col = "grey30", log = "y",
     ylim = range(c(d$ratio[ok], d$ratio_expected[ok])),
     xlab = "insertion angle (degrees from vertical)",
     ylab = "curvature / axis span  (log scale)",
     main = "Circumference-to-depth ratio")
lines(d$angle[ok], d$ratio_expected[ok], lty = 2, col = "coral3", lwd = 2)
abline(v = 45, col = "grey70", lty = 3)
legend("topright", bty = "n", lwd = c(2, 2), lty = c(1, 2),
       col = c("grey30", "coral3"), legend = c("measured", "expected"))
par(op)

# --- 3. depth down the rows, one line per angle ------------------------------
# This is the curvature wave on its own: one full turn of the tube, from the
# tube top (row 1, shallowest) round to the bottom and back.
ramp <- hcl.colors(length(angles), "Viridis")
op <- par(mar = c(4, 4, 2.5, 1))
plot(NA, xlim = c(1, nr), ylim = range(unlist(prof), na.rm = TRUE),
     xlab = "row (around the tube circumference)", ylab = "depth (cm)",
     main = "Depth down the rows, by insertion angle")
for (i in seq_along(prof))
  if (!all(is.na(prof[[i]]))) lines(seq_len(nr), prof[[i]], col = ramp[i], lwd = 1.6)
key <- which(angles %% 15 == 0)
legend("topright", bty = "n", lwd = 2, col = ramp[key], ncol = 2,
       legend = sprintf("angle = %d", angles[key]))
par(op)

# --- 4. does the angle change the geometry? ---------------------------------------------------
# A ratio that moves with the angle means depth is genuinely being traded
# between the tube's length and its circumference.
spread <- diff(range(d$ratio[is.finite(d$ratio)]))
cat(sprintf("\nmeasured ratio varies by %.3g over angle 0-90\n", spread))
cat(if (spread < 1e-6)
      "-> flat: the angle only rescales the map, which is wrong.\n"
    else
      "-> the ratio moves with the angle, as tube geometry requires.\n")
