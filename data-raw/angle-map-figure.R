# The root angle panel, on its own.
#
# The six-panel figure cut down to panel 5: the D8 direction class of every
# skeleton pixel, over a black root field, with a zoom inset and a colour
# wheel. Everything the other five panels needed is gone -- the RGB scan, the
# diameter fill, the depth bands, the skeleton classes as a drawn layer.
#
# What stays is the depth map, because the directions are derived from it, and
# the skeleton, because they are drawn on it.

library(Rootopia)
library(terra)

# ---- paths ----
seg_path <- "C:/Users/jocu0013/Desktop/Oulanka/Scans_Blended_Segmented/Oulanka2025_02/Split_top_Oulanka2025_Session02_T067.png/Split_top_Oulanka2025_Session02_T067.png.segmentation.png"
out_path <- "C:/Users/jocu0013/Desktop/Oulanka/Presentation/angle_T067.png"

# ---- settings: scan and depth map ----
dpi                 <- 300
tube_thicc          <- 7
angle_from_vertical <- 45        # degrees from vertical; create_depthmap() reads `tilt` as 90 - tilt
start_soil          <- 0
center_offset       <- 0.465
depth_range_cm      <- c(0, 40)  # cut the scan to this depth range (NULL = keep everything)
crop_top_frac       <- 0.10      # also cut this fraction of the scan length at the shallow end (tape, mask); 0 = none

# ---- settings: colours ----
# All eight D8 classes get a colour, because upward classes do occur: the tube
# curvature term puts depth on the circumference as well as on the length, so a
# step up the page can still reach a deeper pixel. 90 right, 180 down = deeper,
# 270 left. Index k + 1 is the class k * 45, so the colours stay bound to the
# same directions from tube to tube and two figures compare directly.
# Which classes a scan actually contains is reported below and drives the wheel.
ang_all  <- seq(0, 315, by = 45)
ang_cols <- c("#FFFFFF",   #   0  up
              "#FF8A3D",   #  45
              "#FFD60A",   #  90  right
              "#4D7CFF",   # 135
              "#3DDC84",   # 180  down, deeper
              "#E0218A",   # 225
              "#00E5FF",   # 270  left
              "#B388FF")   # 315
mask_col <- "grey85"                            # masked (foreign object) areas

# ---- settings: zoom inset ----
zoom_on   <- TRUE
zoom_ctr  <- c(0.81, 0.39)   # c(across, down) as fractions of the panel, or "densest" = window with the most skeleton pixels
zoom_w    <- 0.07            # window width as a fraction of the panel width (square window); magnification = inset_w / zoom_w
inset_w   <- 0.65            # inset width as a fraction of the panel width
inset_pos <- "topright"      # "topright", "topleft", "bottomright", "bottomleft"
zoom_col  <- "firebrick4"

# ---- settings: output size and text ----
agg        <- 2              # source pixels per output pixel (1 = full resolution, huge file)
res_dpi    <- 600 / agg      # keeps the physical figure size and text proportions fixed for any agg
title_cex  <- 3
legend_cex <- 2              # text size of the legend that sits on top of the panel
# title_cex and legend_cex are the values the six-panel figure used, and there
# they were multiplied by the base cex R sets for a multi-figure layout: 0.66
# once a layout has more than two rows or columns. A single panel gets no such
# shrink, so applying it here keeps the text the size it is in that figure.
# Set it to 1 to see the unshrunk sizes.
panel_cex  <- 0.66

# cache files for the slow steps; the settings they depend on are part of the file name
crop_tag   <- if (is.null(depth_range_cm)) "all" else paste0(depth_range_cm[1], "-", depth_range_cm[2], "cm")
cache_path <- function(what, ...) sub("\\.png$", paste0("_", paste(c(what, crop_tag, ...), collapse = "_"), ".tif"), out_path)

# ---- load ----
seg_img <- load_flexible_image(seg_path, output_format = "spatrast", scale = "none")

# ---- mask: red pixels in the segmentation = foreign objects (1 = masked) ----
seg_max <- max(global(seg_img, "max", na.rm = TRUE)[[1]])
red     <- seg_img[[1]] / seg_max
green   <- seg_img[[2]] / seg_max
mask    <- (red > 0.5 & green < 0.5) * 1

# ---- depth map (on the full scan, so depths stay correct) ----
depth_map <- create_depthmap(img = seg_img, mask = mask,
                             sinoid = TRUE, tube_thicc = tube_thicc,
                             angle_from_vertical = angle_from_vertical,
                             dpi = dpi, start_soil = start_soil, center_offset = center_offset)
ext(depth_map) <- ext(seg_img)

# ---- cut the scan to the chosen depth range (depth runs along the columns) ----
if (!is.null(depth_range_cm) || crop_top_frac > 0) {
  dm   <- as.matrix(depth_map, wide = TRUE)
  # keep only columns that lie completely inside the range, so the cut is a straight rectangular edge
  # (depth varies along the rows because of the tube curvature, so some deeper pixels are lost at the edge)
  keep <- seq_len(ncol(dm))
  if (!is.null(depth_range_cm)) {
    col_max <- apply(dm, 2, max, na.rm = TRUE)
    col_min <- apply(dm, 2, min, na.rm = TRUE)
    keep <- which(col_max <= depth_range_cm[2] & (depth_range_cm[1] <= 0 | col_min >= depth_range_cm[1]))
  }
  # cut the shallow end (tape, mask): the first columns, i.e. the top of the rotated figure
  keep <- keep[keep > min(keep) + round(crop_top_frac * length(keep))]
  crop_e <- ext(xmin(seg_img) + (min(keep) - 1) * xres(seg_img),
                xmin(seg_img) +  max(keep)      * xres(seg_img),
                ymin(seg_img), ymax(seg_img))
  seg_img   <- crop(seg_img,   crop_e)
  depth_map <- crop(depth_map, crop_e)
  mask      <- crop(mask,      crop_e)
  crop_tag <- paste0(crop_tag, "_cols", min(keep), "-", max(keep))   # cached files from another crop are never reused
  rm(dm)
}
for (nm in c("seg_img", "depth_map", "mask")) {     # plain cell coordinates from here on
  r <- get(nm); ext(r) <- c(0, ncol(r), 0, nrow(r)); assign(nm, r)
}

# ---- root area (0/1) ----
root_bin <- seg_img[[2]]
root_bin[is.na(root_bin)] <- 0
root_bin <- (root_bin / global(root_bin, "max", na.rm = TRUE)[[1]] > 0.5) * 1

# ---- skeleton: full resolution (1 px wide), cached ----
if (file.exists(cache_path("skeleton"))) {
  skel <- rast(cache_path("skeleton"))
} else {
  skel <- skeletonize_image(root_bin, verbose = FALSE)
  writeRaster(skel, cache_path("skeleton"), overwrite = TRUE)
}
skel <- (skel[[1]] > 0) * 1
skel[is.na(skel)] <- 0
ext(skel) <- ext(seg_img)

# ---- root angle on the skeleton (D8 classes from deep_drive), cached ----
# NOTE: deep_drive() returns its map already transposed + flipped, i.e. turned the same way as the
# figure, so it is added AFTER the rotation below and must not be rotated again
if (file.exists(cache_path("angle", dpi, tube_thicc, angle_from_vertical, center_offset, start_soil))) {
  ang_full <- rast(cache_path("angle", dpi, tube_thicc, angle_from_vertical, center_offset, start_soil))
} else {
  ang_full <- deep_drive(DepthMap = depth_map, RootMap = skel, return = "all")$angle_map
  writeRaster(ang_full, cache_path("angle", dpi, tube_thicc, angle_from_vertical, center_offset, start_soil), overwrite = TRUE)
}
stopifnot(all(dim(ang_full)[1:2] == rev(dim(seg_img)[1:2])))     # width and height are swapped
ext(ang_full) <- c(0, ncol(ang_full), 0, nrow(ang_full))

# ---- full-resolution layers ----
# root carries no colour of its own here: it is drawn flat black, so only its
# shape and the mask it leaves behind reach the panel.
skel_cls <- deepcopy(root_bin)           # 0 = background, 1 = root area, 2 = skeleton
skel_cls[skel == 1] <- 2

is_masked <- mask == 1
is_masked[is.na(is_masked)] <- 0

full <- list(root  = root_bin,
             skel  = skel_cls,
             depth = depth_map,
             ismk  = is_masked)
for (nm in c("root", "skel")) full[[nm]][is_masked] <- NA

# ---- rotate 90 degrees clockwise ----
full <- lapply(full, function(r) {
  r <- flip(t(r), direction = "horizontal")
  ext(r) <- c(0, ncol(r), 0, nrow(r))
  r
})

# the angle map is already in the rotated frame: add it now, keep it on skeleton pixels, apply the mask
full$ang <- ang_full
full$ang[full$ismk == 1]  <- NA
full$ang[full$skel != 2]  <- NA

# ---- shrink to output resolution so the plot never skips pixels ----
panel <- list(root = aggregate(full$root, fact = agg, fun = "max",   na.rm = TRUE),
              ang  = aggregate(full$ang,  fact = agg, fun = "modal", na.rm = TRUE))
ang_found <- sort(unique(values(panel$ang, mat = FALSE)))
message("angle classes found: ", paste(ang_found, collapse = " "))

# ---- shared plot settings ----
# The terra:: prefix on the plot calls is not decoration: plot() is an S4
# generic, and anything nearer on the search path than package:terra -- a stray
# object in the workspace, a package attached later -- sends a SpatRaster to
# graphics::plot.default, which drops colNA, range, legend and maxcell and then
# fails in plot.xy. Naming the method makes the script independent of that.
black_ramp <- c("black", "black")                                  # plain black background
leg_bg     <- adjustcolor("black", 0.65)                           # backdrop of the on-panel legend
ang_brk    <- seq(-22.5, 337.5, by = 45)
mar_p      <- c(0.5, 0.5, 4, 0.5)

# ---- inset geometry (panel coordinates: x = 0..W, y = 0..H) ----
W <- ncol(full$root); H <- nrow(full$root)
zw  <- zoom_w * W
iw  <- inset_w * W
pad <- 0.03 * W
ix  <- if (grepl("right", inset_pos)) c(W - pad - iw, W - pad) else c(pad, pad + iw)
iy  <- if (grepl("top",   inset_pos)) c(H - pad - iw, H - pad) else c(pad, pad + iw)
inset_e <- ext(ix[1], ix[2], iy[1], iy[2])

if (identical(zoom_ctr, "densest")) {
  # count skeleton pixels in coarse blocks, sum them over a window-sized neighborhood,
  # and take the best window that does not sit under the inset
  step  <- 10
  dens  <- aggregate((full$skel == 2) * 1, fact = step, fun = "sum", na.rm = TRUE)
  wk    <- 2 * round(zw / step / 2) + 1                       # window size in blocks (odd)
  dens  <- focal(dens, w = wk, fun = "sum", na.rm = TRUE)
  xy    <- xyFromCell(dens, 1:ncell(dens))
  v     <- values(dens, mat = FALSE)
  under <- xy[, 1] > xmin(inset_e) - zw & xy[, 1] < xmax(inset_e) + zw &
           xy[, 2] > ymin(inset_e) - zw & xy[, 2] < ymax(inset_e) + zw
  edge  <- xy[, 1] < zw / 2 | xy[, 1] > W - zw / 2 | xy[, 2] < zw / 2 | xy[, 2] > H - zw / 2
  v[under | edge] <- NA
  best  <- xy[which.max(v), ]
  zoom_ctr <- c(best[1] / W, 1 - best[2] / H)
}
# the square zoom window, centred on zoom_ctr (moved inwards if it would stick out of the panel)
h   <- zw / 2
wx  <- min(max(zoom_ctr[1] * W,       h), W - h)
wy  <- min(max((1 - zoom_ctr[2]) * H, h), H - h)
win <- ext(wx - h, wx + h, wy - h, wy + h)

# the depth range the window covers, for the caption
d_rng <- range(values(crop(full$depth, win), mat = FALSE), na.rm = TRUE)
message(sprintf("zoom window depth range: %.1f-%.1f cm", d_rng[1], d_rng[2]))

# draws the zoom window of a full-resolution layer into the inset box
# new = TRUE starts the inset (black backdrop); new = FALSE adds another layer on top
draw_inset <- function(r, win, new = TRUE, ...) {
  if (!zoom_on) return(invisible(NULL))
  if (new) rect(xmin(inset_e), ymin(inset_e), xmax(inset_e), ymax(inset_e), col = "black", border = NA)
  r_win <- crop(r, win)
  ext(r_win) <- inset_e
  terra::plot(r_win, add = TRUE, legend = FALSE, axes = FALSE, maxcell = ncell(r_win), ...)
  rect(xmin(inset_e), ymin(inset_e), xmax(inset_e), ymax(inset_e), border = zoom_col, lwd = 4)
  rect(xmin(win),     ymin(win),     xmax(win),     ymax(win),     border = zoom_col, lwd = 4)
}

# ---- figure size in pixels ----
# One margin line, in pixels, at the cex the panel is drawn with -- mar is
# counted in lines of text, so the top margin shrinks with the text.
png(tempfile(fileext = ".png"), width = 200, height = 200, res = res_dpi)
line_px <- par("csi") * res_dpi * panel_cex
invisible(dev.off())
panel_w <- ncol(panel$root)
panel_h <- nrow(panel$root)
fig_w   <- round(panel_w + line_px)
fig_h   <- round(panel_h + 4.5 * line_px)

# ---- plot ----
png(out_path, width = fig_w, height = fig_h, res = res_dpi)
par(mar = mar_p, cex = panel_cex)

# root angle (D8 direction, drawn on the skeleton only)
terra::plot(panel$root, col = black_ramp, colNA = mask_col, range = c(0, 1),
     legend = FALSE, axes = FALSE, maxcell = ncell(panel$root), mar = mar_p)
terra::plot(panel$ang, breaks = ang_brk, col = ang_cols, legend = FALSE, axes = FALSE,
     maxcell = ncell(panel$ang), add = TRUE)
draw_inset(full$root, win = win, col = black_ramp, colNA = mask_col, range = c(0, 1))
draw_inset(full$ang, win = win, new = FALSE, breaks = ang_brk, col = ang_cols)
# colour wheel, bottom left: 0 deg = up, clockwise (90 = right, 180 = down = deeper)
# the directions that never occur (315, 0, 45) are left out, i.e. transparent
wc <- c(0.33 * W, 0.33 * W); r_in <- 0.10 * W; r_out <- 0.21 * W
rect(0.02 * W, 0.02 * W, 0.64 * W, 0.64 * W + 0.09 * W, col = leg_bg, border = NA)
for (k in ang_found / 45) {
  tt <- seq((k * 45 - 22.5) * pi / 180, (k * 45 + 22.5) * pi / 180, length.out = 20)
  polygon(wc[1] + c(r_in * sin(tt), rev(r_out * sin(tt))),
          wc[2] + c(r_in * cos(tt), rev(r_out * cos(tt))), col = ang_cols[k + 1], border = NA)
  text(wc[1] + 0.275 * W * sin(k * pi / 4), wc[2] + 0.275 * W * cos(k * pi / 4),
       paste0(k * 45, "\u00b0"), col = "white", cex = legend_cex * 0.6)
}
text(0.05 * W, 0.64 * W + 0.045 * W, "Angle (D8)", adj = c(0, 0.5), col = "white", cex = legend_cex)
mtext("Angle", side = 3, line = 1, cex = title_cex)

dev.off()
