# Replace the angle-panels chunk of the overview vignette.
#
# Context-free on purpose: it finds the chunk by its name and swaps the whole
# thing, so it does not care what is in there now. Run it once from the package
# root, then delete it.
#
#   source("data-raw/apply-angle-chunk.R")

p <- "vignettes/Rootopia_Overview_vignette.Rmd"
stopifnot(file.exists(p))
L <- readLines(p, warn = FALSE, encoding = "UTF-8")

i <- grep("^```\\{r angle-panels", L)
stopifnot(length(i) == 1L)
j <- i + which(L[(i + 1L):length(L)] == "```")[1]
cat("replacing lines", i, "to", j, "\n")

new <- c(
  '```{r angle-panels, fig.width = 4.8, fig.height = 8.6, dpi = 150, dev = "png", out.width = "72%", fig.alt = "The D8 direction class of every skeleton pixel over the black root field, for the shallow part of the tube with its topmost rows cropped away. One window is outlined and magnified as an inset, and a color wheel gives the directions that occur."}',
  'dd_all <- deep_drive(DepthMap = depth_map, RootMap = skl, return = "all")',
  '',
  '# deep_drive() returns its map already transposed and flipped: turned the same',
  '# way rot() turns every other panel. So it is used as it comes and NOT rotated',
  '# again, while the root field underneath is the rotated one the workflow figure',
  '# already built. Rotating both is how this panel ends up mirrored.',
  'ang_full <- dd_all$angle_map',
  'terra::ext(ang_full) <- terra::ext(full$root)',
  'ang_full <- terra::ifel(full$skel == 2, ang_full, NA)   # center line only',
  '',
  '# Index k + 1 is the class k * 45, so a direction keeps its color from tube to',
  '# tube and two figures compare. 45 and 315 are left out of the wheel: a root',
  '# cannot run up and outward at once, so where they occur at all it is noise.',
  'ang_brk  <- seq(-22.5, 337.5, by = 45)',
  'ang_cols <- c("#FFFFFF", "#FF8A3D", "#FFD60A", "#4D7CFF",',
  '              "#3DDC84", "#E0218A", "#00E5FF", "#B388FF")',
  'ang_skip <- c(45, 315)',
  '',
  '# What is shown. Depth runs down the panel, so the shallow end is the top.',
  '# The first rows are tape and mask rather than soil, and the tube is four',
  '# times longer than it is wide: drawn whole it is a sliver, because terra',
  '# draws at asp = 1 and the page can only give it height. Showing part of the',
  '# length is what buys the scale back.',
  'top_cut   <- 0.05    # drop this much of the length at the shallow end',
  'show_frac <- 0.50    # how much of the length to draw',
  'win_ctr   <- c(0.80, 0.42)   # window center, as fractions of what is shown',
  'win_frac  <- 0.17            # window width, as a fraction of the panel width',
  '',
  'y_hi <- H * (1 - top_cut)',
  'y_lo <- max(0, y_hi - show_frac * H)',
  'e    <- terra::ext(0, W, y_lo, y_hi)',
  '',
  'hw <- win_frac * W / 2',
  'wx <- min(max(win_ctr[1] * W, hw), W - hw)',
  'wy <- min(max(y_hi - win_ctr[2] * (y_hi - y_lo), y_lo + hw), y_hi - hw)',
  'win <- terra::ext(wx - hw, wx + hw, wy - hw, wy + hw)',
  '',
  'a_ang     <- terra::aggregate(ang_full, fact = 4, fun = "modal", na.rm = TRUE)',
  'ang_found <- setdiff(sort(unique(terra::values(terra::crop(a_ang, e),',
  '                                               mat = FALSE))), ang_skip)',
  '',
  '# The same two calls the workflow figure makes: the root raster flat black so',
  '# the field is black and colNA paints only the masked tape, then the',
  '# directions over it as the only color in the panel.',
  'op <- graphics::par(mar = mar_p)',
  'terra::plot(terra::crop(panel$root, e), col = black_ramp, colNA = mask_col,',
  '            range = c(0, 1), legend = FALSE, axes = FALSE, maxcell = Inf,',
  '            mar = mar_p)',
  'terra::plot(terra::crop(a_ang, e), breaks = ang_brk, col = ang_cols,',
  '            legend = FALSE, axes = FALSE, maxcell = Inf, add = TRUE)',
  '',
  '# The inset, bottom right of what is shown. inset_e from the helper chunk is',
  '# in whole-frame coordinates, so the box is built against this crop instead.',
  'iw <- 0.52 * W',
  'ie <- terra::ext(e[2] - 0.04 * W - iw, e[2] - 0.04 * W,',
  '                 e[3] + 0.04 * W, e[3] + 0.04 * W + iw)',
  'graphics::rect(ie[1], ie[3], ie[2], ie[4], col = "black", border = NA)',
  'for (lyr in list(list(r = full$root, args = list(col = black_ramp,',
  '                                                 colNA = mask_col,',
  '                                                 range = c(0, 1))),',
  '                 list(r = ang_full,  args = list(breaks = ang_brk,',
  '                                                 col = ang_cols)))) {',
  '  z <- terra::crop(lyr$r, win); terra::ext(z) <- ie',
  '  do.call(terra::plot, c(list(z, add = TRUE, legend = FALSE, axes = FALSE,',
  '                              maxcell = Inf), lyr$args))',
  '}',
  'graphics::rect(ie[1], ie[3], ie[2], ie[4], border = zoom_col, lwd = zoom_lwd)',
  'graphics::rect(win[1], win[3], win[2], win[4], border = zoom_col,',
  '               lwd = zoom_lwd)',
  '',
  '# A wheel rather than a list of numbers: 0 degrees is up and the angle runs',
  '# clockwise, so 180 is straight down the depth axis and 90 and 270 run around',
  '# the circumference. Sectors with no wedge are directions that do not occur.',
  'sz <- 0.20 * W',
  'wc <- c(1.8 * sz, e[4] - 1.8 * sz)',
  'graphics::rect(wc[1] - 1.7 * sz, wc[2] - 2.3 * sz,',
  '               wc[1] + 1.7 * sz, wc[2] + 1.7 * sz, col = leg_bg, border = NA)',
  'for (q in ang_found / 45) {',
  '  tt <- seq((q * 45 - 22.5) * pi / 180, (q * 45 + 22.5) * pi / 180,',
  '            length.out = 24)',
  '  graphics::polygon(wc[1] + c(0.45 * sz * sin(tt), rev(sz * sin(tt))),',
  '                    wc[2] + c(0.45 * sz * cos(tt), rev(sz * cos(tt))),',
  '                    col = ang_cols[q + 1], border = NA)',
  '  graphics::text(wc[1] + 1.35 * sz * sin(q * pi / 4),',
  '                 wc[2] + 1.35 * sz * cos(q * pi / 4),',
  '                 paste0(q * 45, "°"), col = "white", cex = 0.85)',
  '}',
  'graphics::text(wc[1], wc[2] - 1.95 * sz, "Angle (D8)", adj = c(0.5, 0.5),',
  '               col = "white", cex = 0.95)',
  'graphics::mtext("Root direction, shallow end", side = 3, line = 0.4)',
  'graphics::par(op)',
  '```'
)

writeLines(c(if (i > 1L) L[seq_len(i - 1L)], new,
             if (j < length(L)) L[(j + 1L):length(L)]), p,
           useBytes = TRUE)
cat("done:", length(new), "lines written into the chunk\n")
