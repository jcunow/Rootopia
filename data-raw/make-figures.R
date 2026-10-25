# Figures the Overview article cannot compute from the bundled data, written
# from full-resolution originals instead.
#
# The example data in data/ is downsampled and re-encoded to keep the package
# installable, which is fine for showing what a function returns and not fine
# for comparing two skeletons pixel by pixel. Point the three paths below at
# one scan's originals and run this from the package root; the article picks the
# output up automatically and falls back to the bundled data when it is absent.
#
# The panel figure at the end of Part I comes from overview_figure.R, not here.

seg_path <- ""   # the segmentation, e.g. ..._T067.png.segmentation.png
skl_path <- ""   # the segmenter's own skeleton. Left empty, the script looks
                 # for one beside the segmentation, which is where RootDetector
                 # puts it.
dpi      <- 300  # the resolution these scans were taken at

library(Rootopia)
library(terra)

# A device left open by a failed run swallows every later plot in the session,
# including a vignette's, so start clean and close on the way out whatever
# happens.
grDevices::graphics.off()

to_png <- function(file, width, height, expr) {
  grDevices::png(file, width = width, height = height, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
  invisible(file)
}

out <- "vignettes/figures"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# The skeleton comparison needs real scans; the synthetics below do not, so an
# empty path skips that one figure rather than stopping the script.
if (nzchar(seg_path) && !nzchar(skl_path)) {
  near <- list.files(dirname(seg_path), pattern = "skel", ignore.case = TRUE,
                     full.names = TRUE)
  if (length(near) == 1L) {
    skl_path <- near
    message("using the skeleton found beside the segmentation: ",
            basename(skl_path))
  } else if (length(near) > 1L) {
    stop("Several skeleton-looking files beside the segmentation; set ",
         "skl_path explicitly:\n  ", paste(near, collapse = "\n  "),
         call. = FALSE)
  }
}

have_scans <- nzchar(seg_path) && nzchar(skl_path)
if (have_scans)
  for (f in c(seg_path, skl_path))
    if (!file.exists(f)) stop("Not found: ", f, call. = FALSE)

if (have_scans) {
  mask     <- load_flexible_image(seg_path, output_format = "spatrast",
                                  scale = "binary", select_layer = 2)
  supplied <- load_flexible_image(skl_path, output_format = "spatrast",
                                  scale = "binary", select_layer = 2)
  thinned  <- skeletonize_image(mask, verbose = FALSE)

  ext(supplied) <- ext(mask)
  ext(thinned)  <- ext(mask)

  # 0 background, 1 both, 2 supplied only, 3 thinned only
  cmp <- mask * 0
  cmp[supplied == 1] <- 2
  cmp[thinned  == 1] <- 3
  cmp[(supplied == 1) & (thinned == 1)] <- 1

  to_png(file.path(out, "skeleton_compare.png"), 2400, 1500, {
  zoom_plot(cmp, frac = 0.08, center = c(0.8, 0.35), maxcell = Inf,
            overview = FALSE,
            breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5),
            col = c("white", "grey30", "#2455A4", "#E0218A"),
            main = "")
  legend("bottomleft", bty = "n", border = NA,
         fill = c("grey30", "#2455A4", "#E0218A"),
         legend = c("both skeletons", "only the supplied one",
                    "only the thinned one"))
  })

  cat(sprintf("length (cm)  supplied %.1f  thinned %.1f\n",
              root_length(supplied, unit = "cm", dpi = dpi, show_messages = FALSE),
              root_length(thinned,  unit = "cm", dpi = dpi, show_messages = FALSE)))
  message("wrote ", file.path(out, "skeleton_compare.png"))
} else {
  message("seg_path / skl_path are empty: skipping skeleton_compare.png.")
}

# ---------------------------------------------------------------------------
# Synthetic examples. These need no scans, so they run whatever the paths
# above are set to.
# ---------------------------------------------------------------------------

bres <- function(r0, c0, r1, c1) {
  dr <- abs(r1 - r0); dc <- abs(c1 - c0)
  sr <- sign(r1 - r0); sc <- sign(c1 - c0)
  err <- dr - dc; out <- matrix(NA_integer_, dr + dc + 2L, 2L); k <- 0L
  repeat {
    k <- k + 1L; out[k, ] <- c(r0, c0)
    if (r0 == r1 && c0 == c1) break
    e2 <- 2L * err
    if (e2 > -dc) { err <- err - dc; r0 <- r0 + sr }
    if (e2 <  dr) { err <- err + dr; c0 <- c0 + sc }
  }
  out[seq_len(k), , drop = FALSE]
}

paint <- function(size, strokes) {
  m <- matrix(0L, size, size)
  for (st in strokes) m[bres(st[1], st[2], st[3], st[4])] <- 1L
  m
}

# Colour each traced segment separately, so what the tracer calls one root is
# visible. Segments are rasterised through the package's own overlay writer.
seg_raster <- function(et, mask) {
  f <- tempfile(fileext = ".png")
  render_order_overlay(attr(et, "segments"), seq_len(nrow(et)),
                       attr(et, "dims"), f, mask = mask,
                       value_label = "segment")
  a <- png::readPNG(f)
  if (length(dim(a)) == 2L) a <- array(a, c(dim(a), 3L))
  r <- terra::rast(a[, , 1:3, drop = FALSE] * 255)
  terra::RGB(r) <- 1:3
  r
}

# --- 1. what the crossing rule does, case by case -------------------------
# Two rules act at a node, and they are different rules. A node with exactly
# four arms goes to resolve_crossings(), which joins the pair whose directions
# are most nearly opposite. A node with three arms is a junction, and there
# assign_root_order() picks which child continues the parent by straightness
# AND similar diameter together -- diam_weight sets the balance. So the strokes
# below have width: on one-pixel lines the diameter half of that rule has
# nothing to work with.

stamp <- function(m, px, radius) {
  rr <- seq(-radius, radius)
  off <- as.matrix(expand.grid(dr = rr, dc = rr))
  off <- off[rowSums(off^2) <= radius^2, , drop = FALSE]
  for (k in seq_len(nrow(px))) {
    p <- cbind(px[k, 1] + off[, 1], px[k, 2] + off[, 2])
    p <- p[p[, 1] >= 1 & p[, 1] <= nrow(m) & p[, 2] >= 1 & p[, 2] <= ncol(m), ,
           drop = FALSE]
    m[p] <- 1L
  }
  m
}

# Each case is a pair: the one-pixel centre lines the tracer follows, and the
# filled roots it reads diameters from. Passing the filled image as `skel`
# traces almost nothing, which is what the pipeline warns about.
paint <- function(size, strokes) {
  thin <- matrix(0L, size, size)
  thick <- matrix(0L, size, size)
  for (st in strokes) {
    px <- bres(st[1], st[2], st[3], st[4])
    thin[px] <- 1L
    thick <- stamp(thick, px, st[5])
  }
  list(skel = thin, mask = thick)
}

# One colour per traced segment, read straight off the segment list: no PNG
# round-trip, so no duplicated labels and no mismatch between the drawn frame
# and the frame the tracer cropped to.
seg_map <- function(et) {
  d <- attr(et, "dims")
  m <- matrix(NA_integer_, d[1], d[2])
  sg <- attr(et, "segments")
  for (i in seq_along(sg)) m[sg[[i]]$coords] <- i
  r <- terra::rast(m); terra::ext(r) <- c(0, d[2], 0, d[1])
  r
}

# The crossing and junction rules have their own script now,
# data-raw/node-gallery.R, which covers three to six arms and prints the
# scores each decision turns on. What stays here is the order-scheme figure.

# --- 2. one fanning tree under all four order schemes ---------------------
# A trunk with laterals leaving at six heights, two of them carrying their own
# pair of children, and four stubs. Nine branch points, eleven tips, three
# levels deep: enough that the schemes cannot agree by accident.
br <- paint(161, list(
  c(150, 81,  40,  81, 4),   # trunk, thickest
  c(136, 81, 110,  55, 3),   # first order
  c(128, 81, 102, 107, 3),
  c(110, 55,  92,  37, 2),   # second order, off a 3-way junction
  c(110, 55,  98,  67, 2),
  c(102,107,  84, 125, 2),
  c(102,107,  90,  95, 2),
  c(120, 81, 114,  87, 1),   # stubs, at four heights
  c( 96, 81,  90,  75, 1),
  c( 70, 81,  66,  85, 1),
  c( 84, 81,  62,  59, 3),
  c( 58, 81,  40,  99, 3)))

# Drawn trunk-down; turn it the other way up.
br$skel <- br$skel[nrow(br$skel):1, ncol(br$skel):1]
br$mask <- br$mask[nrow(br$mask):1, ncol(br$mask):1]

br_r <- terra::rast(br$mask); terra::ext(br_r) <- c(0, 161, 0, 161)
schemes <- c("order_from_tip", "strahler_order",
             "order_from_base", "root_max_order_from_tip")

max_ord  <- max(vapply(schemes, function(sc)
  max(branch_order_map(br$skel, br$mask, order = sc, unit = "px",
                       template = br_r, verbose = FALSE)$edges[[sc]],
      na.rm = TRUE), numeric(1)))
ord_cols <- hcl.colors(max_ord, "Zissou1")

to_png(file.path(out, "synthetic_orders.png"), 1500, 1500, {
op <- par(mfrow = c(2, 2), mar = c(0, 0, 2.6, 0))
for (sc in schemes) {
  res <- branch_order_map(br$skel, br$mask, order = sc, unit = "px",
                          template = br_r, verbose = FALSE)
  ord <- sort(unique(res$edges[[sc]]))
  # One colour per order value, shared by every panel, so order 2 is the same
  # colour whichever scheme produced it.
  terra::plot(res$class_map, breaks = c(0.5, seq_len(max_ord) + 0.5),
              col = ord_cols, legend = FALSE, axes = FALSE, maxcell = Inf,
              main = sc)
  mtext(paste("orders present:", paste(ord, collapse = "  ")),
        side = 1, line = -1.2, cex = 0.75)
  if (sc == schemes[1])
    legend("topleft", bty = "n", fill = ord_cols, border = NA, cex = 0.8,
           legend = paste("order", seq_len(max_ord)))
  cat(sprintf("%-24s %s\n", sc, paste(ord, collapse = " ")))
}
par(op)
})

message("wrote ", file.path(out, "synthetic_crossings.png"), " and ",
        file.path(out, "synthetic_orders.png"))
