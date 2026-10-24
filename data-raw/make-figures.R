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

seg_path <- ""   # the segmenter's mask,     e.g. ..._T067.png.segmentation.png
skl_path <- ""   # the segmenter's skeleton, if it writes one
dpi      <- 300  # the resolution these scans were taken at

if (!nzchar(seg_path) || !nzchar(skl_path))
  stop("Fill in seg_path and skl_path at the top of this script.",
       call. = FALSE)
for (f in c(seg_path, skl_path))
  if (!file.exists(f)) stop("Not found: ", f, call. = FALSE)

library(Rootopia)
library(terra)

out <- "vignettes/figures"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

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

png(file.path(out, "skeleton_compare.png"),
    width = 2400, height = 1500, res = 200)
zoom_plot(cmp, frac = 0.08, center = c(0.8, 0.35), maxcell = Inf,
          overview = FALSE,
          breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5),
          col = c("white", "grey30", "#2455A4", "#E0218A"),
          main = "")
legend("bottomleft", bty = "n", border = NA,
       fill = c("grey30", "#2455A4", "#E0218A"),
       legend = c("both skeletons", "only the supplied one",
                  "only the thinned one"))
dev.off()

cat(sprintf("length (cm)  supplied %.1f  thinned %.1f\n",
            root_length(supplied, unit = "cm", dpi = dpi, show_messages = FALSE),
            root_length(thinned,  unit = "cm", dpi = dpi, show_messages = FALSE)))
message("wrote ", file.path(out, "skeleton_compare.png"))


# ---------------------------------------------------------------------------
# Synthetic examples. These need no scans, so they run whatever the paths
# above are set to. Check both figures before they go into the article.
# ---------------------------------------------------------------------------

draw <- function(m, title) {
  r <- terra::rast(m)
  terra::ext(r) <- c(0, ncol(m), 0, nrow(m))
  terra::plot(r, col = c("white", "grey20"), legend = FALSE, axes = FALSE,
              maxcell = Inf, main = title)
  invisible(r)
}

# --- 1. four node situations, one per panel -------------------------------
# Each is drawn so the arms meet at exactly one pixel, which is what the
# tracer's crossing rule needs. A junction has three arms, a crossing four.
node <- list()

m <- matrix(0L, 61, 61); m[31, 6:56] <- 1L; m[6:31, 31] <- 1L
node[["T junction (3 arms)"]] <- m

m <- matrix(0L, 61, 61); m[31, 6:56] <- 1L; m[6:56, 31] <- 1L
node[["X crossing (4 arms)"]] <- m

m <- matrix(0L, 61, 61); m[31, 6:56] <- 1L
for (k in -24:24) m[31 - k, 31 + k] <- 1L
node[["oblique crossing"]] <- m

m <- matrix(0L, 61, 61); m[31, 6:56] <- 1L; m[6:31, 26] <- 1L; m[31:56, 36] <- 1L
node[["two junctions, no crossing"]] <- m

png(file.path(out, "synthetic_nodes.png"), width = 2000, height = 650, res = 150)
op <- par(mfrow = c(1, 4), mar = c(0, 0, 2.4, 0))
for (nm in names(node)) {
  draw(node[[nm]], nm)
  et <- root_graph_pipeline(node[[nm]], node[[nm]], verbose = FALSE)
  mtext(sprintf("%d segment(s), %d tip(s), %d branch point(s)",
                nrow(et), sum(et$n_tips), sum(et$n_branch_points)),
        side = 1, line = -1.2, cex = 0.7)
}
par(op); dev.off()

for (nm in names(node)) {
  cat("\n---", nm, "---\n")
  print(root_graph_pipeline(node[[nm]], node[[nm]], verbose = TRUE)[
    , c("from", "to", "length", "n_tips", "n_branch_points")])
}

# --- 2. one branching system under all four order schemes -----------------
# A main axis, three first-order laterals, and two second-order laterals on
# one of them, so the schemes cannot agree by accident.
br <- matrix(0L, 121, 121)
br[61, 11:111] <- 1L                       # main axis
br[31:61, 36] <- 1L                        # lateral 1
br[61:91, 61] <- 1L                        # lateral 2
br[31:61, 86] <- 1L                        # lateral 3
br[31, 36:51] <- 1L                        # sub-lateral on lateral 1
br[16:31, 46] <- 1L                        # sub-sub-lateral

schemes <- c("order_from_tip", "strahler_order",
             "order_from_base", "root_max_order_from_tip")

png(file.path(out, "synthetic_orders.png"), width = 2000, height = 620, res = 150)
op <- par(mfrow = c(1, 4), mar = c(0, 0, 2.4, 0))
for (sc in schemes) {
  res <- branch_order_map(br, br, order = sc, unit = "px", verbose = FALSE)
  terra::plot(res$class_map, col = hcl.colors(6, "Zissou1"), legend = FALSE,
              axes = FALSE, maxcell = Inf, main = sc)
  mtext(paste("orders:", paste(sort(unique(res$edges[[sc]])), collapse = " ")),
        side = 1, line = -1.2, cex = 0.7)
}
par(op); dev.off()

message("wrote ", file.path(out, "synthetic_nodes.png"), " and ",
        file.path(out, "synthetic_orders.png"))
