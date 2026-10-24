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
