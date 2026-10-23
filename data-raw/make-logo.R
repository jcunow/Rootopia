# Regenerate the logo files the package references, from the master artwork.
#
# Two derivatives are needed, and neither can be dropped:
#   man/figures/logo.png  the README points at this path, and pkgdown finds it
#                         there on its own with no _pkgdown.yml entry
#   vignettes/logo.png    a vignette is built in vignettes/ and installed into
#                         doc/, so it cannot reach man/figures/
#
# The master lives here and nowhere else. data-raw/ is in .Rbuildignore, so the
# full-resolution file stays in git without being shipped in the tarball.

master <- list.files("data-raw/logo", full.names = TRUE,
                     pattern = "\\.(png|svg|jpe?g)$", ignore.case = TRUE)
if (length(master) != 1L)
  stop("Put exactly one image in data-raw/logo/; found ", length(master), ".",
       call. = FALSE)

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)

# 240 px wide is what pkgdown's navbar uses and keeps the tarball small.
if (requireNamespace("magick", quietly = TRUE)) {
  img <- magick::image_read(master)
  magick::image_write(magick::image_resize(img, "240"), "man/figures/logo.png")
} else {
  if (tolower(tools::file_ext(master)) != "png")
    stop("magick is needed to convert ", basename(master), " to PNG.",
         call. = FALSE)
  message("magick not installed; copying the master at its own size.")
  file.copy(master, "man/figures/logo.png", overwrite = TRUE)
}

file.copy("man/figures/logo.png", "vignettes/logo.png", overwrite = TRUE)
message("wrote man/figures/logo.png and vignettes/logo.png from ", master)
