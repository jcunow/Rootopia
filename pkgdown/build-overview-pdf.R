# Render the Overview article to PDF and place it where pkgdown will publish it.
#
# pkgdown copies everything in pkgdown/assets/ to the site root, so the article
# can link to ../Rootopia_Overview.pdf. Run this before pkgdown::build_site().
#
# Needs a LaTeX installation. tinytex::install_tinytex() provides one if you
# have none, and tinytex::reinstall_tinytex() refreshes one that reports a new
# TeX Live release.
#
# This renders with the default pdflatex and the default serif. A specific
# font such as Georgia needs latex_engine: xelatex in the article's
# pdf_document block plus mainfont at its YAML top level -- mainfont is a
# pandoc variable, not a pdf_document() argument. xelatex routes through
# xdvipdfmx, which is where an out-of-date TinyTeX tends to fail; for a serif
# close to Georgia without leaving pdflatex, add to the article's YAML:
#
#   header-includes:
#     - \usepackage{XCharter}

# Run this from the package root.
src <- "vignettes/Rootopia_Overview_vignette.Rmd"
if (!file.exists(src))
  stop("Run this from the package root; ", src, " is not here.", call. = FALSE)

# rmarkdown resolves output_dir against the input file's directory in some
# versions, so hand it an absolute path.
dir.create("pkgdown/assets", showWarnings = FALSE, recursive = TRUE)
out <- normalizePath("pkgdown/assets")

if (!requireNamespace("rmarkdown", quietly = TRUE))
  stop("rmarkdown is required to render the PDF.")
if (!nzchar(Sys.which("pdflatex")) && !requireNamespace("tinytex", quietly = TRUE))
  stop("No LaTeX found. Install one, or run tinytex::install_tinytex().")

rmarkdown::render(
  input         = src,
  output_format = "pdf_document",
  output_file   = "Rootopia_Overview.pdf",
  output_dir    = out,
  quiet         = FALSE
)

message("Wrote ", file.path(out, "Rootopia_Overview.pdf"))
