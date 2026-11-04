# Builds data/TurnoverDPC_data.rda from a RootDetector turnover (DPC) image.
#
# Copy the PNG to data-raw/TurnoverDPC.png (or point `src` at it) and run this
# from the package root. The array is stored like the other bundled images:
# rows x columns x 3 layers (R, G, B), values 0-255, alpha dropped.

src <- "data-raw/TurnoverDPC.png"

x <- png::readPNG(src)                 # rows x cols x channels, values 0-1
x <- x[, , 1:3]                        # drop the alpha channel if present
x <- round(x * 255)                    # 0-255, like the other bundled arrays
TurnoverDPC_data <- x

# What the docs in R/data.R have to state: size, and which colors occur.
print(dim(TurnoverDPC_data))
print(table(R = x[, , 1] > 0, G = x[, , 2] > 0, B = x[, , 3] > 0))

save(TurnoverDPC_data, file = "data/TurnoverDPC_data.rda", compress = "xz")
print(file.size("data/TurnoverDPC_data.rda") / 1e6)   # MB
