# Smooth the edges of a root binary image

Smooth the edges of a root binary image

## Usage

``` r
smooth_root_edges(img, kernel_shape = "disk", kernel_size = 3, iterations = 1)
```

## Arguments

- img:

  A binary image of root systems (imager cimg object)

- kernel_shape:

  Shape of the kernel: "square", "diamond", or "disk"

- kernel_size:

  Size of the kernel for morphological operations (odd integer). Use
  larger kernels for higher image resolution

- iterations:

  Number of iterations for the smoothing process

## Value

A smoothed binary image

## Examples

``` r

data("seg_Oulanka2023_Session01_T067")
img <- seg_Oulanka2023_Session01_T067
# Try different kernel shapes
smoothed_square <- smooth_root_edges(img, kernel_shape = "square", kernel_size = 3)
smoothed_diamond <- smooth_root_edges(img, kernel_shape = "diamond", kernel_size = 3)
smoothed_disk <- smooth_root_edges(img, kernel_shape = "disk", kernel_size = 3)
plot(smoothed_disk)

```
