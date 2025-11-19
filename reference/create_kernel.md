# Create a kernel with specified shape for morphological operations

Create a kernel with specified shape for morphological operations

## Usage

``` r
create_kernel(shape = "disk", size = 3)
```

## Arguments

- shape:

  Shape of the kernel: "square", "diamond", or "disk"

- size:

  Size of the kernel (odd integer)

## Value

A kernel as an imager cimg object
