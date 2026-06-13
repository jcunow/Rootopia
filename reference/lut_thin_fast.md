# Zhang–Suen thinning using lookup table (LUT implementation)

Performs iterative skeletonization of a binary raster using a
lookup-table encoding of 3x3 neighbourhood configurations.

## Usage

``` r
lut_thin_fast(img, max_iter = 200L, verbose = FALSE)
```

## Arguments

- img:

  Binary SpatRaster (single layer only)

- max_iter:

  Maximum number of thinning iterations

- verbose:

  Logical. If TRUE, prints iteration count and pixel reduction

## Value

Binary SpatRaster representing skeletonized image

## Details

The algorithm operates as follows:

1\. Raster values are extracted once into a numeric vector. 2. Each
iteration reconstructs a matrix representation of the image. 3. A
zero-padded border is added around the image. 4. For each pixel equal to
1: - The 3x3 neighbourhood is extracted - A weighted sum (mask encoding)
produces a unique neighbourhood code 5. The code is used as an index
into a 256-entry lookup table: - LUT value 1 or 3 → pixel removed in
first sub-step - LUT value 2 or 3 → pixel removed in second sub-step 6.
Pixels are updated in-place in the vector representation 7. Iteration
stops when no pixels are removed or max_iter is reached

Important properties of implementation: - Raster is only written back
once at the end - All intermediate steps operate on in-memory
matrix/vector data - Neighborhood encoding is performed via explicit
loops over image pixels
