# Minirhizotron Scans Analysis with RootScanR in R

## Minirhizotron Scans Analysis with RootScanR

### Introduction

This vignette demonstrates how to analyze minirhizotron scans using the
**RootScanR** package. Minirhizotrons are transparent tubes inserted
into soil to observe root growth via scanning or imaging technology.
RootScanR is designed for comprehensive analysis of these root system
images, providing tools for image preparation, depth mapping, and
extraction of meaningful root traits. To extract root traits, segmenting
the image outside R is needed (see step 2).

### Key Functionalities

A typical workflow to analyze minirhizotron scans could look something
like the following. Step 2 and step 3 depend

1.  **Image Stitching**: Combine sequential minirhizotron images to
    create complete tube scans
2.  **AI-Based Image Segmentation**: Work with images segmented by
    dedicated root detection software
3.  **In-situ Calibration**: Needed to accurately map image pixels to
    soil depth. An alternative is shown here but requires special
    attention
4.  **Depthmap Creation**: Generate depth information for each pixel in
    the root scans
5.  **Root Trait Extraction**: Calculate important root metrics (length,
    density, branching, etc.)
6.  **Spatial Root Features** Extracting spatial-explicit root scape
    features e.g., entropy, root instance count, fractal dimension
7.  **Soil Characterization**: Analyze soil/peat properties from scan
    images

### Installation

``` r
# Install the package from GitHub
# install.packages("remotes")
# remotes::install_github("jcunow/RootScanR")

# Load the package
library(RootScanR)
library(tidyverse)  # For data manipulation and visualization
library(terra)      # For raster image handling
```

### Example Workflow

#### 1. Image Stitching

**Purpose**: Combine multiple sequential images from longer
minirhizotron tubes into a single composite image to enable proper depth
alignment and remove overlap between images.

``` r
# Load example image files (replace with your actual file paths)
image_files <- c("path/to/image1.png", "path/to/image2.png", "path/to/image3.png")

# Stitch images together (removes redundant overlap)
stitched_image <- RootScanR::stitch_sequential_images(image_files, 
                                                      overlap_px = 200, 
                                                      method = "crosscorr",
                                                      side1 = "bottom")  


# Visualize the stitched image
plot(stitched_image, main = "Stitched Minirhizotron Image")
```

#### 2. Image Segmentation Recommendations

**Purpose**: Separate roots from background in the minirhizotron images.

RootScanR works with segmented images but does not perform segmentation
itself. We recommend using specialized tools:

1.  **RootDetector**
    ([GitHub](https://github.com/ExPlEcoGreifswald/RootDetector))
    - Returns segmented and skeletonized images with different
      information channels
    - Can distinguish between tape, roots, and background
    - Supports root tracking for production, decay, and turnover
      analysis
    - Reference: Peters et al. (2023), *Scientific Reports* **13**(1)
2.  **RootPainter** ([GitHub](https://github.com/Abe404/root_painter))
    - Returns segmented images
    - Uses interactive deep learning for improved segmentation
    - Reference: Smith et al. (2022), *New Phytologist* **236**(2)

**Tip**: For images with high root density, consider splitting your
stitched images before segmentation and rejoining afterwards. **Tip**:
Models can be sensitive to resolution. Try

#### 3. Loading Images

**Purpose**: Import your segmented and original RGB images for analysis.

``` r
# Load a segmented image (replace with your file path). The function accepts multiple different input formats.
segmented_image <- load_flexible_image("path/to/segmented_image.tif", output_format = "spatrast", normalize = F, select.layer = NULL, binarize = FALSE)

# Load corresponding RGB image for visual reference
rgb_image <- load_flexible_image("path/to/rgb_image.tif", output_format = "spatrast")

# Display the images
plot(rgb_image, main = "Original RGB Image")
plot(segmented_image, main = "Segmented Image")
```

#### 4. Tube Calibration

**Purpose**: Establish the relationship between image pixels and actual
soil depth.

In-situ calibration is highly recommended for accurate depth
attribution. If calibration data is unavailable, RootScanR provides
estimation functions using tape cover as proxy:

``` r
# Estimate soil surface position based on tape cover. Here, 0.5 cm of tape beyond the true soil surface is assumed to block out light.
soil_surface <- estimate_soil_surface(rgb_image, dpi = 300, tape.overlap = 0.5)
print(paste0("Estimated soil surface begins at column: ", soil_surface$soil0))

# Estimate rotation center. This only works for asymmetric taping around the tube with more on the uside. 
rotation_center <- estimate_rotation_center(rgb_image)
print(paste0("Estimated rotation center at row: ", rotation_center))
```

#### 5. Depthmap Creation

**Purpose**: Create a map that assigns soil depth values to each pixel
in the image.

``` r
# Create a mask to exclude tape from soil analysis
# Assuming RootDetector format where red channel includes tape
tape_mask <- (segmented_image[[1]] - segmented_image[[2]]) / 255
tape_mask <- terra::t(tape_mask)

# Calculate normalized rotation center position
center_offset <- rotation_center / dim(segmented_image)[1]

# Create the depthmap
depth_map <- create_depthmap(
  image = segmented_image, 
  sinusoidal = TRUE,    # Account for tube curvature
  mask = tape_mask,
  soil_start = soil_surface,
  center_offset = center_offset
)

# Visualize the depthmap
plot(depth_map, main = "Depth Map")
```

#### 6. Depth Binning

**Purpose**: Group depths into discrete intervals for statistical
analysis and visualization.

``` r
# Create depth bins (5cm intervals)
binned_map <- bin_depths(depth_map, interval = 5, method = "rounding")

# Correct extent for plotting
terra::ext(binned_map) <- c(0, dim(binned_map)[2], 0, dim(binned_map)[1])

# Visualize binned depthmap
transposed_map <- terra::t(binned_map)
terra::plot(transposed_map, main = "Binned Depth Map (5cm intervals)")
```

#### 7. Root Skeletonization

**Purpose**: Create a single-pixel-wide representation of root structure
for morphological analysis.

``` r
# Create root skeleton using Medial Axis Transform (MAT) method
# Select layer 2 which typically contains root information in RootDetector output
skeleton <- skeletonize_image(
  segmented_image, 
  method = "MAT", 
  select.layer = 2,
  verbose = TRUE
)

# Visualize skeleton
plot(skeleton, main = "Root Skeleton")
```

#### 8. Root Trait Extraction

**Purpose**: Calculate quantitative root system traits for ecological
analysis.

``` r
# Detect endpoints and branching points
skeleton_points <- detect_skeleton_points(skeleton)
root_tips <- px.sum(skeleton_points$endpoints)
branching_points <- px.sum(skeleton_points$branching_points)
print(paste0("Number of root tips: ", root_tips))
print(paste0("Number of branching points: ", branching_points))

# Calculate root length using Kimura's method
root_length <- RootLength(skeleton)
print(paste0("Total root length (cm): ", root_length))

# Calculate root density
# Create binary mask of non-root area
void_area <- (prepare_binary_image(segmented_image) - 1) * -1
root_pixels <- px.sum(segmented_image, layer = 2)
void_pixels <- px.sum(void_area, layer = 2)
root_density <- root_pixels / (root_pixels + void_pixels)
print(paste0("Root density: ", round(root_density, 3)))

# Calculate maximum rooting depth ("deep drive")
deep_drive <- deep_drive(
  RootMap = segmented_image, 
  DepthMap = depth_map
)
print(paste0("Maximum rooting depth (cm): ", deep_drive))

# Calculate landscape metrics for root architecture
root_metrics <- root_scape_metrics(
  img = segmented_image, 
  indexD = 5, 
  metrics = c("lsm_c_ca", "lsm_l_ent", "lsm_c_pd", "lsm_c_enn_mn")
)
print("Root architecture metrics:")
print(root_metrics)
```

#### 9. Soil/Peat Characterization

**Purpose**: Analyze soil properties visible in the minirhizotron
images.

``` r
# Analyze overall tube color
tube_color <- tube_coloration(rgb_image)
print("Overall tube color metrics:")
print(tube_color)

# Create buffer around roots to analyze soil without root & rhizosphere influence
root_buffer <- create_root_buffer(segmented_image, width = 2, halo.only = FALSE)

# Create soil-only image by removing root areas
soil_image <- rgb_image - root_buffer

# Analyze soil color
soil_color <- tube_color(soil_image)
print("Soil color metrics:")
print(soil_color)

# Extract soil texture information using GLCM (Gray-Level Co-occurrence Matrix)
gray_image <- rgb2gray(soil_image)
soil_texture <- analyze_soil_texture(
  gray_image,
  gray_levels = 7,
  window_size = c(3, 3),
  statistics = c("second_moment", "homogeneity"),
  directions = list(c(0, 1), c(1, 1), c(1, 0), c(1, -1))
)
print("Soil texture metrics:")
print(soil_texture)
```

#### 10. Depth-Based Analysis

**Purpose**: Analyze how root traits change with soil depth.

``` r
## create analytical zones; indexD indicates which depth to look at.
root_zone = zoning(img = root_map, 
                   depth = 10,
                   mode = "depth",
                   depth_map =  binning(depthmap = depth_map, 5) )

# go through each slice
bin_width = 2
rl.df  = data.frame(depth = seq(0,30,bin_width), rl = NA,rl.density = NA, rootpx = NA, voidpx = NA)
root_map =  load_flexible_image(root_map,binarize = TRUE, output_format = "spatrast")
for(i in rl.df$depth){
  # analytical unit
  root_zone = zoning(img = root_map , 
                     mode = "depth",
                     depth = i,
                     depth_map =  binning(depthmap = depth_map, bin_width) )
  # root length
  try(rl.df$rl[rl.df$depth == i] <- as.numeric(root_length(root_zone,dpi = 150)))
  
    # available root pixels
  try(rl.df$rootpx[rl.df$depth == i] <-  as.numeric(px.sum(root_zone)))
    # available soil & tape pixels (inverse roots)
    voidpx <- abs(root_zone -1)
  try(rl.df$voidpx[rl.df$depth == i] <- as.numeric(px.sum(voidpx)))
    # calculate rootlength density
  try(rl.df$rl.density[rl.df$depth == i] <- as.numeric(rl.df$rl[rl.df$depth == i] / (rl.df$rootpx[rl.df$depth == i] + rl.df$voidpx[rl.df$depth == i])))
  
}


# Visualize results
ggplot(rl.df, aes(x = depth, y = scale(rl.density))) +
  geom_line() +
  geom_point() +
  geom_line(aes(y = scale(rootpx)),color = "steelblue4") +
  geom_point(aes(y = scale(rootpx)),color="steelblue4") +
    geom_line(aes(y = scale(voidpx+rootpx)/10),color="darkorange") +
  geom_point(aes(y = scale(voidpx+rootpx)/10),color="darkorange") +
  theme_minimal() +
  labs(
    title = "Root Length by Soil Depth",
    x = "Soil Depth (cm)",
    y = "standard deviation"
  )
```

### Conclusion

This vignette demonstrates the core workflows for analyzing
minirhizotron images with RootScanR. For more detailed information about
specific functions and additional features, please refer to the function
documentation and the package GitHub repository at
<https://github.com/jcunow/RootScanR>.
