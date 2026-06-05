# Analyzing Root Traits from Flatbed Scans with RootScanR in R

## Analyzing Root Systems from Flatbed Scans

### Introduction

###### Flatbed Root Scan Analysis: Trait Extraction and Quantification

Flatbed scanning allows for the capture of root morphology in a 2D
plane. RootScanR provides a suite of tools for processing flatbed images
of roots, including trait extraction, root length measurement, and
morphology analysis. This vignette will guide you through extracting
root traits from flatbed scans to help quantify root systems for various
ecological studies.

### Installation

``` r

# Install the package from GitHub
# install.packages("remotes")
# remotes::install_github("jcunow/RootScanR")

# Load the package
library(RootScanR)
library(tidyverse)  # For data manipulation and visualization
```

### Key Functionalities for Flatbed Scan Analysis

When working with flatbed scans, RootScanR offers these primary
functions:

1.  **Image Pre-processing**: Convert and prepare images for analysis
2.  **Skeletonization**: Create single-pixel-wide representation of root
    structures
3.  **Root Trait Extraction**: Calculate key morphological traits
    (length, diameter, density)
4.  **Architectural Analysis**: Assess root system complexity and
    structure

### Example Workflow

#### 1. Loading Images

**Purpose**: Import flatbed scan images for analysis.

``` r

# Load a flatbed scan image (replace with your actual file path)
flatbed_image <- load_flexible_image(
  "path/to/flatbed_scan.tiff", 
  verbose = TRUE
)

# Display the image
plot(flatbed_image, main = "Original Flatbed Scan")
```

#### 2. Image Pre-processing

**Purpose**: Prepare the image for analysis by converting to binary
format if needed.

``` r
# Convert to binary image if not already binary
# (skip this step if your image is already properly segmented)
binary_image <- RootScanR::image_threshold
  flatbed_image,
  threshold = 0.4,  # or specify a threshold value
  select.layer = 2,  
  method = "global",
  binary_01 = TRUE
)

# Display binary image
plot(binary_image, main = "Binary Root Image")
```

#### 3. Skeletonization

**Purpose**: Create a single-pixel-wide representation of the root
system for morphological analysis.

``` r

# Skeletonize the root image
skeleton <- skeletonize_roots(
  binary_image,
  method = "MAT",           # Medial Axis Transform
  output_format = "spatrast" # Return as a SpatRaster object
)

# Visualize the skeleton
plot(skeleton, main = "Root Skeleton")
```

#### 4. Root Trait Extraction

**Purpose**: Calculate quantitative metrics that describe root system
morphology.

##### 4.1 Branching Structure

**Purpose**: Identify key structural elements of the root system.

``` r

# Detect root tips and branching points
skeleton_features <- detect_skeleton_points(skeleton)

# Count root tips (endpoints)
root_tips <- count_pixels(skeleton_features$endpoints)
print(paste0("Number of root tips: ", root_tips))

# Count branching points
branching_points <- count_pixels(skeleton_features$branching_points)
print(paste0("Number of branching points: ", branching_points))

# Calculate branching frequency
branching_frequency <- branching_points / calculate_root_length(skeleton) * 100
print(paste0("Branching frequency (per 100 cm): ", round(branching_frequency, 2)))
```

##### 4.2 Root Length

**Purpose**: Calculate the total length of all roots in the image.

``` r

# Calculate root length using Kimura's method
root_length <- root_length(skeleton)
print(paste0("Total root length (cm): ", round(root_length, 2)))
```

##### 4.3 Root Coverage and Density Conversion

**Purpose**: Assess how much of the scan area is occupied by roots.

``` r

# Create background (non-root) mask
background <- (binary_image - 1) * -1

# Count root and background pixels
root_pixels <- count_pixels(binary_image)
background_pixels <- count_pixels(background)

# Calculate root density (proportion of area covered by roots)
root_density <- root_pixels / (root_pixels + background_pixels)
print(paste0("Root density: ", round(root_density, 4)))

# Calculate root length density (proportion of area covered by roots)
rootlength_density <- root_length / (root_pixels + background_pixels)
print(paste0("Root length density: ", round(rootlength_density, 4)))
```

##### 4.4 Root Diameter

**Purpose**: Estimate the thickness of roots at different positions.

``` r

# Calculate root diameters
diameter_map <- root_diameter(
  root_map, 
  select.layer = 2,
  diagnostics = T,
  skeleton_method = "MAT"  # Distance transform method
)

# Get summary statistics of root diameters
diameter_stats = data.frame(diameter_map$diameters) %>%
  mutate(diameters =lyr.1) %>%
  summarise(
    mean = mean( diameters, na.rm = TRUE),
    sd = sd(diameters, na.rm = TRUE),
    median = median(diameters)
  )

# or
library(psych)
diameter_stats <- describe(diameter_map$diameters)
print("Root diameter statistics (mm):")
print(diameter_stats)

# Visualize diameter distribution
h <- hist(
  diameter_map$diameters, 
  main = "Root Diameter Distribution",
  xlab = "Diameter (mm)", 
  plot = FALSE
)


plot(h$mids, h$counts, type = "h", lwd = 10, col = "lightblue",
     log = "y", xlab = "Log-scaled Data", ylab = "Counts",
     main = "Histogram with Log-Scaled Y-Axis")
```

#### 5. Visualization and Export

**Purpose**: Create visual representations of results and export data
for further analysis.

``` r

library(grid)
# overlay plotting
g1 = ggplotGrob(flatbed_image)
g2 = ggplotGrob(skeleton)

grid.draw(g1)
grid.draw(g2)

# Export results to CSV
export_results <- data.frame(
  sample_id = "example_scan",
  total_length_cm = root_length,
  root_tips = root_tips,
  branching_points = branching_points,
  branching_frequency = branching_frequency,
  mean_diameter_mm = diameter_stats$mean,
  root_density = root_density
)

write.csv(export_results, "flatbed_scan_results.csv", row.names = FALSE)
```

###### 

For more detailed information about specific functions and additional
features, please refer to the function documentation and the package
GitHub repository at <https://github.com/jcunow/RootScanR>.
