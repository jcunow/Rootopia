# Rotation Bias and Circadian Analysis with RootScanR

## Rotation Bias Analysis

### Introduction

This vignette demonstrates how to analyze rotation and methodological
bias of a minirhizotron setup using the **RootScanR** package. In
minirhizotron studies, roots may exhibit preferential growth areas on
the tube surface depending on the tube size and tube insertion angle.
This might also affect how roots behave at the surface itself by e.g.,
changing root growth angles.

### Installation

``` r
# Install the package from GitHub
# install.packages("remotes")
# remotes::install_github("jcunow/RootScanR")

# Load the package
library(RootScanR)
library(ggplot2)  # For enhanced plotting
```

### Key Functionalities

The rotation bias and circadian analysis module includes:

1.  **Sine Curve Fitting**: Model circadian rhythms in root growth data
2.  **Statistical Testing**: Assess the significance of observed
    rhythmic patterns
3.  **Bias Quantification**: Calculate metrics for directional growth
    preferences
4.  **Visualization Tools**: Create informative plots of temporal
    patterns

### Example Workflow

#### 1. Preparing Time-Series Data

**Purpose**: Organize your time-series observations of root growth for
analysis.

``` r
# Example: Loading time-series data
# In a real scenario, you would load your own data from files
# The data should contain timestamps and root measurements

# For demonstration, we'll generate synthetic data
set.seed(32608)  # For reproducibility
period <- 1     # cycle length
n <- 150         # Number of observations
timestamps <- runif(n, 0, 1*period)  # Random sampling times across 3 days
amplitude <- 2   # Strength of the rhythm
phase <- 6       # Time offset (hours)
baseline <- 3    # Average growth rate
noise <- 0.8       # Random variation

# Simulate root growth with circadian pattern
growth_rates <- amplitude * sin(2 * pi / period * (timestamps + phase)) + 
                baseline + rnorm(n, 0, noise)

# Create a data frame for analysis
root_ts_data <- data.frame(
  time = timestamps,
  growth_rate = growth_rates
)

# Visualize the raw time-series data
ggplot(root_ts_data, aes(x = time, y = growth_rate)) +
  geom_point() +
  labs(
    title = "Root Growth Rate Time Series",
    x = "Time (hours)",
    y = "Growth Rate (mm/hour)"
  ) +
  theme_minimal()
```

#### 2. Fitting Circadian Models

**Purpose**: Detect and characterize circadian rhythms in root growth
patterns.

``` r
root_ts_data = root_ts_data %>% arrange(root_ts_data$time)
# Fit a sine curve to the time-series data
fit <- fit_sine_curve(
  tt = root_ts_data$time,
  yy = root_ts_data$growth_rate,
  parStart = list(amp = 2, phase = 6, offset = 3, period = 1)
)

# Display model parameters
print("sine curve fit:")
print(fit)

# Plot the data with fitted curve
plot(
  root_ts_data$time, 
  root_ts_data$growth_rate,
  title = "Circadian Rhythm in Root Growth"
)
lines(root_ts_data$time, fit$predicted, col = "coral" )
```

#### 3. Statistical Testing for Rhythmicity

**Purpose**: Determine if the observed patterns are statistically
significant circadian rhythms.

``` r
# Perform statistical tests to validate circadian rhythms (all identical - check)

# 1. Wald Test (tests if amplitude is significantly different from zero)
wald_result <- rhythmicity_test(
  method = "Wald",
  tt = root_ts_data$time,
  yy = root_ts_data$growth_rate,
  period = 24
)
print("Wald test for circadian rhythmicity:")
print(wald_result)

# 2. Likelihood Ratio Test (compares rhythmic vs. constant models)
lr_result <- rhythmicity_test(
  method = "LR",
  tt = root_ts_data$time,
  yy = root_ts_data$growth_rate,
  period = 24
)
print("Likelihood Ratio test for circadian rhythmicity:")
print(lr_result)

# 3. F Test (tests overall model significance)
f_result <- rhythmicity_test(
  method = "F",
  tt = root_ts_data$time,
  yy = root_ts_data$growth_rate,
  period = 24
)
print("F test for circadian rhythmicity:")
print(f_result)
```

#### 4. Period Analysis

**Purpose**: Identify the dominant period in root growth patterns when
the period is unknown.

``` r
# Scan for the most likely period in the data
period_scan <- analyze_circadian_period(
  time = root_ts_data$time,
  values = root_ts_data$growth_rate,
  period_range = c(12, 36),  # Test periods between 12 and 36 hours
  step = 0.5                 # Step size in hours
)

# Plot period scan results
plot_period_scan(
  period_scan,
  title = "Period Analysis of Root Growth"
)

# Identify the best-fitting period
best_period <- period_scan$period[which.min(period_scan$residual_sum_squares)]
print(paste0("Best-fitting period: ", round(best_period, 2), " hours"))
```

#### 5. Rotation Bias Analysis

**Purpose**: Quantify directional preferences in root growth around
minirhizotron tubes.

``` r
# For rotation bias, we need data that includes angular positions
# In a minirhizotron, this would be the angular position around the tube

# For demonstration, we'll generate synthetic rotation data
set.seed(42)
n_roots <- 100
angular_positions <- runif(n_roots, 0, 360)  # Angles in degrees

# Simulate a bias toward the bottom of the tube (180 degrees)
bias_center <- 180
bias_strength <- 0.5
density <- bias_strength * dnorm(angular_positions, bias_center, 45) + 
           (1 - bias_strength) * dunif(angular_positions, 0, 360)
           
# Normalize to create a probability of root presence
probability <- density / sum(density) * n_roots
root_counts <- rpois(n_roots, probability)

# Create a data frame for analysis
rotation_data <- data.frame(
  angle = angular_positions,
  root_count = root_counts
)

# Analyze rotation bias
rotation_bias <- analyze_rotation_bias(
  angles = rotation_data$angle,
  values = rotation_data$root_count,
  bins = 24  # Number of angular bins
)

print("Rotation bias metrics:")
print(rotation_bias$metrics)

# Visualize rotation bias
plot_rotation_bias(
  rotation_bias,
  title = "Root Growth Rotation Bias"
)
```

#### 6. Combined Temporal and Directional Analysis

**Purpose**: Investigate if rotation bias changes with circadian
rhythms.

``` r
# For this analysis, we need data that combines time and angular position
# Generate synthetic data for demonstration

# Create time points across 3 days
time_points <- seq(0, 72, by = 3)
n_times <- length(time_points)

# Create angular positions
angle_points <- seq(0, 345, by = 15)
n_angles <- length(angle_points)

# Create a matrix of root counts with time-dependent rotation bias
root_matrix <- matrix(0, nrow = n_times, ncol = n_angles)
colnames(root_matrix) <- angle_points
rownames(root_matrix) <- time_points

# Simulate time-varying bias
for (i in 1:n_times) {
  time <- time_points[i]
  
  # Bias center shifts over time (simulating heliotropism)
  bias_center <- (180 + 45 * sin(2 * pi * time / 24)) %% 360
  
  # Generate root counts with the time-dependent bias
  for (j in 1:n_angles) {
    angle <- angle_points[j]
    distance <- min(abs(angle - bias_center), 360 - abs(angle - bias_center))
    root_matrix[i, j] <- rpois(1, lambda = 5 * exp(-distance^2/1000) + 1)
  }
}

# Analyze time-dependent rotation bias
temporal_bias <- analyze_temporal_rotation_bias(
  root_matrix = root_matrix,
  time_points = time_points,
  angle_points = angle_points,
  period = 24
)

print("Temporal rotation bias results:")
print(temporal_bias$summary)

# Visualize time-dependent rotation bias
plot_temporal_bias(
  temporal_bias,
  title = "Time-Dependent Root Growth Direction"
)
```

### Advanced Applications

#### Comparing Rotation Bias Between Treatments

**Purpose**: Determine if different experimental treatments affect root
directional preferences.

``` r
# Simulate data for two treatments
set.seed(123)

# Treatment 1: Strong downward bias
angles_t1 <- rnorm(100, mean = 180, sd = 30) %% 360
counts_t1 <- rpois(100, lambda = 5)

# Treatment 2: Weaker, more diffuse bias
angles_t2 <- rnorm(100, mean = 180, sd = 60) %% 360
counts_t2 <- rpois(100, lambda = 5)

# Analyze rotation bias for each treatment
bias_t1 <- analyze_rotation_bias(angles_t1, counts_t1, bins = 24)
bias_t2 <- analyze_rotation_bias(angles_t2, counts_t2, bins = 24)

# Compare bias strength statistically
bias_comparison <- compare_rotation_bias(
  bias_t1, 
  bias_t2,
  method = "rayleigh"
)

print("Bias comparison results:")
print(bias_comparison)

# Visualize the comparison
plot_bias_comparison(
  list(Treatment1 = bias_t1, Treatment2 = bias_t2),
  title = "Rotation Bias Comparison Between Treatments"
)
```

### Conclusion

This vignette demonstrates how to analyze rotation bias and circadian
rhythms in root growth using the RootScanR package. These analyses can
reveal important aspects of root growth strategies and responses to
environmental cues or methodological quirks.

Key takeaways: 1. Circadian rhythms can be detected and characterized
using sine curve fitting 2. Statistical tests help validate the
significance of observed rhythms 3. Rotation bias analysis quantifies
directional preferences in root growth 4. Combined temporal and
directional analysis can reveal complex growth patterns

For more detailed information about specific functions and additional
features, please refer to the function documentation and the package
GitHub repository at <https://github.com/jcunow/RootScanR>.
