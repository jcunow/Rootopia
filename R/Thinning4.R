
##################################
#' Validate Image Input Parameters
#'
#' Internal function to validate input parameters for image processing functions
#' Used for image skeletonizing
#' @param img The input image to validate
#' @param allow_empty Logical, whether to allow empty images
#' @param min_dim Minimum required dimensions
#' @param require_binary Logical, whether to require binary values
#' @param select.layer Layer to validate for multi-layer images
#'
#' @return List containing validated and processed image data
#' @keywords internal
validate_image_input <- function(img,
                                 allow_empty = FALSE,
                                 min_dim = c(3,3),
                                 require_binary = TRUE,
                                 select.layer = NULL) {
  # Check if input is NULL or missing
  if (is.null(img)) {
    stop("Input image cannot be NULL")
  }
  
  # Check if dimensions are sufficient
  if (any(dim(img)[1:2] < min_dim)) {
    stop(sprintf("Image dimensions must be at least %dx%d", min_dim[1], min_dim[2]))
  }
  
  # Check for empty images if not allowed
  if (!allow_empty && all(terra::values(img) == 0)) {
    stop("Input image cannot be empty (all zeros)")
  }
  
  # Check for single-value images
  if (length(unique(as.vector(img))) == 1) {
    warning("Input image contains only a single value")
  }
  
  # Check for NaN or Inf values
  if (any(is.nan(terra::values(img))) || any(is.infinite(terra::values(img)))) {
    stop("Image contains NaN or Infinite values")
  }
  
  # Handle numerical precision and force binary values if required
  if (require_binary) {
    if (!all(terra::values(img) %in% c(0,1))) {
      if (all(abs(img - round(img)) < 1e-10)) {
        warning("Non-integer values detected, rounding to nearest integer")
        img <- round(img)
      }
      if (!all(img %in% c(0,1))) {
        stop("Image must contain only binary values (0 or 1)")
      }
    }
  }
  
  
  return(list(
    img = img,
    dims = dim(img),
    is_empty = all(img == 0),
    unique_values = unique(as.vector(img))
  ))
}



#' Thin Binary Image using Zhang-Suen Algorithm (Internal)
#'
#' This internal function performs image thinning using the Zhang-Suen thinning algorithm.
#' It reduces binary images to their skeleton while preserving the structure and connectivity of the foreground pixels.
#'
#' @param img A matrix, data frame, or `SpatRaster` object representing the binary image to be thinned.
#' @param verbose Logical. If `TRUE`, prints diagnostic information such as iteration progress and pixel removal counts. Default is `TRUE`.
#' @param select.layer Integer indicating the layer to use if `img` is a multi-layer `SpatRaster`. Default is 2.
#'
#' @details
#' - The function first prepares the image using the \code{\link{load_flexible_image}} function, ensuring binary matrix format.
#' - Thinning is performed iteratively in two subiterations per cycle:
#'   1. Identifying pixels to be removed based on Zhang-Suen conditions (first subiteration).
#'   2. Refining removal decisions in the second subiteration.
#' - The algorithm continues until no pixels are removed in an iteration or a maximum number of iterations is reached (default: 1000).
#'
#' @return A binary matrix representing the thinned image (skeleton).
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Example usage
#' raster <- terra::rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0,0), nrow = 3))
#' thinned_image <- thin_image_zhangsuen(raster, verbose = TRUE, select.layer = NULL)
#' }
thin_image_zhangsuen <- function(img, verbose = TRUE, select.layer = 2) {
  
  img <- load_flexible_image(
    img, 
    select.layer = select.layer, 
    output_format = "matrix", 
    normalize = TRUE, 
    binarize = TRUE
  )
  
  if (verbose) {
    cat("Image dimensions:", nrow(img), "x", ncol(img), "\n")
    cat("Initial foreground pixels:", sum(img == 1), "\n")
  }
  
  img[abs(img - 1) < 1e-10] <- 1
  img[abs(img) < 1e-10] <- 0
  
  if (sum(img == 1) == 0) stop("No foreground pixels found after type conversion")
  
  count_transitions <- function(p) {
    neighbors <- c(p[2:9], p[2])
    sum(neighbors[1:8] == 0 & neighbors[2:9] == 1)
  }
  
  get_neighbors <- function(img, i, j) {
    p <- rep(0, 9)
    p[1] <- img[i, j]
    if (i > 1) {
      p[2] <- img[i-1, j]
      if (j < ncol(img)) p[3] <- img[i-1, j+1]
    }
    if (j < ncol(img)) p[4] <- img[i, j+1]
    if (i < nrow(img)) {
      if (j < ncol(img)) p[5] <- img[i+1, j+1]
      p[6] <- img[i+1, j]
      if (j > 1) p[7] <- img[i+1, j-1]
    }
    if (j > 1) {
      p[8] <- img[i, j-1]
      if (i > 1) p[9] <- img[i-1, j-1]
    }
    return(p)
  }
  
  max_iterations <- 1000
  iteration_count <- 0
  any_changes <- TRUE
  
  while (any_changes && iteration_count < max_iterations) {
    iteration_count <- iteration_count + 1
    any_changes <- FALSE
    pixels_removed <- 0
    
    for (subiter in 1:2) {
      to_delete <- matrix(FALSE, nrow=nrow(img), ncol=ncol(img))
      for (i in 2:(nrow(img)-1)) {
        for (j in 2:(ncol(img)-1)) {
          if (img[i, j] == 1) {
            p <- get_neighbors(img, i, j)
            B <- sum(p[-1])
            A <- count_transitions(p)
            
            if (B >= 2 && B <= 6 && A == 1) {
              if (subiter == 1) {
                if (p[2] * p[4] * p[6] == 0 && p[4] * p[6] * p[8] == 0) {
                  to_delete[i, j] <- TRUE
                  any_changes <- TRUE
                }
              } else {
                if (p[2] * p[4] * p[8] == 0 && p[2] * p[6] * p[8] == 0) {
                  to_delete[i, j] <- TRUE
                  any_changes <- TRUE
                }
              }
            }
          }
        }
      }
      pixels_removed <- pixels_removed + sum(to_delete)
      img[to_delete] <- 0
    }
    
    if (verbose && pixels_removed > 0) {
      cat("Iteration", iteration_count, ": Removed", pixels_removed, "pixels\n")
    }
  }
  
  if (verbose) {
    cat("Final foreground pixels:", sum(img == 1), "\n")
    cat("Total iterations:", iteration_count, "\n")
  }
  
  return(img)
}


#' Enhanced Guo-Hall Thinning Algorithm
#' Thin Binary Image using Guo-Hall Algorithm (Internal)
#'
#' This internal function applies the Guo-Hall thinning algorithm to reduce binary images to their skeletons while preserving connectivity and structure.
#'
#' @param img A matrix, data frame, or `SpatRaster` object representing the binary image to be thinned.
#' @param verbose Logical. If `TRUE`, outputs diagnostic information such as image dimensions, pixel removal counts, and iteration progress. Default is `FALSE`.
#' @param select.layer Integer indicating the layer to use if `img` is a multi-layer `SpatRaster`. Default is 2.
#'
#' @details
#' - The input image is first processed using \code{\link{load_flexible_image}} to ensure it is a binary matrix.
#' - Thinning is performed in an iterative process consisting of two subiterations per cycle:
#'   1. In the first subiteration, pixels are marked for removal based on specific Guo-Hall conditions.
#'   2. In the second subiteration, a different set of conditions is applied to mark additional pixels for removal.
#' - The process continues until no pixels are removed in an iteration or the maximum number of iterations (default: 1000) is reached.
#' - The Guo-Hall algorithm ensures that the skeleton of the image is preserved.
#'
#' @return A binary matrix representing the thinned image (skeleton).
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Example usage
#' raster <- terra::rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0,0), nrow = 3))
#' thinned_image <- thin_image_guohall(raster, verbose = TRUE)
#' }
thin_image_guohall <- function(img, verbose = FALSE, select.layer = 2) {
  
  img <- load_flexible_image(
    img, 
    select.layer = select.layer, 
    output_format = "matrix", 
    normalize = TRUE, 
    binarize = TRUE
  )
  
  if (verbose) {
    cat("Image dimensions:", nrow(img), "x", ncol(img), "\n")
    cat("Initial foreground pixels:", sum(img == 1), "\n")
  }
  
  img[abs(img - 1) < 1e-10] <- 1
  img[abs(img) < 1e-10] <- 0
  
  get_neighbors <- function(img, i, j) {
    p <- rep(0, 9)
    p[1] <- img[i, j]
    if (i > 1) {
      p[2] <- img[i-1, j]
      if (j < ncol(img)) p[3] <- img[i-1, j+1]
    }
    if (j < ncol(img)) p[4] <- img[i, j+1]
    if (i < nrow(img)) {
      if (j < ncol(img)) p[5] <- img[i+1, j+1]
      p[6] <- img[i+1, j]
      if (j > 1) p[7] <- img[i+1, j-1]
    }
    if (j > 1) {
      p[8] <- img[i, j-1]
      if (i > 1) p[9] <- img[i-1, j-1]
    }
    return(p)
  }
  
  check_pixel <- function(p) {
    C1 <- function(p) {
      s <- sum(p[2:9])
      s >= 2 && s <= 6
    }
    
    C2 <- function(p) {
      neighbors <- c(p[2:9], p[2])
      sum(neighbors[1:8] == 0 & neighbors[2:9] == 1) == 1
    }
    
    C3_C4_first <- function(p) {
      ((p[2] * p[4] * p[6] == 0) || (p[4] * p[6] * p[8] == 0)) &&
        ((p[2] * p[4] * p[8] == 0) || (p[2] * p[6] * p[8] == 0))
    }
    
    C3_C4_second <- function(p) {
      ((p[2] * p[4] * p[6] == 0) || (p[2] * p[4] * p[8] == 0)) &&
        ((p[2] * p[6] * p[8] == 0) || (p[4] * p[6] * p[8] == 0))
    }
    
    list(
      first = C1(p) && C2(p) && C3_C4_first(p),
      second = C1(p) && C2(p) && C3_C4_second(p)
    )
  }
  
  max_iterations <- 1000
  iteration_count <- 0
  any_changes <- TRUE
  
  while (any_changes && iteration_count < max_iterations) {
    iteration_count <- iteration_count + 1
    any_changes <- FALSE
    pixels_removed <- 0
    
    for (subiter in c("first", "second")) {
      to_delete <- matrix(FALSE, nrow=nrow(img), ncol=ncol(img))
      for (i in 2:(nrow(img)-1)) {
        for (j in 2:(ncol(img)-1)) {
          if (img[i, j] == 1) {
            p <- get_neighbors(img, i, j)
            if (check_pixel(p)[[subiter]]) {
              to_delete[i, j] <- TRUE
              any_changes <- TRUE
            }
          }
        }
      }
      pixels_removed <- pixels_removed + sum(to_delete)
      img[to_delete] <- 0
    }
    
    if (verbose && pixels_removed > 0) {
      cat("Iteration", iteration_count, ": Removed", pixels_removed, "pixels\n")
    }
  }
  
  if (verbose) {
    cat("Final foreground pixels:", sum(img == 1), "\n")
    cat("Total iterations:", iteration_count, "\n")
  }
  
  return(img)
}


#' Medial Axis Transform (Internal)
#'
#' This internal function computes the medial axis transform of a binary image, identifying the set of skeleton points equidistant to the object's boundaries.
#'
#' @param img A matrix, data frame, or `SpatRaster` object representing the binary image for transformation.
#' @param verbose Logical. If `TRUE`, outputs diagnostic information such as image dimensions, progress of computation, and final skeleton size. Default is `FALSE`.
#' @param select.layer Integer indicating the layer to use if `img` is a multi-layer `SpatRaster`. Default is 2.
#'
#' @details
#' - The input image is first processed using \code{\link{load_flexible_image}} to ensure it is a binary matrix.
#' - The algorithm proceeds through the following steps:
#'   1. **Distance Transform**: Computes the distance of each foreground pixel to the nearest background pixel using a two-pass algorithm.
#'   2. **Local Maxima Detection**: Identifies local maxima in the distance transform to mark potential skeleton points.
#'   3. **Skeleton Refinement**: Ensures connectivity by connecting skeleton points within an 8-neighborhood.
#' - The result is a binary image representing the medial axis of the input object.
#'
#' @return A binary matrix where `1` represents skeleton pixels and `0` represents the background.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Example usage
#' raster <- terra::rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0,1), nrow = 3))
#' skeleton <- medial_axis_transform(raster, verbose = TRUE)
#' }
medial_axis_transform <- function(img, verbose = FALSE, select.layer = NULL) {
  
  
  # Flexible input processing with error handling
  img <- tryCatch({
    result <- load_flexible_image(img, output_format = "spatrast",
                                  normalize = TRUE,binarize = T, select.layer = select.layer)
    matrix(as.numeric(result), nrow = nrow(result))
  }, error = function(e) {
    stop("Failed to load or convert image: ", e$message)
  })
  


  # Check for numerical stability in distance calculations
  check_numerical_stability <- function(dist_map) {
    if (any(is.nan(dist_map)) || any(is.infinite(dist_map))) {
      stop("Numerical instability detected in distance transform")
    }
    invisible(TRUE)
  }
  
  # Enhanced compute_distance_transform with error checking
  compute_distance_transform <- function(binary_img) {
    if (!all(binary_img %in% c(0,1))) {
      stop("Distance transform requires binary input")
    }
    
    dist_map <- matrix(0, nrow = nrow(binary_img), ncol = ncol(binary_img))
    
    tryCatch({
      # First pass - forward scan
      for(i in 2:nrow(binary_img)) {
        for(j in 2:ncol(binary_img)) {
          if(binary_img[i,j] == 1) {
            dist_map[i,j] <- min(
              dist_map[i-1,j],
              dist_map[i,j-1],
              dist_map[i-1,j-1]
            ) + 1
          }
        }
      }
      
      # Check stability after first pass
      check_numerical_stability(dist_map)
      
      # Second pass - backward scan
      for(i in (nrow(binary_img)-1):1) {
        for(j in (ncol(binary_img)-1):1) {
          if(binary_img[i,j] == 1) {
            dist_map[i,j] <- min(
              dist_map[i,j],
              dist_map[i+1,j] + 1,
              dist_map[i,j+1] + 1,
              dist_map[i+1,j+1] + 1
            )
          }
        }
      }
      
      # Final stability check
      check_numerical_stability(dist_map)
      
    }, error = function(e) {
      stop("Error in distance transform computation: ", e$message)
    })
    
    return(dist_map)
  }
  
  # Enhanced find_local_maxima with boundary checking
  find_local_maxima <- function(dist_map) {
    maxima <- matrix(FALSE, nrow = nrow(dist_map), ncol = ncol(dist_map))
    
    tryCatch({
      for(i in 2:(nrow(dist_map)-1)) {
        for(j in 2:(ncol(dist_map)-1)) {
          if(dist_map[i,j] > 0) {
            # Safe neighborhood extraction with bounds checking
            i_range <- max(1, i-1):min(nrow(dist_map), i+1)
            j_range <- max(1, j-1):min(ncol(dist_map), j+1)
            neighborhood <- dist_map[i_range, j_range]
            
            if(dist_map[i,j] >= max(neighborhood)) {
              maxima[i,j] <- TRUE
            }
          }
        }
      }
    }, error = function(e) {
      stop("Error in local maxima detection: ", e$message)
    })
    
    return(maxima)
  }
  
  # Main processing with error handling
  tryCatch({
    dist_transform <- compute_distance_transform(img)
    if(verbose) cat("\nDistance transform computed\n")
    
    skeleton <- terra::t(find_local_maxima(dist_transform))
    result <- matrix(0, nrow = nrow(img), ncol = ncol(img))
    result[skeleton] <- 1
    
    # Validate final result
    if(all(result == 0)) {
      warning("No skeleton points detected in the result")
    }
    
  }, error = function(e) {
    stop("Medial axis transform failed: ", e$message)
  })
  
  return(result)
}

#' Detect Skeleton Points: Branching Points and Endpoints
#'
#' Identifies the branching points and endpoints of a skeletonized binary image.
#'
#' @param img A matrix, data frame, or `SpatRaster` object representing the skeletonized binary image.
#' @param select.layer Integer. Specifies which layer to use if the input is a multi-band image. Default is `2`.
#'
#' @details
#' This function detects key points in a skeletonized binary image:
#' \itemize{
#'   \item \strong{Endpoints}: Pixels with exactly one neighbor in the skeleton.
#'   \item \strong{Branching Points}: Pixels with more than two neighbors in the skeleton.
#' }
#'
#' The function uses a 3x3 neighborhood kernel to count the number of neighbors for each foreground pixel (\code{1}) in the image. Based on the neighbor count, points are classified as endpoints or branching points.
#'
#' The input image should be skeletonized (thin and connected) before using this function. If not already binary, the input image will be binarized internally.
#'
#' @return A named list containing two `SpatRaster` objects:
#' \itemize{
#'   \item \code{endpoints}: A binary raster where endpoints are marked as \code{1}.
#'   \item \code{branching_points}: A binary raster where branching points are marked as \code{1}.
#' }
#'
#' @examples
#' \dontrun{
#' library(terra)
#'
#' # Example skeletonized image
#' skeleton <- rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0), nrow = 4))
#'
#' # Detect endpoints and branching points
#' points <- detect_skeleton_points(skeleton)
#'
#' # Access results
#' endpoints <- points$endpoints
#' branching_points <- points$branching_points
#' }
#'
#' @seealso \code{\link{skeletonize_image}}, \code{\link{thin_image_zhangsuen}}, \code{\link{thin_image_guohall}}
#' @export
detect_skeleton_points <- function(img, select.layer = 2) {
  # Input validation
  tryCatch({
    validated <- validate_image_input(
      img = img,
      allow_empty = FALSE,
      min_dim = c(3,3),
      require_binary = TRUE,
      select.layer = select.layer
    )
    img <- validated$img
  }, error = function(e) {
    stop("Input validation failed: ", e$message)
  })
  
  # Flexible input processing with error handling
  img <- tryCatch({
    result <- load_flexible_image(img, select.layer = select.layer,
                                  output_format = "spatrast", normalize = TRUE)
    matrix(as.numeric(result), nrow = nrow(result))
  }, error = function(e) {
    stop("Failed to load image: ", e$message)
  })
  
  # Enhanced neighbor counting with bounds checking
  count_neighbors <- function(img) {
    if (!all(img %in% c(0,1))) {
      stop("Input must be binary for neighbor counting")
    }
    
    # Define the kernel
    kernel <- matrix(1, nrow = 3, ncol = 3)
    kernel[2, 2] <- 0
    
    # Safe padding
    padded_img <- tryCatch({
      result <- matrix(0, nrow = nrow(img) + 2, ncol = ncol(img) + 2)
      result[2:(nrow(img) + 1), 2:(ncol(img) + 1)] <- img
      result
    }, error = function(e) {
      stop("Failed to pad image: ", e$message)
    })
    
    # Initialize result matrix
    neighbor_count <- matrix(0, nrow = nrow(img), ncol = ncol(img))
    
    # Safe convolution with error checking
    tryCatch({
      for (i in 2:(nrow(padded_img) - 1)) {
        for (j in 2:(ncol(padded_img) - 1)) {
          region <- padded_img[(i - 1):(i + 1), (j - 1):(j + 1)]
          neighbor_count[i - 1, j - 1] <- sum(region * kernel)
        }
      }
    }, error = function(e) {
      stop("Error in neighbor counting: ", e$message)
    })
    
    return(neighbor_count)
  }
  
  # Main processing with validation
  tryCatch({
    neighbor_count <- count_neighbors(img)
    
    endpoints <- (img == 1) & (neighbor_count == 1)
    branching_points <- (img == 1) & (neighbor_count > 2)
    
    # Validate results
    if(sum(endpoints) == 0 && sum(branching_points) == 0) {
      warning("No endpoints or branching points detected")
    }
    
    # Convert to SpatRaster with error handling
    endpoints_rast <- tryCatch({
      terra::rast(endpoints)
    }, error = function(e) {
      stop("Failed to convert endpoints to SpatRaster: ", e$message)
    })
    
    branching_points_rast <- tryCatch({
      terra::rast(branching_points)
    }, error = function(e) {
      stop("Failed to convert branching points to SpatRaster: ", e$message)
    })
    
    return(list(
      endpoints = endpoints_rast,
      branching_points = branching_points_rast
    ))
    
  }, error = function(e) {
    stop("Failed to detect skeleton points: ", e$message)
  })
}


#' Skeletonization Wrapper Function
#'
#' This function serves as a wrapper for applying different skeletonization methods to a binary image, including the Zhang-Suen, Guo-Hall, and Medial Axis Transform (MAT) algorithms.
#'
#' @param img A matrix, data frame, or `SpatRaster` object representing the binary image to be skeletonized.
#' @param methods A character vector specifying the skeletonization methods to apply. Valid options are \code{"ZhangSuen"}, \code{"GuoHall"}, and \code{"MAT"}. Defaults to all three methods.
#' @param verbose Logical. If \code{TRUE}, displays progress and diagnostic messages during processing. Defaults to \code{TRUE}.
#' @param select.layer Integer specifying the layer to use if \code{img} is a multi-layer `SpatRaster`. Defaults to 2.
#'
#' @details
#' This function allows for flexible and streamlined skeletonization of binary images using one or more supported algorithms:
#' \itemize{
#'   \item \code{"ZhangSuen"}: Implements the Zhang-Suen thinning algorithm.
#'   \item \code{"GuoHall"}: Implements the Guo-Hall thinning algorithm.
#'   \item \code{"MAT"}: Computes the Medial Axis Transform to extract the skeleton.
#' }
#'
#' The function processes the input image with the specified methods and returns the results. If multiple methods are chosen, the results are returned as a named list, with each element corresponding to a method.
#'
#' @return If a single method is selected, the function returns a `SpatRaster` object representing the skeletonized image. If multiple methods are selected, a named list of `SpatRaster` objects is returned.
#'
#' @examples
#' 
#' # Load a binary image as a SpatRaster
#' binary_image <- terra::rast(matrix(c(0, 1, 1, 0, 0, 1, 1, 0,0), nrow = 3))
#'
#' # Apply all skeletonization methods
#' skeletons <- skeletonize_image(binary_image, verbose = TRUE)
#' 
#'
#' @seealso \code{\link{thin_image_zhangsuen}}, \code{\link{thin_image_guohall}}, \code{\link{medial_axis_transform}}
#' @export
skeletonize_image <- function(img, methods = c("ZhangSuen", "GuoHall", "MAT", "SteepestAscend"), verbose = TRUE, select.layer = NULL) {
  # Ensure methods are valid
  valid_methods <- c("ZhangSuen", "GuoHall", "MAT", "SteepestAscend")
  methods <- intersect(methods, valid_methods)
  if (length(methods) == 0) {
    stop("No valid methods specified. Choose from: 'ZhangSuen', 'GuoHall', 'MAT', 'SteepestAscend'.")
  }
  
  
  # Process each method
  results <- list()
  for (method in methods) {
    if (verbose) cat("\nApplying method:", method, "\n")
    result <- switch(
      method,
      "ZhangSuen" =
        thin_image_zhangsuen(img, verbose = verbose, select.layer = NULL),
      "GuoHall" = thin_image_guohall(img, verbose = verbose, select.layer = NULL),
      "MAT" = medial_axis_transform(img, verbose = verbose, select.layer = NULL),
      "SteepestAscend" = thin_image_steepest_ascend(img, verbose = verbose, select.layer = select.layer),
      stop(paste("Unsupported method:", method))
    )
    results[[method]] <- result
  }
  
  # Return results as a list
  if (length(results) == 1) {
    return(results[[1]])  # Single method result
  } else {
    return(results)  # Multiple method results
  }
}



#' Steepest Ascend Skeletonization (Internal)
#'
#' Skeletonize a binary image using a ridge-based steepest ascend approximation on the distance transform.
#'
#' @param img A matrix, data frame, or `SpatRaster` object representing the binary image.
#' @param verbose Logical. If `TRUE`, prints diagnostic information. Default is `FALSE`.
#' @param select.layer Integer indicating the layer to use if `img` is a multi-layer `SpatRaster`. Default is `NULL`.
#'
#' @return A binary matrix representing the skeleton obtained from local maxima on the distance transform.
#' @keywords internal
thin_image_steepest_ascend <- function(img, verbose = FALSE, select.layer = NULL) {
  
  
  # Convert to binary matrix with error handling
  img <- tryCatch({
    result <- load_flexible_image(img, select.layer = select.layer,
                                  output_format = "matrix", normalize = TRUE, binarize = TRUE)
    result
  }, error = function(e) {
    stop("Failed to load image: ", e$message)
  })
  
  # Compute distance transform and identify ridge lines
  dt <- tryCatch({
    im <- imager::as.cimg(t(img))  # transpose to match cimg's (x,y) convention
    distmap <- imager::distance_transform(im, value = 0)
    local_maxima <- distmap == imager::isoblur(distmap, sigma = 1)
    ridge <- as.matrix(local_maxima)
    t(ridge) * img  # convert back and mask with original image
  }, error = function(e) {
    stop("Ridge detection failed: ", e$message)
  })
  
  if (verbose) {
    cat("Steepest ascend ridge-based skeletonization complete. Skeleton pixels:", sum(dt), "\n")
  }
  
  return(dt)
}
