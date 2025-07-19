

#' Smooth the edges of a root binary image
#'
#' @param img A binary image of root systems (imager cimg object)
#' @param kernel_shape Shape of the kernel: "square", "diamond", or "disk"
#' @param kernel_size Size of the kernel for morphological operations (odd integer). Use larger kernels for higher image resolution
#' @param iterations Number of iterations for the smoothing process
#' @return A smoothed binary image
#' @export
#' @keywords internal
#'
#' @importFrom imager dilate erode as.cimg is.cimg grayscale
#' @examples
#' 
#' data("seg_Oulanka2023_Session01_T067")
#' img <- seg_Oulanka2023_Session01_T067
#' # Try different kernel shapes
#' smoothed_square <- smooth_root_edges(img, kernel_shape = "square", kernel_size = 3)
#' smoothed_diamond <- smooth_root_edges(img, kernel_shape = "diamond", kernel_size = 3)
#' smoothed_disk <- smooth_root_edges(img, kernel_shape = "disk", kernel_size = 3)
#' plot(smoothed_disk)
#' 
smooth_root_edges <- function(img, kernel_shape = "disk", kernel_size = 3, iterations = 1) {

  # Try to convert to cimg
  img <- load_flexible_image(img, output_format = "cimg",binarize = T)
    

  # Ensure image is grayscale (3D array) if it's color (4D array)
  if (imager::spectrum(img) > 1) {
    img <- imager::grayscale(img)
  }
  
  # Binarize the image if needed
  #img <- img > 0.5
  
  # Create kernel based on selected shape
  kern <- create_kernel(shape = kernel_shape, size = kernel_size)
  
  # Apply closing operation (dilation followed by erosion) multiple times
  smoothed_img <- img
  for (i in 1:iterations) {
    smoothed_img <- imager::dilate(smoothed_img, kern)
    smoothed_img <- imager::erode(smoothed_img, kern)
  }
  
  return(smoothed_img)
}

#' Create a kernel with specified shape for morphological operations
#'
#' @param shape Shape of the kernel: "square", "diamond", or "disk"
#' @param size Size of the kernel (odd integer)
#' @return A kernel as an imager cimg object
#' @keywords internal
create_kernel <- function(shape = "disk", size = 3) {
  
  # Ensure methods are valid
  valid_shapes <- c( "square", "diamond", "disk")
  methods <- intersect(shape, valid_shapes)
  if (length(methods) == 0) {
    stop("No valid shape specified. Choose from: 'square', 'diamond', 'disk'.")
  }
  # Ensure odd size for symmetry
  if (size %% 2 == 0) size <- size + 1
  
  # Create empty kernel matrix
  kern_matrix <- matrix(0, size, size)
  center <- floor(size/2) + 1
  radius <- floor(size/2)
  
  # Fill kernel based on shape
  for (i in 1:size) {
    for (j in 1:size) {
      if (shape == "square") {
        # Square: all positions filled
        kern_matrix[i, j] <- 1
      } else if (shape == "diamond") {
        # Diamond: Manhattan distance <= radius
        if (abs(i - center) + abs(j - center) <= radius) {
          kern_matrix[i, j] <- 1
        }
      } else if (shape == "disk") {
        # Disk/Circle: Euclidean distance <= radius
        if (sqrt((i - center)^2 + (j - center)^2) <= radius) {
          kern_matrix[i, j] <- 1
        }
      } else if (shape == "cross") {
        # Cross shape
        if (i == center || j == center) {
          kern_matrix[i, j] <- 1
        }
      } else if (shape == "x") {
        # X shape (diagonals)
        if (abs(i - center) == abs(j - center)) {
          kern_matrix[i, j] <- 1
        }
      }
    }
  }
  
  # Convert matrix to cimg object
  return(imager::as.cimg(kern_matrix))
}








#################3 real function
#' Fill holes in binary images
#'
#' Identifies internal black regions ("holes") in a binary image and fills them by setting their pixel values to 1. Holes are defined as black areas (value = 0) completely surrounded by white (value = 1), i.e., not connected to the image border.
#'
#' @param img A `cimg` object representing a binary image (values 0 and 1).
#' @param max_size Optional maximum size (in pixels) of holes to fill. If `NULL`, all holes are filled.
#'
#' @return A `cimg` object with holes filled (as 1s).
#' @keywords internal
fill_holes <- function(img, max_size = NULL) {
  
  # 1. Use image directly - holes (0s) are already "foreground" for labeling
  # We want to label the black regions (holes and background)
  holes_and_bg <- 1 - img
  
  # 2. Label all connected components in holes and background
  lbl <- imager::label(holes_and_bg)
  
  # 3. Identify labels that touch the border — these are NOT holes
  border_labels <- unique(c(
    lbl[1, , 1, 1],
    lbl[dim(lbl)[1], , 1, 1],
    lbl[, 1, 1, 1],
    lbl[, dim(lbl)[2], 1, 1]
  ))
  
  # Remove 0 from border labels (background)
  border_labels <- border_labels[border_labels > 0]
  
  # 4. Keep only internal labels (potential holes)
  internal_mask <- !(lbl %in% border_labels)
  internal_labels <- lbl * as.numeric(internal_mask)
  
  # 5. Get hole regions and filter by size
  hole_ids <- unique(as.numeric(internal_labels[internal_labels > 0]))
  
  if (length(hole_ids) == 0) {
    return(img)  # No holes to fill
  }
  
  # Create a copy of the image to modify
  filled_img <- img
  
  for (id in hole_ids) {
    # Find all pixels belonging to this hole
    hole_coords <- which(internal_labels == id, arr.ind = TRUE)
    
    # Calculate the actual size of the hole
    # Count pixels that are currently 0 (black holes) in original image
    hole_size <- sum(img[hole_coords] == 0)
    
    # Only fill holes that are smaller than max_size (if specified)
    if (is.null(max_size) || hole_size <= max_size) {
      # Fill the hole by setting all hole pixels to 1 (white)
      for (i in 1:nrow(hole_coords)) {
        coord <- hole_coords[i, ]
        # Only fill if it's currently a hole (value 0)
        if (img[coord[1], coord[2], 1, 1] == 0) {
          filled_img[coord[1], coord[2], 1, 1] <- 1
        }
      }
    }
  }
  
  return(filled_img)
}

#' Remove small white artifacts from binary images
#'
#' Identifies small internal white regions (objects) in a binary image and removes them by setting their pixel values to 0. Artifacts are defined as white areas (value = 1) not connected to the image border.
#'
#' @param img A `cimg` object representing a binary image (values 0 and 1).
#' @param max_size Optional maximum size (in pixels) of white objects to remove. If `NULL`, all isolated objects are removed.
#'
#' @return A `cimg` object with small white artifacts removed (set to 0).
#' @keywords internal
remove_small_objects <- function(img, max_size = NULL) {
  
  # 1. Label white objects directly (no inversion needed)
  lbl <- imager::label(img)
  
  # 2. Identify labels that touch the border — these are large/connected objects to keep
  border_labels <- unique(c(
    lbl[1, , 1, 1],
    lbl[dim(lbl)[1], , 1, 1],
    lbl[, 1, 1, 1],
    lbl[, dim(lbl)[2], 1, 1]
  ))
  
  # Remove 0 from border labels (background)
  border_labels <- border_labels[border_labels > 0]
  
  # 3. Keep only internal labels (potential small artifacts)
  internal_mask <- !(lbl %in% border_labels)
  internal_labels <- lbl * as.numeric(internal_mask)
  
  # 4. Get artifact regions and filter by size
  artifact_ids <- unique(as.numeric(internal_labels[internal_labels > 0]))
  
  if (length(artifact_ids) == 0) {
    return(img)  # No artifacts to remove
  }
  
  # Create a copy of the image to modify
  cleaned_img <- img
  
  for (id in artifact_ids) {
    # Find all pixels belonging to this artifact
    artifact_coords <- which(internal_labels == id, arr.ind = TRUE)
    
    # Calculate the actual size of the artifact
    # Count pixels that are currently 1 (white objects) in original image
    artifact_size <- sum(img[artifact_coords] == 1)
    
    # Only remove artifacts that are smaller than max_size (if specified)
    if (is.null(max_size) || artifact_size <= max_size) {
      # Remove the artifact by setting all object pixels to 0 (black)
      for (i in 1:nrow(artifact_coords)) {
        coord <- artifact_coords[i, ]
        # Only remove if it's currently an object (value 1)
        if (img[coord[1], coord[2], 1, 1] == 1) {
          cleaned_img[coord[1], coord[2], 1, 1] <- 0
        }
      }
    }
  }
  
  return(cleaned_img)
}



#' Report sizes of holes and white artifacts in binary images
#'
#' Analyzes a binary image and prints the sizes of internal black holes (0-valued regions enclosed by 1s) and isolated white artifacts (1-valued regions not connected to the border). This function is useful for diagnosing what would be affected by `fill_holes()` and `remove_small_objects()` operations.
#'
#' Holes are detected by inverting the image and applying connected component labeling. Any region that does not touch the image border is considered a candidate hole or artifact.
#'
#' @param img A `cimg` object representing a binary image with values 0 and 1.
#'
#' @return None (invisible `NULL`). The function prints human-readable summaries to the console.
#'
#' @details
#' - Holes are black regions (0-valued pixels) completely enclosed by white (1-valued) areas.
#' - Artifacts are small white regions (1-valued pixels) that are not connected to the image border.
#' - Pixel counts are printed for each detected region.
#'
#' @examples
#' #' # Create a complex test image with holes and artifacts
#' 
#'   img <- imager::as.cimg(matrix(0, 150, 150))  # Start with black background
#'
#'   # Create multiple white objects with black holes
#'   img[20:50, 20:50] <- 1       # White square 1
#'   img[30:35, 30:35] <- 0       # Small black hole in square 1
#'
#'   img[70:120, 70:120] <- 1     # White square 2
#'   img[80:85, 80:85] <- 0       # Small black hole 1 in square 2
#'   img[100:115, 100:115] <- 0   # Large black hole 2 in square 2
#'   
#' # Add small artifacts (1-pixel specks)
#' img[10, 140] <- 1
#' img[145, 15] <- 1
#'
#' # Add a 2×2 speck
#' img[130:131, 40:41] <- 1
#'
#' # Add an irregular blob
#' img[100:102, 10] <- 1
#' img[101:102, 11] <- 1
#' img[101, 12] <- 1
#'
#'   # Create a white ring (donut shape)
#'   center_x <- 40
#'   center_y <- 100
#'for (i in 1:150) {
#'  for (j in 1:150) {
#'    dist <- sqrt((i - center_x)^2 + (j - center_y)^2)
#'    if (dist <= 20 && dist >= 10) {
#'      img[i, j,,] <- 1  
#'  }}}
#'  
#' report_image_components(img)
#'
#' @keywords internal
report_image_components <- function(img) {
  cat("=== HOLES (black regions inside objects) ===\n")
  
  # Report holes
  holes_and_bg <- 1 - img
  lbl_holes <- imager::label(holes_and_bg)
  
  border_labels_holes <- unique(c(
    lbl_holes[1, , 1, 1],
    lbl_holes[dim(lbl_holes)[1], , 1, 1],
    lbl_holes[, 1, 1, 1],
    lbl_holes[, dim(lbl_holes)[2], 1, 1]
  ))
  
  border_labels_holes <- border_labels_holes[border_labels_holes > 0]
  
  internal_mask_holes <- !(lbl_holes %in% border_labels_holes)
  internal_labels_holes <- lbl_holes * as.numeric(internal_mask_holes)
  hole_ids <- unique(as.numeric(internal_labels_holes[internal_labels_holes > 0]))
  
  if (length(hole_ids) == 0) {
    cat("No holes found\n")
  } else {
    hole_sizes <- sapply(hole_ids, function(id) {
      hole_coords <- which(internal_labels_holes == id, arr.ind = TRUE)
      sum(img[hole_coords] == 0)
    })
    
    for (i in seq_along(hole_ids)) {
      cat(sprintf("Hole %d: %d pixels\n", hole_ids[i], hole_sizes[i]))
    }
  }
  
  cat("\n=== ARTIFACTS (small isolated white objects) ===\n")
  
  # Report artifacts
  lbl_objects <- imager::label(img)
  
  border_labels_objects <- unique(c(
    lbl_objects[1, , 1, 1],
    lbl_objects[dim(lbl_objects)[1], , 1, 1],
    lbl_objects[, 1, 1, 1],
    lbl_objects[, dim(lbl_objects)[2], 1, 1]
  ))
  
  border_labels_objects <- border_labels_objects[border_labels_objects > 0]
  
  internal_mask_objects <- !(lbl_objects %in% border_labels_objects)
  internal_labels_objects <- lbl_objects * as.numeric(internal_mask_objects)
  artifact_ids <- unique(as.numeric(internal_labels_objects[internal_labels_objects > 0]))
  
  if (length(artifact_ids) == 0) {
    cat("No isolated artifacts found\n")
  } else {
    artifact_sizes <- sapply(artifact_ids, function(id) {
      artifact_coords <- which(internal_labels_objects == id, arr.ind = TRUE)
      sum(img[artifact_coords] == 1)
    })
    
    for (i in seq_along(artifact_ids)) {
      cat(sprintf("Artifact %d: %d pixels\n", artifact_ids[i], artifact_sizes[i]))
    }
  }
}










#' Clean binary images by filling holes, removing small artifacts, and optionally smoothing edges
#'
#' This function performs a comprehensive cleaning operation on a binary image by:
#' 1. Filling internal black holes (regions of 0 surrounded by 1s),
#' 2. Removing small internal white artifacts (1s not connected to the image border), and
#' 3. Optionally applying edge smoothing to refine boundaries of root structures or other objects.
#'
#' Holes and artifacts are detected using connected component labeling. Objects touching the image border are preserved and not modified. Pixel connectivity is assumed to be 4-connected.
#'
#' @param img A `cimg` object representing a binary image with pixel values 0 and 1.
#' @param max_hole_size Maximum size (in pixels) of black holes to fill. If `NULL`, all holes are filled.
#' @param max_artifact_size Maximum size (in pixels) of white artifacts to remove. If `NULL`, all isolated objects are removed.
#' @param edge_smooth Logical; if `TRUE`, applies morphological smoothing to object edges.
#' @param kernel_shape Shape of the morphological kernel used for smoothing. One of `"disk"` (default), `"square"`, etc., depending on your implementation of `smooth_root_edges()`.
#' @param kernel_size Size of the structuring element used for edge smoothing.
#' @param iterations Number of times the smoothing operation is applied.
#' @param report Logical; if `TRUE`, the function returns a list with the cleaned image and a printed report on hole and artifact sizes. Defaults to `FALSE`.
#'
#' @return A cleaned `cimg` object. If `report = TRUE`, returns a list with two elements: the cleaned image and the printed summary.
#'
#' @examples
#' # Create a complex test image with holes and artifacts
#' 
#'   img <- imager::as.cimg(matrix(0, 150, 150))  # Start with black background
#'
#'   # Create multiple white objects with black holes
#'   img[20:50, 20:50] <- 1       # White square 1
#'   img[30:35, 30:35] <- 0       # Small black hole in square 1
#'
#'   img[70:120, 70:120] <- 1     # White square 2
#'   img[80:85, 80:85] <- 0       # Small black hole 1 in square 2
#'   img[100:115, 100:115] <- 0   # Large black hole 2 in square 2
#'   
#' # Add small artifacts (1-pixel specks)
#' img[10, 140] <- 1
#' img[145, 15] <- 1
#'
#' # Add a 2×2 speck
#' img[130:131, 40:41] <- 1
#'
#' # Add an irregular blob
#' img[100:102, 10] <- 1
#' img[101:102, 11] <- 1
#' img[101, 12] <- 1
#'
#'   # Create a white ring (donut shape)
#'   center_x <- 40
#'   center_y <- 100
#'for (i in 1:150) {
#'  for (j in 1:150) {
#'    dist <- sqrt((i - center_x)^2 + (j - center_y)^2)
#'    if (dist <= 20 && dist >= 10) {
#'      img[i, j,,] <- 1  
#'  }}}
#'
#'
#' # Clean with various thresholds
#' cleaned1 <- clean_image(img, max_hole_size = 50, max_artifact_size = 10)
#' cleaned2 <- clean_image(img, max_hole_size = 20, max_artifact_size = 30)
#' cleaned3 <- clean_image(img, max_hole_size = 30, max_artifact_size = 20, 
#'                         edge_smooth = TRUE, kernel_size = 3)
#'
#' # Plot results
#' par(mfrow = c(2, 2))
#' plot(img, main = "Original")
#' plot(cleaned1, main = "Fill ≤50, Remove ≤10")
#' plot(cleaned2, main = "Fill ≤20, Remove ≤30")
#' plot(cleaned3, main = "Fill ≤30, Remove ≤20 + Smooth")
#' par(mfrow = c(1, 1))
#' @export
clean_image <- function(img,
                        max_hole_size = NULL,
                        max_artifact_size = NULL,
                        edge_smooth = TRUE,
                        kernel_shape = "disk",
                        kernel_size = 3,
                        iterations = 1,
                        report = FALSE) {
  # First fill holes
  img_filled <- fill_holes(img, max_hole_size)
  
  # Then remove small artifacts
  img_cleaned <- remove_small_objects(img_filled, max_artifact_size)
  
  # Optionally smooth edges
  if (edge_smooth) {
    img_smooth <- smooth_root_edges(img_cleaned, kernel_shape = kernel_shape,
                                    kernel_size = kernel_size, iterations = iterations)
  } else {
    img_smooth <- img_cleaned
  }
  
  if (report) {
    report_image_components(img)
    return(list(image = img_smooth))
  } else {
    return(img_smooth)
  }
}









