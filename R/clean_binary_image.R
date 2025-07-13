

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
#' \dontrun{
#' data("seg_Oulanka2023_Session01_T067")
#' img <- seg_ seg_Oulanka2023_Session01_T067
#' # Try different kernel shapes
#' smoothed_square <- smooth_root_edges(img, kernel_shape = "square", kernel_size = 3)
#' smoothed_diamond <- smooth_root_edges(img, kernel_shape = "diamond", kernel_size = 3)
#' smoothed_disk <- smooth_root_edges(img, kernel_shape = "disk", kernel_size = 3)
#' plot(smoothed_disk)
#' }
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







#' Fill internal holes in root binary images while preserving edges
#'
#' @param img A binary image of root systems (imager cimg object)
#' @param min_hole_size Minimum size of holes to fill (in pixels)
#' @param max_hole_size Maximum size of holes to fill (in pixels, Inf means no limit)
#' @return A processed binary image with internal holes filled
#' @export
#' @keywords internal
#'
#' @importFrom imager as.cimg is.cimg fill spectrum
#' @examples
#' \dontrun{
#' data("seg_Oulanka2023_Session01_T067")
#' img <- seg_Oulanka2023_Session01_T067
#' filled <- fill_root_holes(img, min_hole_size = 1, max_hole_size = 100)
#' plot(filled)
#' }
fill_root_holes <- function(img, min_hole_size = 1, max_hole_size = Inf) {
  
  img <- load_flexible_image(img, output_format = "cimg", binarize = T)
  
  # Ensure image is grayscale if it's color
  if (imager::spectrum(img) > 1) {
    img <- imager::grayscale(img)
  }
  
  
  # Create a copy for hole filling
  filled_img <- img
  
  # First, label all connected components in the background (inverse of the image)
  # We need to invert the image to identify holes
  inv_img <- !img
  
  # Label connected regions in the inverted image (holes)
  labeled <- imager::label(inv_img)
  
  # Get unique labels (excluding 0, which is the outer background)
  unique_labels <- unique(as.vector(labeled))
  unique_labels <- unique_labels[unique_labels != 0]
  
  # Process each labeled region (potential hole)
  for (lbl in unique_labels) {
    # Create a mask for this hole
    hole_mask <- labeled == lbl
    
    # Calculate hole size
    hole_size <- sum(hole_mask)
    
    # Check if the hole touches the image border (not a true hole)
    border_touch <- any(hole_mask[1,]) || any(hole_mask[nrow(hole_mask),]) || 
      any(hole_mask[,1]) || any(hole_mask[,ncol(hole_mask)])
    
    # Fill the hole if it's within size range and doesn't touch border
    if (!border_touch && hole_size >= min_hole_size && hole_size <= max_hole_size) {
      filled_img[hole_mask] <- TRUE
    }
  }
  
  return(filled_img)
}





#' Remove small foreground objects (e.g., dirt) from root binary images
#'
#' @param img A binary image of root systems (imager cimg object)
#' @param min_size Minimum size in pixels for foreground objects to keep
#' @param connectivity Type of connectivity for connected components: 4 or 8 (default: 8)
#' @return A processed binary image with small foreground objects removed (0=background, 1=foreground)
#' @export
#' @keywords internal
#'
#' @importFrom imager as.cimg is.cimg label spectrum
#' @examples
#' \dontrun{
#' data("seg_Oulanka2023_Session01_T067")
#' img <-  seg_Oulanka2023_Session01_T067
#' cleaned <- remove_small_foreground_objects(img, min_size = 100)
#' plot(cleaned)
#' }
remove_small_foreground_objects <- function(img, min_size = 100, connectivity = 8) {
  
  binary_img <- load_flexible_image(img, output_format = "cimg", binarize = TRUE)
  
  # Ensure image is grayscale if it's color
  if (imager::spectrum(binary_img) > 1) {
    binary_img <- imager::grayscale(binary_img)
  }
  
  
  # Label connected components in the FOREGROUND ONLY
  # Use 4-connectivity or 8-connectivity based on parameter
  labeled <- if (connectivity == 4) {
    imager::label(binary_img, high_connectivity = FALSE)
  } else {
    imager::label(binary_img, high_connectivity = TRUE)
  }
  
  # Start with original image
  cleaned_img <- binary_img
  
  # Get unique labels (excluding 0, which is the background)
  unique_labels <- unique(as.vector(labeled))
  unique_labels <- unique_labels[unique_labels != 0]
  
  # Process each foreground component
  for (lbl in unique_labels) {
    # Create mask for this object
    obj_mask <- labeled == lbl
    
    # Calculate object size
    obj_size <- sum(obj_mask)
    
    # If the foreground object is too small, convert it to background
    if (obj_size < min_size) {
      cleaned_img[obj_mask] <- 0  # Set to 0 (background)
    }
    # Large foreground objects remain unchanged (1)
  }
  
  cleaned_img <-  load_flexible_image(cleaned_img, output_format = "spatrast", binarize = TRUE)
  
  return(cleaned_img)
}






#' Process root images by filling holes and optionally smoothing edges
#'
#' @param img A binary image of root systems (imager cimg object)
#' @param remove_small_particles Boolean indicating whether only big objects are retained
#' @param min_size Minimum particle size (in pixel). Removing all patches below the minimal size
#' @param connectivity Determines which connection type is used to assess whether pixel belong to the same patch: 4 or 8   
#' @param fill_holes Boolean indicating whether to fill holes
#' @param min_hole_size Minimum size of holes to fill (in pixels)
#' @param max_hole_size Maximum size of holes to fill (in pixels, Inf means no limit)
#' @param smooth_edges Boolean indicating whether to smooth edges
#' @param kernel_shape Shape of the kernel for edge smoothing: "square", "diamond", or "disk"
#' @param kernel_size Size of the kernel for edge smoothing (odd integer)
#' @param iterations Number of iterations for the smoothing process
#' @return A processed binary image
#' @export
#' 
#' @seealso \code{\link{smooth_root_edges}}, \code{\link{fill_root_holes}}, \code{\link{remove_small_foreground_objects}}
#'
#' @examples
#' data("seg_Oulanka2023_Session01_T067")
#' img <-  seg_Oulanka2023_Session01_T067
#' # Fill holes only
#' result1 <- clean_root_image(img, fill_holes = TRUE, smooth_edges = FALSE)
#' # Fill holes and smooth edges
#' result2 <- clean_root_image(img, fill_holes = TRUE, smooth_edges = TRUE)
#' plot(result2)
clean_root_image <- function(img,
                             remove_small_particles = TRUE, min_size  = 100, connectivity = 8,
                             fill_holes = TRUE, min_hole_size = 1, max_hole_size = Inf, 
                             smooth_edges = FALSE, kernel_shape = "disk", kernel_size = 3, iterations = 1) {
  
  result_img <- img
  
  # Step 0: Remove small particales
  if(remove_small_particles){
    result_img <- remove_small_foreground_objects(result_img, min_size = min_size,connectivity = 8)
  }
  
  # Step 1: Fill holes if requested
  if (fill_holes) {
    result_img <- fill_root_holes(result_img, min_hole_size, max_hole_size)
  }
  
  # Step 2: Smooth edges if requested
  if (smooth_edges) {
    result_img <- smooth_root_edges(result_img, kernel_shape, kernel_size, iterations)
  }
  
  # ensure binary spatrast output
  result_img <-  load_flexible_image(result_img, output_format = "spatrast", binarize = TRUE)
  
  return(result_img)
}

