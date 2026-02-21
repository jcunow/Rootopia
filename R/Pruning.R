#' ####################################
#' 
#' 
#' 
#' skeleton_pruner <- function(
#'     skl,          # SpatRaster skeleton (0/1)
#'     distmap,          # SpatRaster distance transform
#'     branch_fun = NULL  # optional custom branch-extraction function
#' ) {
#'   stopifnot(inherits(skl, "SpatRaster"), inherits(distmap, "SpatRaster"))
#'   
#'   # 1. Cell numbers and row/col indices of skeleton pixels
#'   cellnums <- terra::cells(skl, 1)
#'   rc <- terra::rowColFromCell(skl, cellnums$lyr.1)
#'   
#'   # 2. Branch extraction
#'   if (is.null(branch_fun)) {
#'     # placeholder: each pixel = its own branch
#'     branch_fun <- function(rc_mat) {
#'       lapply(seq_len(nrow(rc_mat)), function(i) rc_mat[i, , drop = FALSE])
#'     }
#'   }
#'   branches <- branch_fun(rc)
#'   
#'   # 3. Branch scores (length × local radius)
#'   dist_mat <- load_flexible_image(distmap, output_format = "spatrast", normalize = F)
#'     score <- vapply(
#'       branches,
#'       function(b) nrow(b) * min(dist_mat[b]),
#'       numeric(1)
#'     )
#'   branch_tbl <- data.frame(
#'     id    = seq_along(branches),
#'     score = score,
#'     len   = vapply(branches, nrow, integer(1))
#'   )
#'   
#'   # 4. Internal pruning function (fast)
#'   prune_at <- function(thresh) {
#'     keep <- which(branch_tbl$score >= thresh)
#'     out  <- matrix(FALSE, terra::nrow(skl), terra::ncol(skl))
#'     for (i in keep) out[branches[[i]]] <- TRUE
#'     terra::rast(out, extent = terra::ext(skl))
#'   }
#'   
#'   # Return an object with:
#'   structure(
#'     list(
#'       branches   = branches,
#'       branch_tbl = branch_tbl,
#'       prune      = prune_at
#'     ),
#'     class = "skeletonPruner"
#'   )
#' }
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' # Precompute once
#' prn <- skeleton_pruner(
#'   skl = skl,
#'   distmap = distmap
#' )
#' 
#' # Adjust pruning threshold any time—no recomputation
#' low    <- prn$prune(1)
#' medium <- prn$prune(1.1)
#' high   <- prn$prune(2)
#' 
#' terra::plot(low)
#' terra::plot(medium)
#' terra::plot(high)
#' 
#' 
#' 
#' #######################
#' #' Skeleton Branch Pruner with Whole-Branch Extraction
#' #'
#' #' Compute branch significance scores for a skeletonized root image
#' #' and interactively prune “weak” branches without recomputing the skeleton.
#' #'
#' #' This function precomputes all connected skeleton branches and their
#' #' salience scores (branch length × minimum local radius from a distance transform).
#' #' The returned object contains both the branch table and a fast
#' #' `prune(thresh)` method to generate pruned skeletons at any threshold.
#' #'
#' #' **Threshold guidance:**  
#' #' Minirhizotron scans typically produce fine, low-contrast skeletons
#' #' where meaningful side branches have low significance scores.
#' #' Flatbed scans produce thicker roots and tolerate higher thresholds.
#' #'
#' #' @param skl      A binary [terra::SpatRaster] skeleton (0/1 or TRUE/FALSE).
#' #' @param distmap  A [terra::SpatRaster] distance transform of the same extent
#' #'                 and resolution as `skl`.
#' #' @param connectivity 4 or 8; pixel connectivity for branch extraction.
#' #'
#' #' @return An object of class `"skeletonPruner"` containing:
#' #'   \itemize{
#' #'     \item \code{branches}: list of branch pixel coordinates
#' #'     \item \code{branch_tbl}: data.frame with columns \code{id}, \code{score}, \code{len}
#' #'     \item \code{prune}: a function \code{prune(thresh)} returning a
#' #'           [terra::SpatRaster] of the skeleton pruned at the given threshold
#' #'   }
#' #' 
#' #' @examples
#' #' \dontrun{
#' #' prn <- skeleton_pruner(skl, distmap, connectivity = 8)
#' #' low    <- prn$prune(1)
#' #' medium <- prn$prune(1.1)
#' #' high   <- prn$prune(2)
#' #' terra::plot(low)
#' #' terra::plot(medium)
#' #' terra::plot(high)
#' #' }
#' #'
#' #' @export
#' skeleton_pruner <- function(skl, distmap, connectivity = 8) {
#'   stopifnot(inherits(skl, "SpatRaster"), inherits(distmap, "SpatRaster"))
#'   stopifnot(connectivity %in% c(4,8))
#'   
#'   # --- Convert rasters to matrices
#'   mat <- as.matrix(skl)
#'   dist_mat <- as.matrix(distmap)
#'   
#'   # --- Define neighbors for DFS/flood-fill
#'   if(connectivity==4){
#'     nbrs <- function(r,c) rbind(c(r-1,c), c(r+1,c), c(r,c-1), c(r,c+1))
#'   } else {
#'     nbrs <- function(r,c) rbind(
#'       c(r-1,c-1), c(r-1,c), c(r-1,c+1),
#'       c(r,  c-1),           c(r,  c+1),
#'       c(r+1,c-1), c(r+1,c), c(r+1,c+1)
#'     )
#'   }
#'   
#'   # --- Flood-fill DFS to extract branches
#'   visited <- matrix(FALSE, nrow(mat), ncol(mat))
#'   dfs <- function(r0,c0){
#'     stack <- list(c(r0,c0))
#'     pixels <- matrix(ncol=2, nrow=0)
#'     while(length(stack) > 0){
#'       xy <- stack[[length(stack)]]
#'       stack <- stack[-length(stack)]
#'       r <- xy[1]; c <- xy[2]
#'       if(r<1 || c<1 || r>nrow(mat) || c>ncol(mat)) next
#'       if(mat[r,c]==0 || visited[r,c]) next
#'       visited[r,c] <<- TRUE
#'       pixels <- rbind(pixels, c(r,c))
#'       for(i in seq_len(nrow(nbrs(r,c)))) stack[[length(stack)+1]] <- nbrs(r,c)[i,]
#'     }
#'     pixels
#'   }
#'   
#'   branches <- list()
#'   for(r in seq_len(nrow(mat))){
#'     for(c in seq_len(ncol(mat))){
#'       if(mat[r,c]==1 && !visited[r,c]){
#'         branches[[length(branches)+1]] <- dfs(r,c)
#'       }
#'     }
#'   }
#'   
#'   # --- Compute branch scores: length × min radius
#'   score <- vapply(branches,
#'                   function(b) nrow(b) * min(dist_mat[cbind(b[,1], b[,2])]),
#'                   numeric(1))
#'   branch_tbl <- data.frame(
#'     id = seq_along(branches),
#'     score = score,
#'     len = vapply(branches, nrow, integer(1))
#'   )
#'   
#'   # --- Internal pruning function
#'   prune_at <- function(thresh){
#'     keep <- which(branch_tbl$score >= thresh)
#'     out <- matrix(FALSE, nrow(mat), ncol(mat))
#'     for(i in keep){
#'       out[cbind(branches[[i]][,1], branches[[i]][,2])] <- TRUE
#'     }
#'     terra::rast(out, extent=terra::ext(skl))
#'   }
#'   
#'   # --- Return object
#'   structure(
#'     list(
#'       branches = branches,
#'       branch_tbl = branch_tbl,
#'       prune = prune_at
#'     ),
#'     class = "skeletonPruner"
#'   )
#' }
#' 
#' 
#' 
#' img1 = load_flexible_image(seg_Oulanka2023_Session01_T067, output_format = "spatrast",select.layer = 2)
#' img2 = load_flexible_image(seg_Oulanka2023_Session01_T067, output_format = "cimg",select.layer = 2)
#' 
#' 
#' 
#'  path = "C:/Users/jocu0013/Desktop/Rjanka diameter testing/test.png"
#'  img = png::readPNG(path)
#'  img1 = load_flexible_image(img, output_format = "spatrast",select.layer = 2, binarize = T)
#'  img2 = load_flexible_image(img, output_format = "cimg",select.layer = 2, binarize = T)
#' 
#' skl = skeletonize_image(img1, method = "MAT")
#' distmap = imager::distance_transform(img2,value = 1)
#' #distmap = load_flexible_image(distmap, output_format = "spatrast", normalize = F)
#' 
#' # Precompute once
#' prn <- skeleton_pruner(load_flexible_image(skl, output_format = "spatrast"), load_flexible_image(distmap, output_format = "spatrast"))
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' ###################
#' 
#' 
#' # Pure R Implementation of Skeleton Pruning using DSE method
#' # No external dependencies required
#' 
#' # Helper function to get 8-connected neighbors
#' get_neighbors_8 <- function(r, c, nrows, ncols) {
#'   neighbors <- list()
#'   for (dr in -1:1) {
#'     for (dc in -1:1) {
#'       if (dr == 0 && dc == 0) next
#'       nr <- r + dr
#'       nc <- c + dc
#'       if (nr >= 1 && nr <= nrows && nc >= 1 && nc <= ncols) {
#'         neighbors[[length(neighbors) + 1]] <- c(nr, nc)
#'       }
#'     }
#'   }
#'   return(neighbors)
#' }
#' 
#' # Create distance transform using simple approach
#' distance_transform <- function(binary_img) {
#'   # Simple distance transform using iterative dilation
#'   dist <- array(0, dim = dim(binary_img))
#'   dist[binary_img == 0] <- -1  # Background
#'   
#'   # Find foreground pixels
#'   fg_pixels <- which(binary_img == 1, arr.ind = TRUE)
#'   if (nrow(fg_pixels) == 0) return(dist)
#'   
#'   # Initialize queue with boundary pixels
#'   queue <- list()
#'   for (i in 1:nrow(fg_pixels)) {
#'     r <- fg_pixels[i, 1]
#'     c <- fg_pixels[i, 2]
#'     neighbors <- get_neighbors_8(r, c, nrow(binary_img), ncol(binary_img))
#'     
#'     # Check if this is a boundary pixel
#'     is_boundary <- FALSE
#'     for (nb in neighbors) {
#'       if (binary_img[nb[1], nb[2]] == 0) {
#'         is_boundary <- TRUE
#'         break
#'       }
#'     }
#'     
#'     if (is_boundary) {
#'       dist[r, c] <- 1
#'       queue[[length(queue) + 1]] <- c(r, c)
#'     } else {
#'       dist[r, c] <- Inf
#'     }
#'   }
#'   
#'   # Propagate distances
#'   while (length(queue) > 0) {
#'     current <- queue[[1]]
#'     queue <- queue[-1]
#'     
#'     r <- current[1]
#'     c <- current[2]
#'     current_dist <- dist[r, c]
#'     
#'     neighbors <- get_neighbors_8(r, c, nrow(binary_img), ncol(binary_img))
#'     for (nb in neighbors) {
#'       nr <- nb[1]
#'       nc <- nb[2]
#'       if (binary_img[nr, nc] == 1 && dist[nr, nc] > current_dist + 1) {
#'         dist[nr, nc] <- current_dist + 1
#'         queue[[length(queue) + 1]] <- c(nr, nc)
#'       }
#'     }
#'   }
#'   
#'   return(dist)
#' }
#' 
#' # Create disk-shaped structuring element
#' make_disk <- function(radius) {
#'   size <- 2 * radius + 1
#'   disk <- array(FALSE, dim = c(size, size))
#'   center <- radius + 1
#'   
#'   for (i in 1:size) {
#'     for (j in 1:size) {
#'       if (sqrt((i - center)^2 + (j - center)^2) <= radius) {
#'         disk[i, j] <- TRUE
#'       }
#'     }
#'   }
#'   return(disk)
#' }
#' 
#' # Morphological dilation at a specific point
#' dilate_at_point <- function(img, r, c, radius) {
#'   if (radius <= 0) return(img)
#'   
#'   disk <- make_disk(radius)
#'   disk_size <- nrow(disk)
#'   offset <- (disk_size - 1) / 2
#'   
#'   for (i in 1:disk_size) {
#'     for (j in 1:disk_size) {
#'       if (disk[i, j]) {
#'         nr <- r + (i - offset - 1)
#'         nc <- c + (j - offset - 1)
#'         if (nr >= 1 && nr <= nrow(img) && nc >= 1 && nc <= ncol(img)) {
#'           img[nr, nc] <- 1
#'         }
#'       }
#'     }
#'   }
#'   return(img)
#' }
#' 
#' # Simple graph representation using adjacency lists
#' create_graph <- function() {
#'   list(
#'     nodes = list(),      # node_id -> list(coords, degree)
#'     edges = list(),      # edge_id -> list(from, to, pts)
#'     node_count = 0,
#'     edge_count = 0
#'   )
#' }
#' 
#' add_node <- function(graph, coords) {
#'   graph$node_count <- graph$node_count + 1
#'   node_id <- graph$node_count
#'   graph$nodes[[node_id]] <- list(coords = coords, degree = 0)
#'   return(list(graph = graph, node_id = node_id))
#' }
#' 
#' add_edge <- function(graph, from_id, to_id, pts) {
#'   graph$edge_count <- graph$edge_count + 1
#'   edge_id <- graph$edge_count
#'   graph$edges[[edge_id]] <- list(from = from_id, to = to_id, pts = pts)
#'   
#'   # Update degrees
#'   graph$nodes[[from_id]]$degree <- graph$nodes[[from_id]]$degree + 1
#'   graph$nodes[[to_id]]$degree <- graph$nodes[[to_id]]$degree + 1
#'   
#'   return(graph)
#' }
#' 
#' remove_node <- function(graph, node_id) {
#'   if (is.null(graph$nodes[[node_id]])) return(graph)
#'   
#'   # Remove all edges connected to this node
#'   edges_to_remove <- c()
#'   for (edge_id in names(graph$edges)) {
#'     edge <- graph$edges[[edge_id]]
#'     if (edge$from == node_id || edge$to == node_id) {
#'       edges_to_remove <- c(edges_to_remove, edge_id)
#'     }
#'   }
#'   
#'   # Remove edges and update degrees
#'   for (edge_id in edges_to_remove) {
#'     edge <- graph$edges[[edge_id]]
#'     other_node <- if (edge$from == node_id) edge$to else edge$from
#'     if (!is.null(graph$nodes[[other_node]])) {
#'       graph$nodes[[other_node]]$degree <- graph$nodes[[other_node]]$degree - 1
#'     }
#'     graph$edges[[edge_id]] <- NULL
#'   }
#'   
#'   # Remove node
#'   graph$nodes[[node_id]] <- NULL
#'   return(graph)
#' }
#' 
#' # Convert skeleton to simple graph
#' skeleton_to_graph <- function(skeleton) {
#'   graph <- create_graph()
#'   
#'   # Find skeleton pixels
#'   skel_pixels <- which(skeleton == 1, arr.ind = TRUE)
#'   if (nrow(skel_pixels) == 0) return(graph)
#'   
#'   # Simple approach: treat each skeleton pixel as a potential node
#'   # and connect to neighboring skeleton pixels
#'   pixel_to_node <- array(0, dim = dim(skeleton))
#'   
#'   # Create nodes for junction points and endpoints
#'   for (i in 1:nrow(skel_pixels)) {
#'     r <- skel_pixels[i, 1]
#'     c <- skel_pixels[i, 2]
#'     
#'     # Count skeleton neighbors
#'     neighbors <- get_neighbors_8(r, c, nrow(skeleton), ncol(skeleton))
#'     skel_neighbor_count <- 0
#'     for (nb in neighbors) {
#'       if (skeleton[nb[1], nb[2]] == 1) {
#'         skel_neighbor_count <- skel_neighbor_count + 1
#'       }
#'     }
#'     
#'     # Junction (degree > 2) or endpoint (degree 1)
#'     if (skel_neighbor_count != 2) {
#'       result <- add_node(graph, c(r, c))
#'       graph <- result$graph
#'       pixel_to_node[r, c] <- result$node_id
#'     }
#'   }
#'   
#'   return(graph)
#' }
#' 
#' # Reconstruct branch using disk structuring elements
#' reconstruct_by_disk <- function(pts, dist_map, output_mask) {
#'   if (nrow(pts) == 0) return(output_mask)
#'   
#'   for (i in 1:nrow(pts)) {
#'     r <- pts[i, 1]
#'     c <- pts[i, 2]
#'     if (r >= 1 && r <= nrow(dist_map) && c >= 1 && c <= ncol(dist_map)) {
#'       radius <- dist_map[r, c]
#'       if (radius > 0) {
#'         output_mask <- dilate_at_point(output_mask, r, c, radius)
#'       }
#'     }
#'   }
#'   return(output_mask)
#' }
#' 
#' # Calculate weight of reconstruction
#' get_weight <- function(original, reconstructed) {
#'   return(sum(reconstructed * original))
#' }
#' 
#' # Remove terminal branches using DSE
#' remove_branch_by_DSE <- function(graph, reconstruction, dist_map, max_weight) {
#'   if (length(graph$nodes) == 0) return(list(graph = graph, reconstruction = reconstruction))
#'   
#'   # Find terminal nodes (degree 1)
#'   terminal_nodes <- c()
#'   for (node_id in names(graph$nodes)) {
#'     if (graph$nodes[[node_id]]$degree == 1) {
#'       terminal_nodes <- c(terminal_nodes, as.numeric(node_id))
#'     }
#'   }
#'   
#'   nodes_to_remove <- c()
#'   
#'   # Check each terminal node
#'   for (node_id in terminal_nodes) {
#'     if (is.null(graph$nodes[[as.character(node_id)]])) next
#'     
#'     node_coords <- graph$nodes[[as.character(node_id)]]$coords
#'     
#'     # Find the edge connected to this terminal node
#'     connected_edge <- NULL
#'     for (edge_id in names(graph$edges)) {
#'       edge <- graph$edges[[edge_id]]
#'       if (edge$from == node_id || edge$to == node_id) {
#'         connected_edge <- edge
#'         break
#'       }
#'     }
#'     
#'     if (is.null(connected_edge)) next
#'     
#'     # Get other node coordinates
#'     other_node_id <- if (connected_edge$from == node_id) connected_edge$to else connected_edge$from
#'     if (is.null(graph$nodes[[as.character(other_node_id)]])) next
#'     
#'     other_coords <- graph$nodes[[as.character(other_node_id)]]$coords
#'     
#'     # Combine all points for reconstruction
#'     all_pts <- rbind(node_coords, other_coords)
#'     if (!is.null(connected_edge$pts) && nrow(connected_edge$pts) > 0) {
#'       all_pts <- rbind(all_pts, connected_edge$pts)
#'     }
#'     
#'     # Create branch reconstruction
#'     branch_reconstruction <- array(0, dim = dim(reconstruction))
#'     branch_reconstruction <- reconstruct_by_disk(all_pts, dist_map, branch_reconstruction)
#'     
#'     # Calculate weight
#'     weight <- get_weight(reconstruction, branch_reconstruction)
#'     
#'     # Remove if weight is below threshold
#'     if (weight < max_weight) {
#'       nodes_to_remove <- c(nodes_to_remove, node_id)
#'       reconstruction <- reconstruction - branch_reconstruction
#'       reconstruction[reconstruction < 0] <- 0
#'     }
#'   }
#'   
#'   # Remove nodes
#'   for (node_id in nodes_to_remove) {
#'     graph <- remove_node(graph, node_id)
#'   }
#'   
#'   return(list(graph = graph, reconstruction = reconstruction))
#' }
#' 
#' # Convert graph back to skeleton image
#' graph_to_image <- function(graph, image_shape) {
#'   mask <- array(FALSE, dim = image_shape)
#'   
#'   # Draw nodes
#'   for (node_id in names(graph$nodes)) {
#'     coords <- graph$nodes[[node_id]]$coords
#'     r <- coords[1]
#'     c <- coords[2]
#'     if (r >= 1 && r <= image_shape[1] && c >= 1 && c <= image_shape[2]) {
#'       mask[r, c] <- TRUE
#'     }
#'   }
#'   
#'   # Draw edges
#'   for (edge_id in names(graph$edges)) {
#'     edge <- graph$edges[[edge_id]]
#'     if (!is.null(edge$pts) && nrow(edge$pts) > 0) {
#'       for (i in 1:nrow(edge$pts)) {
#'         r <- edge$pts[i, 1]
#'         c <- edge$pts[i, 2]
#'         if (r >= 1 && r <= image_shape[1] && c >= 1 && c <= image_shape[2]) {
#'           mask[r, c] <- TRUE
#'         }
#'       }
#'     }
#'   }
#'   
#'   return(mask)
#' }
#' 
#' # Main skeleton pruning function
#' skel_pruning_DSE <- function(skeleton, distance_map = NULL, min_area_px = 100, return_graph = FALSE) {
#'   # """
#'   # Skeleton pruning using DSE method - Pure R implementation
#'   # 
#'   # Args:
#'   #   skeleton: Binary skeleton image (matrix of 0s and 1s)
#'   #   distance_map: Distance transform of original binary image (optional)
#'   #   min_area_px: Minimum branch area threshold
#'   #   return_graph: Whether to return the graph object
#'   # 
#'   # Returns:
#'   #   Pruned skeleton image (and optionally the graph)
#'   # """
#'   
#'   # Create distance map if not provided
#'   if (is.null(distance_map)) {
#'     # Assume skeleton comes from a binary image
#'     distance_map <- distance_transform(skeleton)
#'   }
#'   
#'   # Convert skeleton to graph
#'   graph <- skeleton_to_graph(skeleton)
#'   
#'   if (length(graph$nodes) == 0) {
#'     if (return_graph) {
#'       return(list(skeleton = skeleton, graph = graph))
#'     } else {
#'       return(skeleton)
#'     }
#'   }
#'   
#'   # Initial reconstruction (simplified)
#'   reconstruction <- skeleton
#'   
#'   # Iteratively remove branches
#'   max_iterations <- 10
#'   for (iteration in 1:max_iterations) {
#'     old_node_count <- length(graph$nodes)
#'     
#'     # Remove terminal branches
#'     result <- remove_branch_by_DSE(graph, reconstruction, distance_map, min_area_px)
#'     graph <- result$graph
#'     reconstruction <- result$reconstruction
#'     
#'     # Stop if no changes
#'     if (length(graph$nodes) == old_node_count) break
#'   }
#'   
#'   # Convert back to skeleton image
#'   pruned_skeleton <- graph_to_image(graph, dim(skeleton))
#'   
#'   if (return_graph) {
#'     return(list(skeleton = pruned_skeleton, graph = graph))
#'   } else {
#'     return(pruned_skeleton)
#'   }
#' }
#' 
#' # Example usage and test function
#' test_skeleton_pruning <- function() {
#'   # Create a simple test skeleton with branches
#'   skeleton <- array(0, dim = c(20, 20))
#'   
#'   # Main branch (horizontal line)
#'   skeleton[10, 5:15] <- 1
#'   
#'   # Small branches
#'   skeleton[8:10, 8] <- 1  # Short branch
#'   skeleton[10:12, 12] <- 1  # Short branch
#'   skeleton[7:10, 15] <- 1  # Longer branch
#'   
#'   cat("Original skeleton has", sum(skeleton), "pixels\n")
#'   
#'   # Prune skeleton
#'   pruned <- skel_pruning_DSE(skeleton, min_area_px = 1)
#'   
#'   cat("Pruned skeleton has", sum(pruned), "pixels\n")
#'   
#'   return(list(original = skeleton, pruned = pruned))
#' }
#' 
#' # Run test
#' # result <- test_skeleton_pruning()
#' 
#' 
#' result = skel_pruning_DSE(load_flexible_image( skl, output_format = "array"), min_area_px = 6)
#' plot(load_flexible_image(result, output_format = "spatrast",binarize = T))
#' plot(skl- load_flexible_image(result, output_format = "spatrast", binarize = T))
