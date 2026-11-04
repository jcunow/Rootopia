#' @title Original Minirhizotron Root Scan - Session 3, Tube 67
#' @description
#' Original RGB root scan image from a sedge fen in northern Finland (October 2023).
#' This image represents a composite of multiple stitched scans.
#'
#' @details
#' The image represents a complete tube scan where:
#' \itemize{
#'   \item Columns correspond to tube length
#'   \item Rows correspond to tube rotation
#'   \item RGB channels represent true color information
#' }
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @format A numeric array (rows x columns x layers), values 0-255. Rebuild a
#'   SpatRaster with \code{terra::rast()}. Dimensions:
#'   \itemize{
#'     \item 1160 rows (height)
#'     \item 4899 columns (width)
#'     \item 3 layers (RGB channels)
#'   }
#' @usage data(rgb_Oulanka2023_Session03_T067)
#' @source
#' Original minirhizotron images acquired by Johannes Cunow and Gesche Blume-Werry
#' at Oulanka Research Station, Finland, in October 2023.
#' @examples
#' \dontrun{
#'   data(rgb_Oulanka2023_Session03_T067)
#'   rgb_Oulanka2023_Session03_T067 = terra::rast(rgb_Oulanka2023_Session03_T067)
#'   terra::plotRGB(rgb_Oulanka2023_Session03_T067, maxcell = Inf)
#' }
"rgb_Oulanka2023_Session03_T067"

#' @title Segmented Minirhizotron Root Scan - Session 1, Tube 67
#' @description
#' Segmented root scan image from a sedge fen in northern Finland (June 2023).
#' The image was processed using RootDetector for root segmentation.
#'
#' @details
#' Binary mask, stored as 0 (background) and 255 (foreground):
#' \itemize{
#'   \item Layer 1: foreign objects such as tape
#'   \item Layers 2 and 3 (identical): roots
#' }
#' Spatial dimensions correspond to physical tube measurements:
#' columns = tube length, rows = tube rotation.
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @format A numeric array (rows x columns x layers), values 0-255. Rebuild a
#'   SpatRaster with \code{terra::rast()}. Dimensions:
#'   \itemize{
#'     \item 4900 columns (width)
#'     \item 1144 rows (height)
#'     \item 3 layers (channels)
#'   }
#' @usage data(seg_Oulanka2023_Session01_T067)
#' @source
#' Original minirhizotron images acquired by Johannes Cunow and Gesche Blume-Werry
#' at Oulanka Research Station, Finland, in June 2023.
#' @examples
#' \dontrun{
#'   data(seg_Oulanka2023_Session01_T067)
#'   seg_Oulanka2023_Session01_T067 = terra::rast(seg_Oulanka2023_Session01_T067)
#'   terra::plot(seg_Oulanka2023_Session01_T067, maxcell = Inf)
#' }
"seg_Oulanka2023_Session01_T067"

#' @title Segmented Minirhizotron Root Scan - Session 3, Tube 67
#' @description
#' Segmented root scan image from a sedge fen in northern Finland (October 2023).
#' The image was processed using RootDetector for root segmentation.
#'
#' @details
#' Binary mask, stored as 0 (background) and 255 (foreground):
#' \itemize{
#'   \item Layer 1: foreign objects such as tape
#'   \item Layers 2 and 3 (identical): roots
#' }
#' Spatial dimensions correspond to physical tube measurements:
#' columns = tube length, rows = tube rotation.
#'
#' @format A numeric array (rows x columns x layers), values 0-255. Rebuild a
#'   SpatRaster with \code{terra::rast()}. Dimensions:
#'   \itemize{
#'     \item 4900 columns (width)
#'     \item 1161 rows (height)
#'     \item 3 layers (channels)
#'   }
#' @usage data(seg_Oulanka2023_Session03_T067)
#' @author Johannes Cunow
#' @source
#' Original minirhizotron images acquired by Johannes Cunow
#' and Gesche Blume-Werry at Oulanka Research Station, Finland (October 2023).
#' Images were processed using RootDetector.
#'
#' @references
#' Peters B, Blume-Werry G, Gillert A, et al. (2023)
#' RootDetector: a convolutional neural network for root detection.
#' Scientific Reports 13:1399. https://doi.org/10.1038/s41598-023-28400-x
#' @examples
#' \dontrun{
#'   data(seg_Oulanka2023_Session03_T067)
#'   seg_Oulanka2023_Session03_T067 = terra::rast(seg_Oulanka2023_Session03_T067)
#'   terra::plot(seg_Oulanka2023_Session03_T067, maxcell = Inf)
#' }
"seg_Oulanka2023_Session03_T067"

#' @title Skeletonized Root Scan - Session 1, Tube 67
#' @description
#' Skeletonized root scan image from a sedge fen in northern Finland (June 2023).
#' The image was processed using RootDetector for segmentation and skeletonization.
#'
#' @details
#' Binary mask, stored as 0 (background) and a single non-zero value
#' (foreground):
#' \itemize{
#'   \item Layer 1: foreign objects such as tape, filled
#'   \item Layers 2 and 3 (identical): root skeleton
#' }
#' At the stored resolution the skeleton's one-pixel lines are broken into
#' many short pieces. For traits
#' that depend on connectivity (tips, branch points, branching order),
#' skeletonize the matching segmentation with \code{skeletonize_image()}.
#'
#' @format A numeric array (rows x columns x layers), values 0-255. Rebuild a
#'   SpatRaster with \code{terra::rast()}. Dimensions:
#'   \itemize{
#'     \item 4900 columns (width)
#'     \item 1144 rows (height)
#'     \item 3 layers (channels)
#'   }
#' @usage data(skl_Oulanka2023_Session01_T067)
#' @author Johannes Cunow
#' @source
#' Original minirhizotron images acquired by Johannes Cunow
#' and Gesche Blume-Werry at Oulanka Research Station, Finland (June 2023).
#' Images were processed using RootDetector.
#'
#' @references
#' Peters B, Blume-Werry G, Gillert A, et al. (2023)
#' RootDetector: a convolutional neural network for root detection.
#' Scientific Reports 13:1399. https://doi.org/10.1038/s41598-023-28400-x
#' @examples
#' \dontrun{
#'   data(skl_Oulanka2023_Session01_T067)
#'   skl_Oulanka2023_Session01_T067 = terra::rast(skl_Oulanka2023_Session01_T067)
#'   terra::plot(skl_Oulanka2023_Session01_T067, maxcell = Inf)
#' }
"skl_Oulanka2023_Session01_T067"

#' @title Skeletonized Root Scan - Session 3, Tube 67
#' @description
#' Skeletonized root scan image from a sedge fen in northern Finland (October 2023).
#' The image was processed using RootDetector for segmentation and skeletonization.
#'
#' @details
#' Binary mask, stored as 0 (background) and a single non-zero value
#' (foreground):
#' \itemize{
#'   \item Layer 1: foreign objects such as tape, filled
#'   \item Layers 2 and 3 (identical): root skeleton
#' }
#' At the stored resolution the skeleton's one-pixel lines are broken into
#' many short pieces. For traits
#' that depend on connectivity (tips, branch points, branching order),
#' skeletonize the matching segmentation with \code{skeletonize_image()}.
#'
#' @format A numeric array (rows x columns x layers), values 0-255. Rebuild a
#'   SpatRaster with \code{terra::rast()}. Dimensions:
#'   \itemize{
#'     \item 4900 columns (width)
#'     \item 1161 rows (height)
#'     \item 3 layers (channels)
#'   }
#' @usage data(skl_Oulanka2023_Session03_T067)
#' @author Johannes Cunow
#' @source
#' Original minirhizotron images acquired by Johannes Cunow
#' and Gesche Blume-Werry at Oulanka Research Station, Finland (October 2023).
#' Images were processed using RootDetector.
#'
#' @references
#' Peters B, Blume-Werry G, Gillert A, et al. (2023)
#' RootDetector: a convolutional neural network for root detection.
#' Scientific Reports 13:1399. https://doi.org/10.1038/s41598-023-28400-x
#' @examples
#' \dontrun{
#'   data(skl_Oulanka2023_Session03_T067)
#'   skl_Oulanka2023_Session03_T067 = terra::rast(skl_Oulanka2023_Session03_T067)
#'   terra::plot(skl_Oulanka2023_Session03_T067, maxcell = Inf)
#' }
"skl_Oulanka2023_Session03_T067"

#' @title Root Turnover Analysis Data
#' @description
#' Root turnover analysis from a sedge fen in northern Finland, comparing
#' root presence between two time points (June 2023 and October 2023) using
#' RootDetector root tracking feature.
#'
#' @details
#' An RGB image in which each turnover class has its own color:
#' \tabular{lll}{
#'   \strong{Class} \tab \strong{Color} \tab \strong{R, G, B} \cr
#'   background \tab dark slate \tab 39, 54, 59 \cr
#'   foreign objects (e.g. tape) \tab red \tab 255, 0, 0 \cr
#'   constant (present at both dates) \tab white \tab 255, 255, 255 \cr
#'   production (new at the second date) \tab green \tab 96, 209, 130 \cr
#'   decay (gone at the second date) \tab rose \tab 226, 106, 116
#' }
#' \code{root_turnover(method = "dpc")} reads these classes with its default
#' \code{decay_layer = 1} and \code{product_layer = 2}.
#'
#' @format A numeric array (rows x columns x layers), values 0-255. Rebuild a
#'   SpatRaster with \code{terra::rast()}. Dimensions:
#'   \itemize{
#'     \item 2273 rows (height)
#'     \item 2550 columns (width)
#'     \item 3 layers (red, green, blue)
#'   }
#' @usage data(TurnoverDPC_data)
#' @author Johannes Cunow
#' @source
#' Original minirhizotron images acquired by Johannes Cunow
#' and Gesche Blume-Werry at Oulanka Research Station, Finland (June and
#' October 2023). Images were processed using RootDetector.
#'
#' @references
#' Peters B, Blume-Werry G, Gillert A, et al. (2023)
#' RootDetector: a convolutional neural network for root detection.
#' Scientific Reports 13:1399. https://doi.org/10.1038/s41598-023-28400-x
#' @examples
#' \dontrun{
#'   data(TurnoverDPC_data)
#'   img <- terra::rast(TurnoverDPC_data)
#'   terra::plotRGB(img, maxcell = Inf)
#'   root_turnover(img, method = "dpc")
#' }
"TurnoverDPC_data"

#' @title Example Flatbed Root Scan (segmented)
#' @description
#' A segmented flatbed root scan used to demonstrate the depth-independent
#' flatbed workflow (skeletonization, root length, diameter, branching, and
#' branch order). Roots are foreground, background is 0.
#'
#' @details
#' Flatbed scans capture the full root system in a 2D plane, so traits are
#' computed globally rather than per depth bin. The image is a binary
#' segmentation (root vs. background) stored across three channels; layer 2 is
#' used as the root channel in the vignette.
#'
#' @author Johannes Cunow \email{johannes.cunow@gmail.com}
#' @format A 3-dimensional numeric array (rows x columns x layers), matching the
#'   other bundled datasets. Rebuild a SpatRaster with \code{terra::rast()}.
#'   Dimensions:
#'   \itemize{
#'     \item 2000 rows (height)
#'     \item 2000 columns (width)
#'     \item 3 layers (segmentation channels)
#'   }
#' @usage data(flatbed_scan_example)
#' @source
#' Original flatbed scans from boreal pine forests in northern Finland, acquired by
#' Johannes Cunow and others (field campaign 2022).
#'
#' @references
#' Cunow J, Pijcke F, Olofsson J, Väisänen M, Blume-Werry G (accepted, 2026).
#' Reindeer grazing induces spatial and functional shifts in root systems
#' of boreal pine forests. Oikos. https://doi.org/10.1002/oik.12211
#' 
#' 
#' @examples
#' \dontrun{
#'   data(flatbed_scan_example)
#'   seg <- terra::rast(flatbed_scan_example)
#'   seg <- terra::ifel(seg[[2]] > 0, 1, 0)
#'   terra::plot(seg, maxcell = Inf)
#' }
"flatbed_scan_example"
