
###########
#' Calculate global root production and root turnover from temporal comparison
#'
#' @param im.t1 SpatRaster object for timepoint 1
#' @param im.t2 SpatRaster object for timepoint 2
#' @param method Analysis method: "kimura" or "rootpx"
#' @param unit Unit of root length measurement (only for method = "kimura"). Default: "cm"
#' @param dpi Image resolution (only for method = "kimura"). Default: 300
#' @param select_layer Integer specifying the layer to use in both timesteps if \code{img} is a multi-layer `SpatRaster`. Defaults to 2.
#'
#' @return data.frame containing:
#'   - standingroot_t1: Standing roots at first timepoint
#'   - standingroot_t2: Standing roots at second timepoint
#'   - production: Root production between timepoints
#'   - newroot%per_t1: Percentage of new roots compared to starting conditions
#'   - newroot%per_t2: Percentage of new roots at second timepoint
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'   data(skl_Oulanka2023_Session01_T067)
#'   data(skl_Oulanka2023_Session03_T067)
#'   time1 <- terra::rast(skl_Oulanka2023_Session01_T067)
#'   time2 <- terra::rast(skl_Oulanka2023_Session03_T067)
#'   turnover.values <- turnover_tc(
#'     im.t1 = time1,
#'     im.t2 = time2,
#'     method = "kimura")
#'     }
turnover_tc = function(im.t1, im.t2, method="kimura", unit="cm", dpi=300, select_layer = 2) {

  # Input validation module
  validate_inputs <- function(im.t1,im.t2) {
    # Check for missing required arguments
    if (missing(im.t1) || missing(im.t2)) {
      stop("Both im.t1 and im.t2 must be provided")
    }

    # Validate images
    if (!inherits(im.t1, "SpatRaster") || !inherits(im.t2, "SpatRaster")) {
      stop("Both im.t1 and im.t2 must be SpatRaster objects")
    }

    # Check if images have similar dimensions & amount of cells (within +-5%)
    if (!all(
      c(dim(im.t1) >= dim(im.t2)*0.95 | dim(im.t1) >= dim(im.t2)*1.05) &
      c((terra::ncell(im.t1) / terra::ncell(im.t2)) >=0.95 | (terra::ncell(im.t1) / terra::ncell(im.t2)) <=1.05))
    ) {
      stop("Images must have the similar dimensions within & amount of cells (+-5%) ")
    }

    # Validate method
    if (!method %in% c("kimura", "rootpx")) {
      stop("Method must be either 'kimura' or 'rootpx'")
    }

    # Validate kimura-specific parameters
    if (method == "kimura") {
      if (!is.character(unit)) {
        stop("unit must be a character string")
      }
      if (!is.numeric(dpi) || dpi <= 0) {
        stop("dpi must be a positive numeric value")
      }
    }
  }

  # Edge case handling module
  handle_edge_cases <- function(px1, px2) {
    # Check for zero or negative values
    if (px1 <= 0) {
      warning("Zero or negative root length at timepoint 1")
      return(NULL)
    }

    # Handle division by zero cases
    if (px1 == 0) px.turn1 <- NA else px.turn1 <- round(px.prod / px1, 4)
    if (px2 == 0) px.turn2 <- NA else px.turn2 <- round(px.prod / px2, 4)

    return(list(px.turn1 = px.turn1, px.turn2 = px.turn2))
  }

  # Main execution with error handling
  tryCatch({
    im.t1 <- load_flexible_image(im.t1,select_layer = select_layer, scale = "to_01", output_format = "spatrast")
    im.t2 <- load_flexible_image(im.t2,select_layer = select_layer, scale = "to_01", output_format = "spatrast")
    im.t1 = ceiling(im.t1)
    im.t2 = ceiling(im.t2)
    validate_inputs(im.t1,im.t2)

    if (method == "rootpx") {
      px1 <- count_pixels(im.t1)
      px2 <- count_pixels(im.t2)
    } else {  # method == "kimura"
      px1 <- root_length(im.t1, select_layer = NULL, dpi = dpi, unit = unit)
      px2 <- root_length(im.t2, select_layer = NULL, dpi = dpi, unit = unit)
    }

    px.prod <- px2 - px1

    # Handle edge cases
    results <- handle_edge_cases(px1, px2)
    if (is.null(results)) {
      return(NULL)
    }

    return(data.frame(
      "standingroot_t1" = px1,
      "standingroot_t2" = px2,
      "production" = px.prod,
      "newroot%per_t1" = results$px.turn1,
      "newroot%per_t2" = results$px.turn2
    ))

  }, error = function(e) {
    stop(paste("Error in Turnover.TC:", e$message))
  })
}

#' Extract Root Decay, New Root Production, and No-Change Roots (only 'RootDetector' images)
#'
#' @param img SpatRaster with three layers for production, decay, and stagnation
#' @param product_layer Integer indicating the production layer index (1-3)
#' @param decay_layer Integer indicating the decay & tape layer index (1-3)
#' @param blur_capture Threshold for pixel inclusion (0-1). Default: 0.95
#' @param im_return Logical: return images instead of values? Default: FALSE
#' @param include_virtualroots Logical: consider all roots present at any timepoint? Default: FALSE
#'
#' @return If im_return = FALSE: tibble with pixel sums and ratios
#'         If im_return = TRUE: list of SpatRaster layers for tape, constant, production, and decay
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' data(TurnoverDPC_data)
#' img = terra::rast(TurnoverDPC_data)
#' DPCs = turnover_dpc(img = img, im_return = FALSE)
#' }
turnover_dpc = function(img, product_layer=2, decay_layer=1, blur_capture=0.95,
                        im_return=FALSE, include_virtualroots=FALSE) {
  # Input validation module
  validate_inputs <- function(img) {
    if (missing(img)) {
      stop("Input image must be provided")
    }

    # Validate layer indices
    if (!all(c(product_layer, decay_layer) %in% 1:3)) {
      stop("Layer indices must be between 1 and 3")
    }
    if (product_layer == decay_layer) {
      stop("Product and decay layers must be different")
    }

    # Validate blur_capture
    if (!is.numeric(blur_capture) || blur_capture < 0 || blur_capture > 1) {
      stop("blur_capture must be between 0 and 1")
    }

    # Validate logical parameters
    if (!is.logical(im_return) || !is.logical(include_virtualroots)) {
      stop("im_return and include_virtualroots must be logical values")
    }
  }

  # Edge case handling module
  handle_edge_cases <- function(img) {
    # Check if image has exactly 3 layers
    if (terra::nlyr(img) < 3) {
      stop("Input image must at least contain 3 layers")
    }


    # Check for all-NA layers
    if (any(sapply(1:3, function(i) all(is.na(terra::values(img[[i]])))))) {
      warning("One or more layers contain only NA values")
      return(NULL)
    }

    return(TRUE)
  }

  # Main execution with error handling
  tryCatch({
    validate_inputs(img)

    # Load and validate image
    img <- load_flexible_image(img, output_format = "spatrast", scale = "to_01", select_layer = NULL)
    if (is.null(handle_edge_cases(img))) return(NULL)

    # Calculate indices and layers
    l.indx <- 1:3
    no.change.layer <- which(!l.indx %in% c(product_layer, decay_layer))
    l.pr <- img[[product_layer]]
    l.no <- img[[no.change.layer]]
    l.dc <- img[[decay_layer]]

    # Safe calculations with error checking
    tape <- suppressWarnings({
      temp <- l.no - l.dc
      min_val <- terra::global(temp, "min", na.rm=TRUE)[[1]]
      if (is.na(min_val)) stop("Unable to calculate minimum value for tape layer")
      temp <= min_val * blur_capture
    })

    tape <- tape * terra::global(l.pr, "max", na.rm=TRUE)[[1]]

    # Helper function for layer calculations
    calculate_layer <- function(layer, no_change, tape, stat, blur_capture) {
      result <- layer - no_change - tape
      max_val <- terra::global(result, stat, na.rm=TRUE)[[1]]
      if (is.na(max_val)) stop("Unable to calculate maximum value for layer")
      (result >= max_val * blur_capture) * 1
    }
    # Calculate production, decay, and no-change layers with error checking
    l.pr2 <- calculate_layer(l.pr, l.no, tape, "max", blur_capture)
    l.dc2 <- calculate_layer(l.dc, l.no, tape, "max", blur_capture)
    l.no2 <- l.no >= terra::global(l.no, "max", na.rm=TRUE)[[1]] * blur_capture

    if (im_return) {
      return(list("tape"=tape, "constant"=l.no2, "production"=l.pr2, "decay"=l.dc2))
    } else {

      # Helper function for ratio calculations
      calculate_ratios <- function(tape.px, const.px, prodc.px, decay.px, include_virtualroots) {
        denominator <- prodc.px + const.px + decay.px
        if (denominator == 0) {
          warning("Zero denominator in ratio calculations")
          return(data.frame(
            tape = tape.px,
            constant = const.px,
            production = prodc.px,
            decay = decay.px,
            newgrowth.ratio = NA,
            decay.ratio = NA,
            constant.ratio = NA
          ))
        }

        if (include_virtualroots) {
          newgrowth.ratio <- prodc.px / denominator
          decay.ratio <- decay.px / denominator
        } else {
          newgrowth.ratio <- prodc.px / (prodc.px + const.px)
          decay.ratio <- decay.px / (decay.px + const.px)
        }

        constant.ratio <- const.px / denominator

        data.frame(
          tape = tape.px,
          constant = const.px,
          production = prodc.px,
          decay = decay.px,
          newgrowth.ratio = round(newgrowth.ratio, 4),
          decay.ratio = round(decay.ratio, 4),
          constant.ratio = round(constant.ratio, 4)
        )
      }

      # Calculate pixel sums and ratios
      results <- calculate_ratios(
        tape.px = count_pixels(tape),
        const.px = count_pixels(l.no2),
        prodc.px = count_pixels(l.pr2),
        decay.px = count_pixels(l.dc2),
        include_virtualroots = include_virtualroots
      )
      return(results)
    }

  }, error = function(e) {
    stop(paste("Error in Turnover.DPC:", e$message))
  })
}




#' Unified Root Turnover Analysis
#'
#' Wrapper around the two root-turnover methods. \code{method} selects which
#' one runs:
#' \describe{
#'   \item{\code{"tc"} (Temporal Comparison)}{Compares two timepoint images
#'     (\code{img1}, \code{img2}) and reports standing roots, production, and
#'     new-root percentages. Dispatches to \code{\link{turnover_tc}}. The
#'     \code{tc_method} argument chooses how root amount is measured:
#'     \code{"kimura"} (root length) or \code{"rootpx"} (root pixel count).}
#'   \item{\code{"dpc"} (Decay, Production, Constant)}{Decomposes a single
#'     multi-layer 'RootDetector' image into decayed, newly produced, and
#'     unchanged (constant) root fractions. Dispatches to
#'     \code{\link{turnover_dpc}}; \code{img2} is not used.}
#' }
#'
#' @param img1 Primary SpatRaster input. For \code{method = "tc"} this is the
#'   first timepoint image; for \code{method = "dpc"} this is the multi-layer
#'   DPC image.
#' @param img2 Second timepoint SpatRaster. Required for \code{method = "tc"};
#'   ignored (with a warning) for \code{method = "dpc"}.
#' @param method Which turnover method to run: \code{"tc"} (temporal
#'   comparison of two images) or \code{"dpc"} (Decay/Production/Constant
#'   decomposition of one multi-layer image). Default \code{"tc"}.
#' @param tc_method Measurement sub-method for \code{method = "tc"}:
#'   \code{"kimura"} (root length) or \code{"rootpx"} (root pixel count).
#'   Ignored for \code{method = "dpc"}. Default \code{"kimura"}.
#' @param unit Unit of root length measurement (only for
#'   \code{tc_method = "kimura"}). Default: "cm"
#' @param dpi Image resolution (only for \code{tc_method = "kimura"}). Default: 300
#' @param select_layer Integer or NULL. For \code{method = "tc"} with
#'   multi-layer images, selects which layer to compare. Ignored for
#'   \code{method = "dpc"}.
#' @param product_layer Integer indicating the production layer index for the
#'   DPC method (1-3)
#' @param decay_layer Integer indicating the decay & tape layer index for the
#'   DPC method (1-3)
#' @param blur_capture Threshold for pixel inclusion in the DPC method (0-1).
#'   Default: 0.95
#' @param im_return Logical: return images instead of values for the DPC
#'   method? Default: FALSE
#' @param include_virtualroots Logical: consider all roots present at any
#'   timepoint in the DPC method? Default: FALSE
#'
#' @return Depends on the method:
#'   - \code{"tc"}: data.frame with standing roots, production, and new-root percentages.
#'   - \code{"dpc"}: data.frame of pixel sums and ratios, or (if \code{im_return = TRUE}) a list of SpatRaster layers.
#'
#' @seealso \code{\link{turnover_tc}}, \code{\link{turnover_dpc}}
#'
#' @export
#' @examples
#' # DPC: single multi-layer 'RootDetector' image
#' data(TurnoverDPC_data)
#' img <- terra::rast(TurnoverDPC_data)
#' root_turnover(img, method = "dpc")
#'
#' # TC: two timepoint images compared by root length (kimura)
#' data(skl_Oulanka2023_Session01_T067)
#' data(skl_Oulanka2023_Session03_T067)
#' t1 <- terra::rast(skl_Oulanka2023_Session01_T067)
#' t2 <- terra::rast(skl_Oulanka2023_Session03_T067)
#' root_turnover(t1, t2, method = "tc", tc_method = "kimura")
root_turnover = function(img1, img2 = NULL,
                         method = c("tc", "dpc"),
                         tc_method = c("kimura", "rootpx"),
                         unit = "cm",
                         dpi = 300,
                         select_layer = NULL,
                         product_layer = 2,
                         decay_layer = 1,
                         blur_capture = 0.95,
                         im_return = FALSE,
                         include_virtualroots = FALSE) {

  method    <- match.arg(method)
  tc_method <- match.arg(tc_method)

  # Helper function to handle layer selection (TC method only)
  pick_layer <- function(img, select_layer = NULL) {
    # If select_layer is NULL or img has only one layer, return the image as-is
    if (is.null(select_layer) || terra::nlyr(img) == 1) {
      return(img)
    }

    # Validate select_layer
    if (!is.numeric(select_layer) ||
        select_layer < 1 ||
        select_layer > terra::nlyr(img)) {
      stop("Invalid select_layer. Must be an integer between 1 and ", terra::nlyr(img))
    }

    # Return selected layer
    return(img[[select_layer]])
  }

  # Dispatch on method
  if (method == "dpc") {
    # Decay, Production, Constant -- single multi-layer image
    if (!is.null(img2)) {
      warning("img2 is ignored for method = 'dpc' (uses a single multi-layer image)")
    }
    return(turnover_dpc(
      img = img1,
      product_layer = product_layer,
      decay_layer = decay_layer,
      blur_capture = blur_capture,
      im_return = im_return,
      include_virtualroots = include_virtualroots
    ))
  } else {
    # Temporal Comparison -- two timepoint images
    if (is.null(img2)) {
      stop("method = 'tc' (temporal comparison) requires two images: supply img2")
    }
    img1 <- pick_layer(img1, select_layer)
    img2 <- pick_layer(img2, select_layer)

    return(turnover_tc(
      im.t1 = img1,
      im.t2 = img2,
      method = tc_method,
      unit = unit,
      dpi = dpi,
      select_layer = NULL  # layer selection already handled above
    ))
  }
}


