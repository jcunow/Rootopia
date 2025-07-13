#' Combine multi-dimensional arrays
#'
#' @description
#'Combine multi-dimensional arrays.  This is a
#'generalization of cbind and rbind.  Takes a sequence of
#'vectors, matrices, or arrays and produces a single array of
#'the same or higher dimension.
#'
#' @param ...  Any number of vectors, matrices, arrays, or data frames.
#'The dimensions of all the arrays must match, except on one dimension
#'(specified by \code{along=}).  If these arguments are named, the name
#'will be used for the name of the dimension along which the arrays are
#'joined.  Vectors are treated as having a dim attribute of length one.
#'
#'Alternatively, there can be one (and only one) list argument supplied,
#'whose components are the objects to be bound together.  Names of the
#'list components are treated in the same way as argument names.
#' @param along (optional) The dimension along which to bind the arrays.
#'The default is the last dimension, i.e., the maximum length of the dim
#'attribute of the supplied arrays.  \code{along=} can take any
#'non-negative value up to the minimum length of the dim attribute of
#'supplied arrays plus one.  When \code{along=} has a fractional value, a
#'value less than 1, or a value greater than N (N is the maximum of the
#'lengths of the dim attribute of the objects to be bound together), a new
#'dimension is created in the result.  In these cases, the dimensions of
#'all arguments must be identical.
#' @param rev.along (optional)
#'Alternate way to specify the dimension along which to bind the arrays:
#'  \code{along = N + 1 - rev.along}.  This is provided mainly to allow easy
#'specification of \code{along = N + 1} (by supplying
#'                                       \code{rev.along=0}).  If both \code{along} and \code{rev.along} are
#'supplied, the supplied value of \code{along} is ignored.
#' @param new.names (optional)
#'If new.names is a list, it is the first choice for the
#'dimnames attribute of the result.  It should have the same
#'structure as a dimnames attribute.  If the names for a
#'particular dimension are \code{NULL}, names for this dimension are
#'constructed in other ways.
#'
#'If \code{new.names} is a character vector, it is used for dimension
#'names in the same way as argument names are used.  Zero
#'length ("") names are ignored.
#' @param force.array (optional) If \code{FALSE}, rbind or cbind are
#'called when possible, i.e., when the arguments are all vectors, and
#'along is not 1, or when the arguments are vectors or matrices or data
#'frames and along is 1 or 2.  If rbind or cbind are used, they will
#'preserve the data.frame classes (or any other class that r/cbind
#'preserve).  Otherwise, abind will convert objects to class array.  Thus,
#'to guarantee that an array object is returned, supply the argument
#'\code{force.array=TRUE}.  Note that the use of rbind or cbind introduces
#'some subtle changes in the way default dimension names are constructed:
#'  see the examples below.
#' @param make.names (optional)
#'If \code{TRUE}, the last resort for dimnames for the along
#'dimension will be the deparsed versions of anonymous
#'arguments.  This can result in cumbersome names when
#'arguments are expressions.
#'
#'<p>The default is \code{FALSE}.
#' @param use.anon.names (optional)
#'\code{use.anon.names}
#'is a deprecated synonym for \code{make.names}.
#' @param use.first.dimnames (optional)
#'When dimension names are present on more than one
#'argument, should dimension names for the result be take from
#'the first available (the default is to take them from the
#'                     last available, which is the same behavior as
#'                     \code{rbind} and \code{cbind}.)
#' @param hier.names (optional)
#'If \code{TRUE}, dimension names on the concatenated dimension will be
#'composed of the argument name and the dimension names of the objects
#'being bound.  If a single list argument is supplied, then the names of
#'the components serve as the argument names.  \code{hier.names} can
#'also have values \code{"before"} or \code{"after"}; these determine
#'the order in which the argument name and the dimension name are put
#'together (\code{TRUE} has the same effect as \code{"before"}).
#' @param use.dnns (default \code{FALSE}) Use names on dimensions, e.g.,
#'so that \code{names(dimnames(x))} is non-empty.  When there are
#'multiple possible sources for names of dimnames, the value of
#'\code{use.first.dimnames} determines the result.
#'
#' @details
#' The dimensions of the supplied vectors or arrays do not need
#'to be identical, e.g., arguments can be a mixture of vectors
#'and matrices.  \code{abind} coerces arguments by the addition
#'of one dimension in order to make them consistent with other
#'arguments and \code{along=}.  The extra dimension is
#'added in the place specified by \code{along=}.
#'
#'The default action of abind is to concatenate on the last
#'dimension, rather than increase the number of dimensions.
#'For example, the result of calling abind with vectors is a
#'longer vector (see first example below).  This differs from
#'the action of \code{rbind} and cbind which is to return a matrix when
#'called with vectors.  abind can be made to behave like cbind
#'on vectors by specifying \code{along=2}, and like rbind by
#'specifying \code{along=0}.
#'
#'The dimnames of the returned object are pieced together
#'from the dimnames of the arguments, and the names of the
#'arguments.  Names for each dimension are searched for in the
#'following order: new.names, argument name, dimnames (or
#'names) attribute of last argument, dimnames (or names)
#'attribute of second last argument, etc.  (Supplying the
#'                                          argument \code{use.first.dimnames=TRUE} changes this to
#'                                          cause \code{abind} to use dimnames or names from the
#'                                          first argument first.  The default behavior is the same as
#'                                          for \code{rbind} and \code{cbind}: use dimnames
#'                                          from later arguments.)  If some names are supplied for the
#'along dimension (either as argument names or dimnames in
#'                 arguments), names are constructed for anonymous arguments
#'unless \code{use.anon.names=FALSE}.
#'
#' @author Tony Plate \email{tplate@acm.org} and Richard Heiberger
#' @details
#' sourced from the 'abind' package:  https://doi.org/10.32614/CRAN.package.abind under MIT-license
#'
#'
#' @import utils
#' @return merged multidimensional arrays
#' @keywords internal
abind2 = function (..., along = N, rev.along = NULL, new.names = NULL,
                   force.array = TRUE, make.names = use.anon.names, use.anon.names = FALSE,
                   use.first.dimnames = FALSE, hier.names = FALSE, use.dnns = FALSE)
{
  if (is.character(hier.names))
    hier.names <- match.arg(hier.names, c("before", "after",
                                          "none"))
  else hier.names <- if (hier.names)
    "before"
  else "no"
  arg.list <- list(...)
  if (is.list(arg.list[[1]]) && !is.data.frame(arg.list[[1]])) {
    if (length(arg.list) != 1)
      stop("can only supply one list-valued argument for ...")
    if (make.names)
      stop("cannot have make.names=TRUE with a list argument")
    arg.list <- arg.list[[1]]
    have.list.arg <- TRUE
  }
  else {
    N <- max(1, sapply(list(...), function(x) length(dim(x))))
    have.list.arg <- FALSE
  }
  if (any(discard <- sapply(arg.list, is.null)))
    arg.list <- arg.list[!discard]
  if (length(arg.list) == 0)
    return(NULL)
  N <- max(1, sapply(arg.list, function(x) length(dim(x))))
  if (!is.null(rev.along))
    along <- N + 1 - rev.along
  if (along < 1 || along > N || (along > floor(along) && along <
                                 ceiling(along))) {
    N <- N + 1
    along <- max(1, min(N + 1, ceiling(along)))
  }
  if (length(along) > 1 || along < 1 || along > N + 1)
    stop(paste("\"along\" must specify one dimension of the array,",
               "or interpolate between two dimensions of the array",
               sep = "\n"))
  if (!force.array && N == 2) {
    if (!have.list.arg) {
      if (along == 2)
        return(cbind(...))
      if (along == 1)
        return(rbind(...))
    }
    else {
      if (along == 2)
        return(do.call("cbind", arg.list))
      if (along == 1)
        return(do.call("rbind", arg.list))
    }
  }
  if (along > N || along < 0)
    stop("along must be between 0 and ", N)
  pre <- seq(from = 1, len = along - 1)
  post <- seq(to = N - 1, len = N - along)
  perm <- c(seq(len = N)[-along], along)
  arg.names <- names(arg.list)
  if (is.null(arg.names))
    arg.names <- rep("", length(arg.list))
  if (is.character(new.names)) {
    arg.names[seq(along = new.names)[nchar(new.names) > 0]] <- new.names[nchar(new.names) >
                                                                           0]
    new.names <- NULL
  }
  if (any(arg.names == "")) {
    if (make.names) {
      dot.args <- match.call(expand.dots = FALSE)$...
      if (is.call(dot.args) && identical(dot.args[[1]],
                                         as.name("list")))
        dot.args <- dot.args[-1]
      arg.alt.names <- arg.names
      for (i in seq(along = arg.names)) {
        if (arg.alt.names[i] == "") {
          if (utils::object.size(dot.args[[i]]) < 1000) {
            arg.alt.names[i] <- paste(deparse(dot.args[[i]],
                                              40), collapse = ";")
          }
          else {
            arg.alt.names[i] <- paste("X", i, sep = "")
          }
          arg.names[i] <- arg.alt.names[i]
        }
      }
    }
    else {
      arg.alt.names <- arg.names
      arg.alt.names[arg.names == ""] <- paste("X", seq(along = arg.names),
                                              sep = "")[arg.names == ""]
    }
  }
  else {
    arg.alt.names <- arg.names
  }
  use.along.names <- any(arg.names != "")
  names(arg.list) <- arg.names
  arg.dimnames <- matrix(vector("list", N * length(arg.names)),
                         nrow = N, ncol = length(arg.names))
  dimnames(arg.dimnames) <- list(NULL, arg.names)
  arg.dnns <- matrix(vector("list", N * length(arg.names)),
                     nrow = N, ncol = length(arg.names))
  dimnames(arg.dnns) <- list(NULL, arg.names)
  dimnames.new <- vector("list", N)
  arg.dim <- matrix(integer(1), nrow = N, ncol = length(arg.names))
  for (i in seq(len = length(arg.list))) {
    m <- arg.list[[i]]
    m.changed <- FALSE
    if (is.data.frame(m)) {
      m <- as.matrix(m)
      m.changed <- TRUE
    }
    else if (!is.array(m) && !is.null(m)) {
      if (!is.atomic(m))
        stop("arg '", arg.alt.names[i], "' is non-atomic")
      dn <- names(m)
      m <- as.array(m)
      if (length(dim(m)) == 1 && !is.null(dn))
        dimnames(m) <- list(dn)
      m.changed <- TRUE
    }
    new.dim <- dim(m)
    if (length(new.dim) == N) {
      if (!is.null(dimnames(m))) {
        arg.dimnames[, i] <- dimnames(m)
        if (use.dnns && !is.null(names(dimnames(m))))
          arg.dnns[, i] <- as.list(names(dimnames(m)))
      }
      arg.dim[, i] <- new.dim
    }
    else if (length(new.dim) == N - 1) {
      if (!is.null(dimnames(m))) {
        arg.dimnames[-along, i] <- dimnames(m)
        if (use.dnns && !is.null(names(dimnames(m))))
          arg.dnns[-along, i] <- as.list(names(dimnames(m)))
        dimnames(m) <- NULL
      }
      arg.dim[, i] <- c(new.dim[pre], 1, new.dim[post])
      if (any(perm != seq(along = perm))) {
        dim(m) <- c(new.dim[pre], 1, new.dim[post])
        m.changed <- TRUE
      }
    }
    else {
      stop("'", arg.alt.names[i], "' does not fit: should have `length(dim())'=",
           N, " or ", N - 1)
    }
    if (any(perm != seq(along = perm)))
      arg.list[[i]] <- aperm(m, perm)
    else if (m.changed)
      arg.list[[i]] <- m
  }
  conform.dim <- arg.dim[, 1]
  for (i in seq(len = ncol(arg.dim))) {
    if (any((conform.dim != arg.dim[, i])[-along])) {
      stop("arg '", arg.alt.names[i], "' has dims=", paste(arg.dim[,
                                                                   i], collapse = ", "), "; but need dims=", paste(replace(conform.dim,
                                                                                                                           along, "X"), collapse = ", "))
    }
  }
  if (N > 1)
    for (dd in seq(len = N)[-along]) {
      for (i in (if (use.first.dimnames)
        seq(along = arg.names)
        else rev(seq(along = arg.names)))) {
        if (length(arg.dimnames[[dd, i]]) > 0) {
          dimnames.new[[dd]] <- arg.dimnames[[dd, i]]
          if (use.dnns && !is.null(arg.dnns[[dd, i]]))
            names(dimnames.new)[dd] <- arg.dnns[[dd,
                                                 i]]
          break
        }
      }
    }
  for (i in seq(len = length(arg.names))) {
    if (arg.dim[along, i] > 0) {
      dnm.along <- arg.dimnames[[along, i]]
      if (length(dnm.along) == arg.dim[along, i]) {
        use.along.names <- TRUE
        if (hier.names == "before" && arg.names[i] !=
            "")
          dnm.along <- paste(arg.names[i], dnm.along,
                             sep = ".")
        else if (hier.names == "after" && arg.names[i] !=
                 "")
          dnm.along <- paste(dnm.along, arg.names[i],
                             sep = ".")
      }
      else {
        if (arg.dim[along, i] == 1)
          dnm.along <- arg.names[i]
        else if (arg.names[i] == "")
          dnm.along <- rep("", arg.dim[along, i])
        else dnm.along <- paste(arg.names[i], seq(length = arg.dim[along,
                                                                   i]), sep = "")
      }
      dimnames.new[[along]] <- c(dimnames.new[[along]],
                                 dnm.along)
    }
    if (use.dnns) {
      dnn <- unlist(arg.dnns[along, ])
      if (length(dnn)) {
        if (!use.first.dimnames)
          dnn <- rev(dnn)
        names(dimnames.new)[along] <- dnn[1]
      }
    }
  }
  if (!use.along.names)
    dimnames.new[along] <- list(NULL)
  out <- array(unlist(arg.list, use.names = FALSE), dim = c(arg.dim[-along,
                                                                    1], sum(arg.dim[along, ])), dimnames = dimnames.new[perm])
  if (any(order(perm) != seq(along = perm)))
    out <- aperm(out, order(perm))
  if (!is.null(new.names) && is.list(new.names)) {
    for (dd in seq(len = N)) {
      if (!is.null(new.names[[dd]])) {
        if (length(new.names[[dd]]) == dim(out)[dd])
          dimnames(out)[[dd]] <- new.names[[dd]]
        else if (length(new.names[[dd]]))
          warning(paste("Component ", dd, " of new.names ignored: has length ",
                        length(new.names[[dd]]), ", should be ",
                        dim(out)[dd], sep = ""))
      }
      if (use.dnns && !is.null(names(new.names)) && names(new.names)[dd] !=
          "")
        names(dimnames(out))[dd] <- names(new.names)[dd]
    }
  }
  if (use.dnns && !is.null(names(dimnames(out))) && any(i <- is.na(names(dimnames(out)))))
    names(dimnames(out))[i] <- ""
  out
}




#' Calculate a circular mean to determine average Directionality
#'
#' @param angles Numeric vector of input angles
#' @param input_units Character string specifying input units ("radians" or "degrees")
#' @param output_units Character string specifying output units ("radians" or "degrees")
#'
#' @return Numeric value representing the average angle
#' @export
#'
#' @examples circular_mean(angles = c(360,90,0), input_units = "degrees", output_units = "degrees")
circular_mean <- function(angles, input_units = "degrees", output_units = "degrees") {
  # Input validation
  tryCatch({
    if (missing(angles)) {
      stop("angles parameter is required")
    }

    if (!is.numeric(angles)) {
      stop("angles must be numeric")
    }

    if (length(angles) == 0) {
      stop("angles vector is empty")
    }

    # Validate units parameters
    valid_units <- c("degrees", "radians")
    if (!input_units %in% valid_units) {
      stop("input_units must be 'degrees' or 'radians'")
    }
    if (!output_units %in% valid_units) {
      stop("output_units must be 'degrees' or 'radians'")
    }

    # Handle NA/Inf values
    if (any(is.na(angles)) || any(is.infinite(angles))) {
      warning("Removing NA and infinite values from angles")
      angles <- angles[!is.na(angles) & !is.infinite(angles)]
      if (length(angles) == 0) {
        stop("No valid angles remain after removing NA/infinite values")
      }
    }

    # Convert angles to radians if they are in degrees
    if (input_units == "degrees") {
      angles <- angles * pi / 180
    }

    # Calculate the sine and cosine of the angles
    sin_sum <- sum(sin(angles))
    cos_sum <- sum(cos(angles))

    # Check for zero denominator
    if (abs(sin_sum) < .Machine$double.eps && abs(cos_sum) < .Machine$double.eps) {
      warning("Near-zero sums detected, result may be unstable")
    }

    # Calculate the circular mean
    mean_angle <- atan2(sin_sum, cos_sum)

    # Ensure the mean angle is in the range [0, 2*pi)
    if (mean_angle < 0) {
      mean_angle <- mean_angle + 2 * pi
    }

    # Convert the mean angle to the desired output unit
    if (output_units == "degrees") {
      mean_angle <- mean_angle * 180 / pi
    }

    return(mean_angle)

  }, error = function(e) {
    stop(sprintf("Error in circular_mean: %s", e$message))
  })
}



#' Threshold or deblur an image to binarize features
#'
#' Applies global or adaptive thresholding to a raster image, with optional
#' masking and simple deblurring based on structural layer separation.
#'
#' @param img SpatRaster object
#' @param threshold Numeric value (0–1). For global thresholding, it's a fraction of the global max.
#'                  For adaptive thresholding, it's a fraction of local mean. For deblurring, it's the fraction of max used to recover structure.
#' @param method "global" or "adaptive". Ignored if `deblur = TRUE`.
#' @param window_size Integer (odd), only used for adaptive thresholding.
#' @param select.layer Integer or NULL. Which layer to use for thresholding or deblurring.
#'                     If NULL and multilayer, the mean of all layers is used.
#' @param mask.layer Integer or NULL. If set, used to preserve masked regions or enhance structure.
#' @param binary_01 Logical. If TRUE, binarized output uses 0/1. If FALSE, it retains max value of input.
#' @param deblur Logical. If TRUE, applies deblurring logic instead of standard thresholding.
#'
#' @return SpatRaster object
#' @export
#'
#' @examples
#' img = terra::rast(seg_Oulanka2023_Session03_T067)
#' image_threshold(img, threshold = 0.3, method = "global")
#' image_threshold(img, threshold = 0.9, method = "adaptive", window_size = 15, binary_01 = TRUE)
#' image_threshold(img, threshold = 0.4, select.layer = 2, mask.layer = 1, deblur = TRUE)
image_threshold <- function(img, threshold = 0.4, method = "global", window_size = 15,
                            select.layer = 2, mask.layer = 1, binary_01 = FALSE,
                            deblur = FALSE) {
  tryCatch({
    if (missing(img)) stop("img parameter is required")
    if (!inherits(img, "SpatRaster")) stop("img must be a SpatRaster object")
    if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
      stop("threshold must be a numeric value between 0 and 1")
    }
    
    nlyr <- terra::nlyr(img)
    if (all(is.na(terra::values(img)))) stop("Input raster contains only NA values")
    
    img2 <- img
    
    # ---- DEBLURRING MODE ----
    if (deblur) {
      if (nlyr < 2) stop("Deblurring requires a multi-layer image")
      if (select.layer == mask.layer) stop("select.layer and mask.layer cannot be the same")
      
      mx <- terra::global(img[[select.layer]], "max", na.rm = TRUE)[[1]]
      if (is.na(mx) || mx == 0) stop("Invalid maximum in select.layer")
      
      # Deblurring: isolate structure from blurry background
      img2[[mask.layer]] <- img[[mask.layer]] - img[[select.layer]]
      img2[[mask.layer]] <- (img2[[mask.layer]] >= (mx * threshold)) * mx
      img2[[select.layer]] <- (img[[select.layer]] >= (mx * threshold)) * mx
      
      other.layer <- setdiff(1:nlyr, c(select.layer, mask.layer))
      if (length(other.layer) > 0) {
        img2[[other.layer]] <- (img[[select.layer]] >= (mx * threshold)) * mx
      }
      
      return(img2)
    }
    
    # ---- THRESHOLDING MODE ----
    # Determine processing layer
    if (is.null(select.layer) && nlyr > 1) {
      layer_to_process <- terra::app(img, fun = mean, na.rm = TRUE)
    } else if (!is.null(select.layer) && select.layer <= nlyr) {
      layer_to_process <- img[[select.layer]]
    } else if (nlyr == 1) {
      layer_to_process <- img
    } else {
      stop("Invalid select.layer for given image")
    }
    
    mx <- terra::global(layer_to_process, "max", na.rm = TRUE)[[1]]
    if (is.na(mx) || mx == 0) stop("Invalid maximum in processing layer")
    output_value <- ifelse(binary_01, 1, mx)
    
    if (tolower(method) == "global") {
      threshold_mask <- layer_to_process >= (mx * threshold)
    } else if (tolower(method) == "adaptive") {
      if (!is.numeric(window_size) || window_size < 3 || window_size %% 2 == 0) {
        stop("window_size must be an odd integer >= 3")
      }
      w <- matrix(1, nrow = window_size, ncol = window_size)
      local_mean <- terra::focal(layer_to_process, w = w, fun = mean, na.rm = TRUE)
      threshold_mask <- layer_to_process >= (local_mean * threshold)
    } else {
      stop("method must be 'global' or 'adaptive'")
    }
    
    # Apply threshold to layers
    if (nlyr > 1) {
      for (i in 1:nlyr) {
        if (!is.null(mask.layer) && i == mask.layer && (!is.null(select.layer) && mask.layer != select.layer)) {
          mask_layer_data <- img[[mask.layer]]
          if (tolower(method) == "global") {
            img2[[i]] <- (mask_layer_data >= (mx * threshold)) * output_value
          } else {
            mask_mean <- terra::focal(mask_layer_data, w = w, fun = mean, na.rm = TRUE)
            img2[[i]] <- (mask_layer_data >= (mask_mean * threshold)) * output_value
          }
        } else {
          img2[[i]] <- threshold_mask * output_value
        }
      }
    } else {
      img2 <- threshold_mask * output_value
    }
    
    return(img2)
  }, error = function(e) {
    stop(sprintf("Error in image_threshold: %s", e$message))
  })
}




#' Calculate root accumulation
#'
#' @param x Data frame containing group, depth, and variable columns
#' @param group Character vector specifying grouping variable(s)
#' @param depth Character string specifying depth column name
#' @param variable Character string specifying accumulating values column
#' @param stdrz Character string specifying standardization method
#'
#' @return Numeric vector of accumulated values
#'
#' @examples
#'df = data.frame(depth = c(seq(0,80,20),seq(0,80,20)),
#'                Plot = c(rep("a",5),rep("b",5)), rootpx = c(5,50,20,15,5,10,40,30,10,5) )
#' accum_root = root_accumulation(df,group = "Plot", depth = "depth", variable = "rootpx")
root_accumulation = function(x, group, depth, variable, stdrz = "counts") {
  tryCatch({
    # Input validation
    if (missing(x) || missing(group) || missing(depth) || missing(variable)) {
      stop("All parameters (x, group, depth, variable) are required")
    }

    if (!is.data.frame(x)) {
      stop("x must be a data frame")
    }

    # Validate column names
    if (!group %in% names(x)) {
      stop(sprintf("group column '%s' not found in data frame", group))
    }
    if (!depth %in% names(x)) {
      stop(sprintf("depth column '%s' not found in data frame", depth))
    }
    if (!variable %in% names(x)) {
      stop(sprintf("variable column '%s' not found in data frame", variable))
    }

    # Validate standardization method
    valid_stdrz <- c("counts", "additive", "relative")
    if (!stdrz %in% valid_stdrz) {
      stop(sprintf("Invalid standardization method. Must be one of: %s",
                   paste(valid_stdrz, collapse = ", ")))
    }

    # Check for empty data frame
    if (nrow(x) == 0) {
      stop("Input data frame is empty")
    }

    # Split data by group
    split_df <- split(x, x[,group])

    # Initialize an empty list to store results
    result_list <- list()

    # Loop over each group
    for (grp in names(split_df)) {
      # Sort the data within the group by depth
      sorted_group <- split_df[[grp]][order(split_df[[grp]][[depth]]), ]

      # Handle NA values in variable column
      if (all(is.na(sorted_group[[variable]]))) {
        warning(sprintf("Group %s contains only NA values", grp))
        cs <- rep(NA, nrow(sorted_group))
      } else {
        # Compute cumulative sum of variable
        cs <- cumsum(dplyr::coalesce(sorted_group[[variable]], 0)) +
          sorted_group[[variable]]*0

        # Standardization
        if (stdrz == "additive") {
          mx.roots <- max(cs, na.rm = TRUE)
          if (mx.roots == 0) {
            warning(sprintf("Zero maximum value in group %s", grp))
            cs <- cs
          } else {
            cs <- cs / mx.roots
          }
        } else if (stdrz == "relative") {
          sm.roots <- sum(cs, na.rm = TRUE)
          if (sm.roots == 0) {
            warning(sprintf("Zero sum in group %s", grp))
            cs <- cs
          } else {
            cs <- cs / sm.roots
          }
        }
      }

      sorted_group$cs <- cs
      result_list[[grp]] <- sorted_group
    }

    # Combine the list back into a single data frame
    result_df <- do.call(rbind, result_list)

    # Reorder the result to match the original row order
    result_df <- result_df[order(rownames(result_df)), ]

    return(result_df$cs)

  }, error = function(e) {
    stop(sprintf("Error in root.accumulation: %s", e$message))
  })
}

#' Convert RGB image to grayscale with optimized memory management and parallel processing
#'
#' @param img SpatRaster RGB image
#' @param r Weight for red channel
#' @param g Weight for green channel
#' @param b Weight for blue channel
#'
#' @examples
#' data(seg_Oulanka2023_Session01_T067)
#' img = seg_Oulanka2023_Session01_T067
#' gray.raster = rgb2gray(img)
rgb2gray = function(img, r = 0.21, g = 0.72, b = 0.07) {

  img <- load_flexible_image(img,normalize = FALSE, output_format = "spatrast", select.layer = NULL)
  tryCatch({
    # Input validation
    if (missing(img)) {
      stop("img parameter is required")
    }

    if (!inherits(img, "SpatRaster")) {
      stop("img must be a SpatRaster object")
    }

    # Validate number of layers
    if (terra::nlyr(img) != 3) {
      stop("Input image must have exactly 3 layers (RGB)")
    }

    # Validate weights
    if (!all(is.numeric(c(r, g, b)))) {
      stop("RGB weights must be numeric")
    }

    if (abs((r + g + b) - 1) > .Machine$double.eps * 100) {
      warning("RGB weights do not sum to 1")
    }

    # Check for NA values
    if (all(is.na(terra::values(img)))) {
      stop("Input image contains only NA values")
    }

    # Convert to grayscale
    gray.im <- img[[1]] * r + img[[2]] * g + img[[3]] * b

    # Validate output
    if (all(is.na(terra::values(gray.im)))) {
      warning("Resulting grayscale image contains only NA values")
    }

    return(gray.im)

  }, error = function(e) {
    stop(sprintf("Error in rgb2gray: %s", e$message))
  })
}




#' Find Modal Peaks in a Density Estimation
#'
#' This function identifies peaks and valleys in the density estimation of the input data.
#' It computes the prominence of each peak and allows the user to filter peaks based on a prominence threshold.
#' The function also provides the option to visualize the data as either a density plot, a raw line plot, or no display.
#'
#' @param x A numeric vector of data for which the modal peaks are to be identified.
#' @param prominence_threshold A numeric value specifying the minimum prominence of peaks to be considered.
#'        If `NULL`, no prominence threshold is applied. Default is 0.005.
#' @param display_type A character string indicating the type of output display. Options are:
#'        - `"density"` for a density plot with peaks and valleys marked (default),
#'        - `"raw"` for a raw line plot with peaks and valleys marked,
#'        - `"none"` for no plot, just returning the peak and valley data.
#'
#' @return A list containing:
#'         - `peak_x`: The x-coordinates of the detected peaks.
#'         - `peak_y`: The y-coordinates (density values) of the detected peaks.
#'         - `valley_x`: The x-coordinates of the detected valleys (or NULL if none).
#'         - `valley_y`: The y-coordinates (density values) of the detected valleys (or NULL if none).
#'         - `prominences`: The prominence of each peak (or NULL if none).
#'         - A message indicating the number of modes in the distribution.
#'
#' @examples
#' # Simulate a bimodal distribution
#' set.seed(123)
#' x <- c(rnorm(500, mean = -2, sd = 1), rnorm(500, mean = 3, sd = 1))
#'
#' # Call the function with density plot (default)
#' results_density <- modal_peaks(x, prominence_threshold = 0.005, display_type = "density")
#'
#' # Call the function with raw line plot
#' results_raw <- modal_peaks(x, prominence_threshold = 0.005, display_type = "raw")
#'
#' # Call the function with no display (just returns results)
#' results_none <- modal_peaks(x, prominence_threshold = 0.005, display_type = "none")
modal_peaks <- function(x, prominence_threshold = 0.005, display_type = "density") {

  # Input validation for display_type
  if (!display_type %in% c("density", "raw", "none")) {
    stop("Invalid 'display_type'. Choose 'density', 'raw', or 'none'.")
  }

  # Estimate density
  dens <- stats::density(x)

  # Find local maxima (peaks)
  peaks <- which(diff(sign(diff(dens$y))) == -2) + 1  # Local maxima

  # Find local minima (valleys)
  valleys <- which(diff(sign(diff(dens$y))) == 2) + 1  # Local minima

  # Compute prominence of peaks
  prominences <- numeric(length(peaks))
  for (i in seq_along(peaks)) {
    # Check if there are valleys before and after the peak, handle edge cases
    left_valley <- if (any(valleys < peaks[i])) max(valleys[valleys < peaks[i]], na.rm = TRUE) else NA
    right_valley <- if (any(valleys > peaks[i])) min(valleys[valleys > peaks[i]], na.rm = TRUE) else NA

    peak_height <- dens$y[peaks[i]]

    # If there are no valleys, we set valley_height to 0 or a very low value
    valley_height <- if (!is.na(left_valley) && !is.na(right_valley)) {
      max(dens$y[left_valley], dens$y[right_valley], na.rm = TRUE)
    } else if (!is.na(left_valley)) {
      dens$y[left_valley]
    } else if (!is.na(right_valley)) {
      dens$y[right_valley]
    } else {
      0  # No valleys, assume lowest value
    }

    prominences[i] <- peak_height - valley_height
  }

  # Apply prominence threshold if provided
  if (!is.null(prominence_threshold)) {
    valid_peaks <- which(prominences >= prominence_threshold)
    peaks <- peaks[valid_peaks]
    prominences <- prominences[valid_peaks]
  }

  # Handle cases where no valleys or prominences are detected
  if (length(peaks) <= 1) {
    valleys <- NULL
    prominences <- NULL
    message(sprintf("The distribution is unimodal (%d peak detected).", length(peaks)))
  } else if (length(peaks) == 2) {
    message("The distribution is bimodal (2 peaks detected).")
  } else {
    message(sprintf("The distribution is multimodal (%d peaks detected).", length(peaks)))
  }

  # Display output based on the type selected
  if (display_type == "density") {
    # Plot the density with the peaks and valleys marked
    graphics::plot(dens, main = "Density with Peaks and Valleys")
    graphics::abline(v = dens$x[peaks], col = "red", lty = 2)  # Mark peaks
    if (!is.null(valleys)) {
      graphics::abline(v = dens$x[valleys], col = "blue", lty = 2)  # Mark valleys
    }
    graphics::legend("topright", legend = c("Peaks", if (!is.null(valleys)) "Valleys" else NULL),
           col = c("red", if (!is.null(valleys)) "blue" else NULL), lty = 2, cex = 0.8)
  } else if (display_type == "raw") {
    # Plot raw data as a line with peaks and valleys marked
    graphics::plot(x, type = "l", main = "Raw Data with Peaks and Valleys", xlab = "Index", ylab = "Value")
    graphics::points(peaks, x[peaks], col = "red", pch = 19)  # Mark peaks in raw data
    if (!is.null(valleys)) {
      graphics::points(valleys, x[valleys], col = "blue", pch = 19)  # Mark valleys in raw data
    }
    graphics::legend("topright", legend = c("Peaks", if (!is.null(valleys)) "Valleys" else NULL),
           col = c("red", if (!is.null(valleys)) "blue" else NULL), pch = 19, cex = 0.8, xpd = TRUE)
  } # No display for "none"

  # Return results as a list (peaks, valleys, and their prominences)
  return(list(peak_x = dens$x[peaks], peak_y = dens$y[peaks],
              valley_x = if (!is.null(valleys)) dens$x[valleys] else NULL,
              valley_y = if (!is.null(valleys)) dens$y[valleys] else NULL,
              prominences = prominences))
}
